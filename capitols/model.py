import heapq
import logging
import time

from ortools.sat.python import cp_model

from .grid import monotone_paths

log = logging.getLogger("capitols")


class CapitolsModel:
    """The CP-SAT model of a Problem. One method per constraint family, called in dependency order."""

    def __init__(self, prob, cuts=(), exclude=()):
        self.p = prob
        self.grid = prob.grid
        self.cfg = prob.cfg
        self.clues = prob.clues
        self.bd = prob.bounds
        self.n = prob.grid.n
        self.adj = prob.grid.adj
        self.edges = prob.grid.edges()
        self.m = cp_model.CpModel()
        self._t = prob.t_start
        for i in range(self.n):
            if self.bd.sz_lb[i] > self.bd.sz_ub[i] or self.bd.d_lb[i] > self.bd.d_ub[i]:
                self.m.add_bool_or([])
        self._add_labels()
        self._add_same_edges()
        self._add_connectivity()
        self._add_sizes()
        self._add_symmetries()
        self._add_capitols()
        self._add_distances()
        self._add_clue_distance_cuts()
        self._add_cuts(cuts)
        self._add_exclusions(exclude)
        self._add_search_hints()

    def _stage(self, name):
        proto = self.m.Proto()
        log.info(
            "[build] %-22s done in %6.2fs | vars=%7d constraints=%7d",
            name,
            time.time() - self._t,
            len(proto.variables),
            len(proto.constraints),
        )
        self._t = time.time()

    # -- partition into states -------------------------------------------------------------------

    def _add_labels(self):
        """lab[i] is the lowest-numbered cell of i's state (its root); b[i, r] is the one-hot form."""
        m, n = self.m, self.n
        self.lab = [m.new_int_var(0, i, f"lab{i}") for i in range(n)]
        self.b = {}
        for i in range(n):
            bs = []
            for r in range(i + 1):
                self.b[i, r] = m.new_bool_var("")
                if (i, r) not in self.p.regions.possible_b:
                    m.add(self.b[i, r] == 0)
                bs.append(self.b[i, r])
            m.add_exactly_one(bs)
            m.add(self.lab[i] == sum(r * self.b[i, r] for r in range(i + 1)))
            for r in range(i):
                m.add_implication(self.b[i, r], self.b[r, r])
        self._stage("labels")

    def _channel_same(self, s, u, v):
        """s <-> (u and v share a root), as clauses over b."""
        m, b = self.m, self.b
        for r in range(u + 1):
            m.add_bool_or([s.Not(), b[u, r].Not(), b[v, r]])
            m.add_bool_or([s.Not(), b[v, r].Not(), b[u, r]])
            m.add_bool_or([b[u, r].Not(), b[v, r].Not(), s])
        for r in range(u + 1, v + 1):
            m.add_implication(s, b[v, r].Not())

    def _add_same_edges(self):
        m, lab = self.m, self.lab
        self.same = {}
        for i, j in self.edges:
            s = m.new_bool_var("")
            if not self.cfg.lab_channel:
                m.add(lab[i] == lab[j]).OnlyEnforceIf(s)
                m.add(lab[i] != lab[j]).OnlyEnforceIf(s.Not())
            else:
                self._channel_same(s, i, j)
            self.same[i, j] = self.same[j, i] = s
        self._stage("same-state edges")

    def _add_connectivity(self):
        """Each non-root cell has a same-state neighbour with smaller h, so a state is a tree off its root."""
        m, b, same = self.m, self.b, self.same
        h = [m.new_int_var(0, self.n, "") for _ in range(self.n)]
        for i in range(self.n):
            m.add(h[i] == 0).OnlyEnforceIf(b[i, i])
            m.add(h[i] >= 1).OnlyEnforceIf(b[i, i].Not())
            par = []
            for j in self.adj[i]:
                p = m.new_bool_var("")
                m.add_implication(p, same[i, j])
                m.add(h[j] < h[i]).OnlyEnforceIf(p)
                par.append(p)
            m.add_bool_or(par + [b[i, i]])
        self._stage("connectivity")

    # -- sizes -------------------------------------------------------------------------------------

    def _add_sizes(self):
        """sz is counted at the state's root, then carried along same-state edges (states are connected,
        so each cell reaches its root through such edges). max_size caps every state's size."""
        m, bd, b, same, adj = self.m, self.bd, self.b, self.same, self.adj
        n = self.n
        self.sz = [
            m.new_int_var(bd.sz_lb[i], max(bd.sz_lb[i], bd.sz_ub[i]), "")
            for i in range(n)
        ]
        self.is_sing = [m.new_bool_var(f"sing{i}") for i in range(n)]
        sz, is_sing = self.sz, self.is_sing
        for i in range(n):
            m.add(sz[i] == 1).OnlyEnforceIf(is_sing[i])
            m.add(sz[i] >= 2).OnlyEnforceIf(is_sing[i].Not())
            m.add_implication(is_sing[i], b[i, i])
            for j in adj[i]:
                m.add_implication(is_sing[i], same[i, j].Not())
            m.add_bool_or([is_sing[i]] + [same[i, j] for j in adj[i]])
        if self.cfg.clue_cuts:
            self._add_clue_size_bounds()
        for r in range(n):
            m.add(sz[r] == sum(b[i, r] for i in range(r, n))).OnlyEnforceIf(b[r, r])
        for i, j in self.edges:
            m.add(sz[i] == sz[j]).OnlyEnforceIf(same[i, j])
        self._stage("state sizes")

    def _sz_ge(self, j, k):
        """Cached literal for sz[j] >= k, so paths sharing a cell don't allocate duplicates."""
        if (j, k) not in self._sz_ge_cache:
            bv = self.m.new_bool_var("")
            self.m.add(self.sz[j] >= k).OnlyEnforceIf(bv)
            self.m.add(self.sz[j] <= k - 1).OnlyEnforceIf(bv.Not())
            self._sz_ge_cache[j, k] = bv
        return self._sz_ge_cache[j, k]

    def _add_clue_size_bounds(self):
        """For clues c1, c2 up to 6 steps apart, distances differ by at most the sum of the edge weights
        on any path between them, which bounds the sizes of the cells on it."""
        m, sz, bd, idx = self.m, self.sz, self.bd, self.grid.idx
        self._sz_ge_cache = {}
        clue_list = list(self.clues.items())
        for idx_a in range(len(clue_list)):
            c1, v1 = clue_list[idx_a]
            for idx_b in range(idx_a + 1, len(clue_list)):
                c2, v2 = clue_list[idx_b]
                dist = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
                gap = abs(v1 - v2)
                if not (2 <= dist <= 6 and gap > 1):
                    continue
                u0, uk = idx[c1], idx[c2]
                if dist == 2 and bd.sz_lb[u0] + bd.sz_lb[uk] < gap:
                    m.add(sz[u0] + sz[uk] >= gap)
                sum_sz = (gap + 1) // 2
                min_max_sz = (gap + dist - 1) // dist
                for path in monotone_paths(c1, c2):
                    mid_cells = [idx[step[1]] for step in path[:-1]]
                    if dist == 2:
                        u1 = mid_cells[0]
                        if bd.sz_lb[u0] + bd.sz_lb[u1] < gap:
                            m.add(sz[u0] + sz[u1] >= gap)
                        if bd.sz_lb[u1] + bd.sz_lb[uk] < gap:
                            m.add(sz[u1] + sz[uk] >= gap)
                    if sum_sz > 1 and sum(bd.sz_lb[j] for j in mid_cells) < sum_sz:
                        m.add(sum(sz[j] for j in mid_cells) >= sum_sz)
                    if min_max_sz > 1 and all(
                        bd.sz_lb[j] < min_max_sz for j in mid_cells
                    ):
                        if len(mid_cells) == 1:
                            m.add(sz[mid_cells[0]] >= min_max_sz)
                        else:
                            m.add_bool_or(
                                [self._sz_ge(j, min_max_sz) for j in mid_cells]
                            )

    # -- symmetry ------------------------------------------------------------------------------------

    def _add_symmetries(self):
        """Every transform in the model is an involution (i -> j means j -> i). If i uses a transform,
        j is in its state, and a connected state shares one transform, so j uses it too: i and j share a
        single variable T, and the pair needs only one lab equality."""
        self._add_transform_choice()
        if self.cfg.tclauses:
            self._add_shared_transform_clauses()
        else:
            self._add_shared_transform_ints()
        if self.cfg.edge_images:
            self._add_edge_image_clauses()
        if self.cfg.require_free_r2:
            self._add_require_free_r2()
        if self.cfg.min_max_size:
            self._add_min_max_size(self.cfg.min_max_size)
        self._stage("symmetry")

    def _add_transform_choice(self):
        m, bd, mt, reg = self.m, self.bd, self.p.mt, self.p.regions
        n = self.n
        self.T = {}
        T = self.T
        for gi, img in enumerate(mt.images):
            for i in range(n):
                j = img[i]
                if j is None or i > j:
                    continue
                assert img[j] == i
                best_s = 0
                for r in range(i + 1):
                    region = reg.regions.get((r, gi))
                    if region and i in region and j in region:
                        best_s = max(best_s, reg.max_s[r, gi])
                if best_s < max(bd.sz_lb[i], bd.sz_lb[j], 2 if i < j else 1):
                    continue
                v = m.new_bool_var("")
                T[i, gi] = T[j, gi] = v
                if i < j:
                    m.add(self.lab[i] == self.lab[j]).OnlyEnforceIf(v)
                    m.add_implication(v, self.b[j, j].Not())
        for i in range(n):
            m.add_exactly_one([T[i, gi] for gi in range(len(mt)) if (i, gi) in T])

    def _add_shared_transform_ints(self):
        m, T, same = self.m, self.T, self.same
        G = len(self.p.mt)
        tvar = [m.new_int_var(0, G - 1, "") for _ in range(self.n)]
        for i in range(self.n):
            opts_idx = [(gi, T[i, gi]) for gi in range(G) if (i, gi) in T]
            m.add(tvar[i] == sum(gi * v for gi, v in opts_idx))
        for i, j in self.edges:
            m.add(tvar[i] == tvar[j]).OnlyEnforceIf(same[i, j])

    def _add_shared_transform_clauses(self):
        """Neighbours in one state use the same transform."""
        m, T, same = self.m, self.T, self.same
        for i, j in self.edges:
            for gi in range(len(self.p.mt)):
                in_i = (i, gi) in T
                in_j = (j, gi) in T
                if in_i and in_j:
                    vi, vj = T[i, gi], T[j, gi]
                    if vi is not vj:
                        m.add_bool_or([same[i, j].Not(), vi.Not(), vj])
                        m.add_bool_or([same[i, j].Not(), vj.Not(), vi])
                elif in_i:
                    m.add_bool_or([same[i, j].Not(), T[i, gi].Not()])
                elif in_j:
                    m.add_bool_or([same[i, j].Not(), T[j, gi].Not()])

    def _add_edge_image_clauses(self):
        """If a state uses g, an edge inside it maps to an edge inside it."""
        m, T, same, adj = self.m, self.T, self.same, self.adj
        seen = set()
        for gi, img in enumerate(self.p.mt.images):
            for i in range(self.n):
                if (i, gi) not in T:
                    continue
                u = img[i]
                if u is None or i > u:
                    continue
                var_t = T[i, gi]
                t_id = id(var_t)
                for j in adj[i]:
                    v = img[j]
                    if v is None or (j, gi) not in T:
                        m.add_bool_or([var_t.Not(), same[i, j].Not()])
                    else:
                        e1 = (min(i, j), max(i, j))
                        e2 = (min(u, v), max(u, v))
                        if e1 != e2:
                            key = (t_id, min(e1, e2), max(e1, e2))
                            if key not in seen:
                                seen.add(key)
                                m.add_bool_or([var_t.Not(), same[e1].Not(), same[e2]])
                                m.add_bool_or([var_t.Not(), same[e1], same[e2].Not()])

    def _add_require_free_r2(self):
        mt, T = self.p.mt, self.T
        chosen, seen = [], set()
        for gi, g in enumerate(mt.G):
            if g.kind == "R2" and not mt.fixed_at[gi]:
                for i in range(self.n):
                    if (i, gi) in T and id(T[i, gi]) not in seen:
                        seen.add(id(T[i, gi]))
                        chosen.append(T[i, gi])
        self.m.add_bool_or(chosen)

    def _add_min_max_size(self, M):
        """Some state has at least M cells (use only a proven lower bound, or it is a restriction like
        max_size). Its root has that many cells at or after it, and its transform needs that many cells
        in its valid root-region."""
        m, T, sz, reg = self.m, self.T, self.sz, self.p.regions
        G = len(self.p.mt)
        big = []
        for r in range(self.n - M + 1):
            if not any(
                reg.max_s.get((r, gi), 0) >= M for gi in range(G) if (r, gi) in T
            ):
                continue
            v = m.new_bool_var(f"big{r}")
            m.add_implication(v, self.b[r, r])
            m.add(sz[r] >= M).OnlyEnforceIf(v)
            for gi in range(G):
                if (r, gi) in T and reg.max_s.get((r, gi), 0) < M:
                    m.add_implication(v, T[r, gi].Not())
            big.append(v)
        m.add_bool_or(big)

    # -- capitols ------------------------------------------------------------------------------------

    def _eq(self, i, j):
        """Literal for 'i and j are in the same state', reusing the adjacent-pair variables."""
        if (i, j) in self.same:
            return self.same[i, j]
        key = (min(i, j), max(i, j))
        if key not in self._eqcache:
            u, v = key
            s = self.m.new_bool_var("")
            if not self.cfg.eq_channel:
                self.m.add(self.lab[u] == self.lab[v]).OnlyEnforceIf(s)
                self.m.add(self.lab[u] != self.lab[v]).OnlyEnforceIf(s.Not())
            else:
                self._channel_same(s, u, v)
            self._eqcache[key] = s
        return self._eqcache[key]

    def _add_capitols(self):
        m, bd, adj, T = self.m, self.bd, self.adj, self.T
        mt, finder = self.p.mt, self.p.finder
        n, clue_cuts = self.n, self.cfg.clue_cuts
        self._eqcache = {}
        self.cap = [m.new_bool_var(f"cap{i}") for i in range(n)]
        cap, is_sing, sz = self.cap, self.is_sing, self.sz
        for i in range(n):
            if not bd.cap_possible[i]:
                m.add(cap[i] == 0)
            elif bd.d_ub[i] == 0:
                m.add(cap[i] == 1)
        for i in range(n):
            m.add_implication(is_sing[i], cap[i])
        for i in range(n):
            ors = []
            c_sz_lb = bd.capitol_size_lb(i, adj, clue_cuts)
            for gi, img in enumerate(mt.images):
                if (
                    img[i] != i or (i, gi) not in T
                ):  # only transforms fixing i can make it a capitol
                    continue
                r_cap = finder.max_region(i, gi, forbid_other_fixed=True)
                cap_ub = len(r_cap) if len(r_cap) % 2 == 1 else len(r_cap) - 1
                others = [j for j in mt.fixed_at[gi] if j != i]
                if not bd.cap_possible[i] or cap_ub < c_sz_lb:
                    # i cannot be this transform's capitol, so choosing it needs another fixed cell in the state
                    if not others:
                        m.add(T[i, gi] == 0)
                    else:
                        m.add_bool_or(
                            [T[i, gi].Not()] + [self._eq(i, j) for j in others]
                        )
                    continue
                if (
                    not others
                ):  # nothing else to share a state with: T itself says i is the capitol
                    ors.append(T[i, gi])
                    continue
                a = m.new_bool_var("")
                m.add_implication(a, T[i, gi])
                for j in others:
                    m.add_implication(a, self._eq(i, j).Not())
                m.add_bool_or([a, T[i, gi].Not()] + [self._eq(i, j) for j in others])
                ors.append(a)
            m.add_bool_or(ors).OnlyEnforceIf(cap[i])
            for a in ors:
                m.add_implication(a, cap[i])
            if not ors:
                m.add(cap[i] == 0)
        # A capitol's state is mapped to itself by an involution that fixes only the capitol, so the
        # remaining cells pair up and the size is odd. Implied, but it links cap and sz for propagation.
        half = [m.new_int_var(0, n // 2, "") for _ in range(n)]
        for i in range(n):
            m.add(sz[i] == 2 * half[i] + 1).OnlyEnforceIf(cap[i])
            # a one-cell state is symmetric about itself and so is always its own capitol
            m.add(sz[i] >= 2).OnlyEnforceIf(cap[i].Not())
            if clue_cuts and bd.cap_possible[i]:
                c_sz_lb = bd.capitol_size_lb(i, adj)
                if c_sz_lb > bd.sz_lb[i]:
                    m.add(sz[i] >= c_sz_lb).OnlyEnforceIf(cap[i])
                for nb in adj[i]:
                    if bd.d_lb[nb] > bd.sz_lb[nb]:
                        m.add(sz[nb] >= bd.d_lb[nb]).OnlyEnforceIf(cap[i])
        self._stage("capitols")

    # -- distances -----------------------------------------------------------------------------------

    def _add_distances(self):
        """d is the shortest-path distance to a capitol: d = 0 at capitols, d[i] <= d[j] + w on every
        edge, and every non-capitol has a tight predecessor."""
        m, bd, adj, sz, cap = self.m, self.bd, self.adj, self.sz, self.cap
        n = self.n
        self.d = [
            m.new_int_var(bd.d_lb[i], max(bd.d_lb[i], bd.d_ub[i]), f"d{i}")
            for i in range(n)
        ]
        d = self.d
        self.w = {}
        w = self.w
        for i, j in self.edges:
            w_l = min(bd.sz_lb[i], bd.sz_lb[j])
            w_u = max(w_l, min(bd.sz_ub[i], bd.sz_ub[j]))
            v = m.new_int_var(w_l, w_u, "")
            m.add_min_equality(v, [sz[i], sz[j]])
            w[i, j] = w[j, i] = v
        for i in range(n):
            m.add(d[i] == 0).OnlyEnforceIf(cap[i])
            m.add(d[i] >= 1).OnlyEnforceIf(cap[i].Not())
            is_known_cap = bd.d_ub[i] == 0
            par = []
            for j in adj[i]:
                m.add(d[i] <= d[j] + w[i, j])
                if (
                    not is_known_cap
                    and bd.d_lb[j] + min(bd.sz_lb[i], bd.sz_lb[j]) <= bd.d_ub[i]
                ):
                    p = m.new_bool_var("")
                    m.add(d[i] >= d[j] + w[i, j]).OnlyEnforceIf(p)
                    par.append(p)
            if not is_known_cap:
                m.add_bool_or(par + [cap[i]])
        idx = self.grid.idx
        for cell, v in self.clues.items():
            m.add(d[idx[cell]] == v)

    def _reachable_capitols(self, x, v):
        """Cells that could be the capitol nearest clue cell x (clue value v), by a lower-bound Dijkstra
        that charges each step the smallest edge weight the bounds allow."""
        bd, adj, n = self.bd, self.adj, self.n
        reach = []
        for cand in range(n):
            if not bd.cap_possible[cand]:
                continue
            c_sz = bd.capitol_size_lb(cand, adj)
            if c_sz > bd.sz_ub[cand]:
                continue
            dist = [10**9] * n
            dist[cand] = 0
            pq = [(0, cand)]
            while pq:
                d_u, u = heapq.heappop(pq)
                if d_u > dist[u] or d_u > v:
                    continue
                if u == x:
                    break
                sz_u = c_sz if u == cand else max(2, bd.sz_lb[u])
                for nb in adj[u]:
                    if nb == cand:
                        continue
                    sz_nb = max(2, bd.sz_lb[nb])
                    if u == cand:
                        sz_nb = max(sz_nb, bd.d_lb[nb])
                    w_lb = min(sz_u, sz_nb)
                    if bd.d_lb[nb] > v:
                        continue
                    eff_nd = max(d_u + w_lb, bd.d_lb[nb])
                    if eff_nd <= v and eff_nd < dist[nb]:
                        dist[nb] = eff_nd
                        heapq.heappush(pq, (eff_nd, nb))
            if dist[x] <= v:
                reach.append(cand)
        return reach

    def _add_clue_distance_cuts(self):
        if self.cfg.clue_cuts:
            m, adj, idx = self.m, self.adj, self.grid.idx
            for i in range(self.n):
                for j in adj[i]:
                    m.add(self.d[j] <= 1).OnlyEnforceIf(self.is_sing[i])
            for cell, v in self.clues.items():
                x = idx[cell]
                if v == 1:
                    m.add_bool_or([self.is_sing[j] for j in adj[x]])
                elif v > 1:
                    m.add_bool_or([self.cap[j] for j in self._reachable_capitols(x, v)])
            clue_list = list(self.clues.items())
            for idx_a in range(len(clue_list)):
                c1, v1 = clue_list[idx_a]
                for idx_b in range(idx_a + 1, len(clue_list)):
                    c2, v2 = clue_list[idx_b]
                    dist_manh = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
                    gap = abs(v1 - v2)
                    if 1 <= dist_manh <= 6 and gap > 0:
                        for path in monotone_paths(c1, c2):
                            w_vars = [self.w[idx[p[0]], idx[p[1]]] for p in path]
                            m.add(sum(w_vars) >= gap)
        self._stage("distances + clues")

    # -- refinement and search hints -------------------------------------------------------------------

    def _add_cuts(self, cuts):
        """Forbid 'state == S with capitol at cell' for each failed candidate."""
        idx, adj, cells = self.grid.idx, self.adj, self.grid.cells
        for S, capcell in cuts:
            rep = min(idx[x] for x in S)
            Sset = set(S)
            bnd = set()
            for x in S:
                for j in adj[idx[x]]:
                    if cells[j] not in Sset and j > rep:
                        bnd.add(j)
            self.m.add_bool_or(
                [self.b[idx[x], rep].Not() for x in S]
                + [self.b[j, rep] for j in bnd]
                + [self.cap[idx[capcell]].Not()]
            )
        self._stage("cuts")

    def _add_exclusions(self, exclude):
        """Forbid label vectors of partitions already found, so repeated solves enumerate every solution."""
        for L in exclude:
            self.m.add_bool_or([self.b[i, L[i]].Not() for i in range(self.n)])

    def _add_search_hints(self):
        if self.cfg.big_first and self.clues:
            idx = self.grid.idx
            sorted_clues = sorted(self.clues.items(), key=lambda x: -x[1])
            top_cells = [idx[pos] for pos, v in sorted_clues if v >= 15]
            if top_cells:
                self.m.add_decision_strategy(
                    [self.sz[i] for i in top_cells],
                    cp_model.CHOOSE_FIRST,
                    cp_model.SELECT_MAX_VALUE,
                )
