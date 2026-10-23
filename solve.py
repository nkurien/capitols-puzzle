import argparse
import logging
import sys
import threading
import time

from ortools.sat.python import cp_model

log = logging.getLogger("capitols")
T0 = time.time()

def stage(m, name, t_start):
    p = m.Proto()
    log.info("[build] %-22s done in %6.2fs | vars=%7d constraints=%7d",
             name, time.time() - t_start, len(p.variables), len(p.constraints))
    return time.time()

def transforms(N):
    G = []  # (kind, params, func, fixedcells(list))
    M = 2 * (N - 1)
    for a in range(M + 1):
        G.append(('H', a, lambda r, c, a=a: (a - r, c)))
        G.append(('V', a, lambda r, c, a=a: (r, a - c)))
    for k in range(-(N - 1), N):
        G.append(('D', k, lambda r, c, k=k: (c + k, r - k)))
    for s in range(M + 1):
        G.append(('A', s, lambda r, c, s=s: (s - c, s - r)))
    for a in range(M + 1):
        for b in range(M + 1):
            G.append(('R2', (a, b), lambda r, c, a=a, b=b: (a - r, b - c)))
            if (a - b) % 2 == 0:
                G.append(('R4', (a, b), lambda r, c, a=a, b=b: ((a - b) // 2 + c, (a + b) // 2 - r)))
    return G

def fixed(g, N):
    f = g[2]
    return [(r, c) for r in range(N) for c in range(N) if f(r, c) == (r, c)]

def full_group_ok(S, N, G):
    """Return capitol cell if state S has one, else None."""
    Sset = set(S)
    caps = set()
    any_sym = False
    for g in G:
        f = g[2]
        if all(f(r, c) in Sset for (r, c) in S):
            any_sym = True
            fx = [x for x in S if f(*x) == x]
            if len(fx) != 1:
                return None, True
            caps.add(fx[0])
    if not any_sym:
        return None, False
    if len(caps) != 1:
        return None, True
    return caps.pop(), True

def undominated(images):
    """Indices of transforms worth modelling. A transform is redundant if another one agrees with it
    everywhere it is defined: any state that can use it can use the other, with the same fixed cells
    inside the state. Of exact duplicates the first is kept."""
    index = {}  # (cell, image) -> transforms sending that cell there, so candidates are found directly
    for k, img in enumerate(images):
        for i, x in enumerate(img):
            if x is not None:
                index.setdefault((i, x), []).append(k)
    keep = []
    for k, img in enumerate(images):
        c0 = next((i for i, x in enumerate(img) if x is not None), None)
        if c0 is None:
            continue
        defined = [(i, x) for i, x in enumerate(img) if x is not None]
        for o in index[c0, img[c0]]:
            if o != k and all(images[o][i] == x for i, x in defined):
                same_table = sum(x is not None for x in images[o]) == len(defined)
                if not same_table or o < k:
                    break
        else:
            keep.append(k)
    return keep

def monotone_paths(c1, c2):
    r1, c1_col = c1
    r2, c2_col = c2
    paths = []
    def dfs(curr, path):
        if curr == c2:
            paths.append(path)
            return
        cr, cc = curr
        if cr != r2:
            dr = 1 if r2 > cr else -1
            nxt = (cr + dr, cc)
            dfs(nxt, path + [(curr, nxt)])
        if cc != c2_col:
            dc = 1 if c2_col > cc else -1
            nxt = (cr, cc + dc)
            dfs(nxt, path + [(curr, nxt)])
    dfs(c1, [])
    return paths

def solve(N, clues, maxd=400, cuts=None, exclude=None, verbose=False, workers=8, tlimit=3000, heartbeat=30,
          max_size=None, no_free_r2=False, min_max_size=None, clue_bounds=False,
          tclauses=True, edge_images=True, lab_channel=False, eq_channel=False, clue_cuts=True,
          sound_distance_bounds=True, big_first=False):
    cuts = cuts or []
    exclude = exclude or []
    log.info('[build] starting model: N=%d, %d clues, %d cuts, maxd=%d', N, len(clues), len(cuts), maxd)
    t_build = time.time()
    G_full = transforms(N)
    # A 90-degree rotation about a centre implies the 180-degree rotation about the same centre, and
    # both fix exactly the same cells, so any state that could use R4 can use R2 with an identical
    # capitol outcome. The model therefore skips R4; G_full keeps it for the check in run().
    G = [g for g in G_full if g[0] != 'R4']
    cells = [(r, c) for r in range(N) for c in range(N)]
    idx = {x: i for i, x in enumerate(cells)}
    n = len(cells)
    # image of every cell under every transform as a cell index (None if it leaves the grid), and
    # the cells each transform fixes: computed once instead of per cell inside the model loops
    images = [[idx.get(g[2](r, c)) for (r, c) in cells] for g in G]
    n_before = len(G)
    keep = undominated(images)
    if no_free_r2:
        # search restriction: drop 180-degree rotations that fix no cell. That rules out states whose only
        # symmetry is such a rotation (S/Z-like shapes); anything found is still a genuine solution
        keep = [k for k in keep if G[k][0] != 'R2' or any(images[k][i] == i for i in range(n))]
    G, images = [G[k] for k in keep], [images[k] for k in keep]
    log.info('[build] %d symmetry transforms in the model (%d with R4, %d before dropping redundant ones)',
             len(G), len(G_full), n_before)
    fixed_at = [[i for i in range(n) if img[i] == i] for img in images]
    m = cp_model.CpModel()
    def nbrs(i):
        r, c = cells[i]
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (r + dr, c + dc) in idx:
                yield idx[(r + dr, c + dc)]
    lab = [m.new_int_var(0, i, f'lab{i}') for i in range(n)]
    b = {}
    for i in range(n):
        bs = []
        for r in range(i + 1):
            b[i, r] = m.new_bool_var('')
            bs.append(b[i, r])
        m.add_exactly_one(bs)
        m.add(lab[i] == sum(r * b[i, r] for r in range(i + 1)))
        for r in range(i):
            m.add_implication(b[i, r], b[r, r])
    t = stage(m, 'labels', t_build)
    edges = [(i, j) for i in range(n) for j in nbrs(i) if i < j]
    same = {}
    for (i, j) in edges:
        s = m.new_bool_var('')
        if not lab_channel:
            m.add(lab[i] == lab[j]).OnlyEnforceIf(s)
            m.add(lab[i] != lab[j]).OnlyEnforceIf(s.Not())
        else:
            for r in range(i + 1):
                m.add_bool_or([s.Not(), b[i, r].Not(), b[j, r]])
                m.add_bool_or([s.Not(), b[j, r].Not(), b[i, r]])
                m.add_bool_or([b[i, r].Not(), b[j, r].Not(), s])
            for r in range(i + 1, j + 1):
                m.add_implication(s, b[j, r].Not())
        same[i, j] = same[j, i] = s
    t = stage(m, 'same-state edges', t)
    # connectivity
    h = [m.new_int_var(0, n, '') for i in range(n)]
    for i in range(n):
        m.add(h[i] == 0).OnlyEnforceIf(b[i, i])
        m.add(h[i] >= 1).OnlyEnforceIf(b[i, i].Not())
        par = []
        for j in nbrs(i):
            p = m.new_bool_var('')
            m.add_implication(p, same[i, j])
            m.add(h[j] < h[i]).OnlyEnforceIf(p)
            par.append(p)
        m.add_bool_or(par + [b[i, i]])
    t = stage(m, 'connectivity', t)
    # sizes: counted at the state's root, then carried to every cell along same-state edges
    # (states are connected, so each cell reaches its root through such edges)
    # max_size caps every state's size: any solution found is still genuine, it only narrows the search
    sz = [m.new_int_var(1, max_size or n, '') for i in range(n)]
    is_sing = [m.new_bool_var(f'sing{i}') for i in range(n)]
    for i in range(n):
        m.add(sz[i] == 1).OnlyEnforceIf(is_sing[i])
        m.add(sz[i] >= 2).OnlyEnforceIf(is_sing[i].Not())
        m.add_implication(is_sing[i], b[i, i])
        for j in nbrs(i):
            m.add_implication(is_sing[i], same[min(i, j), max(i, j)].Not())
    if clue_cuts or clue_bounds:
        # size lower bounds read straight off the clues (presolve doesn't find them, as the edge weights
        # start with upper bound n). A cell between neighbours u and v is crossed by a two-edge route
        # costing at most twice its size, so 2*sz >= |d_u - d_v|; a singleton is a capitol, so all its
        # neighbours are within distance 1, and a cell next to a clue of 2 or more is not one
        cv = {idx[x]: v for x, v in clues.items()}
        for x in range(n):
            near = [cv[u] for u in nbrs(x) if u in cv]
            lb = 2 if any(v >= 2 for v in near) else 1
            if near:
                lb = max(lb, (max(near) - min(near) + 1) // 2)
            if lb > 1:
                m.add(sz[x] >= lb)
        if clue_cuts:
            # 1. Adjacent clues: sz[c1] >= gap, sz[c2] >= gap
            for (r, c), v in clues.items():
                for dr, dc in ((1, 0), (0, 1)):
                    if (r + dr, c + dc) in clues:
                        gap = abs(v - clues[(r + dr, c + dc)])
                        if gap > 1:
                            m.add(sz[idx[(r, c)]] >= gap)
                            m.add(sz[idx[(r + dr, c + dc)]] >= gap)

            # 2. Monotone paths of length 2 to 6 between clues:
            clue_list = list(clues.items())
            for idx_a in range(len(clue_list)):
                c1, v1 = clue_list[idx_a]
                for idx_b in range(idx_a + 1, len(clue_list)):
                    c2, v2 = clue_list[idx_b]
                    dist = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
                    gap = abs(v1 - v2)
                    if 2 <= dist <= 6 and gap > 1:
                        sum_sz = (gap + 1) // 2
                        min_max_sz = (gap + dist - 1) // dist
                        for path in monotone_paths(c1, c2):
                            mid_cells = [idx[step[1]] for step in path[:-1]]
                            if sum_sz > 1:
                                m.add(sum(sz[j] for j in mid_cells) >= sum_sz)
                            if min_max_sz > 1:
                                if len(mid_cells) == 1:
                                    m.add(sz[mid_cells[0]] >= min_max_sz)
                                else:
                                    bools = [m.new_bool_var('') for _ in mid_cells]
                                    for j, b_j in zip(mid_cells, bools):
                                        m.add(sz[j] >= min_max_sz).OnlyEnforceIf(b_j)
                                    m.add_bool_or(bools)
    for r in range(n):
        m.add(sz[r] == sum(b[i, r] for i in range(r, n))).OnlyEnforceIf(b[r, r])
    for (i, j) in edges:
        m.add(sz[i] == sz[j]).OnlyEnforceIf(same[i, j])
    t = stage(m, 'state sizes', t)
    # symmetry. Every transform left in G is an involution (i -> j means j -> i). If i uses the
    # transform, j is in its state, and a connected state shares one transform, so j uses it too: i and
    # j therefore share a single variable, and the pair needs only one lab equality.
    T = {}
    for gi, img in enumerate(images):
        for i in range(n):
            j = img[i]
            if j is None or i > j:
                continue
            assert img[j] == i
            v = m.new_bool_var('')
            T[i, gi] = T[j, gi] = v
            if i < j:
                m.add(lab[i] == lab[j]).OnlyEnforceIf(v)
                m.add_implication(v, b[j, j].Not())
    for i in range(n):
        opts = [T[i, gi] for gi in range(len(G)) if (i, gi) in T]
        m.add_exactly_one(opts)

    if not tclauses:
        tvar = [m.new_int_var(0, len(G) - 1, '') for i in range(n)]
        for i in range(n):
            opts_idx = [(gi, T[i, gi]) for gi in range(len(G)) if (i, gi) in T]
            m.add(tvar[i] == sum(gi * v for gi, v in opts_idx))
        for (i, j) in edges:
            m.add(tvar[i] == tvar[j]).OnlyEnforceIf(same[i, j])
    else:
        for (i, j) in edges:
            for gi in range(len(G)):
                in_i = (i, gi) in T
                in_j = (j, gi) in T
                if in_i and in_j:
                    vi = T[i, gi]
                    vj = T[j, gi]
                    if vi is not vj:
                        m.add_bool_or([same[i, j].Not(), vi.Not(), vj])
                        m.add_bool_or([same[i, j].Not(), vj.Not(), vi])
                elif in_i:
                    m.add_bool_or([same[i, j].Not(), T[i, gi].Not()])
                elif in_j:
                    m.add_bool_or([same[i, j].Not(), T[j, gi].Not()])

    if edge_images:
        seen_clauses = set()
        for gi, img in enumerate(images):
            for i in range(n):
                if (i, gi) not in T:
                    continue
                u = img[i]
                if u is None or i > u:
                    continue
                var_t = T[i, gi]
                t_id = id(var_t)
                for j in nbrs(i):
                    v = img[j]
                    if v is None:
                        m.add_bool_or([var_t.Not(), same[min(i, j), max(i, j)].Not()])
                    else:
                        e1 = (min(i, j), max(i, j))
                        e2 = (min(u, v), max(u, v))
                        if e1 != e2:
                            key = (t_id, min(e1, e2), max(e1, e2))
                            if key not in seen_clauses:
                                seen_clauses.add(key)
                                m.add_bool_or([var_t.Not(), same[e1].Not(), same[e2]])
                                m.add_bool_or([var_t.Not(), same[e1], same[e2].Not()])
    if min_max_size:
        # require some state of at least min_max_size cells (use only a proven lower bound, or it's a
        # restriction like max_size). Its root has that many cells at or after it, and its transform
        # needs at least that many cells with in-grid images
        M = min_max_size
        dom_size = [sum(x is not None for x in img) for img in images]
        big = []
        for r in range(n - M + 1):
            v = m.new_bool_var(f'big{r}')
            m.add_implication(v, b[r, r])
            m.add(sz[r] >= M).OnlyEnforceIf(v)
            for gi in range(len(G)):
                if (r, gi) in T and dom_size[gi] < M:
                    m.add_implication(v, T[r, gi].Not())
            big.append(v)
        m.add_bool_or(big)
    t = stage(m, 'symmetry', t)
    # capitols
    eqcache = {}
    def eq(i, j):
        if (i, j) in same:
            return same[i, j]
        key = (min(i, j), max(i, j))
        if key not in eqcache:
            s = m.new_bool_var('')
            if not eq_channel:
                m.add(lab[i] == lab[j]).OnlyEnforceIf(s)
                m.add(lab[i] != lab[j]).OnlyEnforceIf(s.Not())
            else:
                u, v = key
                for r in range(u + 1):
                    m.add_bool_or([s.Not(), b[u, r].Not(), b[v, r]])
                    m.add_bool_or([s.Not(), b[v, r].Not(), b[u, r]])
                    m.add_bool_or([b[u, r].Not(), b[v, r].Not(), s])
                for r in range(u + 1, v + 1):
                    m.add_implication(s, b[v, r].Not())
            eqcache[key] = s
        return eqcache[key]
    cap = [m.new_bool_var(f'cap{i}') for i in range(n)]
    for (r, c), v in clues.items():
        # a clue of 0 sits on a capitol; any positive clue does not
        m.add(cap[idx[(r, c)]] == (1 if v == 0 else 0))
    for i in range(n):
        m.add_implication(is_sing[i], cap[i])
    for i in range(n):
        ors = []
        for gi, img in enumerate(images):
            if img[i] != i:  # only transforms that fix this cell can make it a capitol
                continue
            others = [j for j in fixed_at[gi] if j != i]
            if not others:  # nothing else to be in the same state as: T itself says i is the capitol
                ors.append(T[i, gi])
                continue
            a = m.new_bool_var('')
            m.add_implication(a, T[i, gi])
            for j in others:
                m.add_implication(a, eq(i, j).Not())
            # if T and all others not same -> a
            m.add_bool_or([a, T[i, gi].Not()] + [eq(i, j) for j in others])
            ors.append(a)
        m.add_bool_or(ors).OnlyEnforceIf(cap[i])
        for a in ors:
            m.add_implication(a, cap[i])
    # A capitol's state is mapped to itself by an involution (reflection or R2) that fixes only the
    # capitol, so the remaining cells pair up and the state size is odd. Implied, but it links the
    # capitol and size variables directly for propagation.
    half = [m.new_int_var(0, n // 2, '') for i in range(n)]
    for i in range(n):
        m.add(sz[i] == 2 * half[i] + 1).OnlyEnforceIf(cap[i])
        # a one-cell state is symmetric about itself and so is always its own capitol
        m.add(sz[i] >= 2).OnlyEnforceIf(cap[i].Not())
    t = stage(m, 'capitols', t)
    # distances
    step_cap = max_size or n
    d = []
    d_lb = []
    d_ub = []
    for i in range(n):
        if sound_distance_bounds and clues:
            r, c = cells[i]
            lb_i = max([0] + [v - step_cap * (abs(r - cr) + abs(c - cc)) for (cr, cc), v in clues.items()])
            ub_i = min([maxd] + [v + step_cap * (abs(r - cr) + abs(c - cc)) for (cr, cc), v in clues.items()])
            if (r, c) in clues:
                lb_i = ub_i = clues[(r, c)]
            d_lb.append(lb_i)
            d_ub.append(ub_i)
            d.append(m.new_int_var(lb_i, ub_i, f'd{i}'))
        else:
            lb_i = clues.get(cells[i], 0)
            ub_i = clues.get(cells[i], maxd)
            d_lb.append(lb_i)
            d_ub.append(ub_i)
            d.append(m.new_int_var(lb_i, ub_i, f'd{i}'))
    w = {}
    for (i, j) in edges:
        v = m.new_int_var(1, n, '')
        m.add_min_equality(v, [sz[i], sz[j]])
        w[i, j] = w[j, i] = v
    for i in range(n):
        m.add(d[i] == 0).OnlyEnforceIf(cap[i])
        m.add(d[i] >= 1).OnlyEnforceIf(cap[i].Not())
        is_known_cap = (clues.get(cells[i]) == 0)
        par = []
        for j in nbrs(i):
            m.add(d[i] <= d[j] + w[i, j])
            if not is_known_cap and d_lb[j] < d_ub[i]:
                p = m.new_bool_var('')
                m.add(d[i] >= d[j] + w[i, j]).OnlyEnforceIf(p)
                par.append(p)
        if not is_known_cap:
            m.add_bool_or(par + [cap[i]])
    for (r, c), v in clues.items():
        m.add(d[idx[(r, c)]] == v)
    if clue_cuts:
        for i in range(n):
            for j in nbrs(i):
                m.add(d[j] <= 1).OnlyEnforceIf(is_sing[i])
        for (r, c), v in clues.items():
            if v == 1:
                x = idx[(r, c)]
                m.add_bool_or([is_sing[j] for j in nbrs(x)])
            elif v > 1:
                # Every path of k steps to a capitol costs >= 2k - 1, so k <= (v + 1) // 2
                max_steps = (v + 1) // 2
                reach = [idx[cand] for cand in cells if abs(r - cand[0]) + abs(c - cand[1]) <= max_steps and clues.get(cand, 0) == 0]
                m.add_bool_or([cap[j] for j in reach])
        clue_list = list(clues.items())
        for idx_a in range(len(clue_list)):
            c1, v1 = clue_list[idx_a]
            for idx_b in range(idx_a + 1, len(clue_list)):
                c2, v2 = clue_list[idx_b]
                dist_manh = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
                gap = abs(v1 - v2)
                if 1 <= dist_manh <= 6 and gap > 0:
                    for path in monotone_paths(c1, c2):
                        w_vars = [w[idx[p[0]], idx[p[1]]] for p in path]
                        m.add(sum(w_vars) >= gap)
    t = stage(m, 'distances + clues', t)
    # cuts: forbid (state==S with capitol at cell)
    for S, capcell in cuts:
        rep = min(idx[x] for x in S)
        Sset = set(S)
        bnd = set()
        for x in S:
            for j in nbrs(idx[x]):
                if cells[j] not in Sset and j > rep:
                    bnd.add(j)
        m.add_bool_or([b[idx[x], rep].Not() for x in S] + [b[j, rep] for j in bnd] + [cap[idx[capcell]].Not()])
    t = stage(m, 'cuts', t)
    # exclude: label vectors of partitions already found, so repeated solves enumerate every solution
    for L in exclude:
        m.add_bool_or([b[i, L[i]].Not() for i in range(n)])
    if big_first and clues:
        sorted_clues = sorted(clues.items(), key=lambda x: -x[1])
        top_cells = [idx[pos] for pos, v in sorted_clues if v >= 15]
        if top_cells:
            m.add_decision_strategy([sz[i] for i in top_cells], cp_model.CHOOSE_FIRST, cp_model.SELECT_MAX_VALUE)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = workers
    solver.parameters.max_time_in_seconds = tlimit
    solver.parameters.log_search_progress = verbose
    solver.parameters.log_to_stdout = False  # route through log_callback only, else every line prints twice
    if verbose:
        solver.log_callback = lambda line: log.debug('[cp-sat] %s', line)
    stop = threading.Event()
    t_solve = time.time()
    def beat():
        while not stop.wait(heartbeat):
            log.info('[solve] still searching... %.0fs elapsed (limit %ds)', time.time() - t_solve, tlimit)
    threading.Thread(target=beat, daemon=True).start()
    log.info('[solve] launching CP-SAT with %d workers, time limit %ds', workers, tlimit)
    st = solver.Solve(m)
    stop.set()
    log.info('[solve] finished: status=%s in %.1fs (conflicts=%d, branches=%d)',
             solver.StatusName(st), time.time() - t_solve, solver.NumConflicts(), solver.NumBranches())
    if st == cp_model.UNKNOWN:
        raise TimeoutError('CP-SAT hit the time limit before finding a solution or proving there is none')
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None
    L = [solver.Value(lab[i]) for i in range(n)]
    caps = [cells[i] for i in range(n) if solver.Value(cap[i])]
    D = [solver.Value(d[i]) for i in range(n)]
    return L, caps, D, G_full

def run(N, clues, verbose=False, **kw):
    cuts = []
    it = 0
    while True:
        it += 1
        log.info('=== refinement round %d ===', it)
        res = solve(N, clues, cuts=cuts, verbose=verbose, **kw)
        if res is None:
            log.warning('no solution found (model is infeasible)'); return None
        L, caps, D, G = res
        states = {}
        for i, l in enumerate(L):
            states.setdefault(l, []).append((i // N, i % N))
        log.info('[check] candidate has %d states, %d declared capitols', len(states), len(caps))
        log.info('[check] state sizes: %s', sorted((len(S) for S in states.values()), reverse=True))
        bad = False
        for l, S in states.items():
            cp, hassym = full_group_ok(S, N, G)
            declared = [x for x in caps if x in S]
            assert hassym
            if declared and cp != declared[0]:
                log.info('[check] state rooted at %s (size %d) declared capitol %s but true capitol is %s -> adding cut', S[0], len(S), declared[0], cp)
                cuts.append((S, declared[0])); bad = True
        if not bad:
            log.info('[check] all capitols verified against full symmetry groups')
            return L, caps, D, states
        log.info('[check] %d cut(s) total; re-solving', len(cuts))

def report(N, res):
    L, caps, _, states = res
    sizes = [len(states[L[i]]) for i in range(N * N)]
    for r in range(N):
        print(' '.join(f'{sizes[r * N + c]:3d}' for c in range(N)))
    print()
    for r in range(N):
        print(' '.join(f'{L[r * N + c]:3d}' for c in range(N)))
    print('caps', caps)
    rows = [sum(sizes[r * N:(r + 1) * N]) for r in range(N)]
    print('rows', rows, 'answer', sum(x * x for x in rows))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('which', choices=['ex', 'main'])
    ap.add_argument('--cpsat-log', action='store_true', help='stream the raw CP-SAT search log (very verbose)')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--time', type=int, default=3000, help='time limit per solve, seconds')
    ap.add_argument('--heartbeat', type=int, default=30, help='seconds between progress messages')
    ap.add_argument('--max-size', type=int, default=None, help='cap state sizes to narrow the search')
    ap.add_argument('--no-free-r2', action='store_true', help='drop 180-degree rotations that fix no cell to narrow the search')
    ap.add_argument('--min-max-size', type=int, default=None, help='require some state of at least this size')
    ap.add_argument('--clue-bounds', action='store_true', help='add size lower bounds derived from the clues')
    ap.add_argument('--no-tclauses', action='store_true', help='disable clause-based transform sharing')
    ap.add_argument('--no-edge-images', action='store_true', help='disable edge-image symmetry clauses')
    ap.add_argument('--lab-channel', action='store_true', help='enable label channelling on edges')
    ap.add_argument('--eq-channel', action='store_true', help='enable eq channelling for capitols')
    ap.add_argument('--no-clue-cuts', action='store_true', help='disable clue-derived path sums and singleton cuts')
    ap.add_argument('--no-sound-bounds', action='store_true', help='disable sound distance bounds')
    ap.add_argument('--big-first', action='store_true', help='branch on size of largest clue cells first')
    ap.add_argument('--logfile', default=None, help='also write the log to this file')
    args = ap.parse_args()
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.logfile:
        handlers.append(logging.FileHandler(args.logfile))
    logging.basicConfig(level=logging.DEBUG if args.cpsat_log else logging.INFO,
                        format='%(asctime)s %(message)s', datefmt='%H:%M:%S', handlers=handlers)
    which = args.which
    if which == 'ex':
        clues = {(0,0):12,(0,4):17,(1,0):8,(1,2):6,(1,3):3,(3,1):7,(3,2):6,(3,4):10,(4,0):25,(4,4):0}
        N = 5
    else:
        N = 11
        raw = {(0,1):8,(0,2):4,(0,5):11,(0,9):5,(1,7):11,(2,2):7,(2,5):1,(2,10):14,
               (3,0):28,(3,3):51,(3,6):1,(3,8):6,(4,1):22,(4,8):4,(4,10):1,
               (5,3):15,(5,5):10,(5,7):0,(6,0):11,(6,2):11,(6,9):9,
               (7,2):17,(7,4):14,(7,7):10,(7,10):13,(8,0):30,(8,5):6,(8,8):45,
               (9,3):10,(10,1):26,(10,5):0,(10,8):77,(10,9):61}
        clues = raw
    try:
        res = run(N, clues, verbose=args.cpsat_log, workers=args.workers, tlimit=args.time, heartbeat=args.heartbeat,
                  max_size=args.max_size, no_free_r2=args.no_free_r2, min_max_size=args.min_max_size, clue_bounds=args.clue_bounds,
                  tclauses=not args.no_tclauses, edge_images=not args.no_edge_images,
                  lab_channel=args.lab_channel, eq_channel=args.eq_channel,
                  clue_cuts=not args.no_clue_cuts, sound_distance_bounds=not args.no_sound_bounds,
                  big_first=args.big_first)
    except TimeoutError as e:
        log.warning('%s', e)
        sys.exit(1)
    if res:
        report(N, res)