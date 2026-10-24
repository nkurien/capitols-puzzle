"""Compile-time reasoning: bounds on state sizes and distances that hold in every solution, and which
(root, transform) pairs can possibly describe a state. Nothing here touches CP-SAT, so it can be tested
directly and costs no model variables or constraints."""

from dataclasses import dataclass


@dataclass
class Bounds:
    sz_lb: list[int]
    sz_ub: list[int]
    d_lb: list[int]
    d_ub: list[int]
    cap_possible: list[bool]

    def capitol_size_lb(self, i, adj, clue_cuts=True):
        """Smallest odd size a state can have with its capitol at cell i."""
        lb = self.sz_lb[i]
        if clue_cuts:
            lb = max([lb] + [self.d_lb[nb] for nb in adj[i]])
        return lb + 1 if lb % 2 == 0 else lb


def initial_bounds(grid, clues, cfg):
    n = grid.n
    use_cuts = cfg.clue_cuts or cfg.clue_bounds
    step_cap = cfg.max_size or n
    bd = Bounds([1] * n, [step_cap] * n, [0] * n, [cfg.maxd] * n, [True] * n)
    for i in range(n):
        if cfg.sound_distance_bounds and clues:
            r, c = grid.cells[i]
            bd.d_lb[i] = max(
                [0]
                + [
                    v - step_cap * (abs(r - cr) + abs(c - cc))
                    for (cr, cc), v in clues.items()
                ]
            )
            bd.d_ub[i] = min(
                [cfg.maxd]
                + [
                    v + step_cap * (abs(r - cr) + abs(c - cc))
                    for (cr, cc), v in clues.items()
                ]
            )
        if grid.cells[i] in clues:
            v = clues[grid.cells[i]]
            bd.d_lb[i] = bd.d_ub[i] = v
            if v > 0:
                bd.cap_possible[i] = False
                if use_cuts:
                    bd.sz_lb[i] = 2
    return bd


def tighten_bounds(grid, clues, cfg):
    """Bounds on sizes, distances and capitol positions, iterated to a fixed point (at most 8 rounds)."""
    bd = initial_bounds(grid, clues, cfg)
    if not (cfg.clue_cuts or cfg.clue_bounds) or not clues:
        return bd
    adj = grid.adj
    cv = {grid.idx[x]: v for x, v in clues.items()}
    for x in range(grid.n):
        near = [cv[u] for u in adj[x] if u in cv]
        lb = 2 if any(v >= 2 for v in near) else 1
        if near:
            lb = max(lb, (max(near) - min(near) + 1) // 2)
        bd.sz_lb[x] = max(bd.sz_lb[x], lb)
    if not cfg.clue_cuts:
        return bd
    for _ in range(8):
        changed = _local_rules(bd, adj)
        changed |= _forced_singleton_capitols(bd, adj)
        changed |= _distance_rules(bd, adj, cfg.sound_distance_bounds)
        changed |= _two_step_size_rule(bd, adj)
        if not changed:
            break
    return bd


def _local_rules(bd, adj):
    changed = False
    for i in range(len(adj)):
        if bd.d_lb[i] > 0 and bd.cap_possible[i]:
            bd.cap_possible[i] = False
            changed = True
        if not bd.cap_possible[i] and bd.d_lb[i] < 1:
            bd.d_lb[i] = 1
            changed = True
        if not bd.cap_possible[i] and bd.sz_lb[i] < 2:
            bd.sz_lb[i] = 2
            changed = True
        if bd.sz_lb[i] == 1 and any(bd.d_lb[nb] >= 2 for nb in adj[i]):
            bd.sz_lb[i] = 2
            changed = True
        if bd.sz_ub[i] == 1:
            # a singleton is its own capitol and so one step from each neighbour
            if bd.d_ub[i] > 0:
                bd.d_ub[i] = 0
                changed = True
            for nb in adj[i]:
                if bd.d_ub[nb] > 1:
                    bd.d_ub[nb] = 1
                    changed = True
    return changed


def _forced_singleton_capitols(bd, adj):
    """A cell at distance exactly 1 that is not a singleton needs a singleton-capitol neighbour; if
    only one neighbour can be one, it is forced."""
    changed = False
    for i in range(len(adj)):
        if bd.d_lb[i] == 1 and bd.d_ub[i] == 1 and bd.sz_lb[i] >= 2:
            sing_cands = [
                nb for nb in adj[i] if bd.sz_lb[nb] == 1 and bd.cap_possible[nb]
            ]
            if len(sing_cands) == 1:
                s = sing_cands[0]
                if bd.sz_ub[s] > 1:
                    bd.sz_ub[s] = 1
                    bd.d_ub[s] = 0
                    changed = True
    return changed


def _distance_rules(bd, adj, sound):
    """Along an edge d[i] <= d[j] + min(sz[i], sz[j]), and every non-capitol has a tight predecessor."""
    changed = False
    for i in range(len(adj)):
        for j in adj[i]:
            if sound:
                w_ub = min(bd.sz_ub[i], bd.sz_ub[j])
                if bd.d_ub[i] > bd.d_ub[j] + w_ub:
                    bd.d_ub[i] = bd.d_ub[j] + w_ub
                    changed = True
            req = bd.d_lb[i] - bd.d_ub[j]
            if req > bd.sz_lb[i]:
                bd.sz_lb[i] = req
                changed = True
            if req > bd.sz_lb[j]:
                bd.sz_lb[j] = req
                changed = True
        if sound and not bd.cap_possible[i]:
            min_inc = min(bd.d_lb[j] + min(bd.sz_lb[i], bd.sz_lb[j]) for j in adj[i])
            if min_inc > bd.d_lb[i]:
                bd.d_lb[i] = min_inc
                changed = True
    return changed


def _two_step_size_rule(bd, adj):
    """Two neighbours a, b of m: d[a] - d[b] <= 2 * sz[m]."""
    changed = False
    for m_cell in range(len(adj)):
        for a in adj[m_cell]:
            for b_nb in adj[m_cell]:
                if a != b_nb:
                    req = (bd.d_lb[a] - bd.d_ub[b_nb] + 1) // 2
                    if req > bd.sz_lb[m_cell]:
                        bd.sz_lb[m_cell] = req
                        changed = True
    return changed


class RegionFinder:
    """Cells a symmetric state could contain, given its transform."""

    def __init__(self, grid, mt, bd):
        self.n = grid.n
        self.adj = grid.adj
        self.images = mt.images
        self.fixed_at = mt.fixed_at
        self.bd = bd

    def _forbidden_sole_fixed(self, gi):
        """A transform fixing exactly one cell puts its capitol there, so that cell must be allowed to be one."""
        if len(self.fixed_at[gi]) == 1:
            f = self.fixed_at[gi][0]
            if not self.bd.cap_possible[f]:
                return {f}
        return set()

    def _shrink(self, valid, seed, img):
        """Largest subset of `valid` connected to seed and closed under the transform."""
        adj = self.adj
        while True:
            comp = {seed}
            queue = [seed]
            while queue:
                curr = queue.pop()
                for nb in adj[curr]:
                    if nb in valid and nb not in comp:
                        comp.add(nb)
                        queue.append(nb)
            new_valid = {u for u in comp if img[u] in comp}
            if len(new_valid) == len(valid):
                return new_valid
            valid = new_valid
            if seed not in valid:
                return set()

    def max_region(self, seed, gi, forbid_other_fixed=False):
        """Largest region for a state containing `seed` with transform gi."""
        img = self.images[gi]
        if img[seed] is None:
            return set()
        forbid = set()
        if forbid_other_fixed:
            for f in self.fixed_at[gi]:
                if f != seed:
                    forbid.add(f)
        forbid |= self._forbidden_sole_fixed(gi)
        bd = self.bd
        valid = set()
        for u in range(self.n):
            v = img[u]
            if v is None or u in forbid or v in forbid:
                continue
            if u != seed and (
                bd.sz_ub[u] < bd.sz_lb[seed] or bd.sz_ub[v] < bd.sz_lb[seed]
            ):
                continue
            valid.add(u)
        if seed not in valid:
            return set()
        return self._shrink(valid, seed, img)

    def valid_region_from_root(self, r, gi):
        """Largest region for a state whose lowest-numbered cell is r. Every cell and its image must be
        at or after r, so this is far tighter than max_region."""
        img = self.images[gi]
        if img[r] is None or img[r] < r:
            return set()
        forbid = self._forbidden_sole_fixed(gi)
        bd = self.bd
        valid = set()
        for u in range(r, self.n):
            v = img[u]
            if v is None or v < r or u in forbid or v in forbid:
                continue
            if u != r and (bd.sz_ub[u] < bd.sz_lb[r] or bd.sz_ub[v] < bd.sz_lb[r]):
                continue
            valid.add(u)
        if r not in valid:
            return set()
        return self._shrink(valid, r, img)


@dataclass
class RootRegions:
    regions: dict  # (root, transform) -> set of cells the state could occupy
    max_s: dict  # (root, transform) -> largest state size that is possible
    possible_b: (
        set  # (cell, root) pairs: the cell can belong to the state rooted at root
    )


def root_regions(grid, mt, bd, finder):
    out = RootRegions({}, {}, set())
    fixed_at = mt.fixed_at
    for r in range(grid.n):
        for gi in range(len(mt)):
            reg = finder.valid_region_from_root(r, gi)
            if not reg:
                continue
            max_s = len(reg)
            # without a fixed cell inside the region the cells pair up, so the size is even
            if (
                len(fixed_at[gi]) == 0
                or (len(fixed_at[gi]) == 1 and fixed_at[gi][0] not in reg)
            ) and max_s % 2 == 1:
                max_s -= 1
            if max_s < bd.sz_lb[r]:
                continue
            out.regions[r, gi] = reg
            out.max_s[r, gi] = max_s
            for j in reg:
                req = max(bd.sz_lb[r], bd.sz_lb[j], 2 if j > r else 1)
                if (
                    max_s >= req
                    and bd.sz_ub[j] >= bd.sz_lb[r]
                    and bd.sz_ub[r] >= bd.sz_lb[j]
                ):
                    out.possible_b.add((j, r))
    return out
