"""Brute-force reference for the Capitols rules, written from the problem statement.

Deliberately independent of solve.py (no shared transform tables) so it can cross-check the
CP-SAT model on tiny grids by enumerating every valid partition outright.
"""
import functools
import heapq

# linear parts of the 8 lattice isometries preserving the square grid: (r, c) -> (a*r + b*c, c_*r + d*c)
LINEAR = [(1, 0, 0, 1), (0, -1, 1, 0), (-1, 0, 0, -1), (0, 1, -1, 0),
          (-1, 0, 0, 1), (1, 0, 0, -1), (0, 1, 1, 0), (0, -1, -1, 0)]
IDENTITY = (1, 0, 0, 1)


def apply(M, t, x):
    a, b, c, d = M
    return (a * x[0] + b * x[1] + t[0], c * x[0] + d * x[1] + t[1])


def neighbours(x, N):
    r, c = x
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        if 0 <= r + dr < N and 0 <= c + dc < N:
            yield (r + dr, c + dc)


def is_connected(S, N):
    S = set(S)
    seen, stack = set(), [next(iter(S))]
    while stack:
        x = stack.pop()
        if x not in seen:
            seen.add(x)
            stack += [y for y in neighbours(x, N) if y in S]
    return seen == S


def symmetries(S):
    """Every nontrivial isometry of the plane mapping the cell set S onto itself, as (M, t)."""
    S = frozenset(S)
    s0 = min(S)
    out = []
    for M in LINEAR:
        base = apply(M, (0, 0), s0)
        for s in S:  # the translation is fixed by where s0 lands
            t = (s[0] - base[0], s[1] - base[1])
            if M == IDENTITY and t == (0, 0):
                continue
            if all(apply(M, t, x) in S for x in S):
                out.append((M, t))
    return out


def classify(S, N):
    """None if S is not a state; otherwise (True, capitol_or_None)."""
    if not is_connected(S, N):
        return None
    syms = symmetries(S)
    if not syms:
        return None
    caps = set()
    for M, t in syms:
        fixed = [x for x in S if apply(M, t, x) == x]
        if len(fixed) != 1:
            return (True, None)
        caps.add(fixed[0])
    # symmetries of a finite set share a common fixed point, so this is a single cell
    return (True, caps.pop() if len(caps) == 1 else None)


@functools.cache
def valid_states(N):
    """{frozenset(cells): capitol or None} for every state that fits in the N x N grid."""
    out = {}
    for mask in range(1, 1 << (N * N)):
        S = frozenset((i // N, i % N) for i in range(N * N) if mask >> i & 1)
        info = classify(S, N)
        if info is not None:
            out[S] = info[1]
    return out


def distances(partition, N, capitol_of):
    """Dijkstra from all capitols; a step costs min(size of the two cells' states). None if no capitol.
    capitol_of maps a state (frozenset of cells) to its capitol cell, or None."""
    size = {x: len(S) for S in partition for x in S}
    caps = [c for S in partition if (c := capitol_of(S)) is not None]
    if not caps:
        return None
    dist = {c: 0 for c in caps}
    pq = [(0, c) for c in caps]
    while pq:
        d, x = heapq.heappop(pq)
        if d > dist[x]:
            continue
        for y in neighbours(x, N):
            nd = d + min(size[x], size[y])
            if nd < dist.get(y, 1 << 30):
                dist[y] = nd
                heapq.heappush(pq, (nd, y))
    return dist


def solutions(N, clues):
    """{partition: (capitols, dist)} for every partition of the grid into valid states (with at least
    one capitol) whose capitol distances match the clues.

    A clue of 0 pins a capitol on that cell and a positive clue rules one out; applying that while
    building the partition keeps the search small without changing the answer."""
    states = valid_states(N)
    cells = [(r, c) for r in range(N) for c in range(N)]
    by_min = {x: [S for S in states if min(S) == x] for x in cells}
    must_cap = {x for x, d in clues.items() if d == 0}
    not_cap = {x for x, d in clues.items() if d > 0}
    out = {}

    def ok(S):
        cap = states[S]
        return all(states[S] == x for x in must_cap & S) and not (cap is not None and cap in not_cap)

    def rec(covered, chosen):
        free = next((x for x in cells if x not in covered), None)
        if free is None:
            part = frozenset(chosen)
            dist = distances(part, N, states.__getitem__)
            if dist is not None and all(dist[x] == d for x, d in clues.items()):
                out[part] = (frozenset(states[S] for S in part if states[S] is not None), dist)
            return
        for S in by_min[free]:
            if not (S & covered) and ok(S):
                rec(covered | S, chosen + [S])

    rec(frozenset(), [])
    return out
