"""End-to-end check of the 5x5 example, validated independently of the CP-SAT model."""

import heapq

import pytest

import solve as S

N = 5
CLUES = {
    (0, 0): 12,
    (0, 4): 17,
    (1, 0): 8,
    (1, 2): 6,
    (1, 3): 3,
    (3, 1): 7,
    (3, 2): 6,
    (3, 4): 10,
    (4, 0): 25,
    (4, 4): 0,
}


@pytest.fixture(scope="module")
def result():
    res = S.run(N, CLUES, workers=4, tlimit=120)
    assert res is not None
    return res


def neighbours(r, c):
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        if 0 <= r + dr < N and 0 <= c + dc < N:
            yield r + dr, c + dc


def test_states_partition_the_grid(result):
    *_, states = result
    cells = sorted(x for S_ in states.values() for x in S_)
    assert cells == [(r, c) for r in range(N) for c in range(N)]


def test_states_are_connected(result):
    states = result[3]
    for S_ in states.values():
        S_ = set(S_)
        seen, stack = set(), [next(iter(S_))]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack += [y for y in neighbours(*x) if y in S_]
        assert seen == S_


def test_every_state_is_symmetric_with_a_single_fixed_cell_per_symmetry(result):
    _, caps, _, states = result
    G = S.transforms(N)
    for S_ in states.values():
        cap, has_sym = S.full_group_ok(S_, N, G)
        assert has_sym
        for c in caps:
            if c in S_:
                assert cap == c


def test_at_most_one_capitol_per_state(result):
    _, caps, _, states = result
    for S_ in states.values():
        assert sum(c in S_ for c in caps) <= 1


def test_clues_are_met_by_independent_dijkstra(result):
    """Distance = cheapest path from a capitol; each step costs min(size of the two cells' states)."""
    L, caps, D, states = result
    size = {(r, c): len(states[L[r * N + c]]) for r in range(N) for c in range(N)}
    dist = {c: 0 for c in caps}
    pq = [(0, c) for c in caps]
    while pq:
        d, x = heapq.heappop(pq)
        if d > dist[x]:
            continue
        for y in neighbours(*x):
            nd = d + min(size[x], size[y])
            if nd < dist.get(y, 1 << 30):
                dist[y] = nd
                heapq.heappush(pq, (nd, y))
    for cell, want in CLUES.items():
        assert dist[cell] == want, cell
    for r in range(N):
        for c in range(N):
            assert dist[(r, c)] == D[r * N + c]


def test_known_answer(result):
    L, _, _, states = result
    sizes = [len(states[L[i]]) for i in range(N * N)]
    rows = [sum(sizes[r * N : (r + 1) * N]) for r in range(N)]
    assert sum(x * x for x in rows) == 8045
