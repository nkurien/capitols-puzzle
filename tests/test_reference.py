"""Sanity checks for the brute-force reference itself, so it can be trusted as the oracle."""

import pytest
import reference as R

import solve as S


@pytest.mark.parametrize(
    "cells,expected",
    [
        ([(2, 2)], (True, (2, 2))),  # singleton: its own capitol
        (
            [(1, 0), (1, 1), (1, 2)],
            (True, None),
        ),  # 1x3: no capitol (from the problem text)
        (
            [(1, 2), (2, 1), (2, 2), (2, 3), (3, 2)],
            (True, None),
        ),  # X pentomino: no capitol (problem text)
        ([(0, 0), (0, 1), (1, 0), (1, 1)], (True, None)),  # 2x2: rotation fixes nothing
        (
            [(1, 1), (1, 2), (2, 2), (3, 2), (3, 3)],
            (True, (2, 2)),
        ),  # S shape: 180 deg only
        (
            [(0, 0), (1, 0), (1, 1), (1, 2), (0, 2)],
            (True, (1, 1)),
        ),  # U pentomino: single mirror, one axis cell
        (
            [(0, 0), (0, 1), (0, 2), (1, 0)],
            None,
        ),  # L tetromino: no symmetry, not a state
        ([(0, 0), (2, 2)], None),  # disconnected: not a state
    ],
)
def test_classify_known_shapes(cells, expected):
    assert R.classify(cells, 5) == expected


@pytest.mark.parametrize("N", [3, 4])
def test_reference_agrees_with_solver_helpers_on_every_connected_subset(N):
    """solve.full_group_ok (table of transforms) vs the first-principles symmetry search."""
    G = S.transforms(N)
    checked = 0
    for mask in range(1, 1 << (N * N)):
        cells = [(i // N, i % N) for i in range(N * N) if mask >> i & 1]
        if not R.is_connected(cells, N):
            continue
        cap, has_sym = S.full_group_ok(cells, N, G)
        ref = R.classify(cells, N)
        assert has_sym == (ref is not None), cells
        if ref is not None:
            assert cap == ref[1], cells
        checked += 1
    assert checked > 100


def test_reference_reproduces_the_worked_example():
    """The 5x5 example from the problem statement: states, capitols, distances and the answer."""
    A = frozenset(
        [(0, 0), (1, 0), (2, 0), (2, 1), (3, 1), (3, 2), (3, 3)]
    )  # 7-cell purple state
    B = frozenset([(0, 1), (0, 2), (1, 1), (1, 2)])  # 4-cell blue state
    C = frozenset([(1, 3), (2, 2), (2, 3)])  # 3-cell green state
    D = (
        frozenset((r, c) for r in range(5) for c in range(5)) - A - B - C
    )  # 11-cell red state
    part = [A, B, C, D]
    infos = {X: R.classify(X, 5) for X in part}
    assert None not in infos.values(), "every part of the example must be a valid state"
    caps = {X: info[1] for X, info in infos.items() if info is not None}
    assert caps == {A: (2, 1), B: None, C: (2, 3), D: (4, 4)}
    dist = R.distances(part, 5, caps.__getitem__)
    assert dist is not None, "the example has capitols"
    clues = {
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
    assert all(dist[x] == d for x, d in clues.items())
    size = {x: len(X) for X in part for x in X}
    rows = [sum(size[(r, c)] for c in range(5)) for r in range(5)]
    assert sum(x * x for x in rows) == 8045
