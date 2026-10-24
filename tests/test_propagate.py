"""Compile-time propagation, tested without building or solving a CP-SAT model."""

import pytest

from capitols.config import SolverConfig
from capitols.grid import Grid
from capitols.problem import prepare
from capitols.propagate import tighten_bounds
from capitols.puzzles import EX_5X5, MAIN_11X11


@pytest.fixture(scope="module")
def main_prob():
    return prepare(11, MAIN_11X11, SolverConfig(min_max_size=26))


def test_main_forces_singleton_capitol_at_5_10(main_prob):
    g, bd = main_prob.grid, main_prob.bounds
    i = g.idx[(5, 10)]
    assert (bd.sz_lb[i], bd.sz_ub[i], bd.d_lb[i], bd.d_ub[i]) == (1, 1, 0, 0)
    assert bd.cap_possible[i]


def test_main_right_edge_cascade(main_prob):
    g, bd = main_prob.grid, main_prob.bounds
    # neighbours of a singleton capitol are at distance 1, and the clues beyond them force big states
    assert bd.sz_lb[g.idx[(6, 10)]] >= 12
    assert bd.sz_lb[g.idx[(5, 9)]] >= 8
    assert bd.sz_lb[g.idx[(3, 10)]] >= 7


@pytest.mark.parametrize("N,clues", [(5, EX_5X5), (11, MAIN_11X11)])
def test_bounds_are_consistent(N, clues):
    g = Grid(N)
    bd = tighten_bounds(g, clues, SolverConfig())
    for i in range(g.n):
        assert 1 <= bd.sz_lb[i] <= bd.sz_ub[i]
        assert 0 <= bd.d_lb[i] <= bd.d_ub[i]
        if not bd.cap_possible[i]:
            assert bd.d_lb[i] >= 1 and bd.sz_lb[i] >= 2
    for cell, v in clues.items():
        i = g.idx[cell]
        assert bd.d_lb[i] == bd.d_ub[i] == v


def test_root_regions_respect_the_root_and_transform(main_prob):
    mt, reg, finder = main_prob.mt, main_prob.regions, main_prob.finder
    assert reg.regions
    for (r, gi), region in reg.regions.items():
        img = mt.images[gi]
        assert r in region
        for u in region:
            assert u >= r and img[u] is not None and img[u] >= r
            assert img[u] in region  # closed under the transform
        assert reg.max_s[r, gi] <= len(region)
        assert finder.valid_region_from_root(r, gi) == region


def test_cell_centred_r2_cannot_root_a_larger_state(main_prob):
    # under 180-degree rotation about cell r, every other cell maps below r, so r can only root itself
    mt, reg = main_prob.mt, main_prob.regions
    for gi, g in enumerate(mt.G):
        if g.kind == "R2":
            for r in mt.fixed_at[gi]:
                assert reg.max_s.get((r, gi), 1) <= 1
