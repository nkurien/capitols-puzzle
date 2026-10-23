"""Solution-set equivalence: the CP-SAT model must accept exactly the partitions the brute-force
reference accepts (same capitols, same distances), not just find *a* valid one.

Any model change that silently loses (or invents) a solution fails here. Run the slower 4x4 cases
with `pytest -m slow`.
"""

import pytest
import reference as R

import solve as S


def cp_solutions(N, clues):
    """Enumerate every solution of the CP-SAT model by re-solving with each found partition excluded."""
    seen, found = [], {}
    while True:
        res = S.run(N, clues, workers=1, tlimit=120, exclude=seen)
        if res is None:
            return found
        L, caps, D, states = res
        seen.append(L)
        part = frozenset(frozenset(cells) for cells in states.values())
        assert part not in found, "solver returned an excluded partition"
        found[part] = (frozenset(caps), {(i // N, i % N): D[i] for i in range(N * N)})


CASES_3 = [
    {(0, 0): 0, (1, 0): 1, (1, 2): 3, (2, 0): 0, (2, 1): 0, (2, 2): 1},
    {(1, 0): 3, (2, 0): 3, (2, 2): 3},
    {(0, 1): 1, (0, 2): 0, (1, 0): 1, (1, 2): 1, (2, 1): 2, (2, 2): 3},
    {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0, (1, 2): 1, (2, 1): 1},
]
CASES_4 = [
    {(0, 0): 0, (1, 2): 4, (2, 0): 0, (2, 1): 1, (2, 2): 3, (3, 2): 2},
    {(1, 3): 7, (2, 0): 1, (2, 2): 4, (2, 3): 7, (3, 0): 1, (3, 3): 4},
]
SLOW_CASES_4 = [
    {(0, 2): 5, (1, 1): 3, (2, 0): 0, (2, 2): 0, (3, 1): 1, (3, 2): 3},
    {(0, 0): 0, (1, 2): 5, (1, 3): 3, (3, 3): 4},
]
UNSAT_3 = {
    (0, 0): 0,
    (0, 1): 7,
}  # a 3x3 grid cannot cost 7 to reach a capitol from the cell next to one


def check(N, clues, expect_solutions=True):
    ref = R.solutions(N, clues)
    assert bool(ref) == expect_solutions, "test case is not what it claims to be"
    got = cp_solutions(N, clues)
    assert set(got) == set(ref), (
        f"{len(set(ref) - set(got))} partition(s) lost, {len(set(got) - set(ref))} invented"
    )
    for part, (caps, dist) in ref.items():
        assert got[part][0] == caps, "capitols differ"
        assert got[part][1] == dist, "distances differ"


@pytest.mark.parametrize("clues", CASES_3, ids=lambda c: f"{len(c)}clues")
def test_3x3_solution_sets_match_reference(clues):
    check(3, clues)


def test_3x3_unsatisfiable_clues_have_no_solutions():
    check(3, UNSAT_3, expect_solutions=False)


@pytest.mark.parametrize("clues", CASES_4, ids=lambda c: f"{len(c)}clues")
def test_4x4_solution_sets_match_reference(clues):
    check(4, clues)


@pytest.mark.slow
@pytest.mark.parametrize("clues", SLOW_CASES_4, ids=lambda c: f"{len(c)}clues")
def test_4x4_solution_sets_match_reference_slow(clues):
    check(4, clues)


def test_flags_equivalence_on_3x3():
    clues = CASES_3[0]
    ref = R.solutions(3, clues)
    for flag_kwargs in [
        {
            "tclauses": True,
            "edge_images": True,
            "lab_channel": True,
            "eq_channel": True,
            "clue_cuts": True,
        },
        {
            "tclauses": False,
            "edge_images": False,
            "lab_channel": False,
            "eq_channel": False,
            "clue_cuts": False,
        },
    ]:
        seen = []
        found = set()
        while True:
            res = S.run(3, clues, workers=1, tlimit=30, exclude=seen, **flag_kwargs)
            if res is None:
                break
            L, _caps, _D, states = res
            seen.append(L)
            found.add(frozenset(frozenset(c) for c in states.values()))
        assert found == set(ref)
