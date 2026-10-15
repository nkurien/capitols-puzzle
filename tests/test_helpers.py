import pytest

import solve as S


@pytest.mark.parametrize("N,expected", [(5, 158), (11, 746)])
def test_transform_count(N, expected):
    assert len(S.transforms(N)) == expected


@pytest.mark.parametrize("N", [3, 5, 6])
def test_transforms_are_injective_on_their_domain(N):
    cells = [(r, c) for r in range(N) for c in range(N)]
    for kind, params, f in [(g[0], g[1], g[2]) for g in S.transforms(N)]:
        imgs = [f(r, c) for r, c in cells]
        inside = [x for x in imgs if 0 <= x[0] < N and 0 <= x[1] < N]
        assert len(inside) == len(set(inside)), (kind, params)


def test_fixed_cells_of_known_transforms():
    N = 5
    G = {(g[0], g[1]): g for g in S.transforms(N)}
    assert S.fixed(G['H', 4], N) == [(2, c) for c in range(N)]   # reflect across middle row
    assert S.fixed(G['V', 4], N) == [(r, 2) for r in range(N)]   # reflect across middle column
    assert S.fixed(G['D', 0], N) == [(i, i) for i in range(N)]   # main diagonal
    assert S.fixed(G['R2', (4, 4)], N) == [(2, 2)]               # 180 deg about centre
    assert S.fixed(G['R2', (3, 4)], N) == []                     # centre falls between cells
    assert S.fixed(G['R4', (4, 4)], N) == [(2, 2)]               # 90 deg about centre


@pytest.mark.parametrize("cells,cap", [
    ([(2, 2)], (2, 2)),                                        # single cell
    ([(1, 1), (1, 2), (2, 2), (3, 2), (3, 3)], (2, 2)),        # S shape: only 180 deg symmetry
    ([(1, 2), (2, 1), (2, 2), (2, 3), (3, 2)], None),          # plus: mirrors fix 3 cells
    ([(0, 0), (0, 1), (0, 2), (1, 1)], None),                  # T: mirror symmetry fixes 2 cells
    ([(1, 0), (1, 1), (1, 2)], None),                          # line: reflections fix several cells
    ([(0, 0), (0, 1), (1, 0), (1, 1)], None),                  # 2x2: rotations fix no cell
    ([(0, 0), (0, 1), (0, 2), (1, 0)], None),                  # L tetromino: no symmetry at all
])
def test_full_group_ok(cells, cap):
    got, _ = S.full_group_ok(cells, 5, S.transforms(5))
    assert got == cap


def test_full_group_ok_reports_no_symmetry():
    _, has_sym = S.full_group_ok([(0, 0), (0, 1), (0, 2), (1, 0)], 5, S.transforms(5))
    assert has_sym is False


@pytest.mark.parametrize("N", [3, 4, 5, 6, 11])
def test_undominated_matches_brute_force(N):
    """The model drops a transform only if another agrees with it everywhere it is defined."""
    G = [g for g in S.transforms(N) if g[0] != 'R4']
    cells = [(r, c) for r in range(N) for c in range(N)]
    idx = {x: i for i, x in enumerate(cells)}
    images = [[idx.get(g[2](r, c)) for (r, c) in cells] for g in G]
    tables = [{(i, x) for i, x in enumerate(img) if x is not None} for img in images]
    dominated = {
        k for k in range(len(G))
        if any(o != k and (tables[k] < tables[o] or (tables[k] == tables[o] and o < k)) for o in range(len(G)))
    }
    assert set(range(len(G))) - set(S.undominated(images)) == dominated
    assert dominated, "expected some redundant transforms (corner and edge-centred rotations)"
