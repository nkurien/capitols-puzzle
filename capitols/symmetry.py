from collections.abc import Callable
from typing import NamedTuple


class Transform(NamedTuple):
    kind: str  # H, V, D, A (reflections), R2, R4 (rotations)
    params: int | tuple[int, int]  # axis, or centre in doubled coordinates
    fn: Callable  # (row, col) -> (row, col), possibly off the grid


def transforms(N):
    G = []
    M = 2 * (N - 1)
    for a in range(M + 1):
        G.append(Transform("H", a, lambda r, c, a=a: (a - r, c)))
        G.append(Transform("V", a, lambda r, c, a=a: (r, a - c)))
    for k in range(-(N - 1), N):
        G.append(Transform("D", k, lambda r, c, k=k: (c + k, r - k)))
    for s in range(M + 1):
        G.append(Transform("A", s, lambda r, c, s=s: (s - c, s - r)))
    for a in range(M + 1):
        for b in range(M + 1):
            G.append(Transform("R2", (a, b), lambda r, c, a=a, b=b: (a - r, b - c)))
            if (a - b) % 2 == 0:
                G.append(
                    Transform(
                        "R4",
                        (a, b),
                        lambda r, c, a=a, b=b: ((a - b) // 2 + c, (a + b) // 2 - r),
                    )
                )
    return G


def fixed(g, N):
    f = g[2]
    return [(r, c) for r in range(N) for c in range(N) if f(r, c) == (r, c)]


def full_group_ok(S, N, G):
    """Return (capitol cell or None, whether S has any symmetry) for state S."""
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


class ModelTransforms:
    """The transforms the model uses, with each one's image of every cell as a cell index (None if it
    leaves the grid) and the cells it fixes, computed once instead of per cell inside model loops."""

    def __init__(self, grid, no_free_r2=False):
        N = grid.N
        self.full = transforms(N)
        # A 90-degree rotation about a centre implies the 180-degree rotation about the same centre, and
        # both fix exactly the same cells, so any state that could use R4 can use R2 with an identical
        # capitol outcome. The model therefore skips R4; `full` keeps it for the check in run().
        G = [g for g in self.full if g.kind != "R4"]
        images = [[grid.idx.get(g.fn(r, c)) for (r, c) in grid.cells] for g in G]
        self.n_before = len(G)
        keep = undominated(images)
        if no_free_r2:
            # search restriction: drop 180-degree rotations that fix no cell. That rules out states whose
            # only symmetry is such a rotation (S/Z-like shapes)
            keep = [
                k
                for k in keep
                if G[k].kind != "R2" or any(images[k][i] == i for i in range(grid.n))
            ]
        self.G = [G[k] for k in keep]
        self.images = [images[k] for k in keep]
        self.fixed_at = [
            [i for i in range(grid.n) if img[i] == i] for img in self.images
        ]

    def __len__(self):
        return len(self.G)
