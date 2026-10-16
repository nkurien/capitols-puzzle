# capitols-puzzle

A CP-SAT (OR-Tools) solver for the Capitols puzzle. The full statement and the 5×5 worked example
are in [`problem.md`](problem.md). This README covers how the puzzle is modelled and how the model
is tested.

## The problem

Partition an N×N grid into *states*: connected sets of cells that are symmetric under some rotation
or reflection. If every nontrivial symmetry of a state fixes exactly one of its cells, that cell is the
state's *capitol*. A 1×3 bar, the X pentomino and a 2×2 block are all valid states with no capitol.

Moving between two orthogonally adjacent cells costs the size of the state you are in, or the smaller
of the two sizes when crossing a border. A clue gives the cheapest cost from that cell to any capitol,
so a clue of `0` marks a capitol. The answer is the sum over rows of the squared row sum, where each
cell is labelled with the size of its state (8,045 for the 5×5 example).

## Approach

The puzzle is encoded as a constraint satisfaction problem and solved with CP-SAT. The model is a
relaxation of the real rules: it enforces connectivity, sizes, distances and symmetry directly, but
checks each state against only one symmetry it chose. The full capitol rule (every symmetry fixes
exactly one cell) is checked afterwards on candidate solutions, and violations are fed back as cuts.

### Model

For an N×N grid with `n = N²` cells indexed row-major:

| Part | Variables | Constraints |
|---|---|---|
| State labels | `lab[i] ∈ [0, i]`, one-hot `b[i, r]` | A cell's label is the index of its state's lowest-numbered cell (the root), so each partition has exactly one encoding. `b[i, r]` implies `b[r, r]`: a label is only used if its root exists. |
| Same-state edges | `same[i, j]` for adjacent cells | `same[i, j]` holds exactly when `lab[i] = lab[j]`. |
| Connectivity | `h[i]` and a parent literal per neighbour | Every non-root cell has a same-state neighbour with a smaller `h`, and the root has `h = 0`. Each state is therefore a tree hanging off its root. |
| Sizes | `sz[i]` | Counted once at the root (`sz[r] = Σ b[i, r]`) and copied along same-state edges. |
| Symmetry | `T[i, g]`, `tvar[i]` | Each cell picks exactly one transform `g`. Adjacent cells in the same state must pick the same one, and the image `g(i)` must lie in the same state. The state is therefore mapped onto itself by `g`. |
| Capitols | `cap[i]` | `cap[i]` holds when the chosen `g` fixes `i` and no other cell fixed by `g` is in the same state. A clue of `0` forces `cap`; a positive clue forbids it. |
| Distances | `d[i]`, edge weight `w[i, j] = min(sz[i], sz[j])` | `d = 0` at capitols, `d[i] ≤ d[j] + w` on every edge, and every non-capitol has a tight parent (`d[i] = d[j] + w`). Together these make `d` the shortest-path distance to the nearest capitol. Clue cells are fixed to their value, and `d` is bounded by `maxd` (default 400). |

The candidate transforms are reflections across rows, columns and both diagonals (`H`, `V`, `D`, `A`)
and 180° and 90° rotations about every half-integer centre (`R2`, `R4`). They are generated per
grid size; for 11×11 that is 746 transforms, of which 441 remain in the model after the reductions below.

### Lazy capitol verification

Because the model only knows the one transform each state chose, it can declare a capitol that is not
one, for example when the state has a second symmetry that fixes several cells. `run()` loops:

1. Solve the model.
2. For each state with a declared capitol, compute its full symmetry group with `full_group_ok`,
   using all transforms including `R4`.
3. If the true capitol differs from the declared one, add a cut that forbids exactly that state with
   that capitol, and solve again.

A cut is a single clause over the state's membership literals and its boundary cells.

### Reductions

Each of these shrinks or tightens the model without changing which partitions, capitols and distances
it accepts:

- **Drop `R4`.** A 90° rotation about a centre implies the 180° rotation about the same centre, and
  both fix the same cells. `R4` is still used in the verification step.
- **Drop dominated transforms.** If transform `o` agrees with `k` everywhere `k` is defined, any
  state that can use `k` lies inside `k`'s domain and can use `o` instead, with the same fixed cells
  inside the state. At 11×11 this removes 84 of 525 transforms.
- **One variable per symmetric pair.** Every remaining transform is an involution, so `g(i) = j` implies
  `g(j) = i`. If `i` uses `g`, then `j` is in its state, and a connected state shares a single
  transform, so `T[i, g]` and `T[j, g]` are the same variable and the pair needs one `lab` equality.
- **Root exclusion.** If `T[i, g]` holds with `g(i) = j > i`, then `j` shares a state with a
  lower-numbered cell and cannot be that state's root.
- **Capitols have odd size.** The chosen transform is an involution that fixes only the capitol, so
  the other cells pair up. A one-cell state is always its own capitol, so a non-capitol has size ≥ 2.
- **Sizes by equalities** instead of an `element` lookup, which propagates better.

Transforms that fix no cell cannot be dropped: states such as a 2×2 block are symmetric only under
such transforms.

## Testing

- `tests/reference.py` is a brute-force oracle written from the problem statement. It finds
  symmetries from first principles (every lattice isometry that maps a cell set onto itself),
  enumerates all valid partitions of a small grid, and computes capitols and distances with its own
  Dijkstra. It shares no transform tables with `solve.py`.
- `tests/test_equivalence.py` enumerates every solution the CP-SAT model accepts, by re-solving with
  each found partition excluded, and requires the same partitions, capitols and distances as the
  oracle on 3×3 and 4×4 grids, including an unsatisfiable clue set. A model change that loses or
  invents a solution fails these tests.
- The equivalence tests were checked by breaking the model on purpose: removing the diagonal
  reflections, or the 180° rotations that fix no cell, makes them fail.
- Smaller checks: the oracle reproduces the 5×5 worked example; the dominated-transform helper is
  compared against a brute-force check on several grid sizes; the solver's 5×5 result is checked with an
  independent Dijkstra.

The oracle enumerates connected cell sets, which grows too quickly to use beyond 4×4. Equivalence is
therefore tested on small grids only.

## Benchmarking

`bench.py` builds random solvable puzzles by generating a random valid partition into symmetric
states, computing its true distances, and keeping a subset as clues. Every solver result is re-checked
against the oracle's rules.

CP-SAT's multi-threaded search is not deterministic, and the same puzzle can vary several-fold in
run time. For comparing model changes, run one thread per puzzle and several puzzles at once
(`--workers 1 --jobs 8`), run versions one after another, and avoid other heavy jobs on the machine.

## Setup

    uv venv
    uv pip install --python .venv/bin/python -r requirements-dev.txt

`python3 -m venv .venv && pip install -r requirements.txt` also works. Run everything with the
project interpreter (`.venv/bin/python`) so the right `ortools` is picked up.

## Usage

    .venv/bin/python solve.py ex                         # 5x5 example
    .venv/bin/python solve.py main --logfile main.log    # 11x11 puzzle

Options: `--workers N`, `--time SECONDS` (limit per solve; hitting it raises and exits with code 1),
`--heartbeat SECONDS`, `--logfile PATH`, `--cpsat-log` (raw CP-SAT log).

## Tests and benchmarks

    .venv/bin/python -m pytest                          # default suite
    .venv/bin/python -m pytest -m "slow or not slow"    # adds the larger 4x4 equivalence cases

    .venv/bin/python bench.py -n 6 7 8 --seeds 4 --workers 1 --jobs 8 --time 120

## Layout

| Path | Contents |
|---|---|
| `solve.py` | model builder, solve loop with lazy capitol verification, CLI |
| `problem.md` | problem statement, 11×11 grid, 5×5 worked example |
| `bench.py` | random solvable puzzles and timing harness |
| `tests/reference.py` | brute-force oracle written from the problem statement |
| `tests/test_*.py` | helper checks, oracle sanity, the 5×5 example, solution-set equivalence |
