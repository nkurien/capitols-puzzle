# Capitols Puzzle (`main`) — Solver Diagnostics, Mathematical Deductions & Benchmark Findings

This document summarizes all mathematical deductions, 3,000s run outcomes, 20s CP-SAT micro-diagnostics, benchmark ablations (`bench.py`), and architectural experiments (including Option A's two-stage metric relaxation) on the $11 \times 11$ Capitols puzzle (`main`).

---

## 1. Executive Summary: What Is Faster vs. What Is Slower

### Kept

- Dedupe and filter the `sz[j] >= min_max_sz` path booleans: removes `3,841` redundant variables (`25,121` $\to$ `21,280`) with identical constraints.
- Compile-time domain tightening of `(sz_lb, sz_ub, d_lb, d_ub, cap_possible)`: proves `(5, 10)` is a singleton capitol (`d=0, sz=1`) at build time and tightens `sz_lb` and `d_ub` along the right edge.
- Root-region pruning (`valid_region_from_root(r, gi)`): valid `--min-max-size 26` `(r, gi)` pairs drop from `7,990` to `4,558`, and `1,302` impossible `(r, gi)` root assignments are pruned.
- Non-singleton incident-edge clause (`is_sing[i] \/ OR_j same[i, j]`) plus Dijkstra capitol reachability: 121 short clauses, and smaller capitol candidate sets (e.g. `(6, 0)=11` goes from 33 to 17 candidates).
- `--require-free-r2` (from the Run 2 `INFEASIBLE` proof): one clause requiring at least one fixed-point-free $180^\circ$ rotation `R2`. Specific to `main`.

### Dropped

- `linearization_level=0` (`--lin-level 0`): raised raw 20s conflict rate, but slowed down both the `--max-size 35` infeasibility proof (median `61s` vs `38s` for default LP) and the full `main` solve (`2,491.7s` vs **`1,580.3s`** for default LP). One run each on `main`, so the size of the gap is uncertain. A possible reason (not verified): with LP off, `default_lp`, `quick_restart`, and `probing` lose their `IB`/`MIR`/`ZeroHalf`/`CG` cuts (`1.8M` in the default-LP run), which may weaken what they share with the `no_lp` worker.
- ~60,000 per-`T[i, gi]` parity and mirror-axis clauses (`180k` $\to$ `255k` constraints): `bench.py -n 6 7 8` went from `158.4s` (median `4.4s`) to `202.2s` (median `8.4s`), because unit propagation got slower.
- Option A, two-stage pure metric relaxation (drops `b[i, r]` and `T[i, gi]`): solves in `0.32s`, but only by building small blobs with claimed sizes summing to more than 250. Once $\sum |S| = 121$ is enforced it times out (`>15s`). A mirage.

---

## 2. What the Three 3,000-Second Runs on `main` Proved

Run 1: baseline (`solve.py main`)
- Model: `31,119` vars / `39,435` constraints.
- Outcome: `UNKNOWN` at `3,000s` (`463k` conflicts on Worker 0).
- Takeaway: the integer `tvar` channel without ternary clauses propagated weakly.

Run 2: `--min-max-size 26 --no-free-r2`
- Model: `16,026` vars / `415,491` constraints.
- Outcome: `INFEASIBLE` in `1,185.5s` (`1.56M` conflicts).
- Theorem: every valid solution to `main` contains at least one state whose only valid symmetry is a fixed-point-free $180^\circ$ rotation (`R2` centered on an edge midpoint or vertex, even size, no capitol).

Run 3: `--min-max-size 26 --clue-bounds` (all 441 transforms)
- Model: `25,121` vars / `180,725` constraints.
- Outcome: `UNKNOWN` at `3,000s` (`485k` conflicts on Worker 0).
- Takeaway: deduplicated `edge_images` cut constraints from `415k` to `180k`, but LP workers and unpruned root/transform pairs bottlenecked the search.

---

## 3. Mathematical Deductions & Compile-Time Proofs on `main`

All deductions below are derived strictly from the rules of Capitols and the 33 clues of `main`, and are computed automatically at model-build time.

### 3.1 Forced Singleton Capitol at `(5, 10)` and Right-Edge Cascade
- Any cell with clue `1` (`d = 1`) must have a neighbor at distance `0` (`cap = 1`) with state size `1` (`sz = 1`, singleton capitol).
- Look at clue `(4, 10) = 1`: its orthogonal neighbors are `(3, 10)`, `(4, 9)`, and `(5, 10)`.
  - `(3, 10)` is adjacent to `(2, 10) = 14`, so `sz[(3, 10)] >= 2` (if `(3, 10)` were a singleton capitol, `d[(2, 10)]` would be $\le 1$).
  - `(4, 9)` is adjacent to `(4, 8) = 4`, so `sz[(4, 9)] >= 2`.
  - Therefore, `(5, 10)` is uniquely forced to be a singleton capitol (`sz[(5, 10)] = 1`, `d[(5, 10)] = 0`, `cap[(5, 10)] = 1`)!
- Because `(5, 10)` is a singleton capitol (`d = 0, sz = 1`), its neighbors `(5, 9)` and `(6, 10)` have edge weight $\min(\text{sz}, 1) = 1$ to `(5, 10)`, forcing:
  - $d[(5, 9)] \le 1$ and $d[(6, 10)] \le 1$ (in fact, since neither can be a capitol, $d[(5, 9)] = d[(6, 10)] = 1$).
  - Since `(6, 9) = 9` is adjacent to `(5, 9)` ($d \le 1$) and `(6, 10)` ($d \le 1$):
    $$\min(\text{sz}[(6, 9)], \text{sz}[(5, 9)]) \ge 9 - 1 = 8 \implies \text{sz}[(6, 9)] \ge 8, \quad \text{sz}[(5, 9)] \ge 8$$
  - Since `(7, 10) = 13` is adjacent to `(6, 10)` ($d \le 1$):
    $$\min(\text{sz}[(7, 10)], \text{sz}[(6, 10)]) \ge 13 - 1 = 12 \implies \text{sz}[(7, 10)] \ge 12, \quad \text{sz}[(6, 10)] \ge 12$$
  - Since `(2, 10) = 14` is 2 steps from `(4, 10) = 1` via `(3, 10)`:
    $$2 \cdot \text{sz}[(3, 10)] \ge 14 - 1 = 13 \implies \text{sz}[(3, 10)] \ge 7$$

### 3.2 General Path Vertex-Cover & Adjacent-Pair Inequalities
For any path $u_0 - u_1 - \dots - u_k$ of $k$ edges between two clues $u_0, u_k$ with distance gap $\Delta = |d[u_0] - d[u_k]|$:
$$\Delta \le \sum_{t=0}^{k-1} w(u_t, u_{t+1}) = \sum_{t=0}^{k-1} \min(\text{sz}[u_t], \text{sz}[u_{t+1}])$$
1. Endpoint & Vertex-Cover Bounds (missed by interior-only $2\sum_{t=1}^{k-1} \text{sz}[u_t] \ge \Delta$):
   - For $k = 2$ ($u_0 - u_1 - u_2$), since $w(u_0, u_1) \le \text{sz}[u_0]$ and $w(u_1, u_2) \le \text{sz}[u_2]$:
     $$\text{sz}[u_0] + \text{sz}[u_2] \ge \Delta, \qquad \text{sz}[u_0] + \text{sz}[u_1] \ge \Delta, \qquad \text{sz}[u_1] + \text{sz}[u_2] \ge \Delta$$
     On `main`:
     - `(2, 2) = 7` and `(3, 3) = 51` ($\Delta = 44$): $\text{sz}[(2, 2)] + \text{sz}[(3, 3)] \ge 44$ (at least one is $\ge 22$).
     - `(5, 3) = 15` and `(3, 3) = 51` ($\Delta = 36$): $\text{sz}[(5, 3)] + \text{sz}[(3, 3)] \ge 36$ (at least one is $\ge 18$).
     - `(8, 8) = 45` and `(10, 8) = 77` ($\Delta = 32$): $\text{sz}[(8, 8)] + \text{sz}[(10, 8)] \ge 32$ (at least one is $\ge 16$).
2. Adjacent-Pair Disjunction (`w(u_t, u_{t+1}) >= ceil(gap / k)`):
   - Since the $k$ edges sum to $\ge \Delta$, at least one edge $(u_t, u_{t+1})$ has weight $\ge M = \lceil \Delta / k \rceil$, which requires both endpoints of that edge to have $\text{sz}[u_t] \ge M \land \text{sz}[u_{t+1}] \ge M$ (not merely a single interior vertex).

### 3.3 Root-Region Pruning (`valid_region_from_root(r, gi)`)
Every state $S$ has a unique minimum-index root cell $r = \min(S)$ (`b[r, r] == 1`) and an involution symmetry $g \in G$ (`T[r, gi] == 1`).
- Every cell $u \in S$ must satisfy both $u \ge r$ and $g(u) \ge r$ (since $g(u) \in S$ and $r = \min(S)$).
- Corollary (Cell-Centered $R_2$ Root Theorem): Under a $180^\circ$ rotation $R_2$ fixing cell $r$, every cell $u > r$ maps to $g(u) = 2r - u < r$. Thus no state of size $\ge 2$ can ever have $r$ as its root under the $R_2$ centered at $r$!
- Computing the maximal connected $g$-invariant subgraph $R_{\text{root}}(r, g) \subseteq \{u \ge r : g(u) \ge r\}$ for each `(r, gi)`:
  - Prunes 1,302 impossible `(r, gi)` root-transform pairs via binary clauses `m.add_bool_or([b[r, r].Not(), T[r, gi].Not()])`.
  - Prunes 208 geometrically impossible `b[j, r]` variables to `0`.
  - Cuts the valid `(r, gi)` pairs supporting `--min-max-size 26` nearly in half (`7,990` $\to$ `4,558`).

---

## 4. CP-SAT Bottlenecks & 20-Second Diagnostic Suite Results

### 4.1 Two Code Issues Found in [`solve.py`](solve.py)
1. Uncached Path-Disjunction Boolean Bloat:
   - Out of 1,149 monotone paths of length $2 \le \text{dist} \le 6$ between clues in `main`, 986 paths were already trivially satisfied at level 0 because at least one interior cell already had `sz_lb[j] >= min_max_sz`.
   - However, `solve.py` still allocated fresh anonymous boolean variables `[m.new_bool_var("") for _ in mid_cells]` for all 1,149 paths without caching `(j, min_max_sz)`, creating 3,841 redundant variables.
   - Skipping already-satisfied paths and caching `(j, k)` drops model variables from `25,121` to `21,280`.
2. Missing Converse Non-Singleton Clause:
   - `solve.py` had `is_sing[i] => not same[i, j]` for all neighbors `j`, but omitted the converse:
     ```python
     m.add_bool_or([is_sing[i]] + [same[i, j] for j in nbrs(i)])
     ```
   - Since 119 of 121 cells have `sz_lb[i] >= 2` (`is_sing[i] == 0`), this immediately gives 119 short clauses forcing at least one active incident `same` edge at every non-singleton cell.

### 4.2 LP Worker Starvation & `linearization_level=0`
Note that `solver.NumConflicts()` on timeout only reports Worker 0 (`default_lp`). Inspecting the full per-worker subsolver tables via `--cpsat-log` in 20s runs on `main` revealed:
- CP-SAT's `'max_lp'` worker builds a 234,854-row, 23,396-column LP (`736,897` nonzeros) and spends the entire run inside simplex iterations, producing only 89–100 conflicts in 20 seconds while consuming a full core and memory bandwidth.
- Setting `solver.parameters.linearization_level = 0` (`--lin-level 0`) raised conflict throughput across the portfolio in these 20s runs. Higher throughput did not translate into faster solves: see sections 1 and 7, where `lin_level=0` was slower to a proof and to the solution.


Raw size: `25,121` vars / `180,725` constraints for the baseline, `21,372` / `255,091` for the three enhanced configs. Presolved: `23,966` / `138,440` for the baseline, `23,684` / `140,689` for the enhanced configs.

Per-worker order below: `default_lp` / `no_lp` / `fs_random_no_lp` / `quick_restart` / `quick_restart_no_lp`.

- Baseline (Run 3)
  - Fixed during search: `39`; shared bounds `~3,100`; learned binary clauses `~340k`.
  - Total conflicts (20s, 8 workers): `~21,500`.
  - Per worker: `1,609` / `16,784` / `13,348` / `376` / `~3,500`.
- Enhanced, default LP
  - Fixed: `47`; shared bounds `6,084`; binary clauses `682,958`.
  - Total: `23,100`.
  - Per worker: `674` / `10,023` / `7,978` / `201` / `4,122`.
- Enhanced, `cut_level=0`
  - Fixed: `44`; shared bounds `12,229`; binary clauses `609,628`.
  - Total: `30,984`.
  - Per worker: `2,216` / `10,692` / `12,196` / `796` / `4,955`.
- Enhanced, `lin_level=0`
  - Fixed: `76-115`; shared bounds `7,061`; binary clauses `1,279,786`.
  - Total: `39,108`.
  - Per worker: `2,904` / `13,287` / `11,265` / `5,596` / `5,904`.

---

## 5. What Was Slower or a Mirage (Negative Results & Lessons Learned)

### 5.1 Why Bundling ~60,000 Extra Clauses Slowed Down `bench.py`
In `scratch/test_new_clauses.py`, we tested adding all candidate redundant clauses at once:
- `T[i, gi] => not is_odd[i]` across all 18,600 `(i, gi)` pairs (`~18,600` clauses)
- Mirror-axis intersection `T[i, gi] => OR_{f in fixed_at[gi]} T[f, gi]` (`~6,600` long clauses)
- Extra `min_max_size` implications and $2\times 2$ cycle clauses

While presolve compressed the model well, total raw constraints rose from `180,725` to `255,091`. On `bench.py -n 6 7 8` (15 random 6×6, 7×7, 8×8 instances):
- Committed [`solve.py`](solve.py): `total 158.4s | solved 14/15 | median 4.4s`
- All extra clauses bundled (`build_enhanced`): `total 202.2s | solved 14/15 | median 8.4s`

Lesson: Avoid flooding the model with `O(| cells | * | transforms |)` redundant clauses (`is_odd` per transform, mirror-axis disjunctions). Keep only clauses that remove variables or prune domain/transform pairs at build time.

### 5.2 Subtle Equivalence Bug Caught by `tests/test_equivalence.py`
When pruning impossible capitol transforms `gi` at cell `i` (where the maximum odd-sized region `cap_ub < c_sz_lb`), simply `continue`-ing in the `ors` loop of `cap[i]` omitted the implication `a => cap[i]`. That allowed the solver to pick `T[i, gi] == 1` with `i` as the unique fixed cell of `gi` while setting `cap[i] = 0` (`test_4x4_solution_sets_match_reference` caught this immediately: `0 lost, 2 invented`).
- Fix: When `not cap_possible[i]` or `cap_ub < c_sz_lb`, cell `i` cannot be the unique fixed point of `gi`, so we must enforce `m.add(T[i, gi] == 0)` (if `fixed_at[gi] == [i]`) or `m.add_bool_or([T[i, gi].Not()] + [eq(i, j) for j in others])` (if `others` is non-empty). This also eliminates the auxiliary variable `a` on all 34 non-capitol cells.

### 5.3 Why Option A (Two-Stage Pure-Metric Relaxation) Was a Mirage
We tested Option A (two-stage metric relaxation, Benders cut-set generation, spanning-tree subtree flow, size-multiplicity equations, and root-representative area sums):
1. Pure Metric (`sz, d, cap, same`, without `b[i, r]` and `T[i, gi]`):
   - `vars = 1,825, cons = 8,372` $\to$ Solves in `0.32s`!
   - Even with exact 1-, 2-, and 3-cell boundary rules and Benders cut-set rounds (`sz[i] >= |C| + 1 => OR_{e in bnd(C)} same[e]`), each round took only `0.35s`.
   - Why it was a mirage: Without global area conservation, the solver satisfied the large distance clues (`77, 61, 51, 45, 30, 28, 26`) by assigning `sz = 77`, `sz = 65`, `sz = 51`, `sz = 48`, `sz = 29`, `sz = 27` to separate 4- to 6-cell blobs—claiming a total board area of $77 + 65 + 51 + 48 + 29 + 27 = 297 > 121$ cells!
2. Adding Area Conservation ($\sum_{\text{states}} |S| = 121$) Back to Stage 1:
   - We tested four different ways to enforce that state sizes actually fit inside the 121-cell board without `T[i, gi]`:
     - Root assignment `b[i, r]` (`with_partition=True`): `vars = 10,067, cons = 19,665` $\to$ `UNKNOWN` at `15s`.
     - Spanning-Tree Subtree Flow (`sub[i] == 1 + sum f_tree[k, i]`): `vars = 3,387, cons = 9,203` $\to$ `UNKNOWN` at `15s`.
     - Exact Size-Multiplicity ($\sum_i [\text{sz}[i] == s] = s \cdot \text{num\_states}[s]$): `vars = 16,157, cons = 8,855` $\to$ `UNKNOWN` at `15s`.
     - Root-Representative Sum ($\sum_i \text{rep\_sz}[i] == 121$): `vars = 2,748, cons = 10,319` $\to$ `UNKNOWN` at `15s`.
     - Threshold Count Inequalities ($\text{sz}[i] \ge K \implies \sum_j [\text{sz}[j] \ge K] \ge K$): `UNKNOWN` at `5s`.
   - Why: Once total state area is capped at 121 cells, the clues `(10, 8)=77, (10, 9)=61, (8, 8)=45` and `(3, 3)=51, (3, 0)=28, (4, 1)=22` and `(8, 0)=30, (10, 1)=26` cannot use independent giant states; their shortest-path trees are forced to share and wind through the same 1–2 large states with multi-step internal paths. Finding that joint partition + shortest-path layout *is* the core NP-hard bottleneck of `main`.

---

## 6. Ported Lean Improvements in [`solve.py`](solve.py)

1. Port Only the Lean, Zero-Bloat Improvements into [`solve.py`](solve.py):
   - Deduplicate & filter the `sz[j] >= min_max_sz` path booleans (`-3,870` vars).
   - Add compile-time iterative domain tightening (fixing `(5, 10)` as a singleton capitol and tightening `sz_lb`, `d_lb`, `d_ub`) + vertex-cover path bounds (`sz[u_0] + sz[u_2] >= gap`).
   - Add root-region pruning (`valid_region_from_root(r, gi)`) to prune impossible `b[j, r]` and `T[i, gi]` variables and tighten `--min-max-size`.
   - Add the 121 non-singleton incident-edge clauses `m.add_bool_or([is_sing[i]] + [same[i, j] for j in nbrs(i)])` and Dijkstra capitol reachability.
   - Add CLI flags `--require-free-r2` and `--lin-level` (`--lin-level` is available but `0` turned out slower; see section 7).
2. Verify on `pytest` and `bench.py -n 6 7 8`:
   - Ensure the lean model passes all 46 equivalence tests (`pytest -m "slow or not slow"`) and improves on the `158.4s` (`median 4.4s`) baseline in `bench.py -n 6 7 8`.

---

## 7. Solved `main` Puzzle (`answer = 408951` in `1,580.3s` Default LP vs `2,491.7s` `lin_level=0`)

Running the lean, zero-bloat model with default LP:
```bash
python solve.py main --min-max-size 26 --require-free-r2 --time 2500 --cpsat-log --logfile main_default_lp.log
```
solved `main` in Refinement Round 1 (zero full-group capitol cuts needed) in **`1,580.3s`** (`26.3m`, found by worker `'no_lp'` at `1,579.71s`), beating the `--lin-level 0` run (`2,491.7s` / `41.5m`) by **911 seconds (`36.6%` faster)**.

### Possible Reasons Default LP Beat `--lin-level 0` (`1,580.3s` vs `2,491.7s`)
The subsolver tables from the `1,580.3s` run suggest a reason why `linearization_level=0` may have hurt both the `--max-size 35` infeasibility proof (`61s` vs `38s`) and the full `main` solve (`2,491.7s` vs `1,580.3s`):
- Only `'max_lp'` builds the giant LP (`162,338` rows, `23,526` cols, producing just `666` conflicts and `7` shared bounds in 26 minutes).
- Meanwhile, `'default_lp'` (`1,069` rows), `'quick_restart'` (`1,337` rows), and `'probing'` (`1,459` rows) maintain compact, fast LPs (`7,469` cols), generating **1.80 million cuts** (`IB`, `MIR`, `CG`, `ZERO_HALF`) and exporting **160,898 clauses** and **402 bound tightenings** into the shared repository.
- CP-SAT's portfolio already includes three dedicated no-LP SAT workers (`'no_lp'`, `'fs_random_no_lp'`, `'quick_restart_no_lp'`). In the default-LP run `'no_lp'` imported `69,955` clauses and read `803,711` binary clauses before finding the solution at `1,579.71s`. These are clauses from all workers, not only LP-derived ones, and we did not capture the same tables for the `lin-level 0` run, so the link between the LP cuts and the slowdown is a hypothesis. Evidence for the slowdown itself is one `main` run each (`1,580.3s` vs `2,491.7s`) plus three runs each on the `--max-size 35` proof (medians `38s` vs `61s`), against run-to-run spread of 10-30s on that proof.

### State Sizes Grid
```text
  4   4   4   4   7  10   5   5   5   5   5
  7   7   7   7   7  10  10  10   7   7   7
  7  37  37  37  37  10   1  10   7  12   7
 11  37  37  37  37  10  10  10   7  12   7
 11  11  37  37  37  37   5  10   2  12  12
 11  37  37  37  37  37   5   5   2  12   1
 11  11  11   2  37  37  37   5   5  12  12
 10  10  11   2   4  37  37  37  37  37  12
 10  11  11   4   4  37  37  37  37  12  12
 10  10  11  10   4   3  37  37  37  37  12
 10  10  10  10   3   3  37  37  37  37  12
```

### State Root Labels Grid
```text
  0   0   0   0   4   5   6   6   6   6   6
  4   4   4   4   4   5   5   5  19  19  19
  4  23  23  23  23   5  28   5  19  31  19
 33  23  23  23  23   5   5   5  19  31  19
 33  33  23  23  23  23  50   5  52  31  31
 33  23  23  23  23  23  50  50  52  31  65
 33  33  33  69  23  23  23  50  50  31  31
 77  77  33  69  81  23  23  23  23  23  31
 77  33  33  81  81  23  23  23  23  31  31
 77  77  33  77  81 104  23  23  23  23  31
 77  77  77  77 104 104  23  23  23  23  31
```

### Breakdown of the 16 States & 8 Capitols
- 8 Capitols: `[(1, 2), (1, 9), (2, 6), (5, 7), (5, 10), (6, 1), (6, 5), (10, 5)]`
  - Includes `(5, 10)` (`root=65`, `sz=1`), proved at compile time by our domain-tightening pass, and `(2, 6)` (`root=28`, `sz=1`).
  - The largest state (`root=23` at `(2, 1)`, `sz=37`) is $180^\circ$ rotationally symmetric (`R2`) about its capitol at `(6, 5)`.
  - Two states—`root=31` at `(2, 9)` (`sz=12`) and `root=5` at `(0, 5)` (`sz=10`)—have only a fixed-point-free $180^\circ$ rotation (`R2` with 0 fixed cells), directly confirming the theorem from Run 2 (`--no-free-r2` $\to$ `INFEASIBLE`).
- Row Sums: `[58, 86, 202, 215, 211, 221, 180, 234, 212, 208, 206]`
- Final Answer (sum of squares of row sums): `408951`
