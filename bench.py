"""Benchmark ladder: random solvable Capitols puzzles of increasing size.

Each instance is built from a random valid partition (so a solution is known to exist), then a random
subset of its true distances becomes the clues. Every solver result is re-checked with the independent
brute-force rules in tests/reference.py, so a speedup can't come from an invalid answer.

    python bench.py                      # N=6,7,8 x 5 seeds, 60s limit each
    python bench.py -n 7 8 --seeds 3 --time 120
    python bench.py -n 8 9 --profile big --workers 1 --jobs 8 --time 120 --opt tclauses=False

The 'small' profile (default) only makes states of up to ~14 cells. The 'big' profile first grows 1-3
states of 15-40 cells and uses fewer clues, which is closer to the real 11x11 puzzle (clue 77 sits
three steps from a capitol, forcing a state of 26+ cells).
"""

import argparse
import ast
import concurrent.futures
import functools
import itertools
import random
import statistics
import time

import solve as S
from tests import reference as R

SMALL_TARGETS = [1, 2, 3, 3, 4, 5, 5, 7, 9, 11, 13]
BIG_MIN, BIG_MAX = 15, 40


def grow(u, target, free, G, N, rnd, tries=30):
    """Best symmetric connected state containing u with at most target + 1 cells, using free cells."""
    best = frozenset([u])
    for _ in range(tries):
        g = rnd.choice(G)[2]
        gu = g(*u)
        if gu == u:
            state = {u}
        elif gu in free and gu in R.neighbours(u, N):
            state = {u, gu}
        else:
            continue
        while len(state) < target:
            options = [
                x
                for s in state
                for x in R.neighbours(s, N)
                if x in free
                and x not in state
                and g(*x) in free
                and 0 <= g(*x)[0] < N
                and 0 <= g(*x)[1] < N
            ]
            if not options:
                break
            x = rnd.choice(options)
            state |= {x, g(*x)}
        if len(state) > len(best) and len(state) <= target + 1:
            best = frozenset(state)
    return best


def random_partition(N, rnd, profile="small"):
    """Random partition of the grid into symmetric connected states (singletons fill any gaps)."""
    G = [g for g in S.transforms(N) if g[0] != "R4"]
    cells = [(r, c) for r in range(N) for c in range(N)]
    free = set(cells)
    part = []
    if profile == "big":
        b_min = min(BIG_MIN, max(5, N * N // 4))
        b_max = max(b_min, min(BIG_MAX, N * N // 2))
        for _ in range(rnd.randint(1, 3)):
            target = rnd.randint(b_min, b_max)
            for _ in range(
                20
            ):  # a big state needs a good seed cell and transform; retry a few
                state = grow(
                    rnd.choice(sorted(free)), target, free, G, N, rnd, tries=60
                )
                if len(state) >= b_min:
                    part.append(state)
                    free -= state
                    break
    while free:
        u = rnd.choice(sorted(free))
        best = grow(u, rnd.choice(SMALL_TARGETS), free, G, N, rnd)
        part.append(best)
        free -= best
    return part


@functools.cache
def capitol_of(state, N):
    info = R.classify(state, N)
    assert info is not None, "generated state is not a valid state"
    return info[1]


def make_instance(N, seed, density=0.35, profile="small"):
    rnd = random.Random(seed)
    b_min = min(BIG_MIN, max(5, N * N // 4))
    while True:
        part = random_partition(N, rnd, profile)
        dist = R.distances(part, N, lambda s: capitol_of(s, N))
        if dist is None or len(part) < N or len(part) > N * N // 2:
            continue  # want a few capitols and a few states, not a sea of singletons
        if profile == "big" and max(len(s) for s in part) < b_min:
            continue
        cells = sorted(dist)
        k = max(4, round(density * N * N))
        clues = {c: dist[c] for c in rnd.sample(cells, k)}
        if any(d == 0 for d in clues.values()):
            return clues


def verify(N, clues, res):
    """Independent check of a solver result against the reference rules."""
    _, caps, _, states = res
    part = [frozenset(cells) for cells in states.values()]
    assert sorted(x for s in part for x in s) == [
        (r, c) for r in range(N) for c in range(N)
    ]
    assert all(capitol_of(s, N) is None or capitol_of(s, N) in caps for s in part)
    dist = R.distances(part, N, lambda s: capitol_of(s, N))
    assert dist is not None
    for cell, want in clues.items():
        assert dist[cell] == want, (
            f"clue {cell}={want} but reference distance is {dist[cell]}"
        )


def run_one(job):
    N, seed, density, profile, workers, tlimit, opts = job
    clues = make_instance(N, seed, density, profile)
    t = time.time()
    try:
        res = S.run(N, clues, workers=workers, tlimit=tlimit, **opts)
        if res:
            verify(N, clues, res)
            status = "solved"
        else:
            status = "INFEASIBLE (bug: a solution exists)"
    except TimeoutError:
        status = "timeout"
    return N, seed, len(clues), time.time() - t, status


def parse_opt(s):
    """'key=value' -> (key, python literal value); bare words stay strings."""
    key, _, raw = s.partition("=")
    try:
        return key, ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return key, raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, nargs="+", default=[6, 7, 8])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--first-seed", type=int, default=0)
    ap.add_argument(
        "--repeat", type=int, default=1, help="solve each instance this many times"
    )
    ap.add_argument(
        "--time", type=int, default=60, help="time limit per instance, seconds"
    )
    ap.add_argument(
        "--workers", type=int, default=8, help="CP-SAT threads per instance"
    )
    ap.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="instances solved at once; use with --workers 1 for reproducible, parallel timings",
    )
    ap.add_argument(
        "--profile",
        choices=["small", "big"],
        default="small",
        help="'big' seeds 1-3 states of 15-40 cells, like the real puzzle",
    )
    ap.add_argument(
        "--density",
        type=float,
        default=None,
        help="fraction of cells given as clues (default 0.35 small, 0.27 big)",
    )
    ap.add_argument(
        "--opt",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="extra keyword argument for solve.run, e.g. --opt tclauses=False (repeatable)",
    )
    args = ap.parse_args()
    density = (
        args.density
        if args.density is not None
        else (0.27 if args.profile == "big" else 0.35)
    )
    opts = dict(parse_opt(s) for s in args.opt)
    print(f"profile={args.profile} density={density} opts={opts}", flush=True)
    jobs = [
        (N, seed, density, args.profile, args.workers, args.time, opts)
        for N in args.n
        for seed, _ in itertools.product(
            range(args.first_seed, args.first_seed + args.seeds), range(args.repeat)
        )
    ]
    total = timeouts = 0
    times = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for N, seed, k, dt, status in pool.map(run_one, jobs):
            total += dt
            timeouts += status == "timeout"
            times.append(dt)
            print(f"N={N} seed={seed} clues={k:3d}  {dt:6.1f}s  {status}", flush=True)
    print(
        f"total {total:.1f}s  solved {len(jobs) - timeouts}/{len(jobs)}  "
        f"median {statistics.median(times):.1f}s  ({timeouts} timeouts)"
    )


if __name__ == "__main__":
    main()
