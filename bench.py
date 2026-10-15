"""Benchmark ladder: random solvable Capitols puzzles of increasing size.

Each instance is built from a random valid partition (so a solution is known to exist), then a random
subset of its true distances becomes the clues. Every solver result is re-checked with the independent
brute-force rules in tests/reference.py, so a speedup can't come from an invalid answer.

    python bench.py                      # N=6,7,8 x 5 seeds, 60s limit each
    python bench.py -n 7 8 --seeds 3 --time 120
"""
import argparse
import concurrent.futures
import functools
import itertools
import random
import time

import solve as S
from tests import reference as R


def random_partition(N, rnd):
    """Random partition of the grid into symmetric connected states (singletons fill any gaps)."""
    G = [g for g in S.transforms(N) if g[0] != 'R4']
    cells = [(r, c) for r in range(N) for c in range(N)]
    free = set(cells)
    part = []
    while free:
        u = rnd.choice(sorted(free))
        target = rnd.choice([1, 2, 3, 3, 4, 5, 5, 7, 9, 11, 13])
        best = frozenset([u])
        for _ in range(30):
            g = rnd.choice(G)[2]
            gu = g(*u)
            if gu == u:
                state = {u}
            elif gu in free and gu in R.neighbours(u, N):
                state = {u, gu}
            else:
                continue
            while len(state) < target:
                options = [x for s in state for x in R.neighbours(s, N)
                           if x in free and x not in state and g(*x) in free
                           and 0 <= g(*x)[0] < N and 0 <= g(*x)[1] < N]
                if not options:
                    break
                x = rnd.choice(options)
                state |= {x, g(*x)}
            if len(state) > len(best) and len(state) <= target + 1:
                best = frozenset(state)
        part.append(best)
        free -= best
    return part


@functools.cache
def capitol_of(state, N):
    info = R.classify(state, N)
    assert info is not None, "generated state is not a valid state"
    return info[1]


def make_instance(N, seed, density=0.35):
    rnd = random.Random(seed)
    while True:
        part = random_partition(N, rnd)
        dist = R.distances(part, N, lambda s: capitol_of(s, N))
        if dist is None or len(part) < N or len(part) > N * N // 2:
            continue  # want a few capitols and a few states, not a sea of singletons
        cells = sorted(dist)
        k = max(4, round(density * N * N))
        clues = {c: dist[c] for c in rnd.sample(cells, k)}
        if any(d == 0 for d in clues.values()):
            return clues


def verify(N, clues, res):
    """Independent check of a solver result against the reference rules."""
    _, caps, _, states = res
    part = [frozenset(cells) for cells in states.values()]
    assert sorted(x for s in part for x in s) == [(r, c) for r in range(N) for c in range(N)]
    assert all(capitol_of(s, N) is None or capitol_of(s, N) in caps for s in part)
    dist = R.distances(part, N, lambda s: capitol_of(s, N))
    assert dist is not None
    for cell, want in clues.items():
        assert dist[cell] == want, f"clue {cell}={want} but reference distance is {dist[cell]}"


def run_one(job):
    N, seed, density, workers, tlimit = job
    clues = make_instance(N, seed, density)
    t = time.time()
    try:
        res = S.run(N, clues, workers=workers, tlimit=tlimit)
        if res:
            verify(N, clues, res)
            status = 'solved'
        else:
            status = 'INFEASIBLE (bug: a solution exists)'
    except TimeoutError:
        status = 'timeout'
    return N, seed, len(clues), time.time() - t, status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-n', type=int, nargs='+', default=[6, 7, 8])
    ap.add_argument('--seeds', type=int, default=5)
    ap.add_argument('--first-seed', type=int, default=0)
    ap.add_argument('--repeat', type=int, default=1, help='solve each instance this many times')
    ap.add_argument('--time', type=int, default=60, help='time limit per instance, seconds')
    ap.add_argument('--workers', type=int, default=8, help='CP-SAT threads per instance')
    ap.add_argument('--jobs', type=int, default=1,
                    help='instances solved at once; use with --workers 1 for reproducible, parallel timings')
    ap.add_argument('--density', type=float, default=0.35, help='fraction of cells given as clues')
    args = ap.parse_args()
    jobs = [(N, seed, args.density, args.workers, args.time)
            for N in args.n
            for seed, _ in itertools.product(range(args.first_seed, args.first_seed + args.seeds), range(args.repeat))]
    total = timeouts = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for N, seed, k, dt, status in pool.map(run_one, jobs):
            total += dt
            timeouts += status == 'timeout'
            print(f'N={N} seed={seed} clues={k:3d}  {dt:6.1f}s  {status}', flush=True)
    print(f'total {total:.1f}s  ({timeouts} timeouts)')


if __name__ == '__main__':
    main()
