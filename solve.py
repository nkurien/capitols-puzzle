import sys, time, threading, argparse, logging
from ortools.sat.python import cp_model

log = logging.getLogger("capitols")
T0 = time.time()

def stage(m, name, t_start):
    p = m.Proto()
    log.info("[build] %-22s done in %6.2fs | vars=%7d constraints=%7d",
             name, time.time() - t_start, len(p.variables), len(p.constraints))
    return time.time()

def transforms(N):
    G = []  # (kind, params, func, fixedcells(list))
    M = 2 * (N - 1)
    for a in range(M + 1):
        G.append(('H', a, lambda r, c, a=a: (a - r, c)))
        G.append(('V', a, lambda r, c, a=a: (r, a - c)))
    for k in range(-(N - 1), N):
        G.append(('D', k, lambda r, c, k=k: (c + k, r - k)))
    for s in range(M + 1):
        G.append(('A', s, lambda r, c, s=s: (s - c, s - r)))
    for a in range(M + 1):
        for b in range(M + 1):
            G.append(('R2', (a, b), lambda r, c, a=a, b=b: (a - r, b - c)))
            if (a - b) % 2 == 0:
                G.append(('R4', (a, b), lambda r, c, a=a, b=b: ((a - b) // 2 + c, (a + b) // 2 - r)))
    return G

def fixed(g, N):
    f = g[2]
    return [(r, c) for r in range(N) for c in range(N) if f(r, c) == (r, c)]

def full_group_ok(S, N, G):
    """Return capitol cell if state S has one, else None."""
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

def solve(N, clues, maxd=400, cuts=None, exclude=None, verbose=False, workers=8, tlimit=3000, heartbeat=30):
    cuts = cuts or []
    exclude = exclude or []
    log.info('[build] starting model: N=%d, %d clues, %d cuts, maxd=%d', N, len(clues), len(cuts), maxd)
    t_build = time.time()
    G_full = transforms(N)
    # A 90-degree rotation about a centre implies the 180-degree rotation about the same centre, and
    # both fix exactly the same cells, so any state that could use R4 can use R2 with an identical
    # capitol outcome. The model therefore skips R4; G_full keeps it for the check in run().
    G = [g for g in G_full if g[0] != 'R4']
    log.info('[build] %d candidate symmetry transforms (%d with R4, used only for verification)', len(G), len(G_full))
    cells = [(r, c) for r in range(N) for c in range(N)]
    idx = {x: i for i, x in enumerate(cells)}
    n = len(cells)
    # image of every cell under every transform as a cell index (None if it leaves the grid), and
    # the cells each transform fixes: computed once instead of per cell inside the model loops
    images = [[idx.get(g[2](r, c)) for (r, c) in cells] for g in G]
    fixed_at = [[i for i in range(n) if img[i] == i] for img in images]
    m = cp_model.CpModel()
    def nbrs(i):
        r, c = cells[i]
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (r + dr, c + dc) in idx:
                yield idx[(r + dr, c + dc)]
    lab = [m.NewIntVar(0, i, f'lab{i}') for i in range(n)]
    b = {}
    for i in range(n):
        bs = []
        for r in range(i + 1):
            b[i, r] = m.NewBoolVar('')
            bs.append(b[i, r])
        m.AddExactlyOne(bs)
        m.Add(lab[i] == sum(r * b[i, r] for r in range(i + 1)))
        for r in range(i):
            m.AddImplication(b[i, r], b[r, r])
    t = stage(m, 'labels', t_build)
    edges = [(i, j) for i in range(n) for j in nbrs(i) if i < j]
    same = {}
    for (i, j) in edges:
        s = m.NewBoolVar('')
        m.Add(lab[i] == lab[j]).OnlyEnforceIf(s)
        m.Add(lab[i] != lab[j]).OnlyEnforceIf(s.Not())
        same[i, j] = same[j, i] = s
    t = stage(m, 'same-state edges', t)
    # connectivity
    h = [m.NewIntVar(0, n, '') for i in range(n)]
    for i in range(n):
        m.Add(h[i] == 0).OnlyEnforceIf(b[i, i])
        m.Add(h[i] >= 1).OnlyEnforceIf(b[i, i].Not())
        par = []
        for j in nbrs(i):
            p = m.NewBoolVar('')
            m.AddImplication(p, same[i, j])
            m.Add(h[j] < h[i]).OnlyEnforceIf(p)
            par.append(p)
        m.AddBoolOr(par + [b[i, i]])
    t = stage(m, 'connectivity', t)
    # sizes
    size = [m.NewIntVar(0, n, '') for r in range(n)]
    for r in range(n):
        m.Add(size[r] == sum(b[i, r] for i in range(r, n)))
    sz = [m.NewIntVar(1, n, '') for i in range(n)]
    for i in range(n):
        m.AddElement(lab[i], size[:i + 1], sz[i])
    t = stage(m, 'state sizes', t)
    # symmetry
    T = {}
    tvar = [m.NewIntVar(0, len(G) - 1, '') for i in range(n)]
    for i in range(n):
        opts = []
        for gi, img in enumerate(images):
            j = img[i]
            if j is None:
                continue
            v = m.NewBoolVar('')
            T[i, gi] = v
            opts.append((gi, v))
            if j != i:
                m.Add(lab[j] == lab[i]).OnlyEnforceIf(v)
        m.AddExactlyOne([v for _, v in opts])
        m.Add(tvar[i] == sum(gi * v for gi, v in opts))
    for (i, j) in edges:
        m.Add(tvar[i] == tvar[j]).OnlyEnforceIf(same[i, j])
    t = stage(m, 'symmetry', t)
    # capitols
    eqcache = {}
    def eq(i, j):
        if (i, j) in same:
            return same[i, j]
        key = (min(i, j), max(i, j))
        if key not in eqcache:
            s = m.NewBoolVar('')
            m.Add(lab[i] == lab[j]).OnlyEnforceIf(s)
            m.Add(lab[i] != lab[j]).OnlyEnforceIf(s.Not())
            eqcache[key] = s
        return eqcache[key]
    cap = [m.NewBoolVar(f'cap{i}') for i in range(n)]
    for (r, c), v in clues.items():
        # a clue of 0 sits on a capitol; any positive clue does not
        m.Add(cap[idx[(r, c)]] == (1 if v == 0 else 0))
    for i in range(n):
        ors = []
        for gi, img in enumerate(images):
            if img[i] != i:  # only transforms that fix this cell can make it a capitol
                continue
            others = [j for j in fixed_at[gi] if j != i]
            a = m.NewBoolVar('')
            m.AddImplication(a, T[i, gi])
            for j in others:
                m.AddImplication(a, eq(i, j).Not())
            # if T and all others not same -> a
            m.AddBoolOr([a, T[i, gi].Not()] + [eq(i, j) for j in others])
            ors.append(a)
        m.AddBoolOr(ors).OnlyEnforceIf(cap[i])
        for a in ors:
            m.AddImplication(a, cap[i])
    t = stage(m, 'capitols', t)
    # distances
    d = [m.NewIntVar(0, maxd, f'd{i}') for i in range(n)]
    w = {}
    for (i, j) in edges:
        v = m.NewIntVar(1, n, '')
        m.AddMinEquality(v, [sz[i], sz[j]])
        w[i, j] = w[j, i] = v
    for i in range(n):
        m.Add(d[i] == 0).OnlyEnforceIf(cap[i])
        m.Add(d[i] >= 1).OnlyEnforceIf(cap[i].Not())
        par = []
        for j in nbrs(i):
            m.Add(d[i] <= d[j] + w[i, j])
            p = m.NewBoolVar('')
            m.Add(d[i] == d[j] + w[i, j]).OnlyEnforceIf(p)
            par.append(p)
        m.AddBoolOr(par + [cap[i]])
    for (r, c), v in clues.items():
        m.Add(d[idx[(r, c)]] == v)
    t = stage(m, 'distances + clues', t)
    # cuts: forbid (state==S with capitol at cell)
    for S, capcell in cuts:
        rep = min(idx[x] for x in S)
        Sset = set(S)
        bnd = set()
        for x in S:
            for j in nbrs(idx[x]):
                if cells[j] not in Sset and j > rep:
                    bnd.add(j)
        m.AddBoolOr([b[idx[x], rep].Not() for x in S] + [b[j, rep] for j in bnd] + [cap[idx[capcell]].Not()])
    t = stage(m, 'cuts', t)
    # exclude: label vectors of partitions already found, so repeated solves enumerate every solution
    for L in exclude:
        m.AddBoolOr([b[i, L[i]].Not() for i in range(n)])
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = workers
    solver.parameters.max_time_in_seconds = tlimit
    solver.parameters.log_search_progress = verbose
    solver.parameters.log_to_stdout = False  # route through log_callback only, else every line prints twice
    if verbose:
        solver.log_callback = lambda line: log.debug('[cp-sat] %s', line)
    stop = threading.Event()
    t_solve = time.time()
    def beat():
        while not stop.wait(heartbeat):
            log.info('[solve] still searching... %.0fs elapsed (limit %ds)', time.time() - t_solve, tlimit)
    threading.Thread(target=beat, daemon=True).start()
    log.info('[solve] launching CP-SAT with %d workers, time limit %ds', workers, tlimit)
    st = solver.Solve(m)
    stop.set()
    log.info('[solve] finished: status=%s in %.1fs (conflicts=%d, branches=%d)',
             solver.StatusName(st), time.time() - t_solve, solver.NumConflicts(), solver.NumBranches())
    if st == cp_model.UNKNOWN:
        raise TimeoutError('CP-SAT hit the time limit before finding a solution or proving there is none')
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None
    L = [solver.Value(lab[i]) for i in range(n)]
    caps = [cells[i] for i in range(n) if solver.Value(cap[i])]
    D = [solver.Value(d[i]) for i in range(n)]
    return L, caps, D, G_full

def run(N, clues, verbose=False, **kw):
    cuts = []
    it = 0
    while True:
        it += 1
        log.info('=== refinement round %d ===', it)
        res = solve(N, clues, cuts=cuts, verbose=verbose, **kw)
        if res is None:
            log.warning('no solution found (model is infeasible)'); return None
        L, caps, D, G = res
        states = {}
        for i, l in enumerate(L):
            states.setdefault(l, []).append((i // N, i % N))
        log.info('[check] candidate has %d states, %d declared capitols', len(states), len(caps))
        log.info('[check] state sizes: %s', sorted((len(S) for S in states.values()), reverse=True))
        bad = False
        for l, S in states.items():
            cp, hassym = full_group_ok(S, N, G)
            declared = [x for x in caps if x in S]
            assert hassym
            if declared and cp != declared[0]:
                log.info('[check] state rooted at %s (size %d) declared capitol %s but true capitol is %s -> adding cut', S[0], len(S), declared[0], cp)
                cuts.append((S, declared[0])); bad = True
        if not bad:
            log.info('[check] all capitols verified against full symmetry groups')
            return L, caps, D, states
        log.info('[check] %d cut(s) total; re-solving', len(cuts))

def report(N, res):
    L, caps, D, states = res
    sizes = [len(states[L[i]]) for i in range(N * N)]
    for r in range(N):
        print(' '.join('%3d' % sizes[r * N + c] for c in range(N)))
    print()
    for r in range(N):
        print(' '.join('%3d' % L[r * N + c] for c in range(N)))
    print('caps', caps)
    rows = [sum(sizes[r * N:(r + 1) * N]) for r in range(N)]
    print('rows', rows, 'answer', sum(x * x for x in rows))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('which', choices=['ex', 'main'])
    ap.add_argument('--cpsat-log', action='store_true', help='stream the raw CP-SAT search log (very verbose)')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--time', type=int, default=3000, help='time limit per solve, seconds')
    ap.add_argument('--heartbeat', type=int, default=30, help='seconds between progress messages')
    ap.add_argument('--logfile', default=None, help='also write the log to this file')
    args = ap.parse_args()
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.logfile:
        handlers.append(logging.FileHandler(args.logfile))
    logging.basicConfig(level=logging.DEBUG if args.cpsat_log else logging.INFO,
                        format='%(asctime)s %(message)s', datefmt='%H:%M:%S', handlers=handlers)
    which = args.which
    if which == 'ex':
        clues = {(0,0):12,(0,4):17,(1,0):8,(1,2):6,(1,3):3,(3,1):7,(3,2):6,(3,4):10,(4,0):25,(4,4):0}
        N = 5
    else:
        N = 11
        raw = {(0,1):8,(0,2):4,(0,5):11,(0,9):5,(1,7):11,(2,2):7,(2,5):1,(2,10):14,
               (3,0):28,(3,3):51,(3,6):1,(3,8):6,(4,1):22,(4,8):4,(4,10):1,
               (5,3):15,(5,5):10,(5,7):0,(6,0):11,(6,2):11,(6,9):9,
               (7,2):17,(7,4):14,(7,7):10,(7,10):13,(8,0):30,(8,5):6,(8,8):45,
               (9,3):10,(10,1):26,(10,5):0,(10,8):77,(10,9):61}
        clues = raw
    try:
        res = run(N, clues, verbose=args.cpsat_log, workers=args.workers, tlimit=args.time, heartbeat=args.heartbeat)
    except TimeoutError as e:
        log.warning('%s', e)
        sys.exit(1)
    if res:
        report(N, res)