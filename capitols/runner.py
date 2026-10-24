import logging
import threading
import time

from ortools.sat.python import cp_model

from .config import SolverConfig
from .model import CapitolsModel
from .problem import prepare
from .symmetry import full_group_ok

log = logging.getLogger("capitols")


def _solver_for(cfg):
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = cfg.workers
    solver.parameters.max_time_in_seconds = cfg.tlimit
    if cfg.lin_level is not None:
        solver.parameters.linearization_level = cfg.lin_level
    for name, value in cfg.params.items():
        if isinstance(
            value, (list, tuple)
        ):  # repeated fields (e.g. ignore_subsolvers) can't be assigned
            getattr(solver.parameters, name).extend(value)
        else:
            setattr(solver.parameters, name, value)
    solver.parameters.log_search_progress = cfg.verbose
    # route through log_callback only, else every line prints twice
    solver.parameters.log_to_stdout = False
    if cfg.verbose:
        solver.log_callback = lambda line: log.debug("[cp-sat] %s", line)
    return solver


def _solve_with_heartbeat(solver, model, cfg):
    stop = threading.Event()
    t_solve = time.time()

    def beat():
        while not stop.wait(cfg.heartbeat):
            log.info(
                "[solve] still searching... %.0fs elapsed (limit %ds)",
                time.time() - t_solve,
                cfg.tlimit,
            )

    threading.Thread(target=beat, daemon=True).start()
    log.info(
        "[solve] launching CP-SAT with %d workers, time limit %ds",
        cfg.workers,
        cfg.tlimit,
    )
    st = solver.Solve(model)
    stop.set()
    log.info(
        "[solve] finished: status=%s in %.1fs (conflicts=%d, branches=%d)",
        solver.StatusName(st),
        time.time() - t_solve,
        solver.NumConflicts(),
        solver.NumBranches(),
    )
    return st


def solve(N, clues, cuts=None, exclude=None, **options):
    """Solve once. Returns (labels, capitol cells, distances, all transforms) or None if infeasible."""
    cfg = SolverConfig(**options)
    cuts = cuts or []
    exclude = exclude or []
    log.info(
        "[build] starting model: N=%d, %d clues, %d cuts, maxd=%d",
        N,
        len(clues),
        len(cuts),
        cfg.maxd,
    )
    prob = prepare(N, clues, cfg)
    cm = CapitolsModel(prob, cuts, exclude)
    solver = _solver_for(cfg)
    st = _solve_with_heartbeat(solver, cm.m, cfg)
    if st == cp_model.UNKNOWN:
        raise TimeoutError(
            "CP-SAT hit the time limit before finding a solution or proving there is none"
        )
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None
    n, cells = prob.grid.n, prob.grid.cells
    L = [solver.Value(cm.lab[i]) for i in range(n)]
    caps = [cells[i] for i in range(n) if solver.Value(cm.cap[i])]
    D = [solver.Value(cm.d[i]) for i in range(n)]
    return L, caps, D, prob.mt.full


def run(N, clues, verbose=False, **kw):
    """Solve, check each declared capitol against the state's full symmetry group, and re-solve with a
    cut for every failure until a candidate passes."""
    cuts = []
    it = 0
    while True:
        it += 1
        log.info("=== refinement round %d ===", it)
        res = solve(N, clues, cuts=cuts, verbose=verbose, **kw)
        if res is None:
            log.warning("no solution found (model is infeasible)")
            return None
        L, caps, D, G = res
        states = {}
        for i, l in enumerate(L):
            states.setdefault(l, []).append((i // N, i % N))
        log.info(
            "[check] candidate has %d states, %d declared capitols",
            len(states),
            len(caps),
        )
        log.info(
            "[check] state sizes: %s",
            sorted((len(S) for S in states.values()), reverse=True),
        )
        bad = False
        for l, S in states.items():
            cp, hassym = full_group_ok(S, N, G)
            declared = [x for x in caps if x in S]
            assert hassym
            if declared and cp != declared[0]:
                log.info(
                    "[check] state rooted at %s (size %d) declared capitol %s but true capitol is %s -> adding cut",
                    S[0],
                    len(S),
                    declared[0],
                    cp,
                )
                cuts.append((S, declared[0]))
                bad = True
        if not bad:
            log.info("[check] all capitols verified against full symmetry groups")
            return L, caps, D, states
        log.info("[check] %d cut(s) total; re-solving", len(cuts))


def report(N, res):
    L, caps, _, states = res
    sizes = [len(states[L[i]]) for i in range(N * N)]
    for r in range(N):
        print(" ".join(f"{sizes[r * N + c]:3d}" for c in range(N)))
    print()
    for r in range(N):
        print(" ".join(f"{L[r * N + c]:3d}" for c in range(N)))
    print("caps", caps)
    rows = [sum(sizes[r * N : (r + 1) * N]) for r in range(N)]
    print("rows", rows, "answer", sum(x * x for x in rows))
