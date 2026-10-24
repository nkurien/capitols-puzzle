import argparse
import ast
import logging
import sys

from .puzzles import PUZZLES
from .runner import report, run

log = logging.getLogger("capitols")


def parse_param(s):
    """'name=value' -> (name, python literal); bare words stay strings."""
    name, _, raw = s.partition("=")
    try:
        return name, ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return name, raw


def build_parser():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=sorted(PUZZLES))
    ap.add_argument(
        "--cpsat-log",
        action="store_true",
        help="stream the raw CP-SAT search log (very verbose)",
    )
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--time", type=int, default=3000, help="time limit per solve, seconds"
    )
    ap.add_argument(
        "--heartbeat", type=int, default=30, help="seconds between progress messages"
    )
    ap.add_argument(
        "--max-size",
        type=int,
        default=None,
        help="cap state sizes to narrow the search",
    )
    ap.add_argument(
        "--no-free-r2",
        action="store_true",
        help="drop 180-degree rotations that fix no cell to narrow the search",
    )
    ap.add_argument(
        "--require-free-r2",
        action="store_true",
        help="require at least one state to use a fixed-point-free 180-degree rotation",
    )
    ap.add_argument(
        "--min-max-size",
        type=int,
        default=None,
        help="require some state of at least this size",
    )
    ap.add_argument(
        "--clue-bounds",
        action="store_true",
        help="add size lower bounds derived from the clues",
    )
    ap.add_argument(
        "--no-tclauses",
        action="store_true",
        help="disable clause-based transform sharing",
    )
    ap.add_argument(
        "--no-edge-images",
        action="store_true",
        help="disable edge-image symmetry clauses",
    )
    ap.add_argument(
        "--lab-channel", action="store_true", help="enable label channelling on edges"
    )
    ap.add_argument(
        "--eq-channel", action="store_true", help="enable eq channelling for capitols"
    )
    ap.add_argument(
        "--no-clue-cuts",
        action="store_true",
        help="disable clue-derived path sums and singleton cuts",
    )
    ap.add_argument(
        "--no-sound-bounds", action="store_true", help="disable sound distance bounds"
    )
    ap.add_argument(
        "--big-first",
        action="store_true",
        help="branch on size of largest clue cells first",
    )
    ap.add_argument(
        "--lin-level",
        type=int,
        default=None,
        help="CP-SAT linearization_level parameter (e.g. 0 to disable LP relaxation)",
    )
    ap.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="set a raw CP-SAT parameter, e.g. --param symmetry_level=0 (repeatable)",
    )
    ap.add_argument("--logfile", default=None, help="also write the log to this file")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.logfile:
        handlers.append(logging.FileHandler(args.logfile))
    logging.basicConfig(
        level=logging.DEBUG if args.cpsat_log else logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=handlers,
    )
    N, clues = PUZZLES[args.which]
    try:
        res = run(
            N,
            clues,
            verbose=args.cpsat_log,
            workers=args.workers,
            tlimit=args.time,
            heartbeat=args.heartbeat,
            max_size=args.max_size,
            no_free_r2=args.no_free_r2,
            require_free_r2=args.require_free_r2,
            min_max_size=args.min_max_size,
            clue_bounds=args.clue_bounds,
            tclauses=not args.no_tclauses,
            edge_images=not args.no_edge_images,
            lab_channel=args.lab_channel,
            eq_channel=args.eq_channel,
            clue_cuts=not args.no_clue_cuts,
            sound_distance_bounds=not args.no_sound_bounds,
            big_first=args.big_first,
            lin_level=args.lin_level,
            params=dict(parse_param(s) for s in args.param),
        )
    except TimeoutError as e:
        log.warning("%s", e)
        sys.exit(1)
    if res:
        report(N, res)
