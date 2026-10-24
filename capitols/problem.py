import logging
import time
from dataclasses import dataclass

from .config import SolverConfig
from .grid import Grid
from .propagate import Bounds, RegionFinder, RootRegions, root_regions, tighten_bounds
from .symmetry import ModelTransforms

log = logging.getLogger("capitols")


@dataclass
class Problem:
    """Everything derived from the clues before any CP-SAT variable exists."""

    grid: Grid
    clues: dict
    cfg: SolverConfig
    mt: ModelTransforms
    bounds: Bounds
    finder: RegionFinder
    regions: RootRegions
    t_start: float


def prepare(N, clues, cfg):
    t_start = time.time()
    grid = Grid(N)
    mt = ModelTransforms(grid, no_free_r2=cfg.no_free_r2)
    log.info(
        "[build] %d symmetry transforms in the model (%d with R4, %d before dropping redundant ones)",
        len(mt),
        len(mt.full),
        mt.n_before,
    )
    bounds = tighten_bounds(grid, clues, cfg)
    finder = RegionFinder(grid, mt, bounds)
    regions = root_regions(grid, mt, bounds, finder)
    return Problem(grid, clues, cfg, mt, bounds, finder, regions, t_start)
