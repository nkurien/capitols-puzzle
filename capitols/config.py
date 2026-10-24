from dataclasses import dataclass, field


@dataclass(frozen=True)
class SolverConfig:
    """Everything that tunes the model or the search; none of it changes what counts as a solution
    unless noted."""

    maxd: int = 400
    verbose: bool = False
    workers: int = 8
    tlimit: int = 3000
    heartbeat: int = 30
    # restrictions: anything found is still a genuine solution, the search is only narrowed
    max_size: int | None = None  # cap every state's size
    no_free_r2: bool = False  # drop 180-degree rotations that fix no cell
    require_free_r2: bool = False  # some state must use such a rotation
    min_max_size: int | None = None  # some state must have at least this many cells
    # model variants
    clue_bounds: bool = False
    tclauses: bool = True
    edge_images: bool = True
    lab_channel: bool = False
    eq_channel: bool = False
    clue_cuts: bool = True
    sound_distance_bounds: bool = True
    big_first: bool = False
    lin_level: int | None = None
    # raw CP-SAT parameters, set by name after the fields above (e.g. {"symmetry_level": 0})
    params: dict = field(default_factory=dict)
