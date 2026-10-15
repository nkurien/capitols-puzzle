# capitols-puzzle

CP-SAT (OR-Tools) solver for the Capitols puzzle. Partitions an N×N grid into
connected, symmetric "states", picks each state's capitol, and satisfies the
given distance clues. Candidate solutions are verified against the full
symmetry group, with cuts added and re-solved until they pass.

## Setup

    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt

## Usage

    python solve.py ex      # 5x5 example
    python solve.py main    # 11x11 puzzle

Options: `--workers N`, `--time SECONDS`, `--heartbeat SECONDS`,
`--logfile PATH`, `--cpsat-log` (raw CP-SAT log).

## Tests

    uv pip install -r requirements-dev.txt
    pytest              # ~15s
    pytest -m slow      # larger 4x4 equivalence cases

`tests/reference.py` is a brute-force oracle written from the problem statement (it shares no code
with `solve.py`). `tests/test_equivalence.py` checks that the CP-SAT model accepts exactly the same
set of solutions as the oracle on 3x3 and 4x4 grids, so a model change can't silently lose one.
