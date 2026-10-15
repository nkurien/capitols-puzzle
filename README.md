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
