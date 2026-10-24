"""Command-line entry point and import shim; the solver lives in the `capitols` package."""

from capitols.cli import main
from capitols.runner import report, run, solve
from capitols.symmetry import fixed, full_group_ok, transforms, undominated

__all__ = [
    "fixed",
    "full_group_ok",
    "main",
    "report",
    "run",
    "solve",
    "transforms",
    "undominated",
]

if __name__ == "__main__":
    main()
