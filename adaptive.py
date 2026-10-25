#!/usr/bin/env python3
"""Adaptive-relocation simulation entry point.

Use ``--experiment`` to select the canonical runner for Fig. 5b, 5c, or 5d.
All remaining arguments are passed unchanged to that runner.
"""

from pathlib import Path
import argparse
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
RUNNERS = {
    "fig5b": ROOT / "scripts" / "run_fig5b_exposure.py",
    "fig5c": ROOT / "scripts" / "run_fig5c_threshold.py",
    "fig5d": ROOT / "scripts" / "run_fig5d_phase_space.py",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=sorted(RUNNERS), required=True)
    args, remainder = parser.parse_known_args()
    runner = RUNNERS[args.experiment]
    if not runner.is_file():
        raise FileNotFoundError(f"Missing canonical adaptive runner: {runner}")
    sys.argv = [str(runner), *remainder]
    runpy.run_path(str(runner), run_name="__main__")


if __name__ == "__main__":
    main()
