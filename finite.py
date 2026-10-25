#!/usr/bin/env python3
"""Finite-steering simulation entry point.

This thin entry point delegates to the canonical repository runner
``scripts/run_finite_steering.py``. The runner supports drift, boundary,
global, and exposure modes through its original command-line interface.
"""

from pathlib import Path
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "run_finite_steering.py"

if __name__ == "__main__":
    if not RUNNER.is_file():
        raise FileNotFoundError(f"Missing canonical finite runner: {RUNNER}")
    sys.argv[0] = str(RUNNER)
    runpy.run_path(str(RUNNER), run_name="__main__")
