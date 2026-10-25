#!/usr/bin/env python3
"""Sparse-steering simulation entry point.

This thin entry point delegates to the canonical repository runner
``scripts/run_sparse_drift.py``. It preserves the original command-line
interface and parameter definitions used for the sparse manuscript data.
"""

from pathlib import Path
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "run_sparse_drift.py"

if __name__ == "__main__":
    if not RUNNER.is_file():
        raise FileNotFoundError(f"Missing canonical sparse runner: {RUNNER}")
    sys.argv[0] = str(RUNNER)
    runpy.run_path(str(RUNNER), run_name="__main__")
