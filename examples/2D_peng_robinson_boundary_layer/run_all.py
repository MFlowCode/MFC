#!/usr/bin/env python3
"""Run all three CPU/serial cases and their comparison from any working directory."""

import subprocess
import sys
from pathlib import Path

from setup_case import VARIANTS

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]

for variant in VARIANTS:
    subprocess.run([str(REPOSITORY / "mfc.sh"), "run", str(HERE / variant / "case.py"), "-t", "pre_process", "simulation", "--no-mpi", "--no-gpu", "-j", "4"], cwd=REPOSITORY, check=True)
subprocess.run([sys.executable, str(HERE / "verify.py")], cwd=REPOSITORY, check=True)
