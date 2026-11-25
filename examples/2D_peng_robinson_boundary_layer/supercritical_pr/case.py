#!/usr/bin/env python3
"""Isothermal flat-plate boundary layer; see ../README.md."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from setup_case import case  # noqa: E402

if __name__ == "__main__":
    print(json.dumps(case(Path(__file__).resolve().parent.name)))
