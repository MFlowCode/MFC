#!/usr/bin/env python3

"""Validate case-optimization output: check every value in D/*.dat for NaN/Inf."""

import math
import os
import re
import sys
from pathlib import Path


def extract_doubles(s: str) -> list:
    return [float(e) for e in re.sub(r"[\n\t\s]+", " ", s).strip().split(" ")]


def check_case_dir(case_dir: str) -> bool:
    D_dir = os.path.join(case_dir, "D")
    if not os.path.isdir(D_dir):
        print(f"ERROR: No D/ directory found in {case_dir}")
        return False

    dat_files = list(Path(D_dir).rglob("*.dat"))
    if not dat_files:
        print(f"ERROR: No .dat files found in {D_dir}")
        return False

    ok = True
    total_values = 0
    for filepath in sorted(dat_files):
        with open(filepath, "r") as f:
            content = f.read()

        try:
            doubles = extract_doubles(content)
        except ValueError:
            print(f"ERROR: Failed to parse {filepath} as floating point numbers")
            ok = False
            continue

        total_values += len(doubles)
        for i, val in enumerate(doubles):
            if math.isnan(val) or math.isinf(val):
                print(f"ERROR: {'NaN' if math.isnan(val) else 'Inf'} at index {i} in {filepath}")
                ok = False
                break

    if ok:
        print(f"OK: {len(dat_files)} files, {total_values} values — no NaN/Inf found")
    return ok


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <case_directory>")
        sys.exit(1)

    case_dir = sys.argv[1]
    if os.path.isfile(case_dir):
        case_dir = os.path.dirname(case_dir)

    if not check_case_dir(case_dir):
        sys.exit(1)
