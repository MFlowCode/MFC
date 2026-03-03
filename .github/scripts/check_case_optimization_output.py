#!/usr/bin/env python3

"""Validate case-optimization output: check D/*.dat for NaN/Inf via the packer."""

import math
import sys
import os

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <case_directory>", file=sys.stderr)
    sys.exit(1)

# Allow importing from the repo root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from toolchain.mfc.packer.pack import compile as pack_compile

case_dir = sys.argv[1]
if os.path.isfile(case_dir):
    case_dir = os.path.dirname(case_dir)

pack, err = pack_compile(case_dir)
if err is not None:
    print(f"ERROR: {err}")
    sys.exit(1)

if not pack.entries:
    print(f"ERROR: No data found in {case_dir}/D/")
    sys.exit(1)

if pack.has_bad_values():
    print("ERROR: NaN or Inf detected in output:")
    for name, entry in pack.entries.items():
        for i, val in enumerate(entry.doubles):
            if math.isnan(val) or math.isinf(val):
                label = 'NaN' if math.isnan(val) else 'Inf'
                print(f"  {label} at index {i} in {name}")
                break
    sys.exit(1)

total = sum(len(e.doubles) for e in pack.entries.values())
print(f"OK: {len(pack.entries)} files, {total} values — no NaN/Inf found")
