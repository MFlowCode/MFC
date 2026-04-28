#!/usr/bin/env python3
"""Generate SILO series file for time-series visualization in ParaView/VisIt."""

import json
import sys
from pathlib import Path


def generate_silo_series(case_dir):
    """Generate collection.silo.series file from existing collection files."""
    root_dir = Path(case_dir) / "silo_hdf5" / "root"
    series_file = Path(case_dir) / "silo_hdf5" / "collection.silo.series"

    if not root_dir.exists():
        return

    collection_files = sorted(
        root_dir.glob("collection_*.silo"),
        key=lambda p: int(p.stem.split("_")[1]),
    )

    if not collection_files:
        return

    files = [{"name": f"root/{f.name}", "time": i} for i, f in enumerate(collection_files)]

    series_data = {"file-series-version": "1.0", "files": files}

    with open(series_file, "w") as f:
        json.dump(series_data, f, indent=2)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: generate_silo_series.py <case_dir>")
        sys.exit(1)

    case_dir = sys.argv[1]
    generate_silo_series(case_dir)
