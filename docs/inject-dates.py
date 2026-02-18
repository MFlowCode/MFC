#!/usr/bin/env python3
"""Inject last-updated dates into docs before Doxygen runs.

Usage: python3 inject-dates.py [source_dir]
  source_dir defaults to current directory.
"""
import subprocess
import sys
from pathlib import Path
from datetime import date

src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
docs_dir = src_dir / "docs" / "documentation"

for md_file in sorted(docs_dir.glob("*.md")):
    if "Page last updated:" in md_file.read_text():
        continue

    result = subprocess.run(
        ["git", "log", "-1", "--format=%as", "--", str(md_file)],
        capture_output=True, text=True,
        cwd=str(src_dir),
    )
    page_date = result.stdout.strip() or str(date.today())

    with open(md_file, "a") as f:
        f.write(
            f"\n\n<div style='text-align:center; font-size:0.75rem; "
            f"color:#888; padding:16px 0 0;'>"
            f"Page last updated: {page_date}</div>\n"
        )
