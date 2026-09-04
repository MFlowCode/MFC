#!/usr/bin/env python3
"""Print a bounded summary of a GPU memory fault in a run log.

For callers that are shell scripts. The ROCm debug agent emits tens of
thousands of lines per fault -- one disassembly and register dump repeated per
faulting wave -- and the part worth reading (the faulting kernel, the fault
reason, the stop-PC distribution) is buried in the middle, so `tail` cannot
find it.

Exits 0 having printed a summary, or 1 having printed nothing when the log has
no agent report, which lets the caller fall back to whatever it did before.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "toolchain"))

from mfc.gpu_diagnostics import summarize_rocm_debug_agent  # noqa: E402


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <run log>", file=sys.stderr)
        return 2

    try:
        with open(sys.argv[1], "r", encoding="utf-8", errors="replace") as log:
            summary = summarize_rocm_debug_agent(log.read())
    except OSError as exc:
        print(f"could not read {sys.argv[1]}: {exc}", file=sys.stderr)
        return 1

    if not summary:
        return 1

    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
