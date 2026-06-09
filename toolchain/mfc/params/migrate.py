"""Rewrite numeric values of enumerated parameters in case-file text to their names."""

import re
from typing import Tuple

from .definitions import CONSTRAINTS


def migrate_text(text: str) -> Tuple[str, int]:
    """Replace integer codes with names for all parameters that define names.

    Only rewrites simple `"param": <int>` (or single-quoted) dict entries.
    Returns (new_text, number_of_replacements).

    Note: operates on raw text, so occurrences inside comments or string
    literals are also rewritten; the caller's reload-and-compare check
    guards semantics but not comment text.
    """
    total = 0
    for param, constraint in sorted(CONSTRAINTS.items()):
        names = constraint.get("names")
        if not names:
            continue
        by_value = {v: n for n, v in names.items()}
        pattern = re.compile(r"""(["']""" + re.escape(param) + r"""["']\s*:\s*)(\d+)(\s*[,}\n])""")

        def repl(m, by_value=by_value):
            nonlocal total
            value = int(m.group(2))
            if value not in by_value:
                return m.group(0)
            total += 1
            return f'{m.group(1)}"{by_value[value]}"{m.group(3)}'

        text = pattern.sub(repl, text)
    return text, total
