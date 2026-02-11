"""Check that file paths referenced in documentation still exist."""

import re
import sys
from pathlib import Path

# Docs to scan for file path references
DOCS = [
    "docs/documentation/contributing.md",
    "docs/documentation/gpuParallelization.md",
    "docs/documentation/running.md",
    "docs/documentation/case.md",
    ".github/copilot-instructions.md",
]

# Match backtick-wrapped strings that look like repo-relative file paths
PATH_RE = re.compile(r"`((?:src|toolchain|\.github|docs|examples|tests)/[^`]+)`")

# Skip paths with placeholders, globs, or patterns (not real file paths)
SKIP_RE = re.compile(r"[<>*()\[\]{}]|/\.\.\.|%|\$")


def check_docs(repo_root: Path) -> list[str]:
    errors = []
    for doc in DOCS:
        doc_path = repo_root / doc
        if not doc_path.exists():
            continue
        text = doc_path.read_text(encoding="utf-8")
        for match in PATH_RE.finditer(text):
            path_str = match.group(1)
            if SKIP_RE.search(path_str):
                continue
            # Strip trailing punctuation that may have leaked in
            path_str = path_str.rstrip(".,;:!?")
            if not (repo_root / path_str).exists():
                errors.append(f"  {doc} references '{path_str}' but it does not exist")
    return errors


def main():
    repo_root = Path(__file__).resolve().parents[2]
    errors = check_docs(repo_root)
    if errors:
        print("Doc reference check failed:")
        for e in errors:
            print(e)
        sys.exit(1)


if __name__ == "__main__":
    main()
