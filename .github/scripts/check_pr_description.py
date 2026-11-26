#!/usr/bin/env python3
"""Fail a pull request whose description does not complete the template.

    check_pr_description.py <body file> <changed files list> <head commit>

Run from a clone that holds the pull request's commits. What a change must fill in follows the
table in .github/pull_request_template.md. The description and file names come from the pull
request: they are read as text, and nothing from the pull request is run. This checks that a
report is complete and newer than the last source change, not that its numbers or figures are
true. An author who says they have no GPU access owes no speed results.
"""

import os
import re
import subprocess
import sys

SECTIONS = {1: "Problem", 2: "Tests", 3: "Speed"}
SECTION_2_ROWS = ("Tests", "Golden files", "Not tested")
SPEED_PATHS = ("src/", "cmake/", "toolchain/cmake/")
SPEED_FILES = ("CMakeLists.txt", "toolchain/modules")
SHA = re.compile(r"\b[0-9a-f]{7,40}\b", re.IGNORECASE)
NOT_APPLICABLE = re.compile(r"^[\s*_>-]*n/a\b(.*)$", re.IGNORECASE)
NO_GPU = re.compile(r"^[\s*_>-]*(?:gpu and compiler:)?[\s*_]*(?:n/a\W*)?(?:i have )?no (?:gpu access|access to (?:a )?gpu)\b", re.IGNORECASE)
# A Markdown or HTML image, a video, or a bare link (how GitHub embeds uploaded video).
MEDIA = re.compile(r"!\[[^\]]*\]\([^)]+\)|<(?:img|video)\b|https?://\S+", re.IGNORECASE)


def needs_speed(path):
    return path.startswith(SPEED_PATHS) or path in SPEED_FILES


def is_doc(path):
    return path.startswith("docs/") or path.endswith(".md") or path in ("LICENSE", "CITATION.cff")


def split_sections(body):
    """Section number -> its lines as (text, inside a code fence)."""
    body = re.sub(r"<!--.*?-->", "", body.replace("\r\n", "\n"), flags=re.DOTALL)
    body = re.sub(r"<details>\s*<summary>Commands.*?</details>", "", body, flags=re.DOTALL)
    sections, current, fence = {}, None, None
    for line in body.split("\n"):
        marker = re.match(r"\s*(`{3,}|~{3,})", line)
        if marker and fence is None:
            fence = marker.group(1)
            continue
        if marker and marker.group(1).startswith(fence):
            fence = None
            continue
        if fence is None and line.startswith("## "):
            heading = re.match(r"##\s+([123])\.", line)
            current = sections.setdefault(int(heading.group(1)), []) if heading else None
        elif current is not None:
            current.append((line, fence is not None))
    return sections


def tables(lines):
    """Each table outside code fences, as rows of cells with the header first."""
    found, rows = [], []
    for line, fenced in [*lines, ("", False)]:
        if not fenced and line.strip().startswith("|"):
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
            if not all(re.fullmatch(r":?-+:?", c) for c in cells):
                rows.append(cells)
        elif rows:
            found.append(rows)
            rows = []
    return found


def na_reason(text):
    """None if the text is not N/A, otherwise whatever follows the N/A."""
    match = NOT_APPLICABLE.match(text)
    return None if match is None else match.group(1).strip(" *_:;,.()[]-–—")


def has_reason(reason):
    return re.search(r"[A-Za-z]{3,}", reason) is not None


def cell_problem(value):
    reason = na_reason(value)
    if reason is not None:
        return None if has_reason(reason) else "says N/A without a reason"
    return None if re.search(r"[A-Za-z0-9]", value) else "is blank"


def field(lines, label):
    for line, fenced in lines:
        if not fenced and line.strip().startswith(label + ":"):
            return line.split(":", 1)[1].strip()
    return ""


def line_starting(lines, label):
    """Index of the first line outside code fences that starts with the label, or None."""
    return next((i for i, (line, fenced) in enumerate(lines) if not fenced and line.strip(" *_").lower().startswith(label)), None)


def filled_rows(lines, first_header):
    """The table whose first header starts as given, without its blank rows; each row padded to the header."""
    table = next((t for t in tables(lines) if t[0][0].startswith(first_header)), None)
    if table is None:
        return []
    return [row + [""] * (len(table[0]) - len(row)) for row in table[1:] if any(row)]


def check_problem(lines, required):
    problems = []
    command_at = line_starting(lines, "command that shows it")
    why_at = line_starting(lines, "why the new result is right")
    if command_at is None:
        problems.append("the 'Command that shows it:' line is missing")
    else:
        command = [lines[command_at][0].split(":", 1)[-1]]
        for i, (line, fenced) in enumerate(lines[command_at + 1 :], command_at + 1):
            if i == why_at or (not fenced and line.strip().startswith("|")):
                break
            command.append(line)
        if not re.search(r"[A-Za-z0-9]", " ".join(command)):
            problems.append("no command is given")

    numbers, figures = filled_rows(lines, "Quantity"), filled_rows(lines, "Figure")
    if not numbers and not figures:
        problems.append("fill in the numbers table or the figures table, master next to this PR")
    bad = [p for row in numbers for p in map(cell_problem, row) if p]
    if bad:
        problems.append(f"a table cell {bad[0]}")
    if not all(MEDIA.search(cell) for row in figures for cell in row):
        problems.append("each figure cell needs an image or video, for master and for this PR")

    why = [] if why_at is None else [lines[why_at][0].partition(":")[2], *[line for line, _ in lines[why_at + 1 :]]]
    if not re.search(r"[A-Za-z0-9]", " ".join(why)):
        problems.append("say why the new result is right, after the 'Why the new result is right' line")
    return problems


def check_tests(lines, required):
    problems = []
    rows = [row for table in tables(lines) for row in table]
    for name in SECTION_2_ROWS:
        row = next((r for r in rows if r[0].startswith(name)), None)
        if row is None:
            problems.append(f"the '{name}' row is missing")
            continue
        value = row[1] if len(row) > 1 else ""
        if name == "Tests" and required and na_reason(value) is not None:
            problems.append("this change requires 'Tests', so it cannot be N/A")
        elif cell_problem(value):
            problems.append(f"'{name}' {cell_problem(value)}")
    return problems


def check_speed(lines, required, changed_since):
    if any(NO_GPU.match(line) for line, fenced in lines if not fenced):
        return []
    problems = []
    if not re.search(r"[A-Za-z0-9]", field(lines, "GPU and compiler")):
        problems.append("'GPU and compiler' is blank (with no GPU, write 'No GPU access')")
    commits = {}
    for label in ("Base commit", "PR commit"):
        match = SHA.search(field(lines, label))
        if match:
            commits[label] = match.group(0)
        else:
            problems.append(f"'{label}' needs a commit id")

    table = next((t for t in tables(lines) if len(t[0]) > 1 and t[0][1] == "Grind speedup"), None)
    if any("Pre Process" in t[0] for t in tables(lines)):
        problems.append("the speedup table is from an older bench_diff; merge master and run the benchmarks again")
    elif table is None or len(table) < 2:
        problems.append("the grind speedup table from build/speed.md is missing")
    else:
        for row in table[1:]:
            value = row[1].strip("*` ") if len(row) > 1 else ""
            if not re.fullmatch(r"\d+(\.\d+)?", value):
                problems.append(f"grind speedup for {row[0].strip('`')} is '{value}', not a number")

    sha = commits.get("PR commit")
    if required and sha:
        changed = changed_since(sha)
        if changed is None:
            problems.append(f"PR commit {sha} is not in this pull request's history (rebased or mistyped)")
        else:
            stale = [path for path in changed if needs_speed(path)]
            if stale:
                more = f" and {len(stale) - 1} more" if len(stale) > 1 else ""
                problems.append(f"results are for {sha}, but {stale[0]}{more} changed since; run the benchmarks again")
    return problems


def check(body, files, changed_since):
    """Problems with a description, given the changed files and a lookup of files changed since a commit."""
    speed = any(needs_speed(path) for path in files)
    required = {1: not all(is_doc(path) for path in files), 2: speed, 3: speed}
    checks = {1: check_problem, 2: check_tests, 3: lambda lines, needed: check_speed(lines, needed, changed_since)}
    sections = split_sections(body)
    report = []
    for number, name in SECTIONS.items():
        lines = sections.get(number)
        if lines is None:
            problems = ["the section is missing; use the pull request template"]
        else:
            reasons = [r for text, fenced in lines if not fenced and (r := na_reason(text)) is not None]
            if reasons and not required[number]:
                problems = [] if any(map(has_reason, reasons)) else ["N/A needs a reason"]
            else:
                problems = checks[number](lines, required[number])
                if problems and reasons:
                    problems = ["this change requires this section, so it cannot be N/A"]
                elif problems and not required[number]:
                    problems = [f"{p} (or write N/A with a reason for the section)" for p in problems]
        # Problems quote the description, and a control character could start a workflow command in the log.
        report += [re.sub(r"[\x00-\x1f\x7f]", " ", f"{number}. {name}: {p}") for p in problems]
    return report


def changed_since(sha, head, repo="."):
    """Files changed from sha to head, or None if sha is not an ancestor of head."""
    git = ["git", "-C", str(repo)]
    if subprocess.run([*git, "merge-base", "--is-ancestor", sha, head], capture_output=True, check=False).returncode != 0:
        return None
    diff = subprocess.run([*git, "diff", "--no-renames", "--name-only", sha, head], capture_output=True, text=True, check=True)
    return diff.stdout.splitlines()


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    body_path, files_path, head = argv
    with open(body_path, encoding="utf-8") as f:
        body = f.read()
    with open(files_path, encoding="utf-8") as f:
        files = f.read().splitlines()

    problems = check(body, files, lambda sha: changed_since(sha, head))
    if problems:
        lines = ["PR description incomplete:", "", *[f"- {p}" for p in problems], "", "Edit the description to match .github/pull_request_template.md; this check reruns on each edit."]
    else:
        lines = ["PR description is complete."]
    text = "\n".join(lines)
    print(text)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as f:
            f.write(text + "\n")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
