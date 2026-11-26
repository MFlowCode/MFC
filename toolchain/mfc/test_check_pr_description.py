"""Tests for .github/scripts/check_pr_description.py.

The PR template asks the author to show that a change solves its problem, creates
no new one and does not slow the code. The check fails a description that skips
any of that, or whose benchmarks predate a later source change. It reads the
description only; it runs none of the pull request's code.

Descriptions here are built by filling in the real template, so a template edit
that the check cannot read fails these tests.
"""

import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

from mfc import state
from mfc.bench import diff

REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("check_pr_description", REPO / ".github" / "scripts" / "check_pr_description.py")
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)

TEMPLATE = (REPO / ".github" / "pull_request_template.md").read_text()

SRC = ["src/simulation/m_rhs.fpp"]
TOOLCHAIN = ["toolchain/mfc/bench.py"]
DOCS = ["docs/documentation/testing.md", "README.md"]

SPEED_PLACEHOLDER = "<!-- Paste build/speed.md here -->"
NUMBERS_EMPTY = "| Quantity | master | this PR | expected |\n|---|---|---|---|\n|  |  |  |  |\n"
NUMBERS = NUMBERS_EMPTY.replace("|  |  |  |  |", "| force difference | 6.0e-5 | 0 | 0 |")
FIGURES_EMPTY = "| Figure: master | Figure: this PR |\n|---|---|\n|  |  |\n"
FIGURES = FIGURES_EMPTY.replace("|  |  |", "| ![before](https://github.com/user-attachments/assets/a) | ![after](https://github.com/user-attachments/assets/b) |")
WHY_LABEL = "Why the new result is right (analytic solution, reference data, convergence, ...):"
WHY = WHY_LABEL + " the force on a sphere at rest must vanish.\n"
COMMAND = "```bash\n# two ranks | one rank\n./mfc.sh run case.py -n 2\n```\n"

SPEED = """### Benchmark

Grind time speedups from `../mfc-base/build/bench.yaml` to `build/bench.yaml`; greater than 1 is faster.

| Case | Grind speedup |
| --- | --- |
| `ibm` | 0.99 |
| `igr` | 1.01 |
"""


def fill(text, old, new):
    assert text.count(old) == 1, f"template no longer has exactly one {old!r}"
    return text.replace(old, new)


def complete():
    body = fill(TEMPLATE, "Command that shows it:\n", "Command that shows it:\n\n" + COMMAND)
    body = fill(body, NUMBERS_EMPTY, NUMBERS)
    body = fill(body, WHY_LABEL + "\n", WHY)
    body = fill(body, "| Tests: command, number passed, number failed |  |", "| Tests: command, number passed, number failed | `./mfc.sh test --select-enforce`: 137 passed, 0 failed |")
    body = fill(body, "| Golden files changed: count, largest difference, why |  |", "| Golden files changed: count, largest difference, why | 0 |")
    body = fill(body, "| Not tested |  |", "| Not tested | Cray GPUs |")
    body = fill(body, "GPU and compiler:", "GPU and compiler: MI300A, amdflang 24.4")
    body = fill(body, "Base commit:\n", "Base commit: ea89e86d5\n")
    body = fill(body, "PR commit:\n", "PR commit: `1234abc`\n")
    return fill(body, SPEED_PLACEHOLDER, SPEED)


def without_speed_results():
    body = fill(complete(), "GPU and compiler: MI300A, amdflang 24.4", "GPU and compiler:")
    body = fill(body, "Base commit: ea89e86d5\n", "Base commit:\n")
    body = fill(body, "PR commit: `1234abc`\n", "PR commit:\n")
    return fill(body, SPEED, "")


def not_applicable(*sections):
    body = complete()
    for number in sections:
        heading = next(line for line in body.splitlines() if line.startswith(f"## {number}."))
        start = body.index(heading) + len(heading)
        end = body.find("\n## ", start)
        body = body[:start] + "\n\nN/A: nothing to measure here\n" + (body[end:] if end >= 0 else "")
    return body


def check(body, files=SRC, since=()):
    return checker.check(body, files, lambda sha: None if since is None else list(since))


def assert_one_problem(problems, *words):
    assert len(problems) == 1, problems
    for word in words:
        assert word in problems[0], problems[0]


@pytest.mark.parametrize("files", [SRC, TOOLCHAIN, DOCS])
def test_a_complete_description_passes(files):
    assert check(complete(), files) == []


def test_line_endings_from_the_github_editor_are_accepted():
    assert check(complete().replace("\n", "\r\n")) == []


@pytest.mark.parametrize("files", [SRC, TOOLCHAIN, DOCS])
def test_the_untouched_template_fails_every_section(files):
    problems = check(TEMPLATE, files)

    for number in (1, 2, 3):
        assert any(p.startswith(f"{number}.") for p in problems), problems


def test_a_description_that_ignores_the_template_fails():
    problems = check("Fixes a bug. Tested locally.")

    assert [p.split(":")[0] for p in problems] == ["1. Problem", "2. Tests", "3. Speed"]
    assert all("missing" in p for p in problems)


@pytest.mark.parametrize(
    "old,new,words",
    [
        (NUMBERS, NUMBERS_EMPTY, ("1.", "numbers", "figures")),
        ("| force difference | 6.0e-5 | 0 | 0 |", "| force difference | 6.0e-5 |  | 0 |", ("1.", "blank")),
        (COMMAND, "", ("1.", "command")),
        (WHY, WHY_LABEL + "\n", ("1.", "why")),
        (WHY, "", ("1.", "why")),
        ("| Not tested | Cray GPUs |", "| Not tested |  |", ("2.", "Not tested", "blank")),
        ("| Not tested | Cray GPUs |", "| Not tested | - |", ("2.", "Not tested", "blank")),
        ("| Not tested | Cray GPUs |", "| Not tested | N/A |", ("2.", "Not tested", "reason")),
        ("| Not tested | Cray GPUs |\n", "", ("2.", "Not tested", "missing")),
        ("`./mfc.sh test --select-enforce`: 137 passed, 0 failed", "N/A: trivial change", ("2.", "Tests", "requires")),
        ("GPU and compiler: MI300A, amdflang 24.4", "GPU and compiler:", ("3.", "GPU and compiler")),
        ("Base commit: ea89e86d5", "Base commit:", ("3.", "Base commit")),
        ("PR commit: `1234abc`", "PR commit: latest", ("3.", "PR commit")),
        (SPEED, "", ("3.", "speedup table")),
        ("| `ibm` | 0.99 |", "| `ibm` | N/A |", ("3.", "speedup", "ibm")),
    ],
)
def test_each_missing_piece_is_named(old, new, words):
    assert_one_problem(check(fill(complete(), old, new)), *words)


def test_a_command_on_the_same_line_as_its_label_counts():
    body = fill(complete(), COMMAND, "")
    body = fill(body, "Command that shows it:", "Command that shows it: `./mfc.sh run case.py -n 2`")

    assert check(body) == []


@pytest.mark.parametrize("number", [1, 2, 3])
def test_a_required_section_cannot_be_not_applicable(number):
    assert_one_problem(check(not_applicable(number)), f"{number}.", "requires")


def test_an_na_remark_inside_a_complete_required_section_is_allowed():
    body = fill(complete(), SPEED, SPEED + "\nN/A: the benchmark cases cover this change, so no extra timing.\n")

    assert check(body) == []


def test_other_changes_may_skip_tests_and_speed_with_a_reason():
    assert check(not_applicable(2, 3), TOOLCHAIN) == []
    assert_one_problem(check(not_applicable(1), TOOLCHAIN), "1.", "requires")


def test_docs_may_skip_all_three_with_reasons():
    assert check(not_applicable(1, 2, 3), DOCS) == []


@pytest.mark.parametrize("bare", ["N/A", "n/a.", "**N/A**"])
def test_na_without_a_reason_fails(bare):
    body = fill(not_applicable(3), "N/A: nothing to measure here", bare)

    assert_one_problem(check(body, DOCS), "3.", "reason")


@pytest.mark.parametrize("path", ["src/common/m_helper.fpp", "CMakeLists.txt", "cmake/GPU.cmake", "toolchain/cmake/x.cmake", "toolchain/modules"])
def test_source_and_build_changes_need_speed_results(path):
    assert_one_problem(check(not_applicable(3), ["README.md", path]), "3.", "requires")


def test_results_older_than_a_source_change_fail():
    problems = check(complete(), since=["README.md", "src/simulation/m_weno.fpp"])

    assert_one_problem(problems, "3.", "1234abc", "src/simulation/m_weno.fpp")


def test_later_changes_outside_source_and_build_keep_results_current():
    assert check(complete(), since=["README.md", "toolchain/mfc/bench.py"]) == []


def test_a_commit_outside_the_pull_request_fails():
    assert_one_problem(check(complete(), since=None), "3.", "1234abc", "history")


def _env_without_git():
    # A commit hook's GIT_DIR would send these commands to the repository being committed.
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def _git(repo, *args):
    command = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.com", *args]
    return subprocess.run(command, check=True, capture_output=True, text=True, env=_env_without_git()).stdout.strip()


def _commit(repo, name):
    (repo / name).write_text(name)
    _git(repo, "add", name)
    _git(repo, "commit", "-qm", name)
    return _git(repo, "rev-parse", "HEAD")


def test_git_fixtures_leave_a_commit_hooks_repository_alone(tmp_path, monkeypatch):
    # The pre-commit hook runs these tests with GIT_DIR and GIT_INDEX_FILE naming the repository
    # being committed, and those override `git -C`. A decoy stands in for that repository.
    decoy, work = tmp_path / "decoy", tmp_path / "work"
    subprocess.run(["git", "init", "-q", str(decoy)], check=True, env={k: v for k, v in os.environ.items() if not k.startswith("GIT_")})
    config = (decoy / ".git" / "config").read_text()
    monkeypatch.setenv("GIT_DIR", str(decoy / ".git"))
    monkeypatch.setenv("GIT_INDEX_FILE", str(decoy / ".git" / "index"))
    work.mkdir()

    _git(work, "init", "-q", "-b", "main")
    _commit(work, "a")

    assert (decoy / ".git" / "config").read_text() == config
    assert list((decoy / ".git" / "refs" / "heads").iterdir()) == []
    assert (work / ".git").is_dir()


def test_changed_since_lists_later_files_and_rejects_foreign_commits(tmp_path, monkeypatch):
    for name in set(os.environ) - set(_env_without_git()):
        monkeypatch.delenv(name)
    _git(tmp_path, "init", "-q", "-b", "main")
    first = _commit(tmp_path, "a")
    _git(tmp_path, "checkout", "-qb", "side")
    foreign = _commit(tmp_path, "side-file")
    _git(tmp_path, "checkout", "-q", "main")
    head = _commit(tmp_path, "b")

    assert checker.changed_since(first[:9], head, tmp_path) == ["b"]
    assert checker.changed_since(head, head, tmp_path) == []
    assert checker.changed_since(foreign, head, tmp_path) is None
    assert checker.changed_since("0123456789", head, tmp_path) is None


def test_main_reports_problems_and_exit_status(tmp_path, capsys):
    body, files = tmp_path / "body.md", tmp_path / "files.txt"
    files.write_text("docs/a.md\n")

    body.write_text(not_applicable(1, 2, 3))
    assert checker.main([str(body), str(files), "HEAD"]) == 0

    body.write_text("")
    assert checker.main([str(body), str(files), "HEAD"]) == 1
    assert "1. Problem" in capsys.readouterr().out


@pytest.mark.parametrize(
    "numbers,figures",
    [
        (NUMBERS, FIGURES_EMPTY),
        (NUMBERS, ""),
        (NUMBERS_EMPTY, FIGURES),
        ("", FIGURES),
        (NUMBERS, FIGURES),
    ],
)
def test_numbers_or_figures_or_both_are_enough(numbers, figures):
    body = fill(fill(complete(), NUMBERS, numbers), FIGURES_EMPTY, figures)

    assert check(body) == []


@pytest.mark.parametrize(
    "cell",
    [
        "![after](https://github.com/user-attachments/assets/b)",
        '<img width="400" src="https://github.com/user-attachments/assets/b">',
        '<video src="https://github.com/user-attachments/assets/b"></video>',
        "https://github.com/user-attachments/assets/b",
    ],
)
def test_each_way_github_embeds_media_counts_as_a_figure(cell):
    figures = FIGURES.replace("![after](https://github.com/user-attachments/assets/b)", cell)

    assert check(fill(fill(complete(), NUMBERS, ""), FIGURES_EMPTY, figures)) == []


@pytest.mark.parametrize("cell", ["", "looks better now"])
def test_figures_need_an_image_on_both_sides(cell):
    figures = FIGURES.replace("![after](https://github.com/user-attachments/assets/b)", cell)

    assert_one_problem(check(fill(complete(), FIGURES_EMPTY, figures)), "1.", "figure")


def test_a_partly_filled_numbers_table_fails_even_with_figures():
    body = fill(fill(complete(), FIGURES_EMPTY, FIGURES), "| force difference | 6.0e-5 | 0 | 0 |", "| force difference | 6.0e-5 |  | 0 |")

    assert_one_problem(check(body), "1.", "blank")


def test_neither_numbers_nor_figures_fails():
    assert_one_problem(check(fill(complete(), NUMBERS, "")), "1.", "numbers", "figures")


def test_the_reason_may_start_on_the_next_line():
    assert check(fill(complete(), WHY, WHY_LABEL + "\n\nIt matches the exact solution.\n")) == []


def test_the_reason_does_not_stand_in_for_the_command():
    body = fill(fill(fill(complete(), COMMAND, ""), NUMBERS, ""), FIGURES_EMPTY, "")

    assert any("no command" in p for p in check(body))


@pytest.mark.parametrize(
    "statement",
    [
        "No GPU access",
        "no access to a GPU.",
        "I have no GPU access, sorry.",
        "N/A: no GPU access",
        "**No GPU access**",
        "GPU and compiler: no GPU access",
    ],
)
def test_saying_there_is_no_gpu_replaces_the_speed_results(statement):
    body = without_speed_results()
    if statement.startswith("GPU and compiler:"):
        body = fill(body, "GPU and compiler:", statement)
    else:
        body = fill(body, "GPU and compiler:", statement + "\n\nGPU and compiler:")

    assert check(body) == []


def test_no_gpu_does_not_excuse_the_other_sections():
    body = fill(without_speed_results(), "GPU and compiler:", "No GPU access\n\nGPU and compiler:")
    body = fill(body, "| Not tested | Cray GPUs |", "| Not tested |  |")

    assert_one_problem(check(body), "2.", "Not tested")


def test_missing_speed_results_mention_the_no_gpu_statement():
    problems = check(without_speed_results())

    assert any("No GPU access" in p for p in problems), problems


def test_an_old_bench_diff_table_says_to_merge_master():
    old = "| Case | Pre Process | Simulation | Post Process |\n| --- | --- | --- | --- |\n| `ibm` | Exec: 0.98 | Exec: 1.01 & Grind: 1.00 | Exec: 1.00 |\n"

    assert_one_problem(check(fill(complete(), SPEED, old)), "3.", "merge master")


def test_the_table_bench_diff_writes_is_accepted(tmp_path, monkeypatch):
    results = tmp_path / "results.yaml"
    results.write_text("metadata: {invocation: [bench], lock: {}}\ncases: {ibm: {output_summary: {simulation: {grind: 2.0}}}}\n")
    summary = tmp_path / "speed.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    monkeypatch.setitem(state.gARG, "lhs", str(results))
    monkeypatch.setitem(state.gARG, "rhs", str(results))

    diff()

    assert check(fill(complete(), SPEED, summary.read_text())) == []


def test_description_text_cannot_start_a_workflow_command():
    body = fill(complete(), "| `ibm` | 0.99 |", "| `ibm`\r::error::forged | slow\r::add-mask::x |")

    problems = check(body)

    assert problems
    assert not any(ch in p for p in problems for ch in "\r\n\x1b"), problems
