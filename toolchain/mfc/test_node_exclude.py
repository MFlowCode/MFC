"""Unit tests for .github/scripts/node-exclude.sh.

When the in-allocation preflight finds the node's GPU unusable it records the
node name in the job's output file. The submit wrapper reads it back, adds it to
the sbatch --exclude list and resubmits, so the next attempt lands somewhere
else instead of failing the same way.

This matters because bad nodes are concentrated, not spread: over 2026-08-18..31
a single Phoenix node (atl1-1-03-007-29-0) accounted for 25 of 29 ECC failures
and two nodes for 32 of ~40. Both are currently in a hand-edited --exclude line
that a human had to notice, diagnose, and commit.
"""

import os
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "node-exclude.sh"


def run(*args):
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        capture_output=True,
        text=True,
        env={**os.environ},
        check=False,
    )


def out(*args):
    result = run(*args)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


@pytest.fixture
def job_output(tmp_path):
    def _write(text):
        path = tmp_path / "build-and-test-gpu-acc.out"
        path.write_text(text)
        return str(path)

    return _write


def test_node_from_reads_the_marker_the_preflight_wrote(job_output):
    path = job_output("some build noise\nMFC_FAULT_NODE=atl1-1-03-007-29-0\nmore noise\n")
    assert out("node-from", path) == "atl1-1-03-007-29-0"


def test_node_from_prints_nothing_when_the_job_recorded_no_fault(job_output):
    assert out("node-from", job_output("ordinary output, no marker\n")) == ""


def test_node_from_prints_nothing_when_the_output_file_is_missing(tmp_path):
    assert out("node-from", str(tmp_path / "absent.out")) == ""


def test_node_from_takes_the_last_marker_when_the_file_holds_several(job_output):
    # submit-slurm-job.sh reuses one output path across resubmits, so a stale
    # marker can precede the current one.
    path = job_output("MFC_FAULT_NODE=oldnode\nMFC_FAULT_NODE=atl1-1-03-007-31-0\n")
    assert out("node-from", path) == "atl1-1-03-007-31-0"


def test_merge_adds_the_first_node_to_an_empty_list():
    assert out("merge", "", "atl1-1-03-007-29-0") == "atl1-1-03-007-29-0"


def test_merge_appends_to_an_existing_list():
    assert out("merge", "nodeA", "nodeB") == "nodeA,nodeB"


def test_merge_does_not_repeat_a_node_already_excluded():
    assert out("merge", "nodeA,nodeB", "nodeA") == "nodeA,nodeB"


def test_merge_does_not_match_a_node_name_that_is_only_a_prefix():
    # "atl1-1-03-007-2" must not suppress excluding "atl1-1-03-007-29-0".
    assert out("merge", "atl1-1-03-007-2", "atl1-1-03-007-29-0") == "atl1-1-03-007-2,atl1-1-03-007-29-0"


def test_merge_leaves_the_list_untouched_when_there_is_no_node_to_add():
    assert out("merge", "nodeA", "") == "nodeA"


def test_an_unknown_subcommand_fails_with_a_usage_error():
    result = run("frobnicate", "x")
    assert result.returncode == 2
    assert "usage" in (result.stdout + result.stderr).lower()
