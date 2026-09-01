"""Unit tests for .github/scripts/preflight.sh.

The preflight runs syscheck at the top of the allocation that will execute the
tests, before the expensive work starts.

Placement is the whole point. Over 2026-08-18..31, 48 of 58 measurable jobs
built on one node and tested on another, and in the ECC failures the build node
was healthy every time while the tests landed on a bad one. A probe that runs
after the build validates the wrong machine.

Cost of getting it wrong: the ECC jobs spent a median of 38 minutes between
syscheck being available and the GPU fault being noticed -- 28.5 hours in two
weeks.
"""

import os
import stat
import subprocess
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"
SCRIPT = SCRIPTS / "preflight.sh"

HEALTHY = 0
NODE_FAULT = 77
OUTAGE = 78


@pytest.fixture
def workspace(tmp_path):
    """A workspace with a stubbed mpirun on PATH and a place for syscheck."""
    binz = tmp_path / "bin"
    binz.mkdir()
    mpirun = binz / "mpirun"
    # Passthrough: drop "-np N" and exec the binary, mirroring how build.sh runs
    # syscheck under mpirun.
    mpirun.write_text('#!/bin/bash\nwhile [ "${1:-}" = "-np" ]; do shift 2; done\nexec "$@"\n')
    mpirun.chmod(mpirun.stat().st_mode | stat.S_IEXEC)
    (tmp_path / "state").mkdir()
    return tmp_path


def write_syscheck(workspace, exit_code, message="syscheck says hello"):
    target = workspace / "build" / "install" / "gpu-acc-abc123" / "bin"
    target.mkdir(parents=True, exist_ok=True)
    binary = target / "syscheck"
    binary.write_text(f'#!/bin/bash\necho "{message}"\nexit {exit_code}\n')
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    return binary


def run(workspace, *args):
    env = {
        **os.environ,
        "PATH": f"{workspace / 'bin'}:{os.environ['PATH']}",
        "MFC_CI_STATE_DIR": str(workspace / "state"),
        "SLURMD_NODENAME": "atl1-1-03-007-29-0",
    }
    return subprocess.run(
        ["bash", str(SCRIPT), *(args or ("phoenix", "gpu"))],
        capture_output=True,
        text=True,
        cwd=workspace,
        env=env,
        check=False,
    )


def test_passes_when_syscheck_succeeds_on_this_node(workspace):
    write_syscheck(workspace, 0)
    assert run(workspace).returncode == HEALTHY


def test_reports_a_node_fault_when_syscheck_fails(workspace):
    write_syscheck(workspace, 1)
    assert run(workspace).returncode == NODE_FAULT


def test_names_the_faulted_node_so_the_wrapper_can_exclude_it(workspace):
    write_syscheck(workspace, 1)
    result = run(workspace)
    assert "MFC_FAULT_NODE=atl1-1-03-007-29-0" in result.stdout + result.stderr


def test_shows_syscheck_output_when_it_fails(workspace):
    # The previous post-build validation sent syscheck's output to /dev/null,
    # which left the CI log with no evidence of why the node was rejected.
    write_syscheck(workspace, 1, message="uncorrectable ECC error encountered")
    result = run(workspace)
    assert "uncorrectable ECC error encountered" in result.stdout + result.stderr


def test_passes_when_no_syscheck_binary_was_built(workspace):
    # A missing binary is a build problem, not a bad node; failing here would
    # requeue onto a healthy node and fail identically.
    assert run(workspace).returncode == HEALTHY


def test_skips_when_the_cluster_is_already_known_to_be_down(workspace):
    write_syscheck(workspace, 0)
    subprocess.run(
        ["bash", str(SCRIPTS / "ci-outage.sh"), "mark", "phoenix", "pypi unreachable"],
        env={**os.environ, "MFC_CI_STATE_DIR": str(workspace / "state")},
        capture_output=True,
        check=True,
    )
    assert run(workspace).returncode == OUTAGE


def test_does_not_report_a_node_fault_merely_because_pmix_printed_a_warning(workspace):
    # PMIX_ERR_NO_PERMISSIONS in dstore_base.c is benign noise: it appears in
    # 16% of passing self-hosted jobs and only 9% of failing ones. Gating on it
    # would fail roughly one healthy job in six.
    write_syscheck(workspace, 0, message="PMIX ERROR: PMIX_ERR_NO_PERMISSIONS in file dstore_base.c at line 238")
    assert run(workspace).returncode == HEALTHY
