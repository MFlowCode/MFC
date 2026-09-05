"""Tests for .github/workflows/frontier/build.sh.

Frontier installs its Python dependencies on the login node, in the "Fetch
Dependencies" step, before any SLURM job exists -- so a failed download costs no
allocation and must fail only its own job. It used to be recorded as a
cluster-wide outage, which then skipped every other job on that cluster.

`outage_recorded` stays as an assertion that nothing does that any more.
"""

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

PYPI_FAILURE = "  |-> Failed to fetch: `https://pypi.org/simple/hatch-vcs/`\n" "mfc: ERROR > (venv) Installation failed.\n"


def _exe(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def workspace(tmp_path):
    shutil.copytree(REPO / ".github", tmp_path / ".github")

    def install_mfc(stdout="", rc=0):
        # `. ./mfc.sh load` is sourced, so the stub must return rather than exit
        # for that subcommand or it would terminate build.sh itself.
        _exe(
            tmp_path / "mfc.sh",
            f'#!/bin/bash\nif [ "$1" = "load" ]; then return 0 2>/dev/null || exit 0; fi\n' f"printf '%s' {_q(stdout)}\nexit {rc}\n",
        )

    def run():
        env = {
            **os.environ,
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            "MFC_BUILD_RETRY_DELAY": "0",
        }
        return subprocess.run(
            ["bash", ".github/workflows/frontier/build.sh", "gpu", "acc"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
            timeout=120,
        )

    return tmp_path, install_mfc, run


def _q(s):
    import shlex

    return shlex.quote(s)


def outage_recorded(tmp_path):
    state = tmp_path / "state"
    return state.exists() and any(state.glob("outage-*"))


def test_a_successful_dependency_fetch_records_nothing(workspace):
    tmp_path, install_mfc, run = workspace
    install_mfc()
    assert run().returncode == 0
    assert not outage_recorded(tmp_path)


def test_a_pypi_failure_is_an_ordinary_build_failure(workspace):
    """A flaky download must not red out the rest of the cluster.

    The dependency install happens before any compute is committed -- on the
    login node for Frontier -- so a failed fetch costs no allocation and there
    is nothing to protect the matrix from. Recording it as a cluster-wide outage
    skipped every other job on that cluster, including ones whose tests had
    already passed, and uv already retries internally.
    """
    tmp_path, install_mfc, run = workspace
    install_mfc(stdout=PYPI_FAILURE, rc=1)

    result = run()

    assert result.returncode == 1, "a failed download must fail only its own job"
    assert not outage_recorded(tmp_path), "a failed download must not trip the cluster breaker"


def test_an_ordinary_dependency_failure_is_not_an_outage(workspace):
    tmp_path, install_mfc, run = workspace
    install_mfc(stdout="CMake Error: could not find HDF5\n", rc=1)
    result = run()
    assert not outage_recorded(tmp_path)
    assert result.returncode != 78
