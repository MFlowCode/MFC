"""Tests for .github/workflows/common/build.sh.

Two behaviours are pinned here.

Order: syscheck is a standalone target that links in 5-19 seconds and is already
built second, right after hipfort. Building it first and running it before the
solver build means a node with a dead GPU is rejected in about a minute instead
of after ~40 minutes of compilation. Across the Aug 2026 ECC failures the gap
between syscheck being available and the fault being noticed was a median of 38
minutes, 28.5 hours in total.

Outage: an unreachable pypi.org is not the node's fault and not fixable by
requeuing. The first job to hit it records it so the rest of the matrix skips
rather than each spending ~33 minutes rediscovering it.
"""

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

PYPI_FAILURE = "  x Failed to build `mfc @ file:///work/toolchain`\n" "  |-> Failed to fetch: `https://pypi.org/simple/hatch-vcs/`\n" "mfc: ERROR > (venv) Installation failed.\n"


def _exe(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def workspace(tmp_path):
    shutil.copytree(REPO / ".github", tmp_path / ".github")
    binz = tmp_path / "bin"
    binz.mkdir()
    passthrough = '#!/bin/bash\nwhile [ "${1:0:1}" = "-" ]; do shift; case "$1" in [0-9]*) shift;; esac; done\nexec "$@"\n'
    _exe(binz / "mpirun", passthrough)
    # frontier launches with srun; stub it so the real /usr/bin/srun on this box
    # cannot stand in and try to submit an actual job.
    _exe(binz / "srun", passthrough)
    trace = tmp_path / "trace.log"

    def install_mfc(full_build_stdout="", full_build_rc=0, syscheck_rc=0, probe_build_stdout="", probe_build_rc=0):
        _exe(
            tmp_path / "mfc.sh",
            f"""#!/bin/bash
echo "mfc.sh $*" >> {trace}
for a in "$@"; do
  if [ "$a" = "syscheck" ]; then
    printf '%s' {shlex_quote(probe_build_stdout)}
    if [ {probe_build_rc} -ne 0 ]; then exit {probe_build_rc}; fi
    mkdir -p build/install/gpu/bin
    printf '#!/bin/bash\\necho probe ran\\nexit {syscheck_rc}\\n' > build/install/gpu/bin/syscheck
    chmod +x build/install/gpu/bin/syscheck
    exit 0
  fi
done
printf '%s' {shlex_quote(full_build_stdout)}
exit {full_build_rc}
""",
        )

    def run():
        env = {
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            "MFC_BUILD_RETRY_DELAY": "0",
            "job_device": "gpu",
            "job_interface": "acc",
            "job_shard": "",
            "job_cluster": "frontier",
            "job_variant": "",
            "SLURMD_NODENAME": "frontier1234",
            # These scripts only ever run inside a SLURM allocation
            # (submit-slurm-job.sh is their sole caller), and the probe
            # refuses to judge a node outside one.
            "SLURM_JOB_ID": "123456",
        }
        return subprocess.run(
            ["bash", ".github/workflows/common/build.sh"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
            timeout=120,
        )

    return tmp_path, install_mfc, run, trace


def shlex_quote(s):
    import shlex

    return shlex.quote(s)


def outage_recorded(tmp_path):
    state = tmp_path / "state"
    return state.exists() and any(state.glob("outage-*"))


def test_the_probe_is_built_before_the_solver(workspace):
    tmp_path, install_mfc, run, trace = workspace
    install_mfc()
    assert run().returncode == 0
    calls = trace.read_text().splitlines()
    assert any("syscheck" in c for c in calls), calls
    assert "syscheck" in calls[0], f"probe must be built first, got {calls}"


def test_a_failing_probe_stops_before_the_solver_build(workspace):
    tmp_path, install_mfc, run, trace = workspace
    install_mfc(syscheck_rc=1)
    result = run()
    assert result.returncode == 77
    assert len(trace.read_text().splitlines()) == 1, "solver build must not be attempted"


def test_a_pypi_failure_is_an_ordinary_build_failure(workspace):
    """A flaky download must not red out the rest of the cluster.

    The dependency install happens before any compute is committed -- on the
    login node for Frontier -- so a failed fetch costs no allocation and there
    is nothing to protect the matrix from. Recording it as a cluster-wide outage
    skipped every other job on that cluster, including ones whose tests had
    already passed, and uv already retries internally.
    """
    tmp_path, install_mfc, run, _ = workspace
    install_mfc(full_build_stdout=PYPI_FAILURE, full_build_rc=1)

    result = run()

    assert result.returncode == 1, "a failed download must fail only its own job"
    assert not outage_recorded(tmp_path), "a failed download must not trip the cluster breaker"


def test_an_ordinary_compile_error_is_not_treated_as_an_outage(workspace):
    tmp_path, install_mfc, run, _ = workspace
    install_mfc(full_build_stdout="NVFORTRAN-S-0034-Syntax error\n", full_build_rc=1)
    result = run()
    assert not outage_recorded(tmp_path)
    assert result.returncode not in (77, 78)


def test_the_build_output_is_still_shown_when_it_fails(workspace):
    # The CCE/flang failures that report no compiler diagnostic are only
    # debuggable if whatever the build did print survives into the CI log.
    tmp_path, install_mfc, run, _ = workspace
    install_mfc(full_build_stdout="ftn-2116 ftn: INTERNAL\n", full_build_rc=1)
    result = run()
    assert "ftn-2116" in result.stdout + result.stderr


def test_an_ordinary_probe_build_failure_is_not_an_outage(workspace):
    tmp_path, install_mfc, run, _ = workspace
    install_mfc(probe_build_stdout="NVFORTRAN-S-0034-Syntax error\n", probe_build_rc=1)
    result = run()
    assert not outage_recorded(tmp_path)
    assert result.returncode not in (77, 78)


def test_the_probe_build_output_is_shown_when_it_fails(workspace):
    tmp_path, install_mfc, run, _ = workspace
    install_mfc(probe_build_stdout="ftn-2116 ftn: INTERNAL\n", probe_build_rc=1)
    result = run()
    assert "ftn-2116" in result.stdout + result.stderr
