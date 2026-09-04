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
    """A bare workspace. No launcher is on PATH unless a test installs one.

    An earlier version of this fixture put a fake mpirun on PATH for every test,
    which is exactly why it could not see that preflight ran mpirun on Frontier,
    where the launcher is srun and Cray MPICH ships no mpirun at all.
    """
    (tmp_path / "bin").mkdir()
    (tmp_path / "state").mkdir()

    # A hermetic PATH holding only the utilities the scripts need. This box has
    # a real mpirun *and* a real /usr/bin/srun, either of which would silently
    # stand in for a launcher the test meant to be absent.
    sysbin = tmp_path / "sysbin"
    sysbin.mkdir()
    for tool in ("bash", "find", "head", "tail", "cat", "sed", "grep", "tr", "cut", "date", "mkdir", "mv", "rm", "hostname", "env", "sort", "wc", "dirname", "basename"):
        for root in ("/usr/bin", "/bin"):
            src = Path(root) / tool
            if src.exists():
                (sysbin / tool).symlink_to(src)
                break
    return tmp_path


def install_launcher(workspace, name):
    """A passthrough launcher that records its argv, mirroring mpirun/srun."""
    launcher = workspace / "bin" / name
    # Skips options and their values by looking for the first executable
    # argument, so it does not need to know which flags take a value.
    launcher.write_text("#!/bin/bash\n" f'echo "$@" >> {workspace}/launched.txt\n' 'while [ $# -gt 0 ] && [ ! -x "$1" ]; do shift; done\n' 'exec "$@"\n')
    launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC)
    return launcher


def launched_with(workspace):
    path = workspace / "launched.txt"
    return path.read_text() if path.exists() else ""


def write_syscheck(workspace, exit_code, message="syscheck says hello"):
    target = workspace / "build" / "install" / "gpu-acc-abc123" / "bin"
    target.mkdir(parents=True, exist_ok=True)
    binary = target / "syscheck"
    binary.write_text(f'#!/bin/bash\necho "{message}"\nexit {exit_code}\n')
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    return binary


def run(workspace, *args, **overrides):
    env = {
        **os.environ,
        # A controlled PATH, not the inherited one: this box has a real mpirun,
        # and inheriting it makes "no launcher available" silently untestable.
        "PATH": f"{workspace / 'bin'}:{workspace / 'sysbin'}",
        "MFC_CI_STATE_DIR": str(workspace / "state"),
        "SLURMD_NODENAME": "atl1-1-03-007-29-0",
        # Present by default: these tests describe behaviour inside a job.
        "SLURM_JOB_ID": "123456",
    }
    if overrides.pop("MFC_NO_SLURM", None):
        env.pop("SLURM_JOB_ID", None)
    env.update(overrides)
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


def test_does_not_report_a_node_fault_merely_because_pmix_printed_a_warning(workspace):
    # PMIX_ERR_NO_PERMISSIONS in dstore_base.c is benign noise: it appears in
    # 16% of passing self-hosted jobs and only 9% of failing ones. Gating on it
    # would fail roughly one healthy job in six.
    write_syscheck(workspace, 0, message="PMIX ERROR: PMIX_ERR_NO_PERMISSIONS in file dstore_base.c at line 238")
    assert run(workspace).returncode == HEALTHY


def test_uses_mpirun_on_phoenix(workspace):
    install_launcher(workspace, "mpirun")
    write_syscheck(workspace, 0)
    assert run(workspace, "phoenix", "gpu").returncode == HEALTHY
    assert "syscheck" in launched_with(workspace)


def test_uses_srun_on_frontier(workspace):
    # Frontier and frontier_amd launch every binary with srun
    # (toolchain/templates/frontier.mako, frontier_amd.mako); Cray MPICH ships no
    # mpirun, so running one there is not merely wrong, it cannot work.
    install_launcher(workspace, "srun")
    write_syscheck(workspace, 0)
    assert run(workspace, "frontier", "gpu").returncode == HEALTHY
    assert "syscheck" in launched_with(workspace)


def test_a_missing_launcher_is_not_blamed_on_the_node(workspace):
    # Without this, a launcher absent from PATH returns 127, preflight calls the
    # node bad, and the requeue loop burns three allocations and blacklists three
    # healthy nodes before declaring a cluster-wide problem.
    write_syscheck(workspace, 0)
    assert run(workspace, "frontier", "gpu").returncode == HEALTHY


def test_it_refuses_to_judge_a_node_outside_a_slurm_allocation(workspace):
    # mfc.sh load is used for building on login nodes too (bench.yml and
    # frontier/build.sh both load the GPU module set there). A probe that ran in
    # that context would find no usable GPU, call the login node bad, and have
    # the wrapper exclude it and requeue. Only a real allocation can be judged.
    write_syscheck(workspace, 1)
    env_without_slurm = {"MFC_NO_SLURM": "1"}
    assert run(workspace, "phoenix", "gpu", **env_without_slurm).returncode == HEALTHY


def test_a_gpu_binary_in_a_cpu_allocation_is_not_a_node_fault(workspace):
    """The failure that condemned three healthy Phoenix nodes.

    The case-optimization pre-build is submitted as a cpu allocation but builds
    GPU binaries, so the only syscheck available asserts a device exists and
    exits non-zero on a node that has no GPU by design. That is a mismatch
    between what the job asked for and what it built -- never evidence about the
    node -- so the probe must decline to judge rather than blame it.
    """
    install_launcher(workspace, "mpirun")
    target = workspace / "build" / "install" / "gpu-mp-abc123" / "bin"
    target.mkdir(parents=True)
    probe = target / "syscheck"
    probe.write_text("#!/bin/bash\necho 'num_devices == 0'\nexit 1\n")
    probe.chmod(probe.stat().st_mode | stat.S_IEXEC)

    result = run(workspace, "phoenix", "cpu")
    assert result.returncode == HEALTHY
    assert "MFC_FAULT_NODE" not in result.stdout + result.stderr


def test_a_gpu_binary_in_a_gpu_allocation_is_still_judged(workspace):
    # The guard above must not blunt the actual feature.
    install_launcher(workspace, "mpirun")
    target = workspace / "build" / "install" / "gpu-mp-abc123" / "bin"
    target.mkdir(parents=True)
    probe = target / "syscheck"
    probe.write_text("#!/bin/bash\nexit 1\n")
    probe.chmod(probe.stat().st_mode | stat.S_IEXEC)

    assert run(workspace, "phoenix", "gpu").returncode == NODE_FAULT
