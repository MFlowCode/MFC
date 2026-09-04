"""The preflight must judge the node, not the launcher.

Open MPI's default binding fails on some Phoenix nodes before the binary is
launched at all ("hwloc_set_cpubind returned Error"). Read as a node fault that
excluded three healthy nodes and failed the run two jobs deep.
"""

import os
import stat
import subprocess
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"

LAUNCH_FAILURE = """\
--------------------------------------------------------------------------
Open MPI tried to bind a new process, but something went wrong.  The
process was killed without launching the target application.
  Error message:     hwloc_set_cpubind returned "Error" for bitmap "0"
--------------------------------------------------------------------------
mpirun was unable to start the specified application as it encountered an error:
Error name: The specified application failed to start
"""


def _fake_cluster(tmp_path, mpirun_body):
    """A tree with a syscheck binary and a stubbed mpirun."""
    binz = tmp_path / "bin"
    binz.mkdir()
    mpirun = binz / "mpirun"
    mpirun.write_text(mpirun_body)
    mpirun.chmod(mpirun.stat().st_mode | stat.S_IEXEC)

    install = tmp_path / "build" / "install" / "gpu-acc-test" / "bin"
    install.mkdir(parents=True)
    (install / "syscheck").write_text("#!/bin/bash\nexit 1\n")
    (install / "syscheck").chmod(0o755)
    return binz


def _run(tmp_path, binz):
    return subprocess.run(
        ["bash", str(SCRIPTS / "preflight.sh"), "phoenix", "gpu"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
        timeout=120,
        env={**os.environ, "PATH": f"{binz}:{os.environ['PATH']}", "SLURM_JOB_ID": "1", "SLURMD_NODENAME": "n1"},
    )


def test_a_launch_failure_does_not_condemn_the_node(tmp_path):
    binz = _fake_cluster(tmp_path, f"#!/bin/bash\ncat <<'EOF'\n{LAUNCH_FAILURE}EOF\nexit 1\n")
    result = _run(tmp_path, binz)
    assert result.returncode == 0, f"excluded a node for a launcher failure:\n{result.stdout}"
    assert "not evidence about the node" in result.stdout


def test_a_binary_that_runs_and_fails_still_condemns_the_node(tmp_path):
    binz = _fake_cluster(tmp_path, '#!/bin/bash\nshift $((OPTIND)); echo "GPU init failed"; exit 1\n')
    result = _run(tmp_path, binz)
    assert result.returncode == 77, f"a real node fault must still be caught:\n{result.stdout}"


def test_the_probe_does_not_ask_open_mpi_to_bind(tmp_path):
    binz = _fake_cluster(tmp_path, '#!/bin/bash\necho "ARGS: $*"; exit 0\n')
    result = _run(tmp_path, binz)
    assert "--bind-to none" in result.stdout, "a single-rank probe must not rely on default binding"


def test_a_launcher_that_rejects_our_options_still_probes(tmp_path):
    """Otherwise the preflight switches itself off silently.

    A failed launch is treated as inconclusive, so a launcher that rejects an
    option we pass would fail every probe and quietly stop the preflight ever
    catching a real node fault -- worse than failing loudly.
    """
    binz = _fake_cluster(
        tmp_path,
        "#!/bin/bash\n" 'if [ "$1" = "--bind-to" ]; then echo "mpirun: unrecognized option \'--bind-to\'"; exit 1; fi\n' 'echo "RETRIED-WITHOUT-FLAGS"; exit 1\n',
    )
    result = _run(tmp_path, binz)

    assert "rejected the probe's options" in result.stdout
    assert "RETRIED-WITHOUT-FLAGS" in result.stdout


def test_a_login_node_is_never_judged(tmp_path):
    """`mfc.sh load` is also used for building on login nodes, where there is no
    GPU to probe and srun may not even be permitted."""
    binz = _fake_cluster(tmp_path, "#!/bin/bash\nexit 1\n")
    env = {**os.environ, "PATH": f"{binz}:{os.environ['PATH']}"}
    env.pop("SLURM_JOB_ID", None)
    result = subprocess.run(
        ["bash", str(SCRIPTS / "preflight.sh"), "phoenix", "gpu"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
        timeout=120,
        env=env,
    )
    assert result.returncode == 0
    assert "not inside a SLURM allocation" in result.stdout
