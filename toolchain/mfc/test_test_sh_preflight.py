"""The test allocation must probe its own node before running the suite.

On Frontier the Build and Test steps are two separate SLURM submissions with no
node affinity, and they landed on different nodes in 26 of 29 measurable
post-#1763 jobs. In the Aug 2026 ECC failures the build node was healthy every
time while the tests landed on a bad one, so a probe that only runs during the
build proves nothing about the machine that runs the tests.
"""

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


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
    # job_cluster here is frontier, whose launcher is srun. Without a stub the
    # real /usr/bin/srun on this box would try to submit an actual job.
    _exe(binz / "srun", passthrough)
    _exe(binz / "nvidia-smi", "#!/bin/bash\necho 'GPU 0: fake'\n")
    trace = tmp_path / "trace.log"
    _exe(tmp_path / "mfc.sh", f'#!/bin/bash\necho "mfc.sh $*" >> {trace}\nexit 0\n')

    def install_probe(exit_code):
        target = tmp_path / "build" / "install" / "gpu" / "bin"
        target.mkdir(parents=True, exist_ok=True)
        _exe(target / "syscheck", f'#!/bin/bash\necho "probe ran"\nexit {exit_code}\n')

    def run():
        env = {
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            "job_device": "gpu",
            "job_interface": "acc",
            "job_shard": "",
            "job_cluster": "frontier",
            "job_variant": "",
            "GITHUB_EVENT_NAME": "push",
            "SLURMD_NODENAME": "frontier9999",
            # These scripts only ever run inside a SLURM allocation
            # (submit-slurm-job.sh is their sole caller), and the probe
            # refuses to judge a node outside one.
            "SLURM_JOB_ID": "123456",
        }
        return subprocess.run(
            ["bash", ".github/workflows/common/test.sh"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
            timeout=120,
        )

    return install_probe, run, trace


def test_a_healthy_node_runs_the_suite(workspace):
    install_probe, run, trace = workspace
    install_probe(0)
    assert run().returncode == 0
    assert "test" in trace.read_text()


def test_a_bad_node_is_reported_before_the_suite_starts(workspace):
    install_probe, run, trace = workspace
    install_probe(1)
    result = run()
    assert result.returncode == 77
    assert not trace.exists() or "test" not in trace.read_text()
