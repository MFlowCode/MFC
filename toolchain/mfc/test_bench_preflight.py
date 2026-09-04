"""The benchmark and case-optimization paths must probe their node too.

The probe was originally wired only into common/build.sh and common/test.sh,
which the test.yml `self` job uses. But three of the five entry points that
submit SLURM jobs -- the benchmark job and both case-optimization jobs -- go
straight to submit-slurm-job.sh and never touch those scripts. That left 134 of
540 first-attempt failures on unprobed paths, including every ECC failure
outside the `self` job.

Ordering is load-bearing here, hence the second test. These scripts nuke and
rebuild build/ themselves -- bench.sh does so on Phoenix precisely because "compute
nodes are heterogeneous -> ISA mismatch risk" -- so a probe placed before that
would test a leftover binary from a previous job, quite possibly built for
another microarchitecture. The SIGILL would be reported as a bad node, and the
wrapper would exclude a perfectly healthy one.
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
    trace = tmp_path / "trace.log"

    passthrough = '#!/bin/bash\nwhile [ "${1:0:1}" = "-" ]; do shift; case "$1" in [0-9]*) shift;; esac; done\nexec "$@"\n'
    for name in ("mpirun", "srun"):
        _exe(binz / name, passthrough)
    _exe(binz / "nvidia-smi", "#!/bin/bash\necho 'GPU 0: fake'\n")

    def install_mfc(probe_exit=0):
        _exe(
            tmp_path / "mfc.sh",
            f"""#!/bin/bash
echo "mfc.sh $*" >> {trace}
for a in "$@"; do
  if [ "$a" = "build" ]; then
    mkdir -p build/install/gpu/bin
    printf '#!/bin/bash\\nexit {probe_exit}\\n' > build/install/gpu/bin/syscheck
    chmod +x build/install/gpu/bin/syscheck
  fi
done
exit 0
""",
        )
        # A stale binary from a "previous job" that a too-early probe would find.
        stale = tmp_path / "build" / "install" / "stale-gpu" / "bin"
        stale.mkdir(parents=True, exist_ok=True)
        _exe(stale / "syscheck", "#!/bin/bash\necho 'stale SIGILL'\nexit 132\n")

    def run():
        env = {
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            "MFC_BUILD_RETRY_DELAY": "0",
            "SLURM_JOB_ID": "123456",
            "SLURMD_NODENAME": "frontier4242",
            "job_device": "gpu",
            "job_interface": "acc",
            "job_slug": "bench-gpu-acc",
            "job_cluster": "frontier",
        }
        return subprocess.run(
            ["bash", ".github/workflows/common/bench.sh"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
            timeout=120,
        )

    return tmp_path, install_mfc, run, trace


def test_the_benchmark_job_probes_its_node(workspace):
    _, install_mfc, run, _ = workspace
    install_mfc(probe_exit=0)
    result = run()
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Preflight:" in result.stdout


def test_a_bad_node_stops_the_benchmark_before_it_runs(workspace):
    _, install_mfc, run, trace = workspace
    install_mfc(probe_exit=1)
    result = run()
    assert result.returncode == 77
    assert "mfc.sh bench" not in trace.read_text()


def test_the_probe_runs_after_the_build_not_before_it(workspace):
    # If the probe ran first it would find the stale binary planted above, which
    # exits 132 like a SIGILL, and condemn a healthy node.
    _, install_mfc, run, trace = workspace
    install_mfc(probe_exit=0)
    result = run()
    assert result.returncode == 0, result.stdout + result.stderr
    calls = trace.read_text().splitlines()
    assert any("build" in c for c in calls), calls
    build_at = next(i for i, c in enumerate(calls) if "build" in c)
    probe_at = result.stdout.index("Preflight:")
    build_marker = result.stdout.find("mfc.sh build")
    assert build_at == 0
    if build_marker != -1:
        assert probe_at > build_marker
