"""The case-optimization pre-build must actually probe, in both shard modes.

The probe needs a syscheck binary to exist. In the sharded mode shard 1 builds
one under the script's marker coordination and the others wait for it; in the
unsharded mode nothing has been built when the probe runs, so it has to build
the binary itself.

That asymmetry failed silently the first time: on Phoenix, where
case-optimization is unsharded, the probe reported "no syscheck binary under
build/install; skipping node probe" and did nothing at all. A guard that skips
looks exactly like a guard that passed, which is why this is pinned.
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
    _exe(
        binz / "mpirun",
        '#!/bin/bash\nwhile [ "${1:0:1}" = "-" ]; do shift; case "$1" in [0-9]*) shift;; esac; done\nexec "$@"\n',
    )
    trace = tmp_path / "trace.log"

    def install_mfc(probe_exit=0):
        _exe(
            tmp_path / "mfc.sh",
            f"""#!/bin/bash
echo "mfc.sh $*" >> {trace}
if [ "$1" = "load" ]; then return 0 2>/dev/null || exit 0; fi
for a in "$@"; do
  if [ "$a" = "syscheck" ]; then
    mkdir -p build/install/gpu-acc-abc/bin
    printf '#!/bin/bash\\nexit {probe_exit}\\n' > build/install/gpu-acc-abc/bin/syscheck
    chmod +x build/install/gpu-acc-abc/bin/syscheck
  fi
done
exit 0
""",
        )

    def run(shard="", cluster="phoenix"):
        env = {
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            "SLURM_JOB_ID": "123456",
            "SLURMD_NODENAME": "atl1-1-02-008-1-1",
            "job_cluster": cluster,
            "job_device": "gpu",
            "job_interface": "acc",
            "job_shard": shard,
        }
        return subprocess.run(
            ["bash", ".github/scripts/prebuild-case-optimization.sh"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
            timeout=120,
        )

    return tmp_path, install_mfc, run, trace


def test_the_unsharded_prebuild_really_probes(workspace):
    # Phoenix runs case-optimization unsharded. This is the path that silently
    # skipped.
    _, install_mfc, run, _ = workspace
    install_mfc(probe_exit=0)
    out = run().stdout
    assert "Preflight: probing" in out
    assert "no syscheck binary" not in out


def test_a_bad_node_stops_the_unsharded_prebuild(workspace):
    _, install_mfc, run, trace = workspace
    install_mfc(probe_exit=1)
    result = run()
    assert result.returncode == 77
    assert "Pre-building" not in result.stdout


def test_the_sharded_prebuild_does_not_rebuild_syscheck_concurrently(workspace):
    # Shard 2 waits on shard 1's marker; building syscheck here would race the
    # shared build/install that the coordination exists to serialize.
    tmp_path, install_mfc, run, trace = workspace
    install_mfc(probe_exit=0)
    # frontier_amd, not phoenix: sharding is a frontier_amd mode, and phoenix's
    # clean_build would move build/ aside and take the marker with it.
    (tmp_path / "build").mkdir(exist_ok=True)
    (tmp_path / "build" / ".prebuild-shared-targets-done").touch()
    run(shard="2/3", cluster="frontier_amd")
    builds = [c for c in trace.read_text().splitlines() if "build -t syscheck" in c]
    assert builds == [], f"shard 2 must not build syscheck itself, got {builds}"
