"""The infrastructure exit code must survive the trip back to the submit wrapper.

preflight.sh exits 77 (node fault) inside the SLURM job.
That becomes the job's ExitCode, which monitor_slurm_job.sh reads and
run_monitored_slurm_job.sh relays to the resubmit loop. If either layer flattens
it to 1 -- as both did for every non-zero code before -- the loop sees a
generic failure and the node is never excluded.
"""

import os
import stat
import subprocess
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"


def _exe(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def slurm(tmp_path):
    """Stub SLURM reporting a finished job whose ExitCode the test chooses."""
    binz = tmp_path / "bin"
    binz.mkdir()

    def configure(exit_code, state="COMPLETED"):
        _exe(binz / "squeue", "#!/bin/bash\nexit 0\n")
        _exe(binz / "sacct", f'#!/bin/bash\nfor a in "$@"; do [ "$a" = "--format=ExitCode" ] && {{ echo "{exit_code}"; exit 0; }}; done\necho "{state}"\n')
        _exe(binz / "scontrol", f'#!/bin/bash\necho "ExitCode={exit_code}"\n')
        _exe(binz / "scancel", "#!/bin/bash\nexit 0\n")
        out = tmp_path / "job.out"
        out.write_text("job output\n")
        return out

    return tmp_path, binz, configure


def run_script(tmp_path, binz, name, *args):
    env = {
        **os.environ,
        "PATH": f"{binz}:{os.environ['PATH']}",
        # This script really sleeps between polls; without these the file
        # cost ~36s in every PR's Lint Gate.
        "MFC_MONITOR_POLL_SECONDS": "0",
        "MFC_MONITOR_RECHECK_SECONDS": "0",
    }
    return subprocess.run(
        ["bash", str(SCRIPTS / name), *args],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=env,
        check=False,
        timeout=180,
    )


@pytest.mark.parametrize("job_exit,expected", [("77:0", 77)])
def test_monitor_relays_the_infrastructure_exit_code(slurm, job_exit, expected):
    tmp_path, binz, configure = slurm
    out = configure(job_exit)
    assert run_script(tmp_path, binz, "monitor_slurm_job.sh", "1234", str(out)).returncode == expected


def test_monitor_still_reports_an_ordinary_failure_as_one(slurm):
    tmp_path, binz, configure = slurm
    out = configure("2:0")
    assert run_script(tmp_path, binz, "monitor_slurm_job.sh", "1234", str(out)).returncode == 1


@pytest.mark.parametrize("monitor_exit", [77])
def test_the_runner_relays_the_infrastructure_exit_code(tmp_path, monitor_exit):
    # Stub the inner monitor so this exercises only the relaying layer.
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    for name in ("run_monitored_slurm_job.sh", "monitor_slurm_job.sh"):
        (scripts / name).write_text((SCRIPTS / name).read_text())
    _exe(scripts / "monitor_slurm_job.sh", f"#!/bin/bash\nexit {monitor_exit}\n")

    result = subprocess.run(
        ["bash", str(scripts / "run_monitored_slurm_job.sh"), "1234", str(tmp_path / "job.out")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
        timeout=180,
    )
    assert result.returncode == monitor_exit
