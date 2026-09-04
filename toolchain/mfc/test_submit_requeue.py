"""Integration tests for the resubmit loop in .github/scripts/submit-slurm-job.sh.

The loop already resubmits on preemption (exit 76). These tests cover the node
-fault path (exit 77): the wrapper must add the faulted node to sbatch's
--exclude and try again, so a persistently bad node stops eating jobs without
anyone hand-editing a constant and committing it.

SLURM is stubbed. The scripts are copied into a temp tree so the monitor can be
replaced with one that returns a scripted sequence of exit codes; nothing here
waits on a scheduler.
"""

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"


def _exe(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def rig(tmp_path):
    """A workspace with stubbed SLURM commands and a scriptable monitor."""
    binz = tmp_path / "bin"
    binz.mkdir()
    scripts = tmp_path / "scripts"
    shutil.copytree(SCRIPTS, scripts)
    submissions = tmp_path / "submissions"
    submissions.mkdir()

    # submit-slurm-job.sh inlines the payload script into the sbatch heredoc.
    payload = tmp_path / ".github" / "workflows" / "common"
    payload.mkdir(parents=True)
    (payload / "build-and-test.sh").write_text("echo payload\n")

    # sbatch records each submitted script, then reports a fresh job id.
    _exe(
        binz / "sbatch",
        f"""#!/bin/bash
n=$(ls {submissions} | wc -l)
cat > {submissions}/submission-$n.sh
echo "Submitted batch job 100$n"
""",
    )
    for name in ("squeue", "scancel", "scontrol"):
        _exe(binz / name, "#!/bin/bash\nexit 0\n")
    _exe(binz / "sacct", "#!/bin/bash\necho COMPLETED\n")
    _exe(binz / "sinfo", "#!/bin/bash\necho idle\n")

    # The monitor returns the next code from MONITOR_CODES on each call.
    _exe(
        scripts / "run_monitored_slurm_job.sh",
        f"""#!/bin/bash
n=$(cat {tmp_path}/monitor_calls 2>/dev/null || echo 0)
echo $((n + 1)) > {tmp_path}/monitor_calls
code=$(echo "$MONITOR_CODES" | cut -d, -f$((n + 1)))
# Mimic the preflight's marker landing in the job output file.
if [ "$code" = "77" ]; then
    echo "MFC_FAULT_NODE=badnode-$n" >> "$2"
fi
exit "${{code:-0}}"
""",
    )

    def run(monitor_codes, script="common/build-and-test.sh", **env_extra):
        env = {
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MONITOR_CODES": monitor_codes,
            "GITHUB_EVENT_NAME": "pull_request",
            "MFC_CI_STATE_DIR": str(tmp_path / "state"),
            **env_extra,
        }
        result = subprocess.run(
            ["bash", str(scripts / "submit-slurm-job.sh"), f".github/workflows/{script}", "gpu", "acc", "phoenix"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            check=False,
        )
        result.submissions = sorted(submissions.glob("submission-*.sh"))
        return result

    run.scripts = scripts
    run.state_dir = tmp_path / "state"
    return run


def excludes(submission: Path):
    for line in submission.read_text().splitlines():
        if "--exclude=" in line:
            return line.split("--exclude=", 1)[1].strip().strip('"')
    return ""


def test_a_healthy_job_is_submitted_exactly_once(rig):
    result = rig("0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(result.submissions) == 1


def test_a_node_fault_causes_a_resubmission(rig):
    result = rig("77,0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(result.submissions) == 2


def test_the_resubmission_excludes_the_faulted_node(rig):
    result = rig("77,0")
    assert "badnode-0" in excludes(result.submissions[1])


def test_the_resubmission_keeps_the_nodes_already_excluded(rig):
    result = rig("77,0")
    first = excludes(result.submissions[0])
    assert first, "expected the baseline --exclude list to be non-empty"
    for node in first.split(","):
        assert node in excludes(result.submissions[1])


def test_node_faults_stop_after_the_bounded_number_of_resubmits(rig):
    result = rig("77,77,77,77,77", MFC_MAX_NODE_RESUBMITS="2")
    assert result.returncode != 0
    assert len(result.submissions) == 3  # original + 2 resubmits


def test_a_known_outage_is_not_resubmitted(rig):
    # Requeuing cannot fix pypi.org being unreachable; trying again just burns
    # another allocation.
    result = rig("78")
    assert len(result.submissions) == 1


def test_the_default_bound_is_one_requeue(rig):
    """One requeue, not two.

    A wrong probe costs one node per attempt: the run that condemned three
    healthy Phoenix nodes was a bounded loop doing exactly what it was told.
    Bad nodes are concentrated -- one accounted for 25 of 29 ECC failures -- so a
    single requeue captures nearly all the benefit at half the blast radius.
    """
    result = rig("77,77,77,77")
    assert result.returncode != 0
    assert len(result.submissions) == 2  # original + 1 requeue
