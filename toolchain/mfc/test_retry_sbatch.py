"""retry_sbatch (.github/scripts/retry-sbatch.sh) after an sbatch that timed out.

A timeout does not mean nothing was queued: on Phoenix the first of two shard-1 submissions
timed out yet ran, the retry ran too, and the two copies deleted each other's test inputs.
SLURM is stubbed; `sleep` is a no-op function.
"""

import stat
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "retry-sbatch.sh"


def _exe(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture
def submit(tmp_path):
    binz = tmp_path / "bin"
    binz.mkdir()
    calls = tmp_path / "sbatch_calls"
    # The first sbatch times out (having queued job 500); later ones succeed.
    _exe(
        binz / "sbatch",
        f"""#!/bin/bash
cat > /dev/null
echo x >> {calls}
[ "$(wc -l < {calls})" -eq 1 ] && {{ echo "sbatch: error: Socket timed out on send/recv operation"; exit 1; }}
echo "Submitted batch job 600"
""",
    )

    def run(queued):
        # queued: the squeue line for job 500, or "" when the timed-out submission queued nothing.
        _exe(binz / "squeue", f"#!/bin/bash\nprintf '%s' {queued!r}\n")
        result = subprocess.run(
            ["bash", "-c", f'PATH="{binz}:$PATH"; source "{SCRIPT}"; sleep() {{ :; }}; retry_sbatch "$1"', "_", "#!/bin/bash\n#SBATCH -J shb-test-1-of-2\n"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            check=False,
            timeout=60,
        )
        n = len(calls.read_text().split()) if calls.exists() else 0
        return result, n

    return tmp_path, run


def test_a_timed_out_submission_that_was_queued_is_used_not_submitted_again(submit):
    tmp_path, run = submit
    result, n = run(f"500|{tmp_path}|2999-01-01T00:00:00\n")
    assert result.returncode == 0 and result.stdout.strip() == "500" and n == 1, result.stderr


def test_a_timed_out_submission_that_queued_nothing_is_retried(submit):
    _, run = submit
    result, n = run("")
    assert result.returncode == 0 and result.stdout.strip() == "600" and n == 2, result.stderr


@pytest.mark.parametrize("where", ["other-dir", "older"])
def test_a_job_from_another_directory_or_an_earlier_run_is_not_taken(submit, where):
    tmp_path, run = submit
    d, t = ("/elsewhere", "2999-01-01T00:00:00") if where == "other-dir" else (str(tmp_path), "2000-01-01T00:00:00")
    result, n = run(f"500|{d}|{t}\n")
    assert result.returncode == 0 and result.stdout.strip() == "600" and n == 2, result.stderr
