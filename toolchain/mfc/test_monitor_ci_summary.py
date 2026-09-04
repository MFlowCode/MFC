"""What a failed job puts on the run's summary page.

Learning why a job failed used to mean opening a log of tens of thousands of
lines, and an infrastructure fault looked exactly like a test failure until you
did. These paths were shipped once with only a syntax check behind them, so
they are exercised here for real.
"""

import os
import stat
import subprocess
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"

GPU_FAULT_OUTPUT = """\
[  0%]  Time step 1
=== GPU fault summary (rocm-debug-agent, 14276 lines collapsed) ===
Memory access fault by GPU node-5 on address 0x1557da700000.
faulting kernel(s): s_tvd_rk$m_time_steppers_$ck_L486_6
faulting waves:     125  [MEMORY_VIOLATION x125]
NOTE: waves halt on fault detection.
"""


@pytest.fixture
def slurm(tmp_path):
    """A finished SLURM job whose exit code and output the test chooses."""
    binz = tmp_path / "bin"
    binz.mkdir()

    def configure(exit_code, output):
        def exe(name, body):
            path = binz / name
            path.write_text(body)
            path.chmod(path.stat().st_mode | stat.S_IEXEC)

        exe("squeue", "#!/bin/bash\nexit 0\n")
        exe("sacct", f'#!/bin/bash\nfor a in "$@"; do [ "$a" = "--format=ExitCode" ] && {{ echo "{exit_code}"; exit 0; }}; done\necho FAILED\n')
        exe("scontrol", f'#!/bin/bash\necho "ExitCode={exit_code}"\n')
        exe("scancel", "#!/bin/bash\nexit 0\n")

        out = tmp_path / "job.out"
        out.write_text(output)
        summary = tmp_path / "summary.md"
        summary.write_text("")
        return out, summary

    return tmp_path, binz, configure


def run(tmp_path, binz, out, summary):
    return subprocess.run(
        ["bash", str(SCRIPTS / "monitor_slurm_job.sh"), "1234", str(out)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
        timeout=180,
        env={
            **os.environ,
            "PATH": f"{binz}:{os.environ['PATH']}",
            "MFC_MONITOR_POLL_SECONDS": "0",
            "MFC_MONITOR_RECHECK_SECONDS": "0",
            "GITHUB_STEP_SUMMARY": str(summary),
        },
    )


def test_an_infrastructure_fault_names_the_node_and_says_it_is_not_a_test_failure(slurm):
    tmp_path, binz, configure = slurm
    out, summary = configure("77:0", "some output\nMFC_FAULT_NODE=frontier9999\n")

    assert run(tmp_path, binz, out, summary).returncode == 77

    text = summary.read_text()
    assert "frontier9999" in text
    assert "not a code or test failure" in text


def test_a_gpu_fault_puts_the_faulting_kernel_on_the_summary_page(slurm):
    tmp_path, binz, configure = slurm
    out, summary = configure("1:0", GPU_FAULT_OUTPUT)

    assert run(tmp_path, binz, out, summary).returncode == 1

    text = summary.read_text()
    assert "s_tvd_rk$m_time_steppers_$ck_L486_6" in text, "the whole point is naming the kernel"
    assert "GPU memory fault" in text


def test_a_successful_job_writes_nothing(slurm):
    tmp_path, binz, configure = slurm
    out, summary = configure("0:0", "all good\n")

    assert run(tmp_path, binz, out, summary).returncode == 0
    assert summary.read_text() == ""


def test_the_job_output_is_not_printed_twice(slurm):
    """It used to be: streamed with tail -f, then cat in full at the end.

    Measured at three copies of every line on a GPU job, and 65,000 lines of
    offload diagnostics repeated for a single fault.
    """
    tmp_path, binz, configure = slurm
    lines = 400
    out, summary = configure("1:0", "\n".join(f"job output line {n}" for n in range(lines)) + "\n")

    result = run(tmp_path, binz, out, summary)

    # Volume, not a marker: `tail -f` on a file that is already complete shows
    # only its last lines, whereas in CI it follows the file as it is written.
    # What must hold either way is that the whole file is not emitted again --
    # before the fix, stdout carried a full copy of it.
    assert result.stdout.count("job output line ") < lines, "the whole log was reprinted"
