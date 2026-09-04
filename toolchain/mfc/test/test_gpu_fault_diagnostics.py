"""GPU memory faults should explain themselves the first time.

A fault reaches CI as an address and, unaided, nothing else. Measured against a
deliberate out-of-bounds write at m_time_steppers.fpp:486, on all four GPU
lanes: NVHPC and AFAR name the faulting kernel for free, CCE names nothing and
no CRAY_ACC_* variable helps, and the ROCm debug agent names it everywhere it
is reachable -- at ROCr level, with no recompile.

Nearly every failure in this area was silent, so these assert on content and on
the hand-offs between parts, never on presence alone.
"""

import pathlib
import subprocess
import sys

from mfc.gpu_diagnostics import (
    GPU_FAULT_MARKER,
    ROCM_DEBUG_AGENT,
    fault_diagnostic_env,
    is_gpu_memory_fault,
    rocm_debug_agent_path,
    summarize_rocm_debug_agent,
)

# Real agent output, verbatim, from two ROCm versions. The formats differ in
# ways that silently defeated a parser written against only one: 7.2.0 puts
# kernel_code_entry=/kernargs= between pc= and "(stopped, reason:", and words
# the fault line "OFFLOAD ERROR: memory access fault ... at virtual address".
ROCM_AGENT_FIXTURE_631 = """\
Memory access fault by GPU node-4 (Agent handle: 0x3b24f40) on address 0x7ffb6a0f6000. Reason: Write access to a read-only page.
Disassembly for function s_tvd_rk$m_time_steppers_$ck_L486_6:
    code object: file:///path/gpu-acc-c819d00b45/bin/simulation#offset=4657152&size=10843816
    loaded at: [0x7ff77da00000-0x7ff77feca9f9]
 => 0x7ff77e253430 <+6192>:    s_waitcnt vmcnt(0) lgkmcnt(0)
    0x7ff77e253434 <+6196>:    v_sub_co_u32_e32 v5, vcc, v26, v42
End of disassembly.
wave_0: pc=0x7ff77e253408 (stopped, reason: MEMORY_VIOLATION)

scalar registers:
            s0: d9800000            s1: 80007ffe
wave_1: pc=0x7ff77e253408 (stopped, reason: MEMORY_VIOLATION)

scalar registers:
            s0: d9800000            s1: 80007ffe
wave_2: pc=0x7ff77e2533fc (stopped, reason: MEMORY_VIOLATION)

scalar registers:
            s0: d9800000            s1: 80007ffe
"""

ROCM_AGENT_FIXTURE_720 = """\
OFFLOAD ERROR: memory access fault by GPU 4 (agent 0x55555896f810) at virtual address 0x7ffb61f02000. Reasons: Write access to a read-only page
Disassembly for function __omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486:
    code object: memory://415882#offset=0x7ff736b8c040&size=17250800
    loaded at: [0x7ff734000000-0x7ff736b0a4e8]
 => 0x7ff734dcbf3c <+7484>:    s_waitcnt vmcnt(0) lgkmcnt(0)
    0x7ff734dcbf40 <+7488>:    v_mul_f64 v[12:13], v[52:53], v[16:17]
End of disassembly.
wave_249: pc=0x7ff734dcbf3c (kernel_code_entry=0x7ff736b8c040 <__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486>, kernargs=0x7ffb61f00000) (stopped, reason: MEMORY_VIOLATION)

scalar registers:
            s0: d9800000            s1: 80007ffe
wave_250: pc=0x7ff734dcbf3c (kernel_code_entry=0x7ff736b8c040 <__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486>, kernargs=0x7ffb61f00000) (stopped, reason: MEMORY_VIOLATION)

scalar registers:
            s0: d9800000            s1: 80007ffe
"""


# Detection.


def test_recognises_each_runtime_wording():
    # AMD/HSA (CCE), libomptarget (AFAR), and NVHPC, which words it nothing like
    # the others -- matching only AMD's left 189 Phoenix faults unclassified.
    assert is_gpu_memory_fault("Memory access fault by GPU node-9 on address 0x1")
    assert is_gpu_memory_fault("OFFLOAD ERROR: memory access fault by GPU 4")
    assert is_gpu_memory_fault("Accelerator Fatal Error: ... (CUDA_ERROR_ILLEGAL_ADDRESS)")


def test_does_not_fire_on_ordinary_failures():
    # PMIX noise appears in 16% of *passing* self-hosted jobs, and NVHPC uses
    # "Accelerator Fatal Error" for out-of-memory too -- neither is a fault.
    assert not is_gpu_memory_fault("Test x: Failed to execute MFC.")
    assert not is_gpu_memory_fault("PMIX ERROR: PMIX_ERR_NO_PERMISSIONS in file dstore_base.c")
    assert not is_gpu_memory_fault("Accelerator Fatal Error: call to cuMemAlloc returned error 2: Out of memory")


def test_the_marker_survives_the_hand_off_and_rich():
    """The failure site tags the exception; classify_error reads the tag back.

    The first version wrote "[gpu-memory-fault]" while the reader searched for
    "memory access fault by gpu", so the two never matched and the feature was
    dead while seven source-inspecting tests passed. Square brackets also make
    Rich delete the marker as a style tag before it reaches the log.
    """
    import io

    from rich.console import Console

    raised = f"Test x: Failed to execute MFC {GPU_FAULT_MARKER}."
    assert is_gpu_memory_fault(raised)

    console = Console(file=io.StringIO(), force_terminal=False)
    console.print(raised)
    assert GPU_FAULT_MARKER in console.file.getvalue()


def test_a_gpu_fault_gets_its_own_failure_class():
    # Otherwise detecting it is inert: classify_error bucketed anything saying
    # "failed to execute" as a generic execution failure.
    from mfc.common import MFCException
    from mfc.test.test import classify_error

    assert classify_error(MFCException(f"Test x: Failed to execute MFC {GPU_FAULT_MARKER}.")) == "GPU memory fault"
    assert classify_error(MFCException("Test x: Failed to execute MFC.")) == "execution failed"


# What the run environment carries.


def test_only_the_agent_is_set():
    """Everything else was measured to be worse than nothing.

    CRAY_ACC_DEBUG named the wrong kernel in 81 of 102 traced faults, because
    CCE dispatches async and its trace's tail is whatever ran next.
    OFFLOAD_TRACK_ALLOCATION_TRACES and _NUM_KERNEL_LAUNCH_TRACES instrument
    every allocation and every kernel launch: on an MI210 either one alone
    turned a 5.94 s test into a >400 s timeout.
    """
    env = fault_diagnostic_env({})

    assert "CRAY_ACC_DEBUG" not in env
    assert "OFFLOAD_TRACK_ALLOCATION_TRACES" not in env
    assert "OFFLOAD_TRACK_NUM_KERNEL_LAUNCH_TRACES" not in env


def test_the_env_is_a_copy_and_keeps_what_it_was_given():
    # Cases run in worker threads; mutating a shared environment would leak
    # settings into every concurrent case.
    base = {"PATH": "/usr/bin", "HOME": "/home/x"}
    env = fault_diagnostic_env(base)

    assert env["PATH"] == "/usr/bin" and env["HOME"] == "/home/x"
    assert "HSA_TOOLS_LIB" not in base


def tmp_agent_dir() -> str:
    """A directory laid out like a ROCm install, for the gate to find."""
    import os
    import tempfile

    root = tempfile.mkdtemp()
    os.makedirs(os.path.join(root, "lib"), exist_ok=True)
    open(os.path.join(root, "lib", ROCM_DEBUG_AGENT), "w", encoding="utf-8").close()
    return root


def test_the_agent_gate_re_reads_the_environment(monkeypatch):
    """It must not be captured at import.

    On Frontier the library is on disk the whole time but only reaches
    LD_LIBRARY_PATH once `mfc.sh load` runs, so an import-time gate reports
    "absent" on the one machine this is for -- indistinguishably from Phoenix,
    where it genuinely is missing. Pinned rather than trusting the host, which
    may have a real ROCm install.
    """
    import os

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert rocm_debug_agent_path() is None
    assert "HSA_TOOLS_LIB" not in fault_diagnostic_env({})

    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())
    assert rocm_debug_agent_path() is not None
    assert fault_diagnostic_env({})["HSA_TOOLS_LIB"] == ROCM_DEBUG_AGENT

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", os.path.join(tmp_agent_dir(), "lib"))
    assert rocm_debug_agent_path() is not None


def test_a_developer_debugging_by_hand_is_left_alone(monkeypatch):
    """`mfc.sh test` and `mfc.sh bench` are not only CI entry points.

    The agent is mutually exclusive with a ROCr core dump, so enabling it behind
    someone collecting one gives them "Failed to enable debug interface" and no
    dump, caused by the harness rather than anything they did.
    """
    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())

    assert "HSA_TOOLS_LIB" in fault_diagnostic_env({})
    assert "HSA_TOOLS_LIB" not in fault_diagnostic_env({"HSA_ENABLE_DEBUG": "1"})
    assert fault_diagnostic_env({"HSA_TOOLS_LIB": "libmine.so"})["HSA_TOOLS_LIB"] == "libmine.so"


# Collapsing the agent's output.


def test_both_rocm_formats_are_recognised():
    """A parser written against one version returns '' for the other.

    That happened: 65,210 lines of real 7.2.0 output produced nothing, on the
    lane the summarizer exists to serve, with no error to explain it. Neither
    format may be fixed at the other's expense.
    """
    for name, fixture in (("6.3.1", ROCM_AGENT_FIXTURE_631), ("7.2.0", ROCM_AGENT_FIXTURE_720)):
        summary = summarize_rocm_debug_agent(fixture)
        assert summary, f"ROCm {name} agent output was not recognised"
        assert "memory access fault" in summary.lower()


def test_the_summary_keeps_the_kernel_and_the_whole_pc_histogram():
    """The two things a fixed tail cannot give.

    On the real 14,635-line report the first 80 lines are one wave's registers
    and the last 80 another's, so the kernel name is in neither. And the waves
    stop at several PCs whose modal one is a load while the fault is a write --
    quoting one PC alone names the wrong instruction.
    """
    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_720)

    assert "__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486" in summary
    assert "0x7ff734dcbf3c x2" in summary

    for field in ("faulting kernel(s): ", "faulting waves:     ", "stop PCs:           ", "--- disassembly (1 of "):
        assert field in summary, f"the summarizer no longer emits {field!r}"


def test_every_measured_symbol_form_survives():
    """Three manglings -- one per compiler, not one per offload model.

    CCE emits the same scheme for OpenACC and OpenMP offload, differing only in
    a trailing counter, while AFAR's Flang form is different again. Reading any
    two lanes suggests the offload model decides.
    """
    for symbol in (
        "s_tvd_rk$m_time_steppers_$ck_L486_6",
        "s_tvd_rk$m_time_steppers_$ck_L486_16",
        "__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486",
    ):
        summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_631.replace("s_tvd_rk$m_time_steppers_$ck_L486_6", symbol))
        assert symbol in summary and "486" in summary


def test_output_without_an_agent_report_falls_back():
    assert summarize_rocm_debug_agent("Memory access fault by GPU node-4 on address 0x1") == ""
    assert summarize_rocm_debug_agent("") == ""


def test_a_missing_agent_report_on_a_gpu_fault_is_called_out():
    # The drift above is silent by nature -- raw output where a summary should
    # be, and nothing saying why. If the agent is reachable and the failure IS a
    # GPU fault, an unrecognised report has to say so.
    import inspect

    from mfc.test.test import _handle_case

    src = inspect.getsource(_handle_case)
    assert "rocm_debug_agent_path() is not None" in src
    assert "format has changed" in src


# The other two callers.


def test_the_bench_runner_summarizes_rather_than_tailing(tmp_path):
    """bench.py ran GPU cases with no fault handling at all.

    Padded past log_tail's 60-line window on purpose: with a short fixture the
    tail contains the kernel name and this passes against the old behaviour.
    """
    from mfc.bench import bench_failure_report
    from mfc.common import log_tail

    log = tmp_path / "case.out"
    log.write_text(ROCM_AGENT_FIXTURE_720 + "\n".join(f"    v{n}: 0x0" for n in range(200)), encoding="utf-8")
    assert "_QMm_time_steppersPs_tvd_rk_l486" not in log_tail(str(log)), "fixture too short to distinguish"

    assert "_QMm_time_steppersPs_tvd_rk_l486" in bench_failure_report(str(log))

    plain = tmp_path / "plain.out"
    plain.write_text("ordinary failure\nsomething went wrong\n", encoding="utf-8")
    assert "something went wrong" in bench_failure_report(str(plain))


def test_the_shell_summarizer_reports_absence_by_exit_code(tmp_path):
    # The case-optimization script is shell, and needs to know when to fall back.
    script = pathlib.Path(__file__).resolve().parents[3] / ".github" / "scripts" / "summarize_gpu_fault.py"

    agent = tmp_path / "agent.log"
    agent.write_text(ROCM_AGENT_FIXTURE_720, encoding="utf-8")
    found = subprocess.run([sys.executable, str(script), str(agent)], capture_output=True, text=True, check=False)
    assert found.returncode == 0 and "_QMm_time_steppersPs_tvd_rk_l486" in found.stdout

    plain = tmp_path / "plain.log"
    plain.write_text("ordinary failure\n", encoding="utf-8")
    missing = subprocess.run([sys.executable, str(script), str(plain)], capture_output=True, text=True, check=False)
    assert missing.returncode == 1 and missing.stdout.strip() == ""


def test_restart_cases_carry_the_diagnostics_too():
    # They reach the GPU through run_restart, which took no env at all.
    import inspect

    from mfc.test.case import TestCase

    assert "env" in inspect.signature(TestCase.run_restart).parameters
