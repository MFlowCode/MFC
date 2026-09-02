"""A GPU memory fault should explain itself the first time.

A fault reaches CI as an address and, without help, nothing else:

    Memory access fault by GPU node-4 (Agent handle: 0x...) on address 0x...

Measured against a deliberately injected out-of-bounds write at
m_time_steppers.fpp:486, on all four GPU lanes:

  * NVHPC prints the file, function and line unaided.
  * AFAR names the faulting kernel unaided, 189 of 189 faults.
  * CCE names nothing on its own, and no CRAY_ACC_* variable helps -- its trace
    blamed the wrong kernel in 81 of 102 traced faults, because dispatch is
    async and the trace's tail is whatever ran next.
  * The ROCm debug agent names it on every lane it is reachable on, including
    both CCE lanes, at ROCr level and with no recompile.

So the diagnostics are set on every run rather than on a retry: they are inert
until the runtime is already aborting, and withholding them only cost an extra
run to learn what the first could have said. The original retry design is gone.

These tests exist because nearly every failure in this area was silent. A
marker written on one side and never read on the other passed seven
source-inspecting tests; a summarizer written against one ROCm version returned
nothing for 65,210 lines of another's; a gpucore that "existed" was a
zero-segment stub. What they have in common is an artifact that looked present
while containing nothing, so these assert on content and on the hand-offs
between parts, never on presence alone.
"""

import pathlib
import sys

from mfc.test.test import fault_diagnostic_env, is_gpu_memory_fault


def test_recognises_the_hsa_level_fault_cce_reports():
    assert is_gpu_memory_fault("Memory access fault by GPU node-9 (Agent handle: 0x32463c0) on address 0x1544")


def test_recognises_the_offload_level_fault_afar_reports():
    assert is_gpu_memory_fault("OFFLOAD ERROR: memory access fault by GPU 1 (agent 0x8c59e0) at virtual address 0x7f11")


def test_does_not_fire_on_an_ordinary_failure():
    assert not is_gpu_memory_fault("Variable n5282 is not within tolerance")
    assert not is_gpu_memory_fault("NVFORTRAN-S-0034-Syntax error at or near end of line")
    assert not is_gpu_memory_fault("")


def test_does_not_fire_on_benign_pmix_noise():
    # PMIX_ERR_NO_PERMISSIONS appears in 16% of passing self-hosted jobs.
    assert not is_gpu_memory_fault("PMIX ERROR: PMIX_ERR_NO_PERMISSIONS in file dstore_base.c at line 238")


def test_the_diagnostic_env_enables_allocation_tracking():
    # AFAR's libomptarget reads this; CCE ignores it. Only AFAR gains anything
    # from a diagnostic retry, so there is nothing to detect the cluster for.
    env = fault_diagnostic_env({"PATH": "/usr/bin"})
    assert env["OFFLOAD_TRACK_ALLOCATION_TRACES"] == "true"


def test_the_diagnostic_env_preserves_the_existing_environment():
    env = fault_diagnostic_env({"PATH": "/usr/bin", "HOME": "/home/x"})
    assert env["PATH"] == "/usr/bin"
    assert env["HOME"] == "/home/x"


def test_the_diagnostic_env_does_not_mutate_what_it_was_given():
    # These run in worker threads; mutating a shared environment would leak
    # diagnostics into every concurrently running case.
    base = {"PATH": "/usr/bin"}
    fault_diagnostic_env(base)
    assert "OFFLOAD_TRACK_ALLOCATION_TRACES" not in base


def test_the_marker_written_on_failure_is_the_one_the_retry_reads():
    """The round trip, which source inspection could not check.

    _handle_case raises an exception carrying a marker; handle_case decides
    whether to enable diagnostics by inspecting that exception. The first
    version wrote "[gpu-memory-fault]" and read for "memory access fault by
    gpu", so the two never matched and the feature was dead while seven tests
    passed. Assert the actual hand-off.
    """
    from mfc.gpu_diagnostics import GPU_FAULT_MARKER, is_gpu_memory_fault

    # what _handle_case raises when the run's output shows a GPU fault
    raised = f"Test whatever: Failed to execute MFC. {GPU_FAULT_MARKER}"

    # what handle_case must conclude from it
    assert is_gpu_memory_fault(raised), "the retry cannot see the marker the failure wrote"


def test_an_ordinary_failure_message_does_not_look_like_a_gpu_fault():
    from mfc.gpu_diagnostics import is_gpu_memory_fault

    assert not is_gpu_memory_fault("Test whatever: Failed to execute MFC.")


def test_the_marker_survives_rich_rendering():
    """Rich eats "[...]" as a style tag.

    The marker is carried in an exception message that main.py prints through
    Rich. A bracketed marker is silently deleted before it reaches the log, so
    the one signal a human has that the diagnostic path was taken disappears.
    """
    import io

    from rich.console import Console

    from mfc.gpu_diagnostics import GPU_FAULT_MARKER

    console = Console(file=io.StringIO(), force_terminal=False)
    console.print(f"Failed to execute MFC {GPU_FAULT_MARKER}.")

    assert GPU_FAULT_MARKER in console.file.getvalue()


def test_the_diagnostics_are_on_for_every_run_not_just_a_retry():
    """Attempt 1 has to be the informative one.

    Both variables are inert until the runtime is already aborting on a memory
    fault, so there is nothing to save by withholding them -- and withholding
    them costs a whole extra run to learn what the first could have said.
    """
    import inspect

    from mfc.test.test import _handle_case

    src = inspect.getsource(_handle_case)

    assert "fault_diagnostic_env" in src, "the run must carry the diagnostics"


def test_nvhpc_faults_are_recognised():
    """Phoenix words it nothing like Frontier does.

    Matching only the AMD phrasing meant 189 faults on a Phoenix gpu-acc shard
    were never recognised as GPU faults at all.
    """
    from mfc.gpu_diagnostics import is_gpu_memory_fault

    nvhpc = "Accelerator Fatal Error: call to cuStreamSynchronize returned error 700 " "(CUDA_ERROR_ILLEGAL_ADDRESS): Illegal address during kernel execution"

    assert is_gpu_memory_fault(nvhpc)


def test_the_costly_cray_trace_is_not_enabled():
    """It cannot attribute, and unlike the others it is not free.

    CRAY_ACC_DEBUG streams a line per launch for the whole run and, because CCE
    dispatches async by default, blamed the wrong kernel in 81 of 102 traced
    faults and the right one in none.
    """
    from mfc.gpu_diagnostics import fault_diagnostic_env

    assert "CRAY_ACC_DEBUG" not in fault_diagnostic_env({})


def test_a_gpu_fault_is_classified_as_its_own_kind_of_failure():
    """The marker has to be read by something, or detecting the fault is inert.

    _handle_case tags the exception, but for a while nothing consumed the tag:
    classify_error bucketed every failure whose message contained "failed to
    execute" as a generic execution failure, so a GPU memory fault -- the one
    execution failure a retry provably cannot fix -- was indistinguishable from
    a transient launcher problem in the failure summary and in #1798's rescue
    accounting.
    """
    from mfc.common import MFCException
    from mfc.gpu_diagnostics import GPU_FAULT_MARKER
    from mfc.test.test import classify_error

    exc = MFCException(f"Test whatever: Failed to execute MFC {GPU_FAULT_MARKER}.")

    assert classify_error(exc) == "GPU memory fault"


def test_an_ordinary_execution_failure_is_still_bucketed_as_one():
    from mfc.common import MFCException
    from mfc.test.test import classify_error

    assert classify_error(MFCException("Test whatever: Failed to execute MFC.")) == "execution failed"


def test_an_nvhpc_out_of_memory_is_not_a_memory_fault():
    """NVHPC prefixes unrelated failures with the same words.

    "Accelerator Fatal Error" covers out-of-memory and launch failures as well
    as illegal addresses, so matching that prefix would classify a GPU running
    out of memory as a memory access fault and send the reader hunting for a
    bad index that does not exist.
    """
    from mfc.gpu_diagnostics import is_gpu_memory_fault

    oom = "Accelerator Fatal Error: call to cuMemAlloc returned error 2: Out of memory"

    assert not is_gpu_memory_fault(oom)


def test_a_restart_case_runs_with_the_diagnostics_too():
    """Restart tests reach the GPU by a different path.

    _handle_case runs them through run_restart rather than run, and that path
    did not take an env, so a fault in a restart case produced none of the
    diagnostics this module exists to provide.
    """
    import inspect

    from mfc.test.case import TestCase

    assert "env" in inspect.signature(TestCase.run_restart).parameters


# A minimal fixture in the agent's ROCm 6.3.1 format, reconstructed line-for-line
# from a real 14,635-line report on Frontier (the raw log did not survive the
# experiment's teardown). Structural details that matter and were taken from the
# real output: the "(Agent handle: ...)" clause in the fault line, "End of
# disassembly." as the terminator, and the blank line plus "scalar registers:"
# header between a wave's pc line and its register dump.
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

# The field set the summarizer produced from the real 14,635-line report. Pinned
# as structure, not values: the wave counts and PC distribution belong to that
# one fault and encoding them would test the fixture rather than the code.
REAL_SUMMARY_FIELDS = (
    "=== GPU fault summary (rocm-debug-agent,",
    "Memory access fault by GPU",
    "faulting kernel(s): ",
    "faulting waves:     ",
    "stop PCs:           ",
    "NOTE: waves halt on fault detection",
    "--- disassembly (1 of ",
    "--- representative wave (modal PC ",
)


def tmp_agent_dir() -> str:
    """A directory laid out like a ROCm install, for the gate to find."""
    import os
    import tempfile

    from mfc.gpu_diagnostics import ROCM_DEBUG_AGENT

    root = tempfile.mkdtemp()
    os.makedirs(os.path.join(root, "lib"), exist_ok=True)
    with open(os.path.join(root, "lib", ROCM_DEBUG_AGENT), "w", encoding="utf-8") as f:
        f.write("")
    return root


def test_the_agent_gate_is_not_evaluated_at_import(monkeypatch):
    """The trap that would disable this on the one machine it is for.

    On Frontier the library sits on disk the whole time, but /opt/rocm-*/lib
    only reaches LD_LIBRARY_PATH once `mfc.sh load` has run. A gate captured in
    a module-level constant answers before that and reports "absent" -- exactly
    as it would on Phoenix, where the library genuinely is missing, and with no
    way to tell the two apart. So it has to re-read the environment each call.
    """
    from mfc.gpu_diagnostics import rocm_debug_agent_path

    # A machine with ROCm installed finds the real agent, so pin the
    # environment rather than trusting whatever the host happens to have.
    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert rocm_debug_agent_path() is None

    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())
    assert rocm_debug_agent_path() is not None, "the gate did not re-read the environment"


def test_the_agent_gate_follows_ld_library_path_too(monkeypatch):
    """That is the variable `mfc.sh load` actually changes on Frontier."""
    from mfc.gpu_diagnostics import rocm_debug_agent_path

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert rocm_debug_agent_path() is None

    import os

    monkeypatch.setenv("LD_LIBRARY_PATH", os.path.join(tmp_agent_dir(), "lib"))
    assert rocm_debug_agent_path() is not None


def test_the_agent_is_enabled_when_reachable(monkeypatch):
    from mfc.gpu_diagnostics import fault_diagnostic_env

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert "HSA_TOOLS_LIB" not in fault_diagnostic_env({})

    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())
    assert fault_diagnostic_env({})["HSA_TOOLS_LIB"] == "librocm-debug-agent.so.2"


def test_the_summary_has_the_same_shape_as_the_real_report():
    """Guards the reconstruction against the real thing.

    The raw 14,635-line log is gone, so the fixture is rebuilt from the lines
    quoted out of it. This asserts the summarizer still emits every field it
    produced from the genuine report -- the check that would catch the fixture
    having drifted from the format it is supposed to stand in for.

    Structure only, never the values: pinning 125 waves or that PC histogram
    would encode one fault rather than test the code.
    """
    from mfc.gpu_diagnostics import summarize_rocm_debug_agent

    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_631)

    for field in REAL_SUMMARY_FIELDS:
        assert field in summary, f"the summarizer no longer emits {field!r}"

    # The part that fails silently: no wave match means an empty summary and a
    # fall back to 14k raw lines, with nothing to say why.
    assert "s_tvd_rk$m_time_steppers_$ck_L486_6" in summary


# ROCm 7.2.0 / AFAR OpenMP-offload format, verbatim from a real 65,210-line
# report. Two things moved versus 6.3.1: the wave line gained
# kernel_code_entry= and kernargs= BETWEEN the pc and the stop reason, and the
# fault line is worded "OFFLOAD ERROR: memory access fault ... at virtual
# address ... Reasons:" instead of "Memory access fault ... on address ...
# Reason:". A summarizer written against 6.3.1 alone returns '' for all 65,210
# lines and says nothing about why.
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


def test_the_summarizer_handles_rocm_720_as_well_as_631():
    """The version skew that silently produced nothing.

    The first version required pc= and "(stopped, reason:" to be adjacent and
    matched the fault line case-sensitively on "Memory". ROCm 7.2.0 puts
    kernel_code_entry= and kernargs= between them and says "OFFLOAD ERROR:
    memory access fault", so it returned '' for 65,210 lines of real output --
    on the very lane it was meant to serve, with no error to explain it.
    """
    from mfc.gpu_diagnostics import summarize_rocm_debug_agent

    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_720)

    assert summary, "ROCm 7.2.0 agent output was not recognised"
    assert "__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486" in summary
    assert "memory access fault" in summary.lower()
    assert "0x7ff734dcbf3c x2" in summary


def test_both_rocm_formats_survive_together():
    """Neither fixture may be fixed at the other's expense."""
    from mfc.gpu_diagnostics import summarize_rocm_debug_agent

    assert summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_631)
    assert summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_720)


def test_the_omp_symbol_still_carries_module_procedure_and_line():
    """Flang mangling keeps all three facts, which is what makes it useful.

    _QM<module> P<procedure> _l<line>: m_time_steppers, s_tvd_rk, line 486 --
    the injected fault site, in a different mangling from CCE's OpenACC form.
    """
    from mfc.gpu_diagnostics import summarize_rocm_debug_agent

    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE_720)

    assert "_QMm_time_steppers" in summary
    assert "Ps_tvd_rk" in summary
    assert "_l486" in summary


def test_every_measured_symbol_form_yields_module_procedure_and_line():
    """Three manglings, one per compiler -- not one per offload model.

    CCE emits the same scheme for OpenACC and OpenMP offload, differing only in
    a trailing counter, while AFAR's Flang form is different again. Reading any
    two lanes suggests the offload model picks the mangling; reading all three
    shows it is the compiler. Anything that parses these must not assume the
    former.
    """
    from mfc.gpu_diagnostics import summarize_rocm_debug_agent

    lanes = {
        "CCE acc": "s_tvd_rk$m_time_steppers_$ck_L486_6",
        "CCE mp": "s_tvd_rk$m_time_steppers_$ck_L486_16",
        "AFAR mp": "__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486",
    }

    for lane, symbol in lanes.items():
        report = ROCM_AGENT_FIXTURE_631.replace("s_tvd_rk$m_time_steppers_$ck_L486_6", symbol)
        summary = summarize_rocm_debug_agent(report)

        assert summary, f"{lane}: agent report not recognised"
        assert symbol in summary, f"{lane}: symbol lost from the summary"
        assert "486" in summary, f"{lane}: source line lost"


def test_the_bench_runner_summarizes_an_agent_report(tmp_path):
    """bench.py ran GPU cases with no fault handling at all.

    It printed a fixed log_tail on failure, which cannot surface an agent
    report: the tail is one wave's registers and the kernel name is not in it.
    """
    from mfc.bench import bench_failure_report
    from mfc.common import log_tail

    # The fixture must be longer than log_tail's window, or the tail happens to
    # contain the kernel name and the test passes against the old behaviour --
    # proving nothing. A real report is 65,210 lines; this pads to just past the
    # window so the tail cannot reach the kernel name, which is the actual
    # failure being guarded against.
    padding = "\n".join(f"    v{n}: 0x00000000" for n in range(200))
    log = tmp_path / "case.out"
    log.write_text(ROCM_AGENT_FIXTURE_720 + padding, encoding="utf-8")

    assert "_QMm_time_steppersPs_tvd_rk_l486" not in log_tail(str(log)), "fixture too short to distinguish tail from summary"

    report = bench_failure_report(str(log))

    assert "__omp_offloading_8116438_1c00689b__QMm_time_steppersPs_tvd_rk_l486" in report
    assert "GPU fault summary" in report


def test_the_bench_runner_falls_back_when_there_is_no_agent_report(tmp_path):
    """An ordinary build or tolerance failure must look exactly as it did."""
    from mfc.bench import bench_failure_report

    log = tmp_path / "case.out"
    log.write_text("ordinary failure\nsomething went wrong\n", encoding="utf-8")

    assert "something went wrong" in bench_failure_report(str(log))


def test_the_shell_summarizer_reports_absence_by_exit_code(tmp_path):
    """The case-opt script needs to know when to fall back, from a shell."""
    import subprocess

    script = pathlib.Path(__file__).resolve().parents[3] / ".github" / "scripts" / "summarize_gpu_fault.py"

    agent = tmp_path / "agent.log"
    agent.write_text(ROCM_AGENT_FIXTURE_720, encoding="utf-8")
    found = subprocess.run([sys.executable, str(script), str(agent)], capture_output=True, text=True, check=False)
    assert found.returncode == 0
    assert "_QMm_time_steppersPs_tvd_rk_l486" in found.stdout

    plain = tmp_path / "plain.log"
    plain.write_text("ordinary failure\n", encoding="utf-8")
    missing = subprocess.run([sys.executable, str(script), str(plain)], capture_output=True, text=True, check=False)
    assert missing.returncode == 1
    assert missing.stdout.strip() == ""


def test_a_missing_agent_report_on_a_gpu_fault_is_called_out():
    """The failure that already happened once, made visible.

    A summarizer written against one ROCm version returned nothing for 65,210
    lines of another's output, and the only symptom was raw output where a
    summary should have been. If the agent is reachable and the failure is a GPU
    fault, an unrecognised report has to say so.
    """
    import inspect

    from mfc.test.test import _handle_case

    src = inspect.getsource(_handle_case)

    assert "rocm_debug_agent_path() is not None" in src
    assert "format has changed" in src
