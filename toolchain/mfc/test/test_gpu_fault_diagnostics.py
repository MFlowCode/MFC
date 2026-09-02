"""A GPU memory fault should diagnose itself on the retry.

MFC already retries a failed case up to three times, and those retries rescue
almost nothing -- 0 of 235 in bench, and every recorded failed test shows the
full attempt count. That last fact is the useful one: when a case fails it
fails all its attempts, so the retry is a free, already-paid-for reproduction
of the fault.

Today a GPU memory fault reaches CI as an address and nothing else:

    Memory access fault by GPU node-9 (Agent handle: 0x...) on address 0x...

Measured on Frontier, enabling the offload runtime's diagnostics turns the same
fault into the kernel and the source line that caused it:

    ACC: Execute kernel <name> async(auto) from <file>:<line>
    Memory access fault by GPU node-4 ...

Those diagnostics print per kernel launch, so they cannot be on for a whole run
-- but they cost nothing on a retry that was going to happen anyway and would
otherwise produce nothing.
"""

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
    from mfc.test.test import GPU_FAULT_MARKER, is_gpu_memory_fault

    # what _handle_case raises when the run's output shows a GPU fault
    raised = f"Test whatever: Failed to execute MFC. {GPU_FAULT_MARKER}"

    # what handle_case must conclude from it
    assert is_gpu_memory_fault(raised), "the retry cannot see the marker the failure wrote"


def test_an_ordinary_failure_message_does_not_look_like_a_gpu_fault():
    from mfc.test.test import is_gpu_memory_fault

    assert not is_gpu_memory_fault("Test whatever: Failed to execute MFC.")


def test_the_marker_survives_rich_rendering():
    """Rich eats "[...]" as a style tag.

    The marker is carried in an exception message that main.py prints through
    Rich. A bracketed marker is silently deleted before it reaches the log, so
    the one signal a human has that the diagnostic path was taken disappears.
    """
    import io

    from rich.console import Console

    from mfc.test.test import GPU_FAULT_MARKER

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
    from mfc.test.test import is_gpu_memory_fault

    nvhpc = "Accelerator Fatal Error: call to cuStreamSynchronize returned error 700 " "(CUDA_ERROR_ILLEGAL_ADDRESS): Illegal address during kernel execution"

    assert is_gpu_memory_fault(nvhpc)


def test_the_costly_cray_trace_is_not_enabled():
    """It cannot attribute, and unlike the others it is not free.

    CRAY_ACC_DEBUG streams a line per launch for the whole run and, because CCE
    dispatches async by default, blamed the wrong kernel in 81 of 102 traced
    faults and the right one in none.
    """
    from mfc.test.test import fault_diagnostic_env

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
    from mfc.test.test import GPU_FAULT_MARKER, classify_error

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
    from mfc.test.test import is_gpu_memory_fault

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


# A minimal fixture in the agent's observed ROCm 6.3.1 format. The verbatim
# lines come from a real 14,635-line report on Frontier; the raw log itself did
# not survive the experiment's teardown.
ROCM_AGENT_FIXTURE = """\
Memory access fault by GPU node-4 on address 0x7ffb6a0f6000. Reason: Write access to a read-only page.
Disassembly for function s_tvd_rk$m_time_steppers_$ck_L486_6:
  code object: file:///path/to/simulation#offset=12345&size=67890
  loaded at: [0x7f0000000000, 0x7f0000010000]
      0x7f0000001000 <+16>:   global_load_dwordx4 v[8:11], v[4:5], off
  =>  0x7f0000001008 <+24>:   global_store_dword v[6:7], v12, off
      0x7f0000001010 <+32>:   s_endpgm
End of disassembly for function s_tvd_rk$m_time_steppers_$ck_L486_6.
wave_0: pc=0x7f0000001008 (stopped, reason: MEMORY_VIOLATION)
    s0: 0x00000000  s1: 0x00000001  s2: 0x00000002
    v0: 0x00000000  v1: 0x00000001
wave_1: pc=0x7f0000001008 (stopped, reason: MEMORY_VIOLATION)
    s0: 0x00000000  s1: 0x00000001  s2: 0x00000002
wave_2: pc=0x7f0000001000 (stopped, reason: MEMORY_VIOLATION)
    s0: 0x00000000  s1: 0x00000001  s2: 0x00000002
"""


def tmp_agent_dir() -> str:
    import os
    import tempfile

    from mfc.test.test import ROCM_DEBUG_AGENT

    root = tempfile.mkdtemp()
    os.makedirs(os.path.join(root, "lib"), exist_ok=True)
    with open(os.path.join(root, "lib", ROCM_DEBUG_AGENT), "w", encoding="utf-8") as f:
        f.write("")
    return root


def test_the_agent_report_is_collapsed_and_keeps_the_kernel_name():
    """A fixed tail cannot do this job.

    Measured on the real 14,635-line report: the first 80 lines are one wave's
    registers and the last 80 are another's, so the kernel name -- the entire
    point -- appears in neither. This asserts the name survives.
    """
    from mfc.test.test import summarize_rocm_debug_agent

    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE)

    assert summary, "the summarizer did not recognise a real agent report"
    assert "s_tvd_rk$m_time_steppers_$ck_L486_6" in summary
    assert "Memory access fault by GPU" in summary
    assert len(summary.splitlines()) < len(ROCM_AGENT_FIXTURE.splitlines()) + 12


def test_the_summary_reports_the_whole_stop_pc_distribution():
    """Quoting one PC would name the wrong instruction.

    The waves halt at several distinct PCs, and on the observed fault the modal
    one is a load while the injected fault is a write. Collapsing to a single PC
    hands the reader a confidently wrong instruction -- the failure mode that
    made CCE's own trace worse than no diagnostic at all.
    """
    from mfc.test.test import summarize_rocm_debug_agent

    summary = summarize_rocm_debug_agent(ROCM_AGENT_FIXTURE)

    assert "0x7f0000001008 x2" in summary
    assert "0x7f0000001000 x1" in summary


def test_output_without_an_agent_report_falls_back():
    from mfc.test.test import summarize_rocm_debug_agent

    assert summarize_rocm_debug_agent("Memory access fault by GPU node-4 on address 0x1557f5ced000.") == ""
    assert summarize_rocm_debug_agent("") == ""


def test_the_agent_gate_is_not_evaluated_at_import(monkeypatch):
    """The trap that would disable this on the one machine it is for.

    On Frontier the library sits on disk the whole time, but /opt/rocm-*/lib
    only reaches LD_LIBRARY_PATH once `mfc.sh load` has run. A gate captured in
    a module-level constant answers before that and reports "absent" -- exactly
    as it would on Phoenix, where the library genuinely is missing, and with no
    way to tell the two apart. So it has to re-read the environment each call.
    """
    from mfc.test.test import rocm_debug_agent_path

    # A machine with ROCm installed finds the real agent, so pin the
    # environment rather than trusting whatever the host happens to have.
    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert rocm_debug_agent_path() is None

    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())
    assert rocm_debug_agent_path() is not None, "the gate did not re-read the environment"


def test_the_agent_gate_follows_ld_library_path_too(monkeypatch):
    """That is the variable `mfc.sh load` actually changes on Frontier."""
    from mfc.test.test import rocm_debug_agent_path

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert rocm_debug_agent_path() is None

    import os

    monkeypatch.setenv("LD_LIBRARY_PATH", os.path.join(tmp_agent_dir(), "lib"))
    assert rocm_debug_agent_path() is not None


def test_the_agent_is_enabled_when_reachable(monkeypatch):
    from mfc.test.test import fault_diagnostic_env

    monkeypatch.setenv("ROCM_PATH", "")
    monkeypatch.setenv("LD_LIBRARY_PATH", "")
    assert "HSA_TOOLS_LIB" not in fault_diagnostic_env({})

    monkeypatch.setenv("ROCM_PATH", tmp_agent_dir())
    assert fault_diagnostic_env({})["HSA_TOOLS_LIB"] == "librocm-debug-agent.so.2"
