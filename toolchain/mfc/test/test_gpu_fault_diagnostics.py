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
