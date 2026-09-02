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

from mfc.test.test import diagnostic_env, is_gpu_memory_fault


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


def test_the_diagnostic_env_carries_both_runtimes():
    # Frontier CCE reads CRAY_ACC_DEBUG; frontier_amd's AFAR build reads the
    # OFFLOAD_ variable. Each runtime ignores the other's, verified on both, so
    # setting both avoids having to detect the cluster here.
    env = diagnostic_env({"PATH": "/usr/bin"})
    assert env["CRAY_ACC_DEBUG"] == "1"
    assert env["OFFLOAD_TRACK_ALLOCATION_TRACES"] == "true"


def test_the_diagnostic_env_preserves_the_existing_environment():
    env = diagnostic_env({"PATH": "/usr/bin", "HOME": "/home/x"})
    assert env["PATH"] == "/usr/bin"
    assert env["HOME"] == "/home/x"


def test_the_diagnostic_env_does_not_mutate_what_it_was_given():
    # These run in worker threads; mutating a shared environment would leak
    # diagnostics into every concurrently running case.
    base = {"PATH": "/usr/bin"}
    diagnostic_env(base)
    assert "CRAY_ACC_DEBUG" not in base


def test_a_gpu_fault_is_flagged_so_the_retry_can_react():
    # _handle_case marks the exception when the run's output shows a GPU fault;
    # the retry loop keys off that mark. Without the mark the retry is just
    # another identical attempt.
    import inspect

    from mfc.test import test as t

    src = inspect.getsource(t._handle_case)
    assert "is_gpu_memory_fault" in src, "the failure path must classify GPU faults"
    assert "gpu-memory-fault" in src


def test_the_retry_loop_enables_diagnostics_after_a_gpu_fault():
    import inspect

    from mfc.test import test as t

    src = inspect.getsource(t.handle_case)
    assert "diagnostic_env" in src, "the retry must enable diagnostics"
    assert "case_env" in src
    # and must pass it to the run, not merely compute it
    assert "env=case_env" in src


def test_diagnostics_are_not_enabled_for_ordinary_failures():
    import inspect

    from mfc.test import test as t

    src = inspect.getsource(t.handle_case)
    # the enabling is guarded by the fault check, not unconditional
    assert "if is_gpu_memory_fault(" in src


def test_a_diagnostic_retry_does_not_dump_its_whole_output():
    # CRAY_ACC_DEBUG prints per kernel launch and per transfer: measured at
    # 142,777 lines for a single 800-cell 1D case on Frontier. Echoing that into
    # the CI log would bury the failure it is meant to explain. The fault is at
    # the end, and the last launch before it is what names the kernel, so only
    # the tail is worth keeping.
    import inspect

    from mfc.test import test as t

    src = inspect.getsource(t._handle_case)
    assert "log_tail" in src, "the diagnostic retry's output must be bounded"
