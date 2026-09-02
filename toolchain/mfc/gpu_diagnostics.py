"""Offload-runtime diagnostics for GPU memory faults.

Shared by the test harness, the benchmark runner and the case-optimization CI
script -- all three run GPU cases and all three need the same answer when one
faults. Kept out of test/ because bench.py depending on the test module to
explain a crash would be the wrong way round.
"""

import collections
import os
import re
import typing

# The marker _handle_case attaches to the exception it raises, so that
# classify_error can tell a GPU memory fault from any other execution failure.
# Two constraints, both learned the hard way:
#
#   * It must be one of the signatures below verbatim, because classify_error
#     recognises it by running the same matcher over the message. An earlier
#     version wrote "[gpu-memory-fault]" while the reader searched for "memory
#     access fault by gpu", so the two never matched and the feature was dead
#     while seven source-inspecting tests passed.
#   * No square brackets. main.py renders these messages through Rich, which
#     parses "[...]" as a style tag and deletes it -- which is why a CI log
#     showed a bare "Failed to execute MFC. " with the marker missing.
GPU_FAULT_MARKER = "(memory access fault by GPU)"

GPU_FAULT_SIGNATURES = (
    # AMD/HSA -- Frontier, both CCE and AFAR builds.
    "memory access fault by gpu",
    "offload error: memory access fault",
    # NVHPC -- Phoenix. Worded nothing like the AMD ones, so matching only the
    # above meant 189 faults on a Phoenix gpu-acc shard were never recognised.
    # Only the specific error: NVHPC prefixes unrelated failures with
    # "Accelerator Fatal Error" too, including "call to cuMemAlloc returned
    # error 2: Out of memory", which is not a memory fault and must not be
    # classified as one.
    "cuda_error_illegal_address",
)


def is_gpu_memory_fault(text: str) -> bool:
    """Whether output shows a GPU memory fault, as opposed to any other failure.

    Deliberately narrow. PMIX_ERR_NO_PERMISSIONS and friends appear in 16% of
    *passing* self-hosted jobs, so anything broader would fire constantly.
    """
    lowered = (text or "").lower()

    return any(sig in lowered for sig in GPU_FAULT_SIGNATURES)


def fault_diagnostic_env(base: dict) -> dict:
    """`base` plus the one diagnostic cheap enough to leave on.

    Set on every run rather than on a retry: the ROCm debug agent writes nothing
    until the runtime is already aborting on a memory fault, so a first failure
    is explained without spending a second run to reproduce it.

    Two variables that used to live here were removed after measurement -- see
    below. What is left is the agent, which is what names the faulting kernel.
    """
    env = dict(base)

    # OFFLOAD_TRACK_ALLOCATION_TRACES and OFFLOAD_TRACK_NUM_KERNEL_LAUNCH_TRACES
    # were set here and had to be removed. They instrument every allocation and
    # every kernel launch, so a healthy run pays continuously: on an MI210 with
    # amdflang, test AFBCBDFA takes 5.94 s with neither and times out past 400 s
    # with either one alone -- enough, against the 1-hour test timeout, to turn
    # a fault into a timeout and hide what they exist to explain.
    #
    # They looked free only because the A/B that cleared them ran on CCE, whose
    # offload runtime ignores libomptarget variables entirely. What they added,
    # one line on whether the address was ever a real allocation, the agent's
    # kernel name and source line subsume.

    # Skipped when the caller is already debugging by hand: they chose a tool,
    # or they set HSA_ENABLE_DEBUG to collect a GPU core dump, which the agent
    # is mutually exclusive with. Loading it anyway would leave them with
    # "Failed to enable debug interface" and no dump. An attached rocgdb trips
    # the same path.
    #
    # Cost, Frontier CCE --gpu mp over four interleaved pairs: no effect
    # detected on a healthy run (resolution ~0.8%), no output at all until
    # something faults, and +0.387 s on a faulting run -- 0.011% of the test
    # timeout. ~3-4% on an MI210 (n=2). It does not supersede libomptarget's own
    # report on the AFAR lane; it is exclusive with ROCr core dumps only.
    if "HSA_TOOLS_LIB" not in env and not env.get("HSA_ENABLE_DEBUG") and rocm_debug_agent_path() is not None:
        env["HSA_TOOLS_LIB"] = ROCM_DEBUG_AGENT

    return env


ROCM_DEBUG_AGENT = "librocm-debug-agent.so.2"


def rocm_debug_agent_path() -> typing.Optional[str]:
    """Where the ROCm debug agent lives, or None if it is not reachable.

    MUST be evaluated at call time, never cached at import. On Frontier the
    library is on disk the whole time, but /opt/rocm-*/lib only reaches
    LD_LIBRARY_PATH once `mfc.sh load` runs. A gate evaluated at import decides
    "absent" on the one machine this exists for, and does it indistinguishably
    from the Phoenix case where the library really is missing.

    Probes for the file rather than dlopen'ing it: ctypes.CDLL would load a
    debug agent into the test harness's own process to answer a question about
    the subprocess.
    """
    rocm_path = os.environ.get("ROCM_PATH", "")
    search = [os.path.join(rocm_path, "lib")] if rocm_path else []
    search += os.environ.get("LD_LIBRARY_PATH", "").split(os.pathsep)

    for directory in search:
        if directory and os.path.isfile(os.path.join(directory, ROCM_DEBUG_AGENT)):
            return os.path.join(directory, ROCM_DEBUG_AGENT)

    return None


def summarize_rocm_debug_agent(out: str, max_disasm: int = 14) -> str:
    """Collapse librocm-debug-agent output to a bounded, informative summary.

    The agent repeats an identical disassembly block and a 115-line register
    dump per faulting wave -- 125 waves produced 14,635 lines on a 49x39 case.
    Only the kernel name, fault reason, stop-PC distribution and one
    representative wave carry information; the rest is duplicated.

    A fixed tail cannot substitute. Measured on that log: the first 80 lines are
    one wave's registers and the last 80 are another's, and the kernel name --
    the entire point -- appears in neither. The stop-PC histogram is kept
    because the waves halted at four distinct PCs whose modal one is a load
    while the injected fault is a write, so quoting a single PC without the
    distribution hands the reader the wrong instruction.

    Returns '' when there is no agent report, so callers fall back to the raw
    output.

    The format is NOT stable across ROCm versions, and the failure is silent --
    no wave match means an empty summary and a fallback to tens of thousands of
    raw lines, with nothing saying why. Measured between two versions:

        6.3.1  wave_124: pc=0x7ff77e253408 (stopped, reason: MEMORY_VIOLATION)
        7.2.0  wave_250: pc=0x7ff734dcbf3c (kernel_code_entry=0x... <...>,
                         kernargs=0x...) (stopped, reason: MEMORY_VIOLATION)

        6.3.1  Memory access fault by GPU node-4 (Agent handle: ...) on address
        7.2.0  OFFLOAD ERROR: memory access fault by GPU 4 (agent ...) at ...

    An earlier version required pc= and "(stopped, reason:" to be adjacent and
    matched the fault line case-sensitively on "Memory". It returned nothing at
    all for 65,210 lines of real 7.2.0 output. Hence the tolerant separator, and
    reusing is_gpu_memory_fault rather than hardcoding one version's wording.
    Both formats are pinned by fixtures below.

    Validated against three real reports, not one:

        CCE  acc  ROCm 6.3.1   14,635 lines -> 37
        CCE  mp   ROCm 6.3.1   13,826 lines -> 35
        AFAR mp   ROCm 7.2.0   65,210 lines -> 36

    and output from a run with no agent loaded still yields '', so the fallback
    is intact. The stop-PC histogram earns its place most on the CCE OpenMP
    lane, which halts at seven distinct PCs (62/21/19/10/10/2/1) against four
    for CCE OpenACC and one for AFAR: quoting a single PC would be wrong there
    six times in seven.
    """
    waves = re.findall(r"^wave_\d+: pc=(0x[0-9a-f]+).*?\(stopped, reason: (\w+)\)", out, re.M)
    if not waves:
        return ""

    lines = out.splitlines()
    fault = next((line for line in lines if is_gpu_memory_fault(line)), None)
    kernels = sorted({m.group(1) for m in re.finditer(r"^Disassembly for function (.+):$", out, re.M)})
    pcs = collections.Counter(pc for pc, _ in waves)
    reasons = collections.Counter(reason for _, reason in waves)

    summary = [f"=== GPU fault summary (rocm-debug-agent, {len(lines)} lines collapsed) ==="]
    if fault:
        summary.append(fault.strip())
    summary.append("faulting kernel(s): " + (", ".join(kernels) or "<none reported>"))
    summary.append(f"faulting waves:     {len(waves)}  [" + ", ".join(f"{r} x{n}" for r, n in reasons.most_common()) + "]")
    summary.append("stop PCs:           " + ", ".join(f"{pc} x{n}" for pc, n in pcs.most_common()))
    summary.append("NOTE: waves halt on fault detection, so the PC is near -- not necessarily at -- the offending instruction.")

    disasm_starts = [n for n, line in enumerate(lines) if line.startswith("Disassembly for function")]
    if disasm_starts:
        start = disasm_starts[0]
        end = next((n for n, line in enumerate(lines[start:], start) if line.startswith("End of disassembly")), start + max_disasm)
        summary += ["", f"--- disassembly (1 of {len(disasm_starts)} identical blocks) ---"]
        summary += lines[start : min(end + 1, start + max_disasm)]

    modal_pc = pcs.most_common(1)[0][0]
    start = next((n for n, line in enumerate(lines) if re.match(r"^wave_\d+: pc=" + re.escape(modal_pc) + r"(?![0-9a-f])", line)), None)
    if start is not None:
        summary += ["", f"--- representative wave (modal PC {modal_pc}, {pcs[modal_pc]} of {len(waves)} waves) ---"]
        summary += lines[start : start + max_disasm]
        summary.append(f"    ... (registers for {len(waves) - 1} further waves suppressed)")

    return "\n".join(summary)
