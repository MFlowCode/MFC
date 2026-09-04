"""Unit tests for src/syscheck/syscheck.fpp.

syscheck is CI's proof that a compute node can actually run MFC: it is the first
binary launched inside a SLURM allocation, and its exit code is what tells the
CI wrapper whether the node is healthy or should be excluded and requeued.

That only works if syscheck genuinely exercises the device. A path that merely
queries the runtime ("how many GPUs are there?") passes on a node whose GPU is
dead, which is exactly what happened on Phoenix in Aug 2026: on the OpenACC
build syscheck died at cuCtxCreate with CUDA_ERROR_ECC_UNCORRECTABLE, while on
the OpenMP build it reported PASSED up to 111 times on the same broken nodes
before the solver fell over.

These tests pin the offload down so the OpenMP path cannot silently regress to
a host-side query again.
"""

import re
from pathlib import Path

import fypp

SYSCHECK_FPP = Path(__file__).resolve().parents[2] / "src" / "syscheck" / "syscheck.fpp"


def _expand_fypp() -> str:
    return fypp.Fypp(fypp.FyppOptions()).process_text(SYSCHECK_FPP.read_text())


def _preprocess(text: str, defines: set) -> str:
    """Resolve the #ifdef/#else/#endif nesting the way a Fortran preprocessor would.

    syscheck.fpp only uses plain #ifdef, so a small evaluator is enough and keeps
    the test free of a dependency on an external cpp.
    """
    out, stack = [], []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#ifdef "):
            stack.append(stripped.split(None, 1)[1].strip() in defines)
        elif stripped == "#else":
            stack[-1] = not stack[-1]
        elif stripped == "#endif":
            stack.pop()
        elif all(stack):
            out.append(line)
    assert not stack, "unbalanced #ifdef in syscheck.fpp"
    return "\n".join(out)


def _variant(*defines: str) -> str:
    """The source as the compiler sees it for one GPU interface."""
    return _preprocess(_expand_fypp(), {"MFC_MPI", *defines})


def _omp() -> str:
    return _variant("MFC_OpenMP")


def _acc() -> str:
    return _variant("MFC_OpenACC")


def test_openmp_path_launches_a_kernel_on_the_device():
    # A host-side omp_get_num_devices() query cannot fail on a node whose GPU is
    # unusable; only an actual target region creates a context and runs code there.
    assert re.search(r"!\$omp\s+target\s+teams\s+distribute\s+parallel\s+do", _omp())


def test_openmp_path_copies_the_result_back_from_the_device():
    assert re.search(r"!\$omp\s+target\s+update\s+from\s*\(\s*arr", _omp())


def test_openmp_path_checks_the_values_returned_from_the_device():
    # Copying back without inspecting the values would still pass on a device
    # that returns garbage.
    assert re.search(r"call\s+assert\s*\(.*\barr\b", _omp())


def test_openacc_path_checks_the_values_returned_from_the_device():
    assert re.search(r"call\s+assert\s*\(.*\barr\b", _acc())


def test_openmp_device_index_is_reduced_modulo_the_device_count():
    # mod(rank, nRanks) is always < nRanks and ignores how many GPUs exist, so a
    # 4-rank job on a 2-GPU node selects devices 2 and 3, which are not there.
    assert re.search(r"omp_set_default_device\s*\(\s*mod\s*\(\s*rank\s*,\s*num_devices", _omp())


def test_openacc_device_index_is_reduced_modulo_the_device_count():
    assert re.search(r"acc_set_device_num\s*\(\s*mod\s*\(\s*rank\s*,\s*num_devices", _acc())
