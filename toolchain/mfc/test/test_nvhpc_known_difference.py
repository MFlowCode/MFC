"""The NVHPC < 25.11 GPU skip (#2181) must hit only GPU builds with an old nvfortran."""

import pytest

from mfc.test.test import nvhpc_gpu_known_difference

OLD = "nvfortran 25.5-0 64-bit target on x86-64 Linux -tp haswell"


@pytest.mark.parametrize(
    "uuid, gpu, version, skipped",
    [
        ("124D7972", "acc", OLD, True),
        ("124D7972", "mp", OLD, True),
        ("124D7972", "acc", "nvfortran 24.5-0 64-bit target", True),
        ("124D7972", "acc", "nvfortran 25.11-0 64-bit target", False),
        ("124D7972", "acc", "nvfortran 26.1-0 64-bit target", False),
        ("124D7972", "acc", "nvfortran garbage", False),
        ("124D7972", "acc", "", False),
        ("124D7972", "no", OLD, False),
        ("124D7972", "acc", "GNU Fortran (GCC) 13.2.0", False),
        ("124D7972", "acc", "ftn 25.5", False),
        ("DEADBEEF", "acc", OLD, False),
    ],
)
def test_skip_decision(uuid, gpu, version, skipped):
    assert nvhpc_gpu_known_difference(uuid, gpu, version) is skipped
