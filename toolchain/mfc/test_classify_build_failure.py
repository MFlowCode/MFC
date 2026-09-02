"""Tests for .github/scripts/classify-build-failure.sh.

The classifier decides whether a failed build was a cluster-wide dependency
outage. Getting it wrong in either direction is costly: a missed outage means
every matrix job rediscovers it, and a false positive halts CI on an ordinary
compile error.
"""

import os
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "classify-build-failure.sh"

OUTAGE = 78
NOT_AN_OUTAGE = 0


@pytest.fixture
def classify(tmp_path):
    def _run(log_text):
        log = tmp_path / "build.log"
        log.write_text(log_text)
        return subprocess.run(
            ["bash", str(SCRIPT), str(log), "frontier"],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "MFC_CI_STATE_DIR": str(tmp_path / "state")},
        ).returncode

    return _run


@pytest.mark.parametrize(
    "line",
    [
        "  |-> Failed to fetch: `https://pypi.org/simple/hatch-vcs/`",  # uv, backticked
        "Failed to fetch: https://pypi.org/simple/hatch-vcs/",  # plain, no wrapper
        "  Failed to fetch:   https://pypi.org/simple/build/",  # padded
        "mfc: ERROR > (venv) Installation failed.",
        "mfc: WARNING > (venv) uv install failed; clearing the uv cache",
    ],
)
def test_dependency_outages_are_recognised_in_every_formatting_variant(classify, line):
    assert classify(line + "\n") == OUTAGE


@pytest.mark.parametrize(
    "line",
    [
        "NVFORTRAN-S-0034-Syntax error at or near end of line",
        "ftn-2116 ftn: INTERNAL",
        "clang: error: ld.lld command failed with exit code 1",
        "CMake Error: could not find HDF5",
        "Failed to fetch: https://example.com/not-pypi",
    ],
)
def test_ordinary_build_failures_are_not_outages(classify, line):
    assert classify(line + "\n") == NOT_AN_OUTAGE


def test_an_absent_log_is_never_called_an_outage(tmp_path):
    result = subprocess.run(
        ["bash", str(SCRIPT), str(tmp_path / "nope.log"), "frontier"],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "MFC_CI_STATE_DIR": str(tmp_path / "state")},
    )
    assert result.returncode == NOT_AN_OUTAGE
