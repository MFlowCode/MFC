"""Unit tests for .github/scripts/ci-outage.sh.

When an environment-wide outage hits a cluster -- pypi.org unreachable from a
Frontier login node, say -- every queued CI job rediscovers it independently.
On 2026-08-28 that cost 17 Frontier jobs about 33 minutes each to learn the same
fact. The circuit breaker lets the first job that notices record it on the shared
filesystem so later jobs skip immediately instead of queueing SLURM to fail.

A breaker that cannot reset is worse than none, so the time-to-live behaviour is
pinned here as tightly as the tripping behaviour.
"""

import os
import subprocess
import time
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "ci-outage.sh"

CLEAR = 0
TRIPPED = 1


@pytest.fixture
def state_dir(tmp_path):
    return tmp_path / "state"


def run(state_dir, *args, ttl=None):
    env = {**os.environ, "MFC_CI_STATE_DIR": str(state_dir)}
    if ttl is not None:
        env["MFC_CI_OUTAGE_TTL_SECONDS"] = str(ttl)
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def test_check_reports_clear_when_nothing_has_been_recorded(state_dir):
    assert run(state_dir, "check", "phoenix").returncode == CLEAR


def test_check_trips_after_an_outage_is_marked(state_dir):
    run(state_dir, "mark", "phoenix", "pypi unreachable")
    assert run(state_dir, "check", "phoenix").returncode == TRIPPED


def test_check_reports_the_recorded_reason_so_the_log_explains_the_skip(state_dir):
    run(state_dir, "mark", "phoenix", "pypi unreachable")
    result = run(state_dir, "check", "phoenix")
    assert "pypi unreachable" in result.stdout + result.stderr


def test_an_outage_on_one_cluster_leaves_the_other_alone(state_dir):
    run(state_dir, "mark", "frontier", "pypi unreachable")
    assert run(state_dir, "check", "phoenix").returncode == CLEAR


def test_an_outage_older_than_the_ttl_stops_tripping(state_dir):
    run(state_dir, "mark", "phoenix", "pypi unreachable", ttl=1)
    time.sleep(2)
    assert run(state_dir, "check", "phoenix", ttl=1).returncode == CLEAR


def test_clear_resets_the_breaker(state_dir):
    run(state_dir, "mark", "phoenix", "pypi unreachable")
    run(state_dir, "clear", "phoenix")
    assert run(state_dir, "check", "phoenix").returncode == CLEAR


def test_marking_works_when_the_state_directory_does_not_exist_yet(state_dir):
    # The first job on a fresh runner must not fail just because nothing has
    # created the directory.
    assert not state_dir.exists()
    assert run(state_dir, "mark", "phoenix", "pypi unreachable").returncode == 0


def test_a_corrupt_marker_is_treated_as_clear_rather_than_wedging_ci(state_dir):
    # A truncated or garbled marker must never become an un-clearable breaker.
    run(state_dir, "mark", "phoenix", "pypi unreachable")
    marker = next(state_dir.glob("*phoenix*"))
    marker.write_text("not-a-timestamp\n")
    assert run(state_dir, "check", "phoenix").returncode == CLEAR


def test_check_names_the_marker_file_so_a_human_can_clear_it(state_dir):
    run(state_dir, "mark", "phoenix", "pypi unreachable")
    result = run(state_dir, "check", "phoenix")
    assert str(state_dir) in result.stdout + result.stderr


def test_an_unknown_subcommand_fails_with_a_usage_error(state_dir):
    # Distinct from both CLEAR and TRIPPED so a typo in a workflow can never be
    # silently read as "no outage".
    result = run(state_dir, "frobnicate", "phoenix")
    assert result.returncode == 2
    assert "usage" in (result.stdout + result.stderr).lower()
