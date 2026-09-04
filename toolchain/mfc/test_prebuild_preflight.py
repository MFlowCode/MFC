"""The case-optimization pre-build must NOT probe the node.

Every other SLURM job MFC submits runs binaries built for the device it asked
for. This one does not: test.yml submits it as a *cpu* allocation, because it is
a --dry-run that only builds, while the binaries it produces are GPU builds.

A syscheck built with --gpu asserts omp_get_num_devices() > 0 (or the OpenACC
equivalent) and so exits non-zero on a node that has no GPU by design. Probing
here reads that as a bad node. It did exactly that in CI: three healthy Phoenix
nodes were condemned and two added to --exclude before the wrapper gave up,
which is the precise false positive the whole guard exists to avoid.

The GPU allocation that actually runs these cases is probed instead, in
run_case_optimization.sh, where the allocation and the binary agree.
"""

from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / ".github" / "scripts"

# Scripts that run binaries built for the device their allocation requested.
PROBES = ["run_case_optimization.sh"]
# Scripts whose allocation device deliberately differs from what they build.
DOES_NOT_PROBE = ["prebuild-case-optimization.sh"]


@pytest.mark.parametrize("script", DOES_NOT_PROBE)
def test_a_build_only_allocation_does_not_judge_its_node(script):
    body = (SCRIPTS / script).read_text()
    assert "preflight.sh" not in body, f"{script} runs on a cpu allocation but builds GPU binaries, so a probe " "there always fails and excludes a healthy node"


@pytest.mark.parametrize("script", PROBES)
def test_the_allocation_that_runs_the_cases_does_judge_its_node(script):
    assert "preflight.sh" in (SCRIPTS / script).read_text()


def test_the_reason_is_recorded_where_someone_would_re_add_it():
    # The next person to notice case-opt's pre-build is unprobed should find the
    # reason in the file rather than rediscovering it through a red CI run.
    body = (SCRIPTS / "prebuild-case-optimization.sh").read_text()
    assert "cpu" in body and "probe" in body.lower()
