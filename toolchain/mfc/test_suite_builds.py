"""No test-suite case may need a build that CI never pre-builds.

The check itself lives in `lint_test_suite.py`, next to the other repo-invariant linters
and runnable on its own (`./build/venv/bin/python3 toolchain/mfc/lint_test_suite.py` --
it needs the venv, since materializing the chemistry cases imports cantera). This wrapper
is what puts it in the gates: `./mfc.sh lint` runs the toolchain's pytest suite, and both
`./mfc.sh precheck` and the Lint Gate job run `./mfc.sh lint`. Lint Gate gates the whole
Test Suite workflow, so a violation is caught before a single cluster job is queued --
the point of the exercise, given the failure it replaces took ~1.5 h to surface on one
lane and pointed at srun rather than at the missing build.
"""

import os
from pathlib import Path

import pytest

from .lint_test_suite import check_amr_flags_tested, check_no_case_specific_builds

# Both checks materialize every test-suite case (about five minutes), so the pre-commit precheck skips them
# (MFC_SKIP_SUITE_LINT=1). The Lint Gate job runs them through ./mfc.sh lint.
pytestmark = pytest.mark.skipif(os.environ.get("MFC_SKIP_SUITE_LINT", "") == "1", reason="MFC_SKIP_SUITE_LINT=1 — skipping the test-suite linters")


def test_no_test_suite_case_needs_its_own_build():
    repo_root = Path(__file__).resolve().parents[2]

    errors = check_no_case_specific_builds(repo_root)

    assert not errors, "These test cases need a build no CI lane pre-builds (see toolchain/mfc/lint_test_suite.py for how to fix each kind):\n" + "\n".join(errors)


def test_every_amr_vetted_flag_is_tested_or_recorded_as_untested():
    errors = check_amr_flags_tested(Path(__file__).resolve().parents[2])
    assert not errors, "AMR_VETTED_FLAGS and the AMR tests disagree (toolchain/mfc/case_validator.py):\n" + "\n".join(errors)
