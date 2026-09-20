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

from pathlib import Path

from .lint_test_suite import check_no_case_specific_builds


def test_no_test_suite_case_needs_its_own_build():
    repo_root = Path(__file__).resolve().parents[2]

    errors = check_no_case_specific_builds(repo_root)

    assert not errors, "These test cases need a build no CI lane pre-builds (see toolchain/mfc/lint_test_suite.py for how to fix each kind):\n" + "\n".join(errors)
