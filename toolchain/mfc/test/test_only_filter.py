"""Tests for the labels `--only` matches a case against.

The "Chemistry" label used to exist only in whatever trace a case was given by
hand. That made it a name, not a property, and Frontier AMD's GPU lane depends
on it being a property: there the chemistry binaries are built by a separate
SLURM job selected with `-o Chemistry` (.github/workflows/common/build.sh), and
the test job runs `--no-build`, so a chemistry case the filter misses is never
compiled and dies at run time with

    execve(): build/install/gpu-mp-chem-<hash>/bin/syscheck: No such file

Auto-registered Examples never carry the label by hand, and six of them were
already unlabeled. They survived only because every one of them happens to use
h2o2.yaml, which the labeled cases build anyway. The first Example to bring its
own mechanism (2D -> Example -> ibm_reacting_surface, carbon_gasphase_reduced_gri11)
is the one that turned that coincidence into a red lane.

This is what makes its mechanism get built on that lane, so the live suite does
depend on it. The six survivors depend on it too, in the weaker sense that the
next Example to bring its own mechanism would otherwise fail the same silent
way, two hours into a Frontier job.
"""

import types

from mfc.test.test import case_filter_labels


class _FakeCase:
    """A stand-in for TestCaseBuilder, which is what __filter actually holds.

    Params live behind to_case(), not on the builder, and reaching them is the
    expensive step the include_chemistry flag exists to avoid -- so the stub
    counts the calls rather than exposing .params directly.
    """

    def __init__(self, trace, params=None, uuid="DEADBEEF"):
        self.trace = trace
        self._params = params or {}
        self._uuid = uuid
        self.to_case_calls = 0

    def get_uuid(self):
        return self._uuid

    def to_case(self):
        self.to_case_calls += 1
        return types.SimpleNamespace(params=self._params)


def test_trace_elements_and_uuid_are_labels():
    labels = case_filter_labels(_FakeCase("1D -> Bubbles -> QBMM", uuid="CE9DBA3F"))

    assert labels == {"1D", "Bubbles", "QBMM", "CE9DBA3F"}


def test_a_chemistry_case_is_labeled_chemistry_without_saying_so_in_its_trace():
    case = _FakeCase("2D -> Example -> ibm_reacting_surface", {"chemistry": "T"}, "F52F0D4C")

    assert "Chemistry" in case_filter_labels(case, include_chemistry=True)


def test_a_non_chemistry_case_is_not_labeled_chemistry():
    case = _FakeCase("2D -> Example -> rayleigh_taylor")

    assert "Chemistry" not in case_filter_labels(case, include_chemistry=True)


def test_chemistry_off_is_not_labeled_chemistry():
    case = _FakeCase("2D -> Example -> rayleigh_taylor", {"chemistry": "F"})

    assert "Chemistry" not in case_filter_labels(case, include_chemistry=True)


def test_a_case_that_already_says_chemistry_is_not_double_counted():
    case = _FakeCase("1D -> Chemistry -> Perfect Reactor", {"chemistry": "T"}, "5DCF300C")

    assert case_filter_labels(case, include_chemistry=True) == {"1D", "Chemistry", "Perfect Reactor", "5DCF300C"}


def test_the_case_is_not_built_when_chemistry_was_not_asked_for():
    case = _FakeCase("2D -> Example -> ibm_reacting_surface", {"chemistry": "T"}, "F52F0D4C")

    case_filter_labels(case)

    assert case.to_case_calls == 0
