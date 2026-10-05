"""Tests for the --no-build preflight and --no-chemistry filter.

Under --no-build, build() is a no-op, so a case whose build variant was never
compiled (chemistry is the common one: ./mfc.sh build does not produce it) used
to fail only when it ran, typically at the end of a full suite.
"""

import sys

from mfc.test.test import _drop_chemistry_cases, find_unbuilt, unbuilt_message


class FakeCase:
    def __init__(self, trace, uuid, params, slug):
        self.trace = trace
        self.uuid = uuid
        self.params = params
        self.slug = slug

    def get_uuid(self):
        return self.uuid

    def to_input_file(self):
        return self


class FakeTarget:
    def __init__(self, name, installed_slugs):
        self.name = name
        self.installed_slugs = installed_slugs

    def get_slug(self, case):
        return f"{self.name}-{case.slug}"

    def get_install_binpath(self, case):
        return f"/nonexistent/{self.get_slug(case)}/bin/{self.name}" if case.slug not in self.installed_slugs else sys.executable


PLAIN = FakeCase("1D -> bc=-1", "AAAAAAAA", {}, "plain")
PLAIN2 = FakeCase("1D -> bc=-2", "BBBBBBBB", {"chemistry": "F"}, "plain")
CHEM = FakeCase("1D -> Chemistry -> Perfect Reactor", "CCCCCCCC", {"chemistry": "T"}, "chem")
REACTING_EXAMPLE = FakeCase("Example -> 2D -> shock_flame", "DDDDDDDD", {"chemistry": "T"}, "chem2")


def test_nothing_is_reported_when_every_build_exists():
    codes = [FakeTarget("pre_process", {"plain", "chem"}), FakeTarget("simulation", {"plain", "chem"})]
    assert find_unbuilt([PLAIN, PLAIN2, CHEM], codes) == []


def test_a_missing_chemistry_build_is_reported_once_per_target_with_all_its_cases():
    codes = [FakeTarget("pre_process", {"plain"}), FakeTarget("simulation", {"plain"})]
    unbuilt = find_unbuilt([PLAIN, CHEM, PLAIN2, REACTING_EXAMPLE], codes)

    assert [(e["target"], e["slug"]) for e in unbuilt] == [
        ("pre_process", "pre_process-chem"),
        ("simulation", "simulation-chem"),
        ("pre_process", "pre_process-chem2"),
        ("simulation", "simulation-chem2"),
    ]
    assert all(len(e["cases"]) == 1 for e in unbuilt)


def test_cases_sharing_a_missing_build_are_grouped():
    codes = [FakeTarget("simulation", set())]
    unbuilt = find_unbuilt([PLAIN, PLAIN2], codes)

    assert len(unbuilt) == 1
    assert unbuilt[0]["cases"] == [PLAIN, PLAIN2]


def test_message_suggests_no_chemistry_only_when_every_missing_case_uses_chemistry():
    codes = [FakeTarget("simulation", {"plain"})]
    assert "--no-chemistry" in unbuilt_message(find_unbuilt([PLAIN, CHEM], codes))

    codes = [FakeTarget("simulation", set())]
    assert "--no-chemistry" not in unbuilt_message(find_unbuilt([PLAIN, CHEM], codes))


def test_message_counts_distinct_cases_across_targets():
    codes = [FakeTarget("pre_process", set()), FakeTarget("simulation", set())]
    assert "needed by 2 test case(s)" in unbuilt_message(find_unbuilt([PLAIN, CHEM], codes))


def test_no_chemistry_keys_on_the_parameter_not_the_trace():
    # The reacting example has no "Chemistry" trace element but still needs a chemistry build.
    cases = [PLAIN, CHEM, PLAIN2, REACTING_EXAMPLE]
    builders = ["builder-plain", "builder-chem", "builder-plain2", "builder-example"]
    kept, skipped = _drop_chemistry_cases(builders, cases, ["already-skipped"])

    assert kept == [PLAIN, PLAIN2]
    # Skipped entries stay builders, matching what __filter puts there.
    assert skipped == ["already-skipped", "builder-chem", "builder-example"]


def test_a_non_executable_file_is_not_an_installed_binary(tmp_path):
    binpath = tmp_path / "simulation"
    binpath.write_text("")
    binpath.chmod(0o644)

    class Target:
        name = "simulation"

        def get_slug(self, case):
            return case.slug

        def get_install_binpath(self, _case):
            return str(binpath)

    assert [e["target"] for e in find_unbuilt([PLAIN], [Target()])] == ["simulation"]
    binpath.chmod(0o755)
    assert find_unbuilt([PLAIN], [Target()]) == []
