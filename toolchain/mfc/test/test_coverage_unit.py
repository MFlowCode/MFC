"""
Unit tests for toolchain/mfc/test/coverage.py

Run with:
    python3 -m pytest toolchain/mfc/test/test_coverage_unit.py -v

These tests are fully offline (no build, no git, no gcov binary required).
They use mocks and in-memory data structures to verify logic.
"""
# pylint: disable=protected-access,exec-used,too-few-public-methods,wrong-import-position

import gzip
import importlib.util
import json
import os
import sys
import types
import unittest
from unittest.mock import patch

# ---------------------------------------------------------------------------
# Import the module under test.
# We patch the module-level imports that require the full toolchain.
# ---------------------------------------------------------------------------

# Create minimal stubs for toolchain modules so coverage.py can be imported
# without the full MFC toolchain being on sys.path.
def _make_stub(name):
    mod = types.ModuleType(name)
    sys.modules[name] = mod
    return mod


for _mod_name in [
    "toolchain",
    "toolchain.mfc",
    "toolchain.mfc.printer",
    "toolchain.mfc.common",
    "toolchain.mfc.build",
    "toolchain.mfc.test",
    "toolchain.mfc.test.case",
]:
    if _mod_name not in sys.modules:
        _make_stub(_mod_name)

# Provide the attributes coverage.py needs from its relative imports
_printer_stub = sys.modules.get("toolchain.mfc.printer", _make_stub("toolchain.mfc.printer"))


class _FakeCons:
    def print(self, *args, **kwargs):
        pass  # suppress output during tests


_printer_stub.cons = _FakeCons()

_common_stub = sys.modules.get("toolchain.mfc.common", _make_stub("toolchain.mfc.common"))
_common_stub.MFC_ROOT_DIR = "/fake/repo"


class _FakeMFCException(Exception):
    pass


_common_stub.MFCException = _FakeMFCException

_build_stub = sys.modules.get("toolchain.mfc.build", _make_stub("toolchain.mfc.build"))
_build_stub.PRE_PROCESS = "pre_process"
_build_stub.SIMULATION = "simulation"
_build_stub.POST_PROCESS = "post_process"
_build_stub.SYSCHECK = "syscheck"

_case_stub = sys.modules.get("toolchain.mfc.test.case", _make_stub("toolchain.mfc.test.case"))
_case_stub.input_bubbles_lagrange = lambda case: None

# Load coverage.py by injecting stubs into sys.modules so relative imports resolve.
_COVERAGE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coverage.py")

sys.modules.pop("toolchain.mfc.test.coverage", None)  # reset if already loaded

_spec = importlib.util.spec_from_file_location(
    "toolchain.mfc.test.coverage",
    _COVERAGE_PATH,
    submodule_search_locations=[]
)
_coverage_mod = importlib.util.module_from_spec(_spec)
_coverage_mod.__package__ = "toolchain.mfc.test"

sys.modules["toolchain.mfc.test"] = types.ModuleType("toolchain.mfc.test")
sys.modules["toolchain.mfc.test"].__package__ = "toolchain.mfc.test"

with patch.dict("sys.modules", {
    "toolchain.mfc.printer":   _printer_stub,
    "toolchain.mfc.common":    _common_stub,
    "toolchain.mfc.build":     _build_stub,
    "toolchain.mfc.test.case": _case_stub,
}):
    try:
        _spec.loader.exec_module(_coverage_mod)
    except ImportError:
        pass  # fallback below

# If the importlib approach failed (relative imports unresolvable), fall back to exec.
try:
    _parse_diff_files = _coverage_mod._parse_diff_files
    _parse_gcov_json_output = _coverage_mod._parse_gcov_json_output
    _normalize_cache = _coverage_mod._normalize_cache
    _gcda_path_to_fpp = _coverage_mod._gcda_path_to_fpp
    should_run_all_tests = _coverage_mod.should_run_all_tests
    filter_tests_by_coverage = _coverage_mod.filter_tests_by_coverage
    ALWAYS_RUN_ALL = _coverage_mod.ALWAYS_RUN_ALL
    COVERAGE_CACHE_PATH = _coverage_mod.COVERAGE_CACHE_PATH
except AttributeError:
    _globals = {
        "__name__": "toolchain.mfc.test.coverage",
        "__package__": "toolchain.mfc.test",
        "cons": _printer_stub.cons,
        "common": _common_stub,
        "MFCException": _FakeMFCException,
        "PRE_PROCESS": "pre_process",
        "SIMULATION": "simulation",
        "POST_PROCESS": "post_process",
        "SYSCHECK": "syscheck",
    }
    with open(_COVERAGE_PATH, encoding="utf-8") as _f:
        _src = _f.read()

    _src = (
        _src
        .replace("from ..printer import cons", "cons = _globals['cons']")
        .replace("from .. import common", "")
        .replace("from ..common import MFCException", "MFCException = _globals['MFCException']")
        .replace("from ..build import PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK", "")
        .replace("from .case import input_bubbles_lagrange", "input_bubbles_lagrange = lambda case: None")
    )
    exec(compile(_src, _COVERAGE_PATH, "exec"), _globals)  # noqa: S102

    _parse_diff_files = _globals["_parse_diff_files"]
    _parse_gcov_json_output = _globals["_parse_gcov_json_output"]
    _normalize_cache = _globals["_normalize_cache"]
    _gcda_path_to_fpp = _globals["_gcda_path_to_fpp"]
    should_run_all_tests = _globals["should_run_all_tests"]
    filter_tests_by_coverage = _globals["filter_tests_by_coverage"]
    ALWAYS_RUN_ALL = _globals["ALWAYS_RUN_ALL"]
    COVERAGE_CACHE_PATH = _globals["COVERAGE_CACHE_PATH"]


# ---------------------------------------------------------------------------
# Helper: minimal fake test case
# ---------------------------------------------------------------------------

class FakeCase:
    """Minimal stand-in for TestCase — only get_uuid() is needed."""

    def __init__(self, uuid: str):
        self._uuid = uuid

    def get_uuid(self) -> str:
        return self._uuid


# ===========================================================================
# Group 1: _parse_diff_files — git diff --name-only parsing
# ===========================================================================

class TestParseDiffFiles(unittest.TestCase):

    def test_parse_single_file(self):
        result = _parse_diff_files("src/simulation/m_rhs.fpp\n")
        assert result == {"src/simulation/m_rhs.fpp"}

    def test_parse_multiple_files(self):
        text = "src/simulation/m_rhs.fpp\nsrc/simulation/m_weno.fpp\nREADME.md\n"
        result = _parse_diff_files(text)
        assert result == {
            "src/simulation/m_rhs.fpp",
            "src/simulation/m_weno.fpp",
            "README.md",
        }

    def test_parse_empty(self):
        assert _parse_diff_files("") == set()
        assert _parse_diff_files("\n") == set()

    def test_parse_ignores_blank_lines(self):
        text = "src/simulation/m_rhs.fpp\n\n\nsrc/simulation/m_weno.fpp\n"
        result = _parse_diff_files(text)
        assert result == {"src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"}

    def test_parse_mixed_extensions(self):
        text = "src/simulation/m_rhs.fpp\ntoolchain/mfc/test/cases.py\nCMakeLists.txt\n"
        result = _parse_diff_files(text)
        assert len(result) == 3
        assert "toolchain/mfc/test/cases.py" in result
        assert "CMakeLists.txt" in result


# ===========================================================================
# Group 2: should_run_all_tests — ALWAYS_RUN_ALL detection
# ===========================================================================

class TestShouldRunAllTests(unittest.TestCase):

    def test_parallel_macros_triggers_all(self):
        assert should_run_all_tests(
            {"src/common/include/parallel_macros.fpp"}
        ) is True

    def test_acc_macros_triggers_all(self):
        assert should_run_all_tests(
            {"src/common/include/acc_macros.fpp"}
        ) is True

    def test_omp_macros_triggers_all(self):
        assert should_run_all_tests(
            {"src/common/include/omp_macros.fpp"}
        ) is True

    def test_shared_parallel_macros_triggers_all(self):
        assert should_run_all_tests(
            {"src/common/include/shared_parallel_macros.fpp"}
        ) is True

    def test_macros_fpp_triggers_all(self):
        assert should_run_all_tests(
            {"src/common/include/macros.fpp"}
        ) is True

    def test_cases_py_triggers_all(self):
        assert should_run_all_tests(
            {"toolchain/mfc/test/cases.py"}
        ) is True

    def test_case_py_triggers_all(self):
        assert should_run_all_tests(
            {"toolchain/mfc/test/case.py"}
        ) is True

    def test_definitions_py_triggers_all(self):
        assert should_run_all_tests(
            {"toolchain/mfc/params/definitions.py"}
        ) is True

    def test_input_py_triggers_all(self):
        assert should_run_all_tests(
            {"toolchain/mfc/run/input.py"}
        ) is True

    def test_case_validator_triggers_all(self):
        assert should_run_all_tests(
            {"toolchain/mfc/case_validator.py"}
        ) is True

    def test_cmakelists_triggers_all(self):
        assert should_run_all_tests(
            {"CMakeLists.txt"}
        ) is True

    def test_simulation_module_does_not_trigger_all(self):
        assert should_run_all_tests(
            {"src/simulation/m_rhs.fpp"}
        ) is False

    def test_empty_set_does_not_trigger_all(self):
        assert should_run_all_tests(set()) is False

    def test_mixed_one_trigger_fires_all(self):
        assert should_run_all_tests({
            "src/simulation/m_rhs.fpp",
            "src/common/include/macros.fpp",
        }) is True


# ===========================================================================
# Group 3: filter_tests_by_coverage — core file-level selection logic
# ===========================================================================

class TestFilterTestsByCoverage(unittest.TestCase):

    def test_file_overlap_includes_test(self):
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"]}
        changed = {"src/simulation/m_rhs.fpp"}
        cases = [FakeCase("AAAA0001")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        assert len(to_run) == 1
        assert len(skipped) == 0

    def test_no_file_overlap_skips_test(self):
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp"]}
        changed = {"src/simulation/m_weno.fpp"}
        cases = [FakeCase("AAAA0001")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        assert len(to_run) == 0
        assert len(skipped) == 1

    def test_uuid_not_in_cache_is_conservative(self):
        """Newly added test not in cache -> include it (conservative)."""
        cache = {}
        changed = {"src/simulation/m_rhs.fpp"}
        to_run, _ = filter_tests_by_coverage([FakeCase("NEWTEST1")], cache, changed)
        assert len(to_run) == 1

    def test_no_fpp_changes_skips_all(self):
        """Only non-.fpp files changed -> skip all tests."""
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp"]}
        changed = {"toolchain/setup.py", "README.md"}
        cases = [FakeCase("AAAA0001")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        assert len(to_run) == 0
        assert len(skipped) == 1

    def test_empty_changed_files_skips_all(self):
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp"]}
        changed = set()
        to_run, skipped = filter_tests_by_coverage([FakeCase("AAAA0001")], cache, changed)
        assert len(to_run) == 0
        assert len(skipped) == 1

    def test_multiple_tests_partial_selection(self):
        """Only the test covering the changed file should run."""
        cache = {
            "TEST_A": ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"],
            "TEST_B": ["src/simulation/m_bubbles.fpp"],
            "TEST_C": ["src/simulation/m_rhs.fpp"],
        }
        changed = {"src/simulation/m_bubbles.fpp"}
        cases = [FakeCase("TEST_A"), FakeCase("TEST_B"), FakeCase("TEST_C")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        uuids_run = {c.get_uuid() for c in to_run}
        assert uuids_run == {"TEST_B"}
        assert len(skipped) == 2

    def test_multiple_changed_files_union(self):
        """Changing multiple files includes any test that covers any of them."""
        cache = {
            "TEST_A": ["src/simulation/m_rhs.fpp"],
            "TEST_B": ["src/simulation/m_weno.fpp"],
            "TEST_C": ["src/simulation/m_bubbles.fpp"],
        }
        changed = {"src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"}
        cases = [FakeCase("TEST_A"), FakeCase("TEST_B"), FakeCase("TEST_C")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        uuids_run = {c.get_uuid() for c in to_run}
        assert uuids_run == {"TEST_A", "TEST_B"}
        assert len(skipped) == 1

    def test_test_covering_multiple_files_matched_via_second(self):
        """Test matched because m_weno.fpp (its second covered file) was changed."""
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"]}
        changed = {"src/simulation/m_weno.fpp"}
        to_run, _ = filter_tests_by_coverage([FakeCase("AAAA0001")], cache, changed)
        assert len(to_run) == 1

    def test_empty_cache_runs_all_conservatively(self):
        """Empty coverage cache -> all tests included (conservative)."""
        cache = {}
        changed = {"src/simulation/m_rhs.fpp"}
        cases = [FakeCase("T1"), FakeCase("T2"), FakeCase("T3")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        assert len(to_run) == 3
        assert len(skipped) == 0

    def test_mixed_fpp_and_nonfpp_changes(self):
        """Non-.fpp files in changed set are ignored for matching."""
        cache = {"TEST_A": ["src/simulation/m_rhs.fpp"]}
        changed = {"src/simulation/m_rhs.fpp", "README.md", "toolchain/setup.py"}
        to_run, _ = filter_tests_by_coverage([FakeCase("TEST_A")], cache, changed)
        assert len(to_run) == 1

    def test_incomplete_coverage_included_conservatively(self):
        """Test with no simulation coverage but simulation file changed -> include."""
        cache = {
            "GOOD_T": ["src/simulation/m_rhs.fpp", "src/pre_process/m_start_up.fpp"],
            "BAD_T":  ["src/pre_process/m_start_up.fpp", "src/common/m_helper.fpp"],
        }
        changed = {"src/simulation/m_rhs.fpp"}
        cases = [FakeCase("GOOD_T"), FakeCase("BAD_T")]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        uuids_run = {c.get_uuid() for c in to_run}
        assert "GOOD_T" in uuids_run  # direct file overlap
        assert "BAD_T" in uuids_run   # no sim coverage -> conservative include
        assert len(skipped) == 0

    def test_incomplete_coverage_not_triggered_by_preprocess(self):
        """Test with no sim coverage is NOT auto-included for pre_process changes."""
        cache = {
            "BAD_T": ["src/pre_process/m_start_up.fpp"],
        }
        changed = {"src/pre_process/m_data_output.fpp"}
        to_run, skipped = filter_tests_by_coverage([FakeCase("BAD_T")], cache, changed)
        assert len(to_run) == 0  # no sim change, no overlap -> skip
        assert len(skipped) == 1


# ===========================================================================
# Group 4: Corner cases from design discussion
# ===========================================================================

class TestDesignCornerCases(unittest.TestCase):

    def test_gpu_ifdef_file_still_triggers_if_covered(self):
        """
        GPU-specific code lives in the same .fpp file as CPU code.
        At file level, changing any part of the file triggers tests that cover it.
        """
        cache = {"MUSCL_T": ["src/simulation/m_muscl.fpp"]}
        changed = {"src/simulation/m_muscl.fpp"}
        to_run, _ = filter_tests_by_coverage([FakeCase("MUSCL_T")], cache, changed)
        assert len(to_run) == 1

    def test_macro_file_triggers_all_via_should_run_all(self):
        """parallel_macros.fpp in changed files -> should_run_all_tests() is True."""
        assert should_run_all_tests({"src/common/include/parallel_macros.fpp"}) is True

    def test_new_fpp_file_no_coverage_skips(self):
        """
        Brand new .fpp file has no coverage in cache.
        All tests are skipped (no test covers the new file).
        """
        cache = {"AAAA0001": ["src/simulation/m_rhs.fpp"]}
        changed = {"src/simulation/m_brand_new.fpp"}
        to_run, skipped = filter_tests_by_coverage([FakeCase("AAAA0001")], cache, changed)
        assert len(to_run) == 0
        assert len(skipped) == 1

    def test_non_fpp_always_run_all_detected(self):
        """
        End-to-end: diff lists only cases.py (non-.fpp) ->
        _parse_diff_files includes it -> should_run_all_tests fires.
        """
        files = _parse_diff_files("toolchain/mfc/test/cases.py\n")
        assert should_run_all_tests(files) is True

    def test_niche_feature_pruning(self):
        """
        Niche features: most tests don't cover m_bubbles.fpp.
        Changing it skips tests that don't touch it.
        """
        cache = {
            "BUBBLE1": ["src/simulation/m_bubbles.fpp", "src/simulation/m_rhs.fpp"],
            "BUBBLE2": ["src/simulation/m_bubbles.fpp"],
            "BASIC_1": ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"],
            "BASIC_2": ["src/simulation/m_rhs.fpp"],
            "BASIC_3": ["src/simulation/m_weno.fpp"],
        }
        changed = {"src/simulation/m_bubbles.fpp"}
        cases = [FakeCase(u) for u in ["BUBBLE1", "BUBBLE2", "BASIC_1", "BASIC_2", "BASIC_3"]]
        to_run, skipped = filter_tests_by_coverage(cases, cache, changed)
        uuids_run = {c.get_uuid() for c in to_run}
        assert uuids_run == {"BUBBLE1", "BUBBLE2"}
        assert len(skipped) == 3


# ===========================================================================
# Group 5: _parse_gcov_json_output — gcov JSON parsing (file-level)
# ===========================================================================

class TestParseGcovJsonOutput(unittest.TestCase):

    def _make_gcov_json(self, files_data: list) -> bytes:
        """Build a fake gzip-compressed gcov JSON blob."""
        data = {
            "format_version": "2",
            "gcc_version": "15.2.0",
            "files": files_data,
        }
        return gzip.compress(json.dumps(data).encode())

    def test_returns_set_of_covered_fpp_files(self):
        compressed = self._make_gcov_json([{
            "file": "/repo/src/simulation/m_rhs.fpp",
            "lines": [
                {"line_number": 45, "count": 3},
                {"line_number": 46, "count": 0},
                {"line_number": 47, "count": 1},
            ],
        }])
        result = _parse_gcov_json_output(compressed, "/repo")
        assert result == {"src/simulation/m_rhs.fpp"}

    def test_ignores_file_with_zero_coverage(self):
        compressed = self._make_gcov_json([{
            "file": "/repo/src/simulation/m_rhs.fpp",
            "lines": [
                {"line_number": 10, "count": 0},
                {"line_number": 11, "count": 0},
            ],
        }])
        result = _parse_gcov_json_output(compressed, "/repo")
        assert result == set()

    def test_ignores_f90_files(self):
        """Generated .f90 files must not appear in coverage output."""
        compressed = self._make_gcov_json([
            {
                "file": "/repo/build/fypp/simulation/m_rhs.fpp.f90",
                "lines": [{"line_number": 10, "count": 5}],
            },
            {
                "file": "/repo/src/simulation/m_rhs.fpp",
                "lines": [{"line_number": 45, "count": 1}],
            },
        ])
        result = _parse_gcov_json_output(compressed, "/repo")
        assert result == {"src/simulation/m_rhs.fpp"}

    def test_handles_raw_json_gcov12(self):
        """gcov 12 outputs raw JSON (not gzip). Must parse correctly."""
        data = {
            "format_version": "1",
            "gcc_version": "12.3.0",
            "files": [{
                "file": "/repo/src/simulation/m_rhs.fpp",
                "lines": [{"line_number": 45, "count": 3}],
            }],
        }
        raw = json.dumps(data).encode()
        result = _parse_gcov_json_output(raw, "/repo")
        assert result == {"src/simulation/m_rhs.fpp"}

    def test_handles_invalid_data_gracefully(self):
        result = _parse_gcov_json_output(b"not valid gzip or json", "/repo")
        assert result == set()

    def test_handles_empty_files_list(self):
        compressed = self._make_gcov_json([])
        result = _parse_gcov_json_output(compressed, "/repo")
        assert result == set()

    def test_multiple_fpp_files(self):
        compressed = self._make_gcov_json([
            {
                "file": "/repo/src/simulation/m_rhs.fpp",
                "lines": [{"line_number": 45, "count": 1}],
            },
            {
                "file": "/repo/src/simulation/m_weno.fpp",
                "lines": [{"line_number": 200, "count": 2}],
            },
        ])
        result = _parse_gcov_json_output(compressed, "/repo")
        assert result == {"src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"}


# ===========================================================================
# Group 6: _normalize_cache — old format conversion
# ===========================================================================

class TestNormalizeCache(unittest.TestCase):

    def test_converts_old_line_level_format(self):
        """Old format {uuid: {file: [lines]}} -> new format {uuid: [files]}."""
        old_cache = {
            "TEST_A": {
                "src/simulation/m_rhs.fpp": [45, 46, 47],
                "src/simulation/m_weno.fpp": [100, 200],
            },
            "TEST_B": {
                "src/simulation/m_bubbles.fpp": [10],
            },
            "_meta": {"cases_hash": "abc123"},
        }
        result = _normalize_cache(old_cache)
        assert isinstance(result["TEST_A"], list)
        assert set(result["TEST_A"]) == {"src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"}
        assert result["TEST_B"] == ["src/simulation/m_bubbles.fpp"]
        assert result["_meta"] == {"cases_hash": "abc123"}

    def test_new_format_unchanged(self):
        """New format {uuid: [files]} passes through unchanged."""
        new_cache = {
            "TEST_A": ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"],
            "_meta": {"cases_hash": "abc123"},
        }
        result = _normalize_cache(new_cache)
        assert result["TEST_A"] == ["src/simulation/m_rhs.fpp", "src/simulation/m_weno.fpp"]

    def test_empty_coverage_dict_becomes_empty_list(self):
        """Test with 0 coverage (old format: empty dict) -> empty list."""
        old_cache = {"TEST_A": {}, "_meta": {"cases_hash": "abc"}}
        result = _normalize_cache(old_cache)
        assert result["TEST_A"] == []


# ===========================================================================
# Group 7: Cache path format
# ===========================================================================

class TestCachePath(unittest.TestCase):

    def test_cache_path_is_gzipped(self):
        """Cache file must use .json.gz so it can be committed to the repo."""
        assert str(COVERAGE_CACHE_PATH).endswith(".json.gz")


# ===========================================================================
# Group 8: _gcda_path_to_fpp — .gcda path to .fpp source mapping
# ===========================================================================

class TestGcdaPathToFpp(unittest.TestCase):

    def test_fypp_simulation_file(self):
        path = "build/staging/abc123/CMakeFiles/simulation.dir/fypp/simulation/m_rhs.fpp.f90.gcda"
        assert _gcda_path_to_fpp(path) == "src/simulation/m_rhs.fpp"

    def test_fypp_pre_process_file(self):
        path = "build/staging/abc123/CMakeFiles/pre_process.dir/fypp/pre_process/m_grid.fpp.f90.gcda"
        assert _gcda_path_to_fpp(path) == "src/pre_process/m_grid.fpp"

    def test_fypp_common_file_in_simulation(self):
        """Common .fpp compiled into simulation: path says simulation, not common."""
        path = "build/staging/abc123/CMakeFiles/simulation.dir/fypp/simulation/m_helper.fpp.f90.gcda"
        assert _gcda_path_to_fpp(path) == "src/simulation/m_helper.fpp"

    def test_plain_f90_returns_empty(self):
        """Non-.fpp files (plain .f90) should be excluded."""
        path = "build/staging/abc123/CMakeFiles/simulation.dir/src/common/m_compile_specific.f90.gcda"
        assert _gcda_path_to_fpp(path) == ""

    def test_module_file_returns_empty(self):
        """Generated module files should be excluded."""
        path = "build/staging/abc123/CMakeFiles/simulation.dir/modules/simulation/m_thermochem.f90.gcda"
        assert _gcda_path_to_fpp(path) == ""

    def test_ltrans_file_returns_empty(self):
        """LTO artifacts should be excluded."""
        path = "build/staging/abc123/simulation.ltrans0.ltrans.gcda"
        assert _gcda_path_to_fpp(path) == ""

    def test_syscheck_fpp(self):
        path = "build/staging/abc123/CMakeFiles/syscheck.dir/fypp/syscheck/syscheck.fpp.f90.gcda"
        assert _gcda_path_to_fpp(path) == "src/syscheck/syscheck.fpp"


if __name__ == "__main__":
    unittest.main()
