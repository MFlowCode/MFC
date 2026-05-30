import tempfile
from pathlib import Path

from mfc.test.coverage import is_always_run_all, load_map, param_hash, save_map, select_tests


def test_param_hash_is_order_independent():
    a = param_hash({"m": 100, "weno_order": 5, "bubbles_euler": "T"})
    b = param_hash({"weno_order": 5, "bubbles_euler": "T", "m": 100})
    assert a == b


def test_param_hash_changes_with_value():
    a = param_hash({"weno_order": 5})
    b = param_hash({"weno_order": 3})
    assert a != b


def test_param_hash_is_hex_and_short():
    h = param_hash({"m": 1})
    assert len(h) == 16 and all(c in "0123456789abcdef" for c in h)


def test_param_hash_nested_order_independent():
    a = param_hash({"patch": {"x": 1, "y": 2}})
    b = param_hash({"patch": {"y": 2, "x": 1}})
    assert a == b


def test_save_then_load_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        save_map(p, {"abc": ["src/simulation/m_rhs.fpp"]}, n_tests=1, git_sha="deadbee", gfortran_version="13")
        entries, meta = load_map(p)
        assert entries == {"abc": ["src/simulation/m_rhs.fpp"]}
        assert meta["n_tests"] == 1 and meta["git_sha"] == "deadbee"
        assert "built_at" in meta


def test_load_missing_returns_none():
    assert load_map(Path("/nonexistent/m.json.gz")) == (None, None)


def test_load_corrupt_returns_none():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        p.write_bytes(b"not gzip")
        assert load_map(p) == (None, None)


def test_macro_file_forces_all():
    assert is_always_run_all({"src/common/include/parallel_macros.fpp"})


def test_cmake_forces_all():
    assert is_always_run_all({"CMakeLists.txt"})
    assert is_always_run_all({"toolchain/cmake/foo.cmake"})


def test_param_codegen_forces_all():
    assert is_always_run_all({"toolchain/mfc/params/definitions.py"})


def test_ordinary_common_module_does_not_force_all():
    assert not is_always_run_all({"src/common/m_helper.fpp"})


def test_ordinary_sim_module_does_not_force_all():
    assert not is_always_run_all({"src/simulation/m_rhs.fpp"})


class _Case:
    def __init__(self, ph, params=None):
        self._ph = ph
        self.params = params or {}

    def coverage_key(self):
        return self._ph


def _cases(*phs):
    return [_Case(p) for p in phs]


def test_rung1_no_changed_files_runs_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, None)
    assert len(run) == 2 and skip == [] and reason.startswith("rung1")


def test_rung2_always_run_all():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {"a": [], "b": []}, {"CMakeLists.txt"})
    assert len(run) == 2 and reason.startswith("rung2")


def test_rung3_f90_change_runs_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": []}, {"src/common/m_precision_select.f90"})
    assert len(run) == 1 and reason.startswith("rung3")


def test_rung4_changed_fpp_with_zero_coverage_runs_all():
    cases = _cases("a")
    # m_gpu_only.fpp is covered by no test in the map
    run, skip, reason = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_gpu_only.fpp"})
    assert len(run) == 1 and reason.startswith("rung4")


def test_rung5_unmapped_test_is_included():
    cases = _cases("a", "new")  # 'new' not in map
    run, skip, _ = select_tests(cases, {"a": ["src/simulation/m_rhs.fpp"]}, {"src/simulation/m_rhs.fpp"})
    assert {c.coverage_key() for c in run} == {"a", "new"}


def test_rung6_and_7_overlap_selects_subset():
    cases = _cases("hit", "miss")
    cov = {"hit": ["src/simulation/m_bubbles_EE.fpp"], "miss": ["src/simulation/m_rhs.fpp"]}
    run, skip, _ = select_tests(cases, cov, {"src/simulation/m_bubbles_EE.fpp"})
    assert [c.coverage_key() for c in run] == ["hit"]
    assert [c.coverage_key() for c in skip] == ["miss"]


def test_case_coverage_key_matches_param_hash():
    from mfc.test.case import TestCase

    tc = TestCase("1D -> Foo", {"m": 100, "weno_order": 5})
    assert tc.coverage_key() == param_hash({"m": 100, "weno_order": 5})


def test_case_coverage_key_ignores_trace():
    from mfc.test.case import TestCase

    a = TestCase("1D -> Foo", {"m": 100})
    b = TestCase("totally -> different -> trace", {"m": 100})
    assert a.coverage_key() == b.coverage_key()
