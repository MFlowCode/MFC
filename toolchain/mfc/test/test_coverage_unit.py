import tempfile
from pathlib import Path

from mfc.test.coverage import is_always_run_all, load_map, param_hash, save_map


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
