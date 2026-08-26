import os
import subprocess
import sys
import tempfile
import types as _types
from pathlib import Path
from unittest.mock import patch

from mfc.test.coverage import canonicalize_param_paths, entries_equal, format_summary, get_changed_files, is_always_run_all, load_map, map_health, param_hash, save_map, select_tests


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
    def __init__(self, ph, params=None, canary=False):
        self._ph = ph
        self.params = params or {}
        self.canary = canary

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


def test_canary_always_runs_even_when_its_coverage_misses():
    # 'canary' and 'plain' have identical, non-overlapping coverage; only the canary must
    # run. 'other' covers the changed file so rung 4 (run-all) does not fire and we reach
    # the per-test rungs where the canary bypass is exercised.
    cases = [_Case("canary", canary=True), _Case("plain"), _Case("other")]
    cov = {
        "canary": ["src/simulation/m_viscous.fpp"],
        "plain": ["src/simulation/m_viscous.fpp"],
        "other": ["src/simulation/m_rhs.fpp"],
    }
    run, skip, _ = select_tests(cases, cov, {"src/simulation/m_rhs.fpp"})
    run_keys = {c.coverage_key() for c in run}
    assert "canary" in run_keys and "other" in run_keys
    assert [c.coverage_key() for c in skip] == ["plain"]


def test_case_coverage_key_uses_full_params():
    from mfc.test.case import TestCase

    tc = TestCase("1D -> Foo", {"m": 100, "weno_order": 5})
    assert tc.coverage_key() == param_hash(tc.params)


def test_case_coverage_key_changes_with_params():
    from mfc.test.case import TestCase

    a = TestCase("1D -> Foo", {"weno_order": 5})
    b = TestCase("1D -> Foo", {"weno_order": 3})
    assert a.coverage_key() != b.coverage_key()


def test_case_coverage_key_ignores_trace():
    from mfc.test.case import TestCase

    a = TestCase("1D -> Foo", {"m": 100})
    b = TestCase("totally -> different -> trace", {"m": 100})
    assert a.coverage_key() == b.coverage_key()


def test_changed_files_prefers_explicit_list():
    files = get_changed_files("/repo", "master", explicit="src/a.fpp\nsrc/b.fpp\n")
    assert files == {"src/a.fpp", "src/b.fpp"}


def test_changed_files_deepens_then_recovers():
    state = {"deepened": False}

    def fake_run(cmd, **kw):
        sub = cmd[1] if len(cmd) > 1 else ""
        if sub == "fetch":
            state["deepened"] = True
            return _types.SimpleNamespace(returncode=0, stdout="", stderr="")
        if sub == "merge-base":
            return _types.SimpleNamespace(returncode=0 if state["deepened"] else 1, stdout="base\n", stderr="")
        if sub == "diff":
            return _types.SimpleNamespace(returncode=0, stdout="src/x.fpp\n", stderr="")
        return _types.SimpleNamespace(returncode=0, stdout="", stderr="")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") == {"src/x.fpp"}


def test_changed_files_returns_none_when_unrecoverable():
    def fake_run(cmd, **kw):
        rc = 1 if (len(cmd) > 1 and cmd[1] == "merge-base") else 0
        return _types.SimpleNamespace(returncode=rc, stdout="", stderr="boom")

    with patch("subprocess.run", fake_run):
        assert get_changed_files("/repo", "master") is None


def test_summary_mentions_counts_age_reason():
    s = format_summary(
        ran=47,
        total=610,
        reason="selected 47/610 by coverage overlap",
        meta={"built_at": "2026-05-20T00:00:00+00:00"},
        now="2026-05-29T00:00:00+00:00",
    )
    assert "47/610" in s and "9d" in s and "coverage overlap" in s


def test_summary_handles_missing_meta():
    s = format_summary(
        ran=610,
        total=610,
        reason="rung1: changed-file list unavailable",
        meta=None,
        now="2026-05-29T00:00:00+00:00",
    )
    assert "610/610" in s and "map age unknown" in s


def test_health_ok():
    ok, msg = map_health(
        meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 600},
        current_keys=set(str(i) for i in range(600)),
        mapped_keys=set(str(i) for i in range(580)),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert ok, msg


def test_health_stale_fails():
    ok, msg = map_health(
        meta={"built_at": "2026-05-01T00:00:00+00:00", "n_tests": 600},
        current_keys=set(["a"]),
        mapped_keys=set(["a"]),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert not ok and "stale" in msg.lower()


def test_health_undercoverage_fails():
    ok, msg = map_health(
        meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 10},
        current_keys=set(str(i) for i in range(100)),
        mapped_keys=set(str(i) for i in range(50)),
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
    )
    assert not ok and "coverage" in msg.lower()


def test_builder_has_coverage_key_matching_case():
    from mfc.test.case import TestCaseBuilder

    b = TestCaseBuilder(trace="1D -> Foo", mods={"m": 100, "weno_order": 5}, path="", args=[], ppn=1, functor=None)
    assert b.coverage_key() == b.to_case().coverage_key()


def test_rung5_empty_coverage_is_included():
    # a test whose map entry is [] (uncertain) must be RUN, not skipped, on a .fpp change.
    # "anchor" covers the changed .fpp so rung4 passes and we reach the per-test rungs.
    cases = _cases("hasempty", "anchor")
    cov_map = {
        "hasempty": [],
        "anchor": ["src/simulation/m_rhs.fpp"],
    }
    run, skip, _ = select_tests(cases, cov_map, {"src/simulation/m_rhs.fpp"})
    run_keys = {c.coverage_key() for c in run}
    assert "hasempty" in run_keys and skip == []


def test_changed_files_explicit_space_and_comma_separated():
    from mfc.test.coverage import get_changed_files

    assert get_changed_files("/r", "master", explicit="src/a.fpp src/b.fpp") == {"src/a.fpp", "src/b.fpp"}
    assert get_changed_files("/r", "master", explicit="src/a.fpp,src/b.fpp") == {"src/a.fpp", "src/b.fpp"}


def test_sim_include_fpp_forces_all():
    # gcov can't reliably attribute #:include'd files; any src include change runs all.
    assert is_always_run_all({"src/simulation/include/inline_riemann.fpp"})
    assert is_always_run_all({"src/pre_process/include/2dHardcodedIC.fpp"})


def test_empty_explicit_is_uncertainty_not_skipall():
    # An empty/whitespace --changed-files must NOT become an empty set (skip-all under
    # enforce). It falls through to git detection -> None when that fails -> run all.
    def fail_git(cmd, **kw):
        rc = 1 if (len(cmd) > 1 and cmd[1] == "merge-base") else 0
        return _types.SimpleNamespace(returncode=rc, stdout="", stderr="x")

    with patch("subprocess.run", fail_git):
        assert get_changed_files("/r", "master", explicit="") is None
        assert get_changed_files("/r", "master", explicit="  ,  ") is None


def test_run_and_test_infra_force_all():
    assert is_always_run_all({"toolchain/mfc/run/input.py"})
    assert is_always_run_all({"toolchain/mfc/test/case.py"})
    assert is_always_run_all({"toolchain/mfc/test/test.py"})


def test_cases_py_is_not_always_run():
    assert not is_always_run_all({"toolchain/mfc/test/cases.py"})


def test_cases_py_change_runs_new_tests_not_skipall():
    # cases.py-only change must run the NEW/modified tests (rung 5), not skip everything.
    cases = _cases("mapped", "newtest")  # "newtest" absent from map
    run, skip, _ = select_tests(cases, {"mapped": ["src/simulation/m_rhs.fpp"]}, {"toolchain/mfc/test/cases.py"})
    assert [c.coverage_key() for c in run] == ["newtest"]
    assert [c.coverage_key() for c in skip] == ["mapped"]


def test_unattributable_nonsource_change_runs_all():
    cases = _cases("a")
    for f in ("mfc.sh", "toolchain/pyproject.toml", "tests/ABC12345/golden.txt", ".github/workflows/test.yml"):
        run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, {f})
        assert len(run) == 1 and "run all" in reason, f


def test_docs_only_still_skips_all():
    cases = _cases("a")
    for f in ("README.md", "docs/foo.rst", "LICENSE", ".claude/x.md"):
        run, skip, reason = select_tests(cases, {"a": ["src/x.fpp"]}, {f})
        assert run == [] and "rung7" in reason, f


def test_load_map_rejects_malformed_entry(tmp_path):
    import gzip
    import json

    from mfc.test.coverage import load_map

    p = tmp_path / "m.json.gz"
    with gzip.open(p, "wt") as fh:
        json.dump({"_meta": {"built_at": "x"}, "good": ["a.fpp"], "bad": "not-a-list"}, fh)
    assert load_map(Path(p)) == (None, None)


def test_uppercase_fortran_extension_forces_all():
    cases = _cases("a")
    run, skip, reason = select_tests(cases, {"a": []}, {"src/common/m_x.F90"})
    assert len(run) == 1 and reason.startswith("rung3")


def test_toolchain_py_change_forces_all_except_cases():
    assert is_always_run_all({"toolchain/mfc/case.py"})
    assert is_always_run_all({"toolchain/mfc/build.py"})
    assert is_always_run_all({"toolchain/mfc/common.py"})
    assert not is_always_run_all({"toolchain/mfc/test/cases.py"})


def test_empty_map_with_fpp_change_runs_all_rung4():
    cases = _cases("a", "b")
    run, skip, reason = select_tests(cases, {}, {"src/simulation/m_rhs.fpp"})
    assert len(run) == 2 and reason.startswith("rung4")


# --- Map churn: a rebuilt map must be byte-comparable when nothing really changed ---


def test_coverage_key_is_independent_of_checkout_location():
    """The same test built from two checkouts must hash identically.

    to_case() absolutizes file-valued params, so the STL cases carry the runner's
    checkout path. Hashing that gave every runner a different key for the same test.
    """
    runner_a = "/home/runner/work/MFC/MFC"
    runner_b = "/storage/scratch/job1234/MFC"
    a = param_hash(canonicalize_param_paths({"m": 100, "stl_models(1)%model_filepath": f"{runner_a}/examples/2D_ibm_stl/model.stl"}, runner_a))
    b = param_hash(canonicalize_param_paths({"m": 100, "stl_models(1)%model_filepath": f"{runner_b}/examples/2D_ibm_stl/model.stl"}, runner_b))
    assert a == b


def test_canonicalize_preserves_paths_outside_the_repo():
    params = {"f": "/opt/shared/meshes/mesh.stl"}
    assert canonicalize_param_paths(params, "/home/runner/MFC") == params


def test_canonicalize_leaves_non_path_values_untouched():
    params = {"m": 100, "weno_order": 5, "bubbles_euler": "T"}
    assert canonicalize_param_paths(params, "/home/runner/MFC") == params


def test_canonicalize_relativizes_in_repo_name_that_starts_with_dots():
    """`..foo` is a leading-dots filename, not a parent-directory escape."""
    root = "/home/runner/MFC"
    assert canonicalize_param_paths({"f": f"{root}/..cache/model.stl"}, root) == {"f": "..cache/model.stl"}


def test_canonicalize_preserves_sibling_directory_sharing_a_prefix():
    """The near miss: relpath yields a genuine leading `..`, so the path stays absolute."""
    params = {"f": "/home/runner/MFC-scratch/model.stl"}
    assert canonicalize_param_paths(params, "/home/runner/MFC") == params


def test_canonicalized_key_still_distinguishes_different_files():
    root = "/home/runner/MFC"
    a = param_hash(canonicalize_param_paths({"f": f"{root}/examples/a/model.stl"}, root))
    b = param_hash(canonicalize_param_paths({"f": f"{root}/examples/b/model.stl"}, root))
    assert a != b


def test_save_map_writes_a_deterministic_gzip_header():
    """gzip stamps mtime + source filename by default, so identical maps differed on disk."""
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.json.gz"
        save_map(p, {"abc": ["src/simulation/m_rhs.fpp"]}, n_tests=1, git_sha="deadbee", gfortran_version="13")
        raw = p.read_bytes()
        assert raw[4:8] == b"\x00\x00\x00\x00", "gzip MTIME field must be zeroed"
        assert not raw[3] & 0x08, "gzip FNAME flag must be clear"


def test_entries_equal_ignores_coverage_list_order():
    assert entries_equal({"a": ["y.fpp", "x.fpp"]}, {"a": ["x.fpp", "y.fpp"]})


def test_entries_equal_detects_real_differences():
    assert not entries_equal({"a": ["x.fpp"]}, {"a": ["x.fpp", "y.fpp"]})
    assert not entries_equal({"a": ["x.fpp"]}, {"b": ["x.fpp"]})
    assert not entries_equal(None, {"a": ["x.fpp"]})


def test_health_quiet_repo_is_not_stale_when_map_is_current():
    """An old map is fine if nothing coverage-relevant landed since it was built."""
    ok, msg = map_health(
        meta={"built_at": "2026-05-01T00:00:00+00:00", "n_tests": 600},
        current_keys={"a"},
        mapped_keys={"a"},
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
        verified_after_last_change=True,
    )
    assert ok, msg


def test_health_fails_immediately_when_no_refresh_ran_since_last_source_change():
    """Detects a dead refresh on the next relevant push, not 10 days later."""
    ok, msg = map_health(
        meta={"built_at": "2026-05-28T00:00:00+00:00", "n_tests": 600},
        current_keys={"a"},
        mapped_keys={"a"},
        now="2026-05-29T00:00:00+00:00",
        max_age_days=10,
        min_fraction=0.8,
        verified_after_last_change=False,
    )
    assert not ok and "stale" in msg.lower()


# --- coverage_map_changed.py: the commit guard the refresh workflow branches on ---

CHANGED_SCRIPT = Path(__file__).resolve().parents[3] / ".github" / "scripts" / "coverage_map_changed.py"


def _env_without_git():
    """The environment minus every GIT_* variable.

    `git -C <dir>` changes directory but does NOT override an inherited GIT_DIR or
    GIT_INDEX_FILE. Git exports both when it runs a hook, and MFC's pre-commit hook runs
    precheck, which runs this suite -- so without this scrub the commits below are made
    against the real repository instead of the throwaway one.
    """
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def _repo_with_committed_map(d, entries):
    """A throwaway git repo whose HEAD holds `entries` as the coverage map."""
    repo = Path(d)
    env = _env_without_git()
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(repo)]
    subprocess.run([*git, "init", "-q"], check=True, env=env)
    save_map(repo / "tests" / "coverage_map.json.gz", entries, n_tests=len(entries), git_sha="aaa", gfortran_version="13")
    subprocess.run([*git, "add", "-A"], check=True, env=env)
    subprocess.run([*git, "commit", "-q", "--no-verify", "-m", "map"], check=True, env=env)
    return repo


def _run_guard(repo):
    return subprocess.run([sys.executable, str(CHANGED_SCRIPT)], cwd=repo, capture_output=True, text=True, check=False, env=_env_without_git()).returncode


def test_throwaway_repos_are_isolated_from_an_inherited_git_dir():
    """The helper above must not commit into whatever repo GIT_DIR names.

    Precheck runs this suite from the pre-commit hook, where git exports GIT_DIR and
    GIT_INDEX_FILE. Without the scrub, `git -C tmpdir commit` rewrites the developer's
    checked-out branch -- silently, while every test still reports as passing.
    """
    with patch.dict(os.environ, {"GIT_DIR": "/nonexistent.git", "GIT_INDEX_FILE": "/nonexistent.index"}):
        assert not [k for k in _env_without_git() if k.startswith("GIT_")]
        with tempfile.TemporaryDirectory() as d:
            repo = _repo_with_committed_map(d, {"k1": ["src/simulation/m_rhs.fpp"]})
            # The commit is in the throwaway repo, so it went nowhere near GIT_DIR.
            log = subprocess.run(["git", "-C", str(repo), "log", "--oneline"], capture_output=True, text=True, check=True, env=_env_without_git())
            assert log.stdout.strip().endswith("map")


def test_guard_reports_unchanged_with_a_code_that_is_not_the_crash_code():
    """A rebuild differing only in _meta must skip -- but never via exit 1.

    An uncaught exception exits 1, so if "unchanged" were 1 the workflow could not tell a
    no-op refresh from a broken comparison, and a dead guard would run permanently green.
    """
    with tempfile.TemporaryDirectory() as d:
        repo = _repo_with_committed_map(d, {"k1": ["src/simulation/m_rhs.fpp"]})
        # Same entries, fresh _meta -- exactly what every no-op refresh produces.
        save_map(repo / "tests" / "coverage_map.json.gz", {"k1": ["src/simulation/m_rhs.fpp"]}, n_tests=1, git_sha="bbb", gfortran_version="13")
        rc = _run_guard(repo)
    assert rc == 10
    assert rc not in (0, 1)


def test_guard_reports_changed_when_coverage_actually_moves():
    with tempfile.TemporaryDirectory() as d:
        repo = _repo_with_committed_map(d, {"k1": ["src/simulation/m_rhs.fpp"]})
        save_map(repo / "tests" / "coverage_map.json.gz", {"k1": ["src/simulation/m_rhs.fpp"], "k2": ["src/common/m_eos.fpp"]}, n_tests=2, git_sha="bbb", gfortran_version="13")
        assert _run_guard(repo) == 0


def test_guard_errors_when_the_rebuilt_map_is_missing():
    """The SLURM job died before writing a map: fail the job, do not read it as unchanged."""
    with tempfile.TemporaryDirectory() as d:
        repo = _repo_with_committed_map(d, {"k1": ["src/simulation/m_rhs.fpp"]})
        (repo / "tests" / "coverage_map.json.gz").unlink()
        rc = _run_guard(repo)
    assert rc == 2
    assert rc != 10


# --- check_coverage_map_health.py: the freshness signal the health workflow branches on ---

HEALTH_SCRIPT = Path(__file__).resolve().parents[3] / ".github" / "scripts" / "check_coverage_map_health.py"


def _health_module():
    """Import the health script by path; it is a script, not an installed module."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("check_coverage_map_health", HEALTH_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _repo_with_history(d):
    """A throwaway repo: one commit touching src/**/*.fpp, then one that does not."""
    repo = Path(d)
    env = _env_without_git()
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(repo)]
    subprocess.run([*git, "init", "-q", "-b", "master"], check=True, env=env)
    (repo / "src" / "simulation").mkdir(parents=True)
    (repo / "src" / "simulation" / "m_rhs.fpp").write_text("! v1\n")
    subprocess.run([*git, "add", "-A"], check=True, env=env)
    subprocess.run([*git, "commit", "-q", "--no-verify", "-m", "relevant"], check=True, env=env)
    relevant = subprocess.run([*git, "rev-parse", "HEAD"], capture_output=True, text=True, check=True, env=env).stdout.strip()
    (repo / "README.md").write_text("docs\n")
    subprocess.run([*git, "add", "-A"], check=True, env=env)
    subprocess.run([*git, "commit", "-q", "--no-verify", "-m", "irrelevant"], check=True, env=env)
    later = subprocess.run([*git, "rev-parse", "HEAD"], capture_output=True, text=True, check=True, env=env).stdout.strip()
    return repo, relevant, later


def _set_verified(repo, sha):
    subprocess.run(["git", "-C", str(repo), "update-ref", "refs/coverage-map/verified", sha], check=True, env=_env_without_git())


def test_verified_sha_is_none_when_the_ref_was_never_pushed():
    """A fork, or the window before the first refresh: undeterminable, not broken."""
    health = _health_module()
    with tempfile.TemporaryDirectory() as d:
        repo, _, _ = _repo_with_history(d)
        assert health.verified_sha(cwd=repo) is None
        # None must reach map_health as None, which falls back to the wall-clock rule.
        assert health.verified_after_last_change(health.verified_sha(cwd=repo), cwd=repo) is None


def test_verified_after_last_change_true_when_a_refresh_ran_since_the_change():
    health = _health_module()
    with tempfile.TemporaryDirectory() as d:
        repo, relevant, later = _repo_with_history(d)
        _set_verified(repo, later)
        assert health.verified_sha(cwd=repo) == later
        assert health.verified_after_last_change(later, cwd=repo) is True
        # The relevant commit itself counts: a refresh AT the change is current.
        assert health.verified_after_last_change(relevant, cwd=repo) is True


def test_verified_after_last_change_false_when_the_refresh_predates_the_change():
    """The genuine broken-refresh case this check exists to catch."""
    health = _health_module()
    with tempfile.TemporaryDirectory() as d:
        repo, relevant, _ = _repo_with_history(d)
        env = _env_without_git()
        git = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(repo)]
        before = subprocess.run([*git, "rev-parse", "HEAD~1"], capture_output=True, text=True, check=True, env=env).stdout.strip()
        assert before == relevant
        (repo / "src" / "simulation" / "m_rhs.fpp").write_text("! v2\n")
        subprocess.run([*git, "add", "-A"], check=True, env=env)
        subprocess.run([*git, "commit", "-q", "--no-verify", "-m", "relevant again"], check=True, env=env)
        _set_verified(repo, relevant)
        assert health.verified_after_last_change(relevant, cwd=repo) is False


def test_a_no_op_refresh_still_keeps_the_map_healthy():
    """The regression #1683 introduced: unchanged entries push no commit, so the map's
    _meta.git_sha stays behind while the refresh is working. The ref, not git_sha, is what
    the health check reads -- an old git_sha must not read as STALE."""
    health = _health_module()
    with tempfile.TemporaryDirectory() as d:
        repo, relevant, later = _repo_with_history(d)
        _set_verified(repo, later)
        ok, msg = map_health(
            meta={"built_at": "2026-05-28T00:00:00+00:00", "git_sha": "ancient", "n_tests": 1},
            current_keys={"a"},
            mapped_keys={"a"},
            now="2026-05-29T00:00:00+00:00",
            max_age_days=10,
            min_fraction=0.8,
            verified_after_last_change=health.verified_after_last_change(health.verified_sha(cwd=repo), cwd=repo),
        )
    assert ok, msg
