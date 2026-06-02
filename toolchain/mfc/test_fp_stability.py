"""Unit tests for the pure helpers behind the FP-stability dd_line confirmation
pass (#1) and macro-expansion flagging (#2).

The Verrou subprocess machinery is exercised by the ./mfc.sh fp-stability CI job;
here we test only the pure functions that decide what to instrument and how to
label results, so they can run without Verrou or built binaries.
"""

from mfc.fp_stability_metrics import (
    MIN_SIG_BITS,
    _autodetect_compare,
    _build_source_filter,
    _cancellation_severity,
    _confirm_decision,
    _digits_left,
    _macro_context_in_lines,
    _mark_cancellation,
    _rank_locs,
    _sig_bits,
    _statement_bounds_in_lines,
)

# --- #2: fypp macro-expansion context detection ---


def test_macro_context_none_outside_any_block():
    lines = [
        "subroutine s_foo()\n",
        "  a = b - c\n",
        "end subroutine\n",
    ]
    assert _macro_context_in_lines(lines, 2) is None


def test_macro_context_inside_for_loop_body():
    lines = [
        "#:for i in [1, 2, 3]\n",
        "  q(${i}$) = a - b\n",
        "#:endfor\n",
    ]
    assert _macro_context_in_lines(lines, 2) == "#:for"


def test_macro_context_if_block_is_not_duplicating():
    lines = [
        "#:if FOO\n",
        "  a = b - c\n",
        "#:endif\n",
    ]
    assert _macro_context_in_lines(lines, 2) is None


def test_macro_context_reports_innermost_duplicating_block():
    lines = [
        "#:def MACRO(x)\n",
        "  #:if cond\n",
        "    #:for j in range(3)\n",
        "      y = ${x}$ - z\n",
        "    #:endfor\n",
        "  #:endif\n",
        "#:enddef\n",
    ]
    assert _macro_context_in_lines(lines, 4) == "#:for"


def test_macro_context_balances_closers():
    lines = [
        "#:for i in [1, 2]\n",
        "  a = b - c\n",
        "#:endfor\n",
        "d = e - f\n",
    ]
    # line 4 is after the loop closed -> not in any duplicating block
    assert _macro_context_in_lines(lines, 4) is None


def test_macro_context_def_body_when_no_inner_loop():
    lines = [
        "#:def GEOM(n)\n",
        "  r = x - y\n",
        "#:enddef\n",
    ]
    assert _macro_context_in_lines(lines, 2) == "#:def"


def test_macro_context_block_and_call_are_duplicating():
    assert _macro_context_in_lines(["#:block B\n", "  a = b - c\n", "#:endblock\n"], 2) == "#:block"
    assert _macro_context_in_lines(["#:call M()\n", "  a = b - c\n", "#:endcall\n"], 2) == "#:call"


def test_macro_context_unbalanced_close_is_safe():
    # a stray #:endfor with an empty stack must not crash or misreport
    assert _macro_context_in_lines(["#:endfor\n", "  a = b - c\n"], 2) is None


# --- #1: building the symbol-correct --source filter from --gen-source output ---


def test_build_source_filter_keeps_matching_file_and_line_with_symbol():
    gen = [
        "m_riemann_solvers.fpp\t512\ts_hllc_riemann_solver\n",
        "m_riemann_solvers.fpp\t999\ts_other\n",
    ]
    suspects = [("src/simulation/m_riemann_solvers.fpp", 512, 512)]
    out = _build_source_filter(gen, suspects)
    assert out == ["m_riemann_solvers.fpp\t512\ts_hllc_riemann_solver\n"]


def test_build_source_filter_matches_inclusive_range():
    gen = [
        "m_foo.fpp\t10\tsym\n",
        "m_foo.fpp\t11\tsym\n",
        "m_foo.fpp\t12\tsym\n",
        "m_foo.fpp\t13\tsym\n",
    ]
    suspects = [("m_foo.fpp", 11, 12)]
    out = _build_source_filter(gen, suspects)
    assert out == ["m_foo.fpp\t11\tsym\n", "m_foo.fpp\t12\tsym\n"]


def test_build_source_filter_excludes_other_basenames():
    gen = ["m_bar.fpp\t5\tsym\n"]
    suspects = [("m_foo.fpp", 5, 5)]
    assert _build_source_filter(gen, suspects) == []


def test_build_source_filter_matches_on_basename_not_full_path():
    # gen-source emits a basename; dd_line locs are repo-relative paths.
    gen = ["m_foo.fpp\t5\tsym\n"]
    suspects = [("src/common/m_foo.fpp", 5, 5)]
    assert _build_source_filter(gen, suspects) == ["m_foo.fpp\t5\tsym\n"]


def test_build_source_filter_skips_malformed_lines():
    gen = ["garbage-no-tab\n", "m_foo.fpp\tnotanumber\tsym\n", "m_foo.fpp\t5\tsym\n"]
    suspects = [("m_foo.fpp", 5, 5)]
    assert _build_source_filter(gen, suspects) == ["m_foo.fpp\t5\tsym\n"]


# --- #1: confirmation decision ---


def test_confirm_decision_true_when_suspect_reproduces_deviation():
    # perturbing only the suspect lines yields >= dd_threshold deviation
    assert _confirm_decision(suspect_dev=1e-3, dd_threshold=1e-5) is True


def test_confirm_decision_false_when_suspect_is_inert():
    # suspect lines barely move the result -> attribution not reproduced
    assert _confirm_decision(suspect_dev=1e-9, dd_threshold=1e-5) is False


def test_confirm_decision_none_when_measurement_unavailable():
    assert _confirm_decision(suspect_dev=None, dd_threshold=1e-5) is None


# --- Tier 1: per-line confirmation ranking ---


def test_rank_locs_sorts_by_share_dev_descending():
    locs = [
        {"path": "a.fpp", "start": 1, "end": 1, "share_dev": 0.1},
        {"path": "b.fpp", "start": 2, "end": 2, "share_dev": 0.9},
    ]
    ranked = _rank_locs(locs, total=1.0)
    assert [loc["path"] for loc in ranked] == ["b.fpp", "a.fpp"]


def test_rank_locs_computes_share_as_fraction_of_total():
    locs = [{"path": "a.fpp", "start": 1, "end": 1, "share_dev": 0.25}]
    ranked = _rank_locs(locs, total=0.5)
    assert ranked[0]["share"] == 0.5


def test_rank_locs_share_none_when_total_nonpositive():
    locs = [{"path": "a.fpp", "start": 1, "end": 1, "share_dev": 0.25}]
    ranked = _rank_locs(locs, total=0.0)
    assert ranked[0]["share"] is None


def test_rank_locs_treats_missing_share_dev_as_zero_and_sorts_last():
    locs = [
        {"path": "a.fpp", "start": 1, "end": 1, "share_dev": None},
        {"path": "b.fpp", "start": 2, "end": 2, "share_dev": 0.3},
    ]
    ranked = _rank_locs(locs, total=1.0)
    assert [loc["path"] for loc in ranked] == ["b.fpp", "a.fpp"]


# --- Tier 1b: dd_line x cancellation cross-reference ---


def test_mark_cancellation_flags_loc_on_a_cancellation_line():
    locs = [{"path": "src/common/m_foo.fpp", "start": 10, "end": 12}]
    _mark_cancellation(locs, [("m_foo.fpp", 11)])
    assert locs[0]["cancellation"] is True


def test_mark_cancellation_false_when_no_site_in_range():
    locs = [{"path": "src/common/m_foo.fpp", "start": 10, "end": 12}]
    _mark_cancellation(locs, [("m_foo.fpp", 99)])
    assert locs[0]["cancellation"] is False


def test_mark_cancellation_matches_on_basename_not_full_path():
    locs = [{"path": "src/common/m_foo.fpp", "start": 5, "end": 5}]
    _mark_cancellation(locs, [("/abs/build/m_foo.fpp", 5)])
    assert locs[0]["cancellation"] is True


def test_mark_cancellation_false_for_different_basename():
    locs = [{"path": "m_foo.fpp", "start": 5, "end": 5}]
    _mark_cancellation(locs, [("m_bar.fpp", 5)])
    assert locs[0]["cancellation"] is False


# --- per-site cancellation severity (bits lost), from a threshold sweep ---


def test_cancellation_severity_takes_highest_surviving_threshold():
    level_sites = [
        (10, [("a.fpp", 1), ("b.fpp", 2)]),
        (20, [("a.fpp", 1)]),
        (30, [("a.fpp", 1)]),
    ]
    # a.fpp:1 survives to 30 bits; b.fpp:2 only at 10
    assert _cancellation_severity(level_sites) == {("a.fpp", 1): 30, ("b.fpp", 2): 10}


def test_cancellation_severity_empty():
    assert _cancellation_severity([]) == {}


# --- auto-detect which output files to compare (for a user case) ---


def test_autodetect_compare_picks_cons_at_latest_step():
    fns = [
        "cons.1.00.000000.dat",
        "cons.1.00.000050.dat",
        "cons.2.00.000050.dat",
        "prim.1.00.000050.dat",
    ]
    assert _autodetect_compare(fns) == ["cons.1.00.000050.dat", "cons.2.00.000050.dat"]


def test_autodetect_compare_falls_back_to_prim_when_no_cons():
    fns = ["prim.1.00.000010.dat", "prim.3.00.000010.dat"]
    assert _autodetect_compare(fns) == ["prim.1.00.000010.dat", "prim.3.00.000010.dat"]


def test_autodetect_compare_empty_when_no_field_output():
    assert _autodetect_compare(["indices.dat", "pre_time_data.dat", "foo.txt"]) == []


# --- scale-free pass/fail: significant bits retained ---


def test_sig_bits_relative_deviation():
    # max_dev/ref_scale = 1e-14 -> ~46.5 retained bits
    assert 46 < _sig_bits(1e-14, 1.0) < 47


def test_sig_bits_is_scale_free():
    # same relative deviation -> same bits regardless of absolute magnitude
    assert abs(_sig_bits(1e-9, 1.0) - _sig_bits(1e-4, 1e5)) < 1e-9


def test_sig_bits_zero_deviation_is_full_precision():
    assert _sig_bits(0.0, 1.0) == 53.0


def test_sig_bits_zero_scale_is_safe():
    assert _sig_bits(1e-12, 0.0) == 53.0


def test_sig_bits_deviation_at_scale_is_unstable():
    # deviation as large as the field -> <= 0 retained bits
    assert _sig_bits(1.0, 1.0) <= 0.0


def test_min_sig_bits_is_single_precision_floor():
    assert MIN_SIG_BITS == 24


def test_digits_left_full_and_clamped():
    assert 15.5 < _digits_left(0) < 16.0  # full double ~ 16 sig digits
    assert _digits_left(53) == 0.0
    assert _digits_left(60) == 0.0  # clamp: never negative


# --- Fortran line-continuation handling (correct-line labeling) ---


def test_statement_bounds_single_line():
    lines = ["  a = b - c\n"]
    assert _statement_bounds_in_lines(lines, 1) == (1, 1)


def test_statement_bounds_spans_continuation_from_first_line():
    lines = ["  poly = (s_cb(i+3) - s_cb(i+1)) * &\n", "         (s_cb(i+2) - s_cb(i))\n"]
    assert _statement_bounds_in_lines(lines, 1) == (1, 2)


def test_statement_bounds_from_middle_continuation_line():
    # a hit on the continuation fragment must resolve to the statement start
    lines = ["  x = a + &\n", "      b + &\n", "      c\n"]
    assert _statement_bounds_in_lines(lines, 2) == (1, 3)
    assert _statement_bounds_in_lines(lines, 3) == (1, 3)


def test_statement_bounds_ignores_ampersand_in_trailing_comment_logic():
    # a real continuation '&' before a trailing comment still continues
    lines = ["  x = a & ! note\n", "      + b\n"]
    assert _statement_bounds_in_lines(lines, 1) == (1, 2)


def test_statement_bounds_non_continuation_neighbors():
    lines = ["  x = 1\n", "  y = 2\n", "  z = 3\n"]
    assert _statement_bounds_in_lines(lines, 2) == (2, 2)


def test_statement_bounds_with_leading_ampersand_continuation():
    # the MFC WENO style: line ends with '&' and the next line *starts* with '&'
    lines = ["  beta = x**2 &\n", "       & + eps\n"]
    assert _statement_bounds_in_lines(lines, 1) == (1, 2)
    assert _statement_bounds_in_lines(lines, 2) == (1, 2)


# --- report emitters: must survive blank and populated result dicts (CI-only path) ---


def _emit_to_tmp(results, tmp_path, monkeypatch):
    """Run _emit_github_summary into a temp file under the GitHub-Actions env."""
    from mfc import fp_stability_report as report

    out = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(out))
    monkeypatch.setenv("GITHUB_ACTIONS", "1")
    report._emit_github_summary(results, 5)
    return out.read_text()


def test_emit_summary_survives_blank_result(tmp_path, monkeypatch):
    # the dict produced on the per-case error path must not KeyError the emitter
    from mfc.fp_stability import _blank_result

    text = _emit_to_tmp([_blank_result("x")], tmp_path, monkeypatch)
    assert "0 passed, 1 failed" in text


def test_emit_summary_populated_result(tmp_path, monkeypatch):
    from mfc.fp_stability import _blank_result

    r = _blank_result("demo")
    r.update(
        passed=False,
        max_dev=1e-9,
        sig_bits=30.0,
        float_proxy=1e-6,
        vprec=[(52, 1e-14), (23, float("inf"))],  # exercises the "crash" branch
        dd_line_locs=[{"path": "src/x/m_a.fpp", "start": 5, "end": 5, "macro": "#:for", "share": 0.4, "cancellation": True}],
        dd_line_confirmed=False,
        cancellation_locs=[("src/x/m_a.fpp", 5)],
        cancellation_bits={("src/x/m_a.fpp", 5): 40},
        float_max_locs=[("m_a.fpp", 9)],
    )
    text = _emit_to_tmp([r], tmp_path, monkeypatch)
    assert "💥 crash" in text and "digits lost" in text


def test_emit_annotations_downgrade_unconfirmed(tmp_path, monkeypatch, capsys):
    from mfc import fp_stability_report as report
    from mfc.fp_stability import _blank_result

    monkeypatch.setenv("GITHUB_ACTIONS", "1")
    r = _blank_result("demo")
    r.update(dd_line_locs=[{"path": "src/x/m_a.fpp", "start": 5, "end": 5, "macro": None, "share": 0.9, "cancellation": False}], dd_line_confirmed=False)
    report._emit_github_annotations([r])
    out = capsys.readouterr().out
    assert "::notice" in out and "::warning" not in out  # unconfirmed -> notice, not warning


# --- Verrou discovery: a bare system valgrind must read as "Verrou absent" ---


def test_find_verrou_prefers_verrou_home_candidate(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    vbin = tmp_path / "bin" / "valgrind"
    vbin.parent.mkdir(parents=True)
    vbin.write_text("#!/bin/sh\n")
    vbin.chmod(0o755)
    monkeypatch.setenv("VERROU_HOME", str(tmp_path))
    # The candidate must also verify as Verrou-enabled; stub that so the test
    # exercises precedence, not a real valgrind invocation.
    monkeypatch.setattr(runners, "_has_verrou_tool", lambda _bin, _env=None: True)
    assert runners._find_verrou() == str(vbin)


def test_find_verrou_rejects_broken_verrou_home_tree(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    # A valgrind exists at $VERROU_HOME but does not actually run the verrou tool
    # (broken/stale/non-Verrou): it must read as absent, not be returned.
    vbin = tmp_path / "bin" / "valgrind"
    vbin.parent.mkdir(parents=True)
    vbin.write_text("#!/bin/sh\n")
    vbin.chmod(0o755)
    monkeypatch.setenv("VERROU_HOME", str(tmp_path))
    monkeypatch.setattr(runners, "_has_verrou_tool", lambda _bin, _env=None: False)
    monkeypatch.setattr(runners.shutil, "which", lambda _name: None)
    assert runners._find_verrou() == ""


def test_find_verrou_rejects_non_verrou_path_valgrind(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    # VERROU_HOME has no valgrind; a plain valgrind is on PATH but lacks the tool.
    monkeypatch.setenv("VERROU_HOME", str(tmp_path))
    monkeypatch.setattr(runners.shutil, "which", lambda _name: "/usr/bin/valgrind")
    monkeypatch.setattr(runners, "_has_verrou_tool", lambda _bin: False)
    assert runners._find_verrou() == ""


def test_find_verrou_accepts_verrou_enabled_path_valgrind(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    monkeypatch.setenv("VERROU_HOME", str(tmp_path))
    monkeypatch.setattr(runners.shutil, "which", lambda _name: "/opt/verrou/bin/valgrind")
    monkeypatch.setattr(runners, "_has_verrou_tool", lambda _bin: True)
    assert runners._find_verrou() == "/opt/verrou/bin/valgrind"


def test_has_verrou_tool_reflects_exit_code(monkeypatch):
    from mfc import fp_stability_runners as runners

    class _R:
        def __init__(self, rc):
            self.returncode = rc

    monkeypatch.setattr(runners.subprocess, "run", lambda *a, **k: _R(0))
    assert runners._has_verrou_tool("/any/valgrind") is True
    monkeypatch.setattr(runners.subprocess, "run", lambda *a, **k: _R(1))
    assert runners._has_verrou_tool("/any/valgrind") is False

    def _boom(*a, **k):
        raise OSError("not executable")

    monkeypatch.setattr(runners.subprocess, "run", _boom)
    assert runners._has_verrou_tool("/stale/valgrind") is False


# --- env composition for relocated (prebuilt) Verrou trees ---


def test_verrou_env_sets_valgrind_lib_when_libexec_present(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    (tmp_path / "libexec" / "valgrind").mkdir(parents=True)
    monkeypatch.delenv("VALGRIND_LIB", raising=False)
    env = runners._verrou_env(str(tmp_path / "bin" / "valgrind"))
    assert env["VALGRIND_LIB"] == str(tmp_path / "libexec" / "valgrind")


def test_verrou_env_omits_valgrind_lib_when_libexec_absent(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    monkeypatch.delenv("VALGRIND_LIB", raising=False)
    env = runners._verrou_env(str(tmp_path / "bin" / "valgrind"))
    assert "VALGRIND_LIB" not in env


def test_verrou_env_preserves_user_valgrind_lib(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    (tmp_path / "libexec" / "valgrind").mkdir(parents=True)
    monkeypatch.setenv("VALGRIND_LIB", "/user/chosen/lib")
    env = runners._verrou_env(str(tmp_path / "bin" / "valgrind"))
    assert env["VALGRIND_LIB"] == "/user/chosen/lib"  # not clobbered


def test_dd_env_prepends_pythonpath_and_inherits_valgrind_lib(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    (tmp_path / "libexec" / "valgrind").mkdir(parents=True)
    monkeypatch.delenv("VALGRIND_LIB", raising=False)
    monkeypatch.setenv("PYTHONPATH", "/pre/existing")
    monkeypatch.setattr(runners, "_verrou_pythonpath", lambda _b: "/vg/site-packages/valgrind")
    env = runners._dd_env(str(tmp_path / "bin" / "valgrind"))
    assert env["PYTHONPATH"] == "/vg/site-packages/valgrind:/pre/existing"
    assert env["VALGRIND_LIB"] == str(tmp_path / "libexec" / "valgrind")


def test_dd_env_no_leading_colon_when_pythonpath_empty(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    monkeypatch.delenv("PYTHONPATH", raising=False)
    monkeypatch.setattr(runners, "_verrou_pythonpath", lambda _b: "/vg/valgrind")
    env = runners._dd_env(str(tmp_path / "bin" / "valgrind"))
    assert env["PYTHONPATH"] == "/vg/valgrind"  # no stray leading ':'


# --- auto-install hard-fail guards ---


def test_install_verrou_raises_when_bootstrap_fails(monkeypatch):
    import pytest

    from mfc import fp_stability as fps

    monkeypatch.setattr(fps.subprocess, "run", lambda *a, **k: type("R", (), {"returncode": 1})())
    with pytest.raises(fps.MFCException, match="Verrou install failed"):
        fps._install_verrou()


def test_install_verrou_raises_when_no_binary_appears(monkeypatch):
    import pytest

    from mfc import fp_stability as fps

    monkeypatch.setattr(fps.subprocess, "run", lambda *a, **k: type("R", (), {"returncode": 0})())
    monkeypatch.setattr(fps, "_find_verrou", lambda: "")
    with pytest.raises(fps.MFCException, match="no valgrind binary"):
        fps._install_verrou()
