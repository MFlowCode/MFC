"""Unit tests for the pure helpers behind the FP-stability cancellation pass, its
fypp macro-expansion flagging, scale-free pass/fail, and Verrou discovery/install.

The Verrou subprocess machinery is exercised by the ./mfc.sh fp-stability CI job;
here we test only the pure functions that decide what to instrument and how to
label results, so they can run without Verrou or built binaries. We keep the tests
that pin a real behavioral contract or a subtle edge, not every micro-variation.
"""

from mfc.fp_stability_metrics import (
    _autodetect_compare,
    _cancellation_severity,
    _macro_context_in_lines,
    _sig_bits,
)

# --- fypp macro-expansion context detection (a #:for/#:def line maps to N instances) ---


def test_macro_context_inside_for_loop_body():
    lines = [
        "#:for i in [1, 2, 3]\n",
        "  q(${i}$) = a - b\n",
        "#:endfor\n",
    ]
    assert _macro_context_in_lines(lines, 2) == "#:for"


def test_macro_context_if_block_is_not_duplicating():
    # #:if selects code but does not duplicate it, so it must NOT be flagged.
    lines = [
        "#:if FOO\n",
        "  a = b - c\n",
        "#:endif\n",
    ]
    assert _macro_context_in_lines(lines, 2) is None


def test_macro_context_unbalanced_close_is_safe():
    # a stray #:endfor with an empty stack must not crash or misreport
    assert _macro_context_in_lines(["#:endfor\n", "  a = b - c\n"], 2) is None


# --- per-site cancellation severity (highest bit-threshold a site survives) ---


def test_cancellation_severity_takes_highest_surviving_threshold():
    level_sites = [
        (10, [("a.fpp", 1), ("b.fpp", 2)]),
        (20, [("a.fpp", 1)]),
        (30, [("a.fpp", 1)]),
    ]
    # a.fpp:1 survives to 30 bits; b.fpp:2 only at 10
    assert _cancellation_severity(level_sites) == {("a.fpp", 1): 30, ("b.fpp", 2): 10}


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


# --- scale-free pass/fail: significant bits retained ---


def test_sig_bits_is_scale_free():
    # same relative deviation -> same bits regardless of absolute magnitude
    assert abs(_sig_bits(1e-9, 1.0) - _sig_bits(1e-4, 1e5)) < 1e-9


def test_sig_bits_zero_scale_is_safe():
    # a zero/degenerate field scale must not divide-by-zero; report full precision
    assert _sig_bits(1e-12, 0.0) == 53.0


# --- report emitters: must survive the CI-only path without KeyError / regressions ---


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


def test_emit_annotations_cancellation_notes_fypp_ambiguity(tmp_path, monkeypatch, capsys):
    from mfc import fp_stability_report as report
    from mfc.fp_stability import _blank_result

    monkeypatch.setenv("GITHUB_ACTIONS", "1")
    r = _blank_result("demo")
    r.update(
        cancellation_locs=[("src/x/m_a.fpp", 5)],
        cancellation_bits={("src/x/m_a.fpp", 5): 40},
        cancellation_macro={("src/x/m_a.fpp", 5): "#:for"},
    )
    report._emit_github_annotations([r])
    out = capsys.readouterr().out
    assert "::notice" in out
    assert "multiple instances" in out  # fypp-expanded cancellation site flagged


# --- Verrou discovery: a bare/broken valgrind must read as "Verrou absent" ---


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
    monkeypatch.setattr(runners, "_has_verrou_tool", lambda _bin, _env=None: False)
    assert runners._find_verrou() == ""


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


def test_verrou_env_preserves_user_valgrind_lib(tmp_path, monkeypatch):
    from mfc import fp_stability_runners as runners

    (tmp_path / "libexec" / "valgrind").mkdir(parents=True)
    monkeypatch.setenv("VALGRIND_LIB", "/user/chosen/lib")
    env = runners._verrou_env(str(tmp_path / "bin" / "valgrind"))
    assert env["VALGRIND_LIB"] == "/user/chosen/lib"  # not clobbered


# --- auto-install hard-fail guard (a green bootstrap that produced no binary) ---


def test_install_verrou_raises_when_no_binary_appears(monkeypatch):
    import pytest

    from mfc import fp_stability as fps

    monkeypatch.setattr(fps.subprocess, "run", lambda *a, **k: type("R", (), {"returncode": 0})())
    monkeypatch.setattr(fps, "_find_verrou", lambda: "")
    with pytest.raises(fps.MFCException, match="no valgrind binary"):
        fps._install_verrou()
