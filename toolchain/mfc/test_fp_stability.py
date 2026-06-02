"""Unit tests for the pure helpers behind the FP-stability dd_line confirmation
pass (#1) and macro-expansion flagging (#2).

The Verrou subprocess machinery is exercised by the ./mfc.sh fp-stability CI job;
here we test only the pure functions that decide what to instrument and how to
label results, so they can run without Verrou or built binaries.
"""

from mfc.fp_stability import (
    MIN_SIG_BITS,
    _build_source_filter,
    _cancellation_by_file,
    _cancellation_severity,
    _confirm_decision,
    _macro_context_in_lines,
    _mark_cancellation,
    _rank_locs,
    _sig_bits,
    _stability_pass,
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


# --- cancellation-origin view: where cancellation concentrates ---


def test_cancellation_by_file_counts_and_sorts_by_density():
    locs = [
        ("src/simulation/m_weno.fpp", 10),
        ("m_weno.fpp", 20),
        ("a/m_riemann_solvers.fpp", 5),
    ]
    assert _cancellation_by_file(locs) == [("m_weno.fpp", 2), ("m_riemann_solvers.fpp", 1)]


def test_cancellation_by_file_breaks_ties_by_name():
    locs = [("z.fpp", 1), ("a.fpp", 2)]
    assert _cancellation_by_file(locs) == [("a.fpp", 1), ("z.fpp", 1)]


def test_cancellation_by_file_empty():
    assert _cancellation_by_file([]) == []


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


def test_stability_pass_uses_global_floor():
    # well-conditioned: ~46 bits >= floor
    assert _stability_pass(1e-14, 1.0, MIN_SIG_BITS) is True
    # catastrophic: deviation at field scale -> fails
    assert _stability_pass(0.5, 1.0, MIN_SIG_BITS) is False


def test_min_sig_bits_is_single_precision_floor():
    assert MIN_SIG_BITS == 24


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
