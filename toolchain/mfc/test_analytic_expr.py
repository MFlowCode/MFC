"""Tests for AST-based analytic expression translation."""

import pytest

IC_MAP = {
    "x": "x_cc(i)",
    "y": "y_cc(j)",
    "z": "z_cc(k)",
    "lx": "patch_icpp(1)%length_x",
    "ly": "patch_icpp(1)%length_y",
    "r": "patch_icpp(1)%radius",
    "tau_e": "patch_icpp(1)%tau_e",
}


def test_simple_expression():
    from mfc.analytic_expr import fortranize_expr

    out = fortranize_expr("0.5 + 0.2*sin(2.0*pi*x/lx)", IC_MAP)
    assert out.replace(" ", "") == "0.5+0.2*sin(2.0*pi*x_cc(i)/patch_icpp(1)%length_x)"


def test_parentheses_regenerated_by_precedence():
    from mfc.analytic_expr import fortranize_expr

    out = fortranize_expr("2.0 * (x / lx + y / ly)", IC_MAP)
    assert out.replace(" ", "") == "2.0*(x_cc(i)/patch_icpp(1)%length_x+y_cc(j)/patch_icpp(1)%length_y)"


def test_power_and_unary_minus():
    from mfc.analytic_expr import fortranize_expr

    assert fortranize_expr("x**2", IC_MAP).replace(" ", "") == "x_cc(i)**2"
    assert fortranize_expr("-x + y", IC_MAP).replace(" ", "") == "-x_cc(i)+y_cc(j)"
    assert fortranize_expr("y * -x", IC_MAP).replace(" ", "") == "y_cc(j)*(-x_cc(i))"


def test_underscore_variable_maps_correctly():
    # The legacy regex could never match tau_e (it tokenized through the underscore,
    # corrupting the trailing e); the AST path must handle it.
    from mfc.analytic_expr import fortranize_expr

    out = fortranize_expr("tau_e * 2", IC_MAP)
    assert out.replace(" ", "") == "patch_icpp(1)%tau_e*2"


def test_scientific_notation_is_safe():
    # Legacy regex corrupted the 'e' in 2e-3 inside larger expressions.
    from mfc.analytic_expr import fortranize_expr

    out = fortranize_expr("1.0 + 2e-3*x", IC_MAP)
    assert out.replace(" ", "") == "1.0+0.002*x_cc(i)"


def test_unknown_variable_raises_listing_valid():
    from mfc.analytic_expr import AnalyticExprError, fortranize_expr

    with pytest.raises(AnalyticExprError) as excinfo:
        fortranize_expr("x0 + 1", IC_MAP)
    assert "x0" in str(excinfo.value)
    assert "lx" in str(excinfo.value)  # the message lists the valid names


def test_unknown_function_raises():
    from mfc.analytic_expr import AnalyticExprError, fortranize_expr

    with pytest.raises(AnalyticExprError, match="frobnicate"):
        fortranize_expr("frobnicate(x)", IC_MAP)


def test_syntax_error_raises():
    from mfc.analytic_expr import AnalyticExprError, fortranize_expr

    with pytest.raises(AnalyticExprError, match="invalid"):
        fortranize_expr("0.5 + * 2", IC_MAP)


def test_intrinsics_and_pi_pass_through():
    from mfc.analytic_expr import fortranize_expr

    out = fortranize_expr("exp(-x**2) + cos(pi*y) + min(x, y)", IC_MAP)
    assert out.replace(" ", "") == "exp(-x_cc(i)**2)+cos(pi*y_cc(j))+min(x_cc(i),y_cc(j))"


def test_exponent_literal_gets_wp_suffix():
    from mfc.analytic_expr import fortranize_expr

    assert fortranize_expr("1e20", IC_MAP) == "1e+20_wp"
    assert fortranize_expr("1e20 * x", IC_MAP).replace(" ", "") == "1e+20_wp*x_cc(i)"
