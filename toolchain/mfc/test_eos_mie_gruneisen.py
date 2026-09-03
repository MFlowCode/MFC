"""Manufactured checks for the Mie-Gruneisen backend in m_variables_conversion.

The Fortran is s_eos_coefficients: a linear-Hugoniot reference curve mapped onto MFC's
rho e = Gamma p + Pi form. These tests pin the maths that mapping rests on, without a solver run.
"""

import pytest

from mfc.eos import eos_coefficients, jwl_coefficients, jwl_isentrope, jwl_reference, mg_reference

# Copper-like, no calibrated material: order-of-magnitude values only.
RHO0, C0, S, G0 = 8930.0, 3940.0, 1.49, 2.0
STATES = [(r, p) for r in (7500.0, 8930.0, 9800.0, 11000.0) for p in (1e8, 5e9, 2e10)]


def _fd(f, x, h):
    return (f(x + h) - f(x - h)) / (2.0 * h)


@pytest.mark.parametrize("rho", [7500.0, 8931.0, 9800.0, 11000.0])
def test_reference_curve_derivatives(rho):
    """The analytic dp_ref/drho and de_ref/drho match central differences of the curve itself."""
    _, _, dp, de = mg_reference(rho, RHO0, C0, S)
    h = rho * 1e-6
    dp_fd = _fd(lambda r: mg_reference(r, RHO0, C0, S)[0], rho, h)
    de_fd = _fd(lambda r: mg_reference(r, RHO0, C0, S)[1], rho, h)
    assert dp == pytest.approx(dp_fd, rel=1e-8)
    assert de == pytest.approx(de_fd, rel=1e-8)


@pytest.mark.parametrize("rho,pres", STATES)
def test_gamma_pi_form_is_the_mie_gruneisen_pressure(rho, pres):
    """rho e = Gamma p + Pi inverts to exactly p = p_ref + rho Gamma_G (e - e_ref)."""
    gamma, pi, _ = eos_coefficients(rho, RHO0, C0, S, G0)
    e = (gamma * pres + pi) / rho  # the energy MFC stores for this p
    p_ref, e_ref, _, _ = mg_reference(rho, RHO0, C0, S)
    p_mg = p_ref + rho * G0 * (e - e_ref)
    assert p_mg == pytest.approx(pres, rel=1e-12)


@pytest.mark.parametrize("rho,pres", STATES)
def test_sound_speed_matches_isentrope(rho, pres):
    """c^2 = [((Gamma+1)p + Pi)/rho - dPi/drho - p dGamma/drho]/Gamma equals (dp/drho)_s, integrated numerically."""
    gamma, pi, dpi = eos_coefficients(rho, RHO0, C0, S, G0)
    c2 = (((gamma + 1.0) * pres + pi) / rho - dpi) / gamma  # dGamma/drho = 0 for constant Gamma_G

    # Walk the isentrope: de = p/rho^2 drho, then re-evaluate p at the new state.
    e = (gamma * pres + pi) / rho
    h = rho * 1e-6

    def p_at(r, e_):
        g, pi_, _ = eos_coefficients(r, RHO0, C0, S, G0)
        return (r * e_ - pi_) / g

    c2_fd = (p_at(rho + h, e + pres / rho**2 * h) - p_at(rho - h, e - pres / rho**2 * h)) / (2.0 * h)
    assert c2 > 0.0
    assert c2 == pytest.approx(c2_fd, rel=1e-6)


def test_stiffened_gas_is_the_degenerate_member():
    """p_ref = -gamma pi_inf, e_ref = 0, Gamma_G = gamma - 1 reproduces MFC's stored gammas and pi_infs."""
    gam, pinf = 4.4, 6.0e8  # water-like stiffened gas
    Gamma_G = gam - 1.0
    p_ref, e_ref = -gam * pinf, 0.0
    rho = 1000.0
    Gamma = 1.0 / Gamma_G
    Pi = rho * e_ref - p_ref / Gamma_G
    assert Gamma == pytest.approx(1.0 / (gam - 1.0), rel=1e-15)  # what MFC stores as gammas(i)
    assert Pi == pytest.approx(gam * pinf / (gam - 1.0), rel=1e-15)  # what MFC stores as pi_infs(i)


def test_release_branch_is_c1_at_rho0():
    """Compression and release branches meet with matching p_ref and dp_ref at mu = 0."""
    eps = 1e-9 * RHO0
    p_up, _, dp_up, _ = mg_reference(RHO0 + eps, RHO0, C0, S)
    p_dn, _, dp_dn, _ = mg_reference(RHO0 - eps, RHO0, C0, S)
    assert p_up == pytest.approx(C0**2 * eps, rel=1e-6)  # linear through zero from above ...
    assert p_dn == pytest.approx(-(C0**2) * eps, rel=1e-6)  # ... and from below, same slope
    assert dp_up == pytest.approx(dp_dn, rel=1e-6)
    assert dp_up == pytest.approx(C0**2, rel=1e-6)  # bulk modulus rho0 c0^2 over rho0


# JWL, synthetic TNT-like magnitudes: no calibrated material.
JWL = dict(rho0=1630.0, a=3.712e11, b=3.231e9, r1=4.15, r2=0.95)
OMEGA = 0.30


@pytest.mark.parametrize("rho", [1630.0, 1200.0, 800.0, 400.0])
def test_jwl_curve_derivatives(rho):
    """The analytic dp_ref/drho and de_ref/drho match central differences, and de_ref/drho = p_ref/rho^2."""
    p, _, dp, de = jwl_reference(rho, **JWL)
    h = rho * 1e-6
    assert dp == pytest.approx(_fd(lambda r: jwl_reference(r, **JWL)[0], rho, h), rel=1e-8)
    assert de == pytest.approx(_fd(lambda r: jwl_reference(r, **JWL)[1], rho, h), rel=1e-8)
    assert de == pytest.approx(p / rho**2, rel=1e-14)


@pytest.mark.parametrize("rho", [1630.0, 1200.0, 800.0, 400.0])
def test_jwl_sound_speed_matches_isentrope(rho):
    pres = 1.2 * jwl_reference(rho, **JWL)[0] + 1e8
    gamma, pi, dpi = jwl_coefficients(rho, omega=OMEGA, **JWL)
    c2 = (((gamma + 1.0) * pres + pi) / rho - dpi) / gamma
    e = (gamma * pres + pi) / rho
    h = rho * 1e-6

    def p_at(r, e_):
        g, pi_, _ = jwl_coefficients(r, omega=OMEGA, **JWL)
        return (r * e_ - pi_) / g

    c2_fd = (p_at(rho + h, e + pres / rho**2 * h) - p_at(rho - h, e - pres / rho**2 * h)) / (2.0 * h)
    assert c2 > 0.0
    assert c2 == pytest.approx(c2_fd, rel=1e-6)


def test_jwl_isentrope_is_closed_form():
    """Integrating de = p/rho^2 drho through the Gamma/Pi form reproduces p_ref(V) + C V^-(omega+1)."""
    rho, p0 = JWL["rho0"], 2.0e10
    gamma, pi, _ = jwl_coefficients(rho, omega=OMEGA, **JWL)
    e = (gamma * p0 + pi) / rho
    n, r_end = 20000, 800.0
    dr = (rho - r_end) / n
    for _ in range(n):
        gamma, pi, _ = jwl_coefficients(rho, omega=OMEGA, **JWL)
        p = (rho * e - pi) / gamma
        e -= p / rho**2 * dr
        rho -= dr
    assert p == pytest.approx(jwl_isentrope(rho + dr, JWL["rho0"], JWL["a"], JWL["b"], JWL["r1"], JWL["r2"], OMEGA, JWL["rho0"], p0), rel=1e-4)
