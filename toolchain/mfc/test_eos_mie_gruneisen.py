"""Manufactured checks for the Mie-Gruneisen backend in m_variables_conversion.

The Fortran is s_eos_coefficients: a linear-Hugoniot reference curve mapped onto MFC's
rho e = Gamma p + Pi form. These tests pin the maths that mapping rests on, without a solver run.
"""

import pytest

from mfc.eos import eos_coefficients, isentrope_rk4, jwl_coefficients, jwl_reference, mg_reference, reference_isentrope, reference_temperature, sound_speed, temperature, vinet_reference

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
    gamma, pi, _, _ = eos_coefficients(rho, RHO0, C0, S, G0)
    e = (gamma * pres + pi) / rho  # the energy MFC stores for this p
    p_ref, e_ref, _, _ = mg_reference(rho, RHO0, C0, S)
    p_mg = p_ref + rho * G0 * (e - e_ref)
    assert p_mg == pytest.approx(pres, rel=1e-12)


@pytest.mark.parametrize("rho,pres", STATES)
def test_sound_speed_matches_isentrope(rho, pres):
    """c^2 = [((Gamma+1)p + Pi)/rho - dPi/drho - p dGamma/drho]/Gamma equals (dp/drho)_s, integrated numerically."""
    gamma, pi, dpi, _ = eos_coefficients(rho, RHO0, C0, S, G0)
    c2 = (((gamma + 1.0) * pres + pi) / rho - dpi) / gamma  # dGamma/drho = 0 for constant Gamma_G

    # Walk the isentrope: de = p/rho^2 drho, then re-evaluate p at the new state.
    e = (gamma * pres + pi) / rho
    h = rho * 1e-6

    def p_at(r, e_):
        g, pi_, _, _ = eos_coefficients(r, RHO0, C0, S, G0)
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
    gamma, pi, dpi, _ = jwl_coefficients(rho, omega=OMEGA, **JWL)
    c2 = (((gamma + 1.0) * pres + pi) / rho - dpi) / gamma
    e = (gamma * pres + pi) / rho
    h = rho * 1e-6

    def p_at(r, e_):
        g, pi_, _, _ = jwl_coefficients(r, omega=OMEGA, **JWL)
        return (r * e_ - pi_) / g

    c2_fd = (p_at(rho + h, e + pres / rho**2 * h) - p_at(rho - h, e - pres / rho**2 * h)) / (2.0 * h)
    assert c2 > 0.0
    assert c2 == pytest.approx(c2_fd, rel=1e-6)


def test_jwl_isentrope_is_closed_form():
    """Integrating de = p/rho^2 drho through the Gamma/Pi form reproduces p_ref(V) + C V^-(omega+1)."""
    rho, p0 = JWL["rho0"], 2.0e10
    gamma, pi, _, _ = jwl_coefficients(rho, omega=OMEGA, **JWL)
    e = (gamma * p0 + pi) / rho
    n, r_end = 20000, 800.0
    dr = (rho - r_end) / n
    for _ in range(n):
        gamma, pi, _, _ = jwl_coefficients(rho, omega=OMEGA, **JWL)
        p = (rho * e - pi) / gamma
        e -= p / rho**2 * dr
        rho -= dr
    assert p == pytest.approx(reference_isentrope(lambda r: jwl_reference(r, **JWL), OMEGA, rho + dr, JWL["rho0"], p0), rel=1e-4)


@pytest.mark.parametrize("rho,pres", STATES)
def test_density_dependent_gruneisen_sound_speed(rho, pres):
    """With Gamma_G = Gamma_0 + a mu the -p dGamma/drho term is live; c^2 must still be (dp/drho)_s."""
    a = 0.5
    gamma, pi, dpi, dgamma = eos_coefficients(rho, RHO0, C0, S, G0, a)
    assert dgamma != 0.0
    c2 = sound_speed(rho, pres, gamma, pi, dpi, dgamma) ** 2
    e = (gamma * pres + pi) / rho
    h = rho * 1e-6

    def p_at(r, e_):
        g, pi_, _, _ = eos_coefficients(r, RHO0, C0, S, G0, a)
        return (r * e_ - pi_) / g

    c2_fd = (p_at(rho + h, e + pres / rho**2 * h) - p_at(rho - h, e - pres / rho**2 * h)) / (2.0 * h)
    assert c2 == pytest.approx(c2_fd, rel=1e-6)


@pytest.mark.parametrize("xi", [0.8, 0.95, 1.05, 1.25])
def test_rk4_isentrope_matches_the_closed_form(xi):
    """Eight RK4 steps along dp/drho = c^2 reproduce the exact JWL isentrope to ~1e-8 over a 25% density change."""
    rho0, p0 = JWL["rho0"], 2.0e10
    p_rk4 = isentrope_rk4(lambda r: jwl_coefficients(r, omega=OMEGA, **JWL), rho0, p0, xi * rho0)
    p_exact = reference_isentrope(lambda r: jwl_reference(r, **JWL), OMEGA, xi * rho0, rho0, p0)
    assert p_rk4 == pytest.approx(p_exact, rel=1e-7)


CV, T0 = 1000.0, 300.0


@pytest.mark.parametrize("rho", [1630.0, 1200.0, 800.0])
def test_jwl_reference_temperature_is_closed_form(rho):
    """Along an isentrope reference the ODE collapses to T0 (rho/rho0)^omega."""
    T = reference_temperature(lambda r: jwl_reference(r, **JWL), lambda r: OMEGA, JWL["rho0"], T0, CV, rho)
    assert T == pytest.approx(T0 * (rho / JWL["rho0"]) ** OMEGA, rel=1e-7)


@pytest.mark.parametrize("rho,pres", [(8930.0, 5e9), (9800.0, 2e10), (8500.0, 1e8)])
def test_hugoniot_temperature_satisfies_the_maxwell_relation(rho, pres):
    """(de/dV)_T = T (dp/dT)_V - p for the complete EOS e = e_H + cv (T - T_H), p = p_H + Gamma_G cv (T - T_H)/V."""
    a = 0.5

    def ref(r):
        return mg_reference(r, RHO0, C0, S)

    def gru(r):
        return G0 + a * (r / RHO0 - 1.0)

    T = temperature(rho, pres, ref, gru, RHO0, T0, CV)

    def e_of(r, T_):
        return ref(r)[1] + CV * (T_ - reference_temperature(ref, gru, RHO0, T0, CV, r))

    def p_of(r, T_):
        return ref(r)[0] + gru(r) * r * CV * (T_ - reference_temperature(ref, gru, RHO0, T0, CV, r))

    V, h = 1.0 / rho, 1e-6 / rho
    de_dV = (e_of(1.0 / (V + h), T) - e_of(1.0 / (V - h), T)) / (2.0 * h)
    dp_dT = (p_of(rho, T * (1 + 1e-6)) - p_of(rho, T * (1 - 1e-6))) / (2e-6 * T)
    assert de_dV == pytest.approx(T * dp_dT - p_of(rho, T), rel=1e-5)


@pytest.mark.parametrize("rho", [8931.0, 9800.0, 11000.0])
def test_cubic_hugoniot_derivatives_and_linear_limit(rho):
    """With s2, s3 the Newton curve differentiates consistently; with s2 = s3 = 0 it is the linear closed form."""
    assert mg_reference(rho, RHO0, C0, S, 0.0, 0.0) == mg_reference(rho, RHO0, C0, S)
    s2, s3 = 0.1 / C0, 0.02 / C0**2
    p, e, dp, de = mg_reference(rho, RHO0, C0, S, s2, s3)
    h = rho * 1e-6
    assert dp == pytest.approx(_fd(lambda r: mg_reference(r, RHO0, C0, S, s2, s3)[0], rho, h), rel=1e-7)
    assert de == pytest.approx(_fd(lambda r: mg_reference(r, RHO0, C0, S, s2, s3)[1], rho, h), rel=1e-7)
    # the fit is stiffer than the linear one at the same compression
    assert p > mg_reference(rho, RHO0, C0, S)[0]


@pytest.mark.parametrize("rho", [8930.0, 9500.0, 11000.0, 8000.0])
def test_vinet_curve_derivatives(rho):
    """dp_c/drho matches central differences and the cold curve is an isentrope: de_c/drho = p_c/rho^2."""
    k0, k0p = 1.4e11, 5.0
    p, e, dp, de = vinet_reference(rho, RHO0, k0, k0p)
    h = rho * 1e-6
    assert dp == pytest.approx(_fd(lambda r: vinet_reference(r, RHO0, k0, k0p)[0], rho, h), rel=1e-7)
    assert de == pytest.approx(_fd(lambda r: vinet_reference(r, RHO0, k0, k0p)[1], rho, h), rel=1e-7)
    assert vinet_reference(RHO0, RHO0, k0, k0p)[0] == 0.0
