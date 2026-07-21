"""Manufactured verification for src/common/m_eos_mie_gruneisen.fpp.

The Mie-Gruneisen backend is not yet runtime-selectable, so it has no golden case; these
tests mirror the Fortran formulas line for line and run in CI with the toolchain unit
tests. They check exact reduction to stiffened gas, the analytic sound speed and
thermodynamic derivatives against finite differences of the frozen identity, and the
isentrope self-consistency of the reference curves (Arienti et al., GALCIT FM99-8, 2004;
Menikoff and Plohr, Rev. Mod. Phys. 61, 75, 1989).
"""

from math import exp


# Generic Mie-Gruneisen assembler, mirroring m_eos_mie_gruneisen.fpp
def mg_pressure(gamma_mg, rho, e_int, p_ref, e_ref):
    return p_ref + gamma_mg * rho * (e_int - e_ref)


def mg_internal_energy(gamma_mg, rho, pres, p_ref, e_ref):
    return e_ref + (pres - p_ref) / (gamma_mg * rho)


def mg_dp_drho(gamma_mg, rho, pres, p_ref, dp_ref, de_ref):
    return dp_ref + (pres - p_ref) / rho - gamma_mg * rho * de_ref


def mg_dp_de(gamma_mg, rho):
    return gamma_mg * rho


def mg_sound_speed_sq(gamma_mg, rho, pres, p_ref, dp_ref, de_ref):
    return mg_dp_drho(gamma_mg, rho, pres, p_ref, dp_ref, de_ref) + pres / rho**2 * mg_dp_de(gamma_mg, rho)


def stiffened_reference(gamma_s, pi_inf, rho):
    """Stiffened gas as a Mie-Gruneisen reference curve."""
    return gamma_s - 1.0, -pi_inf, 0.0, pi_inf / rho, -pi_inf / rho**2


def jwl_reference(rho, A, B, R1, R2, omega, rho0):
    """Synthetic (non-calibrated) JWL-form principal isentrope, exercising a density
    dependent reference curve."""
    V = rho0 / rho
    e1, e2 = A * exp(-R1 * V), B * exp(-R2 * V)
    p_ref = e1 + e2
    e_ref = e1 / (rho0 * R1) + e2 / (rho0 * R2)
    dp_ref = (V / rho) * (R1 * e1 + R2 * e2)
    de_ref = p_ref / rho**2
    return omega, p_ref, dp_ref, e_ref, de_ref


# ideal air, water (stiffened), a stiff solid, monatomic: (gamma_s, pi_inf, rho, e)
SG_CASES = [
    (1.4, 0.0, 1.2, 2.5e5),
    (4.4, 6.0e8, 1000.0, 3.0e5),
    (3.0, 1.0e5, 50.0, 2.0e4),
    (1.667, 0.0, 0.1, 2.0e6),
]

JWL_PARAMS = (5.0e11, 8.0e9, 4.5, 1.2, 0.3, 1600.0)
JWL_STATES = [(rho, e) for rho in (900.0, 1600.0, 2400.0) for e in (1.0e6, 4.0e6)]


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


def test_stiffened_gas_reduction():
    """The generic assembler fed the stiffened reference curve must reproduce the
    analytic stiffened-gas pressure, energy round-trip, and sound speed exactly."""
    for gamma_s, pi_inf, rho, e in SG_CASES:
        gamma_mg, p_ref, dp_ref, e_ref, de_ref = stiffened_reference(gamma_s, pi_inf, rho)
        p_mg = mg_pressure(gamma_mg, rho, e, p_ref, e_ref)
        assert rel(p_mg, (gamma_s - 1.0) * rho * e - gamma_s * pi_inf) < 1.0e-13
        assert rel(mg_internal_energy(gamma_mg, rho, p_mg, p_ref, e_ref), e) < 1.0e-13
        c2 = mg_sound_speed_sq(gamma_mg, rho, p_mg, p_ref, dp_ref, de_ref)
        assert rel(c2, gamma_s * (p_mg + pi_inf) / rho) < 1.0e-13


def test_derivatives_and_sound_speed_vs_finite_difference():
    """Analytic dp/drho|_e, dp/de|_rho, and c^2 against central finite differences of
    p(rho, e) for the density-dependent JWL-form reference curve."""
    A, B, R1, R2, omega, rho0 = JWL_PARAMS

    def p_of(r, en):
        g, p_r, _, e_r, _ = jwl_reference(r, A, B, R1, R2, omega, rho0)
        return mg_pressure(g, r, en, p_r, e_r)

    for rho, e in JWL_STATES:
        g, p_ref, dp_ref, e_ref, de_ref = jwl_reference(rho, A, B, R1, R2, omega, rho0)
        pres = mg_pressure(g, rho, e, p_ref, e_ref)

        h = rho * 1.0e-6
        fd_dp_drho = (p_of(rho + h, e) - p_of(rho - h, e)) / (2.0 * h)
        assert rel(mg_dp_drho(g, rho, pres, p_ref, dp_ref, de_ref), fd_dp_drho) < 1.0e-6

        he = e * 1.0e-6
        fd_dp_de = (p_of(rho, e + he) - p_of(rho, e - he)) / (2.0 * he)
        assert rel(mg_dp_de(g, rho), fd_dp_de) < 1.0e-9

        c2_fd = fd_dp_drho + (pres / rho**2) * mg_dp_de(g, rho)
        assert rel(mg_sound_speed_sq(g, rho, pres, p_ref, dp_ref, de_ref), c2_fd) < 1.0e-6


def test_isentrope_self_consistency():
    """Both reference curves must satisfy de_ref/drho = p_ref/rho^2, the density form
    of de = -p dv along an isentrope; this is what makes them valid MG references and
    catches per-mass vs per-initial-volume scaling errors in e_ref."""
    A, B, R1, R2, omega, rho0 = JWL_PARAMS
    for rho in (900.0, 1600.0, 2400.0):
        _, p_ref, _, e_ref, de_ref = jwl_reference(rho, A, B, R1, R2, omega, rho0)
        assert rel(de_ref, p_ref / rho**2) < 1.0e-13
        h = rho * 1.0e-6
        _, _, _, e_hi, _ = jwl_reference(rho + h, A, B, R1, R2, omega, rho0)
        _, _, _, e_lo, _ = jwl_reference(rho - h, A, B, R1, R2, omega, rho0)
        assert rel((e_hi - e_lo) / (2.0 * h), de_ref) < 1.0e-6
    for gamma_s, pi_inf, rho, _ in SG_CASES:
        _, p_ref, _, _, de_ref = stiffened_reference(gamma_s, pi_inf, rho)
        assert rel(de_ref, p_ref / rho**2) < 1.0e-13


def test_admissibility_scan():
    """c^2 > 0 over the synthetic envelope for both reference curves."""
    A, B, R1, R2, omega, rho0 = JWL_PARAMS
    for rho, e in JWL_STATES:
        g, p_ref, dp_ref, e_ref, de_ref = jwl_reference(rho, A, B, R1, R2, omega, rho0)
        pres = mg_pressure(g, rho, e, p_ref, e_ref)
        assert mg_sound_speed_sq(g, rho, pres, p_ref, dp_ref, de_ref) > 0.0
    for gamma_s, pi_inf, rho, e in SG_CASES:
        g, p_ref, dp_ref, e_ref, de_ref = stiffened_reference(gamma_s, pi_inf, rho)
        # admissible states carry thermal energy on top of the reference energy;
        # below e_ref the stiffened gas is in tension past -pi_inf and c^2 < 0 is correct
        pres = mg_pressure(g, rho, e_ref + e, p_ref, e_ref)
        assert mg_sound_speed_sq(g, rho, pres, p_ref, dp_ref, de_ref) > 0.0
