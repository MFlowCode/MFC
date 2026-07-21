#!/usr/bin/env python3
"""Manufactured verification for the generic Mie-Gruneisen EOS core in
src/common/m_eos_mie_gruneisen.fpp.

The Mie-Gruneisen backend is not yet selectable at runtime, so it has no golden
regression case; this standalone check is its verification instead. It mirrors the
Fortran formulas line for line and exercises the two properties that fix the generic
assembler: exact reduction to the existing stiffened-gas EOS, and agreement of the
analytic sound speed with a finite-difference of the frozen identity for a density
dependent reference curve.

Reference: Arienti, Morano and Shepherd, GALCIT FM99-8 (2004); Menikoff and Plohr,
Rev. Mod. Phys. 61, 75 (1989). Run with python3; exits nonzero on any failure.
"""

import sys


# --- generic Mie-Gruneisen assembler (mirrors m_eos_mie_gruneisen.fpp) --------------
def mg_pressure(gamma_mg, rho, e_int, p_ref, e_ref):
    return p_ref + gamma_mg * rho * (e_int - e_ref)


def mg_internal_energy(gamma_mg, rho, pres, p_ref, e_ref):
    return e_ref + (pres - p_ref) / (gamma_mg * rho)


def mg_sound_speed_sq(gamma_mg, rho, pres, p_ref, dp_ref, de_ref):
    return dp_ref - gamma_mg * rho * de_ref + ((1.0 + gamma_mg) * pres - p_ref) / rho


def stiffened_reference(gamma_s, pi_inf, rho):
    """Stiffened gas as a Mie-Gruneisen reference curve."""
    return gamma_s - 1.0, -pi_inf, 0.0, pi_inf / rho, -pi_inf / rho**2


# --- analytic stiffened-gas oracle --------------------------------------------------
def sg_pressure(gamma_s, rho, e_int, pi_inf):
    return (gamma_s - 1.0) * rho * e_int - gamma_s * pi_inf


def sg_sound_speed_sq(gamma_s, rho, pres, pi_inf):
    return gamma_s * (pres + pi_inf) / rho


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


# --- Test A: exact reduction to stiffened gas (ideal gas is the pi_inf = 0 sub-case) -
def test_reduction():
    # (gamma_s, pi_inf, rho, e): ideal air, water (stiffened), a stiff solid, monatomic
    cases = [
        (1.4, 0.0, 1.2, 2.5e5),
        (4.4, 6.0e8, 1000.0, 3.0e5),
        (3.0, 1.0e5, 50.0, 2.0e4),
        (1.667, 0.0, 0.1, 2.0e6),
    ]
    worst = 0.0
    ok = True
    for gamma_s, pi_inf, rho, e in cases:
        gamma_mg, p_ref, dp_ref, e_ref, de_ref = stiffened_reference(gamma_s, pi_inf, rho)
        p_mg = mg_pressure(gamma_mg, rho, e, p_ref, e_ref)
        e_back = mg_internal_energy(gamma_mg, rho, p_mg, p_ref, e_ref)
        c2_mg = mg_sound_speed_sq(gamma_mg, rho, p_mg, p_ref, dp_ref, de_ref)
        errs = (
            rel(p_mg, sg_pressure(gamma_s, rho, e, pi_inf)),
            rel(e_back, e),
            rel(c2_mg, sg_sound_speed_sq(gamma_s, rho, p_mg, pi_inf)),
        )
        worst = max(worst, *errs)
        ok = ok and all(x < 1.0e-13 for x in errs)
    print(f"[A] stiffened-gas reduction: worst rel err = {worst:.2e}  ->  {'PASS' if ok else 'FAIL'}")
    return ok


# --- Test B: analytic c^2 vs finite difference for a JWL-form isentrope reference ----
# Synthetic (non-calibrated) products isentrope, exercising a density-dependent p_ref.
def jwl_reference(rho, A, B, R1, R2, omega, rho0):
    V = rho0 / rho
    e1, e2 = A * pow(2.718281828459045, -R1 * V), B * pow(2.718281828459045, -R2 * V)
    p_ref = e1 + e2
    e_ref = e1 / (rho0 * R1) + e2 / (rho0 * R2)
    dp_ref = (V / rho) * (R1 * e1 + R2 * e2)
    de_ref = p_ref / rho**2
    return omega, p_ref, dp_ref, e_ref, de_ref


def test_soundspeed_fd():
    A, B, R1, R2, omega, rho0 = 5.0e11, 8.0e9, 4.5, 1.2, 0.3, 1600.0
    worst = 0.0
    ok = True
    for rho in (900.0, 1600.0, 2400.0):
        for e in (1.0e6, 4.0e6):

            def p_of(r, en):
                g, p_r, _, e_r, _ = jwl_reference(r, A, B, R1, R2, omega, rho0)
                return mg_pressure(g, r, en, p_r, e_r)

            g, p_r, dp_r, e_r, de_r = jwl_reference(rho, A, B, R1, R2, omega, rho0)
            pres = mg_pressure(g, rho, e, p_r, e_r)
            c2_analytic = mg_sound_speed_sq(g, rho, pres, p_r, dp_r, de_r)

            # frozen identity c^2 = dp/drho|_e + (p/rho^2) dp/de|_rho, dp/de|_rho = omega*rho
            h = rho * 1.0e-6
            dpdrho_e = (p_of(rho + h, e) - p_of(rho - h, e)) / (2.0 * h)
            c2_fd = dpdrho_e + (pres / rho**2) * (omega * rho)

            err = rel(c2_analytic, c2_fd)
            worst = max(worst, err)
            ok = ok and err < 1.0e-6
    print(f"[B] sound speed vs finite difference: worst rel err = {worst:.2e}  ->  {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    passed = test_reduction() & test_soundspeed_fd()
    print("ALL PASS" if passed else "FAILURE")
    sys.exit(0 if passed else 1)
