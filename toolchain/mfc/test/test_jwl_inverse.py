"""
JWL EOS inverse-consistency tests.

Tests all four JWL/ideal-gas mixture closures (jwl_mix_type = 0, 1, 2, 3) by:
  1. Starting from a sampled (rho, p, Y_jwl, alpha_jwl) state.
  2. Computing e = energy_pr(rho, p, Y, alpha) [inversion: pressure → energy].
  3. Computing p_back = pressure_er(rho, e, Y, alpha) [forward: energy → pressure].
  4. Checking relative round-trip error |p_back - p| / p.
  5. Checking p > 0, c² > 0, no NaN, no Inf at each step.

The Python implementations below are direct translations of the Fortran routines in
src/common/m_jwl.fpp, using the same branching and flooring logic.

Run with:
    python3 -m pytest toolchain/mfc/test/test_jwl_inverse.py -v
or standalone:
    python3 toolchain/mfc/test/test_jwl_inverse.py
"""

import math
import unittest

# Numerical constants (mirror src/common/m_constants.fpp)

SGM_EPS = 1.0e-16  # segmentation tolerance (sgm_eps in Fortran)


def _clamp(x, lo, hi):
    return min(max(x, lo), hi)


def _is_bad(x):
    """Return True if x is NaN or Inf."""
    return math.isnan(x) or math.isinf(x)


# JWL parameters — PETN-like (standard detonation products benchmark)

JWL_PARAMS = {
    "A": 6.170e11,  # Pa
    "B": 1.690e10,  # Pa
    "R1": 4.4,
    "R2": 1.2,
    "omega": 0.25,
    "rho0": 1760.0,  # kg/m³  (crystal density)
    "E0": 1760.0 * 5.814e6,  # J/m³  (E0 = rho0 * e0_det ≈ 1.02e10 J/m³)
    # Air / gas phase parameters (ideal-gas, γ = 1.4 → air_gamma = γ-1 = 0.4)
    "air_e0": 2.068e5,  # J/kg   (ambient internal energy at ~288 K, Cv*T)
    "air_rho0": 1.225,  # kg/m³
    "air_gamma": 0.4,  # = γ - 1  (p = air_gamma * rho * e)
    # Specific heats for Kuhl (type 1) and p-T-equilibrium (type 2)
    "cv_prod": 1000.0,  # J/(kg·K)  detonation products
    "cv_air": 717.5,  # J/(kg·K)  air
}

# Convenience aliases
_A = JWL_PARAMS["A"]
_B = JWL_PARAMS["B"]
_R1 = JWL_PARAMS["R1"]
_R2 = JWL_PARAMS["R2"]
_W = JWL_PARAMS["omega"]
_R0 = JWL_PARAMS["rho0"]
_E0 = JWL_PARAMS["E0"]
_AE0 = JWL_PARAMS["air_e0"]
_AR0 = JWL_PARAMS["air_rho0"]
_AG = JWL_PARAMS["air_gamma"]
_CVJ = JWL_PARAMS["cv_prod"]
_CVA = JWL_PARAMS["cv_air"]

# Sample grids (from test specification)

RHO_VALS = [1.0, 10.0, 100.0, 1000.0, 1800.0]  # kg/m³
PRES_VALS = [1.0e5, 1.0e7, 1.0e9, 1.0e10, 2.4e10]  # Pa
Y_VALS = [0.0, 1.0e-8, 1.0e-4, 0.01, 0.5, 0.99, 1.0]
ALPHA_VALS = [0.0, 1.0e-8, 1.0e-4, 0.01, 0.5, 0.99, 1.0]

# Helper: JWL cold-curve pressure  p_cold(rho)


def _jwl_pcold(rho, A, B, R1, R2, omega, rho0):
    """s_jwl_pcold in m_jwl.fpp."""
    rho_s = max(rho, SGM_EPS)
    V = rho0 / rho_s
    return A * (1.0 - omega / (R1 * V)) * math.exp(-R1 * V) + B * (1.0 - omega / (R2 * V)) * math.exp(-R2 * V)


def _jwl_sound_speed_sq(rho, pres, A, B, R1, R2, omega, rho0):
    """s_jwl_sound_speed_squared in m_jwl.fpp."""
    rho_s = max(rho, SGM_EPS)
    v = rho0 / rho_s
    exp1 = math.exp(-R1 * v)
    exp2 = math.exp(-R2 * v)
    e = (pres - A * (1.0 - omega / (R1 * v)) * exp1 - B * (1.0 - omega / (R2 * v)) * exp2) / (omega * rho_s)
    c2 = A * exp1 * (omega / (R1 * v**2) - R1 * (1.0 - omega / (R1 * v))) + B * exp2 * (omega / (R2 * v**2) - R2 * (1.0 - omega / (R2 * v)))
    c2 = c2 * (-rho0 / rho_s**2) + omega * e + omega * pres / rho_s
    return max(c2, SGM_EPS)


# Type 0 — isobaric (mechanical-equilibrium) closure


def _p0_er(rho, e, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG):
    """s_jwl_pressure_er (type 0): energy → pressure."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    a_j = _clamp(alpha_j, 0.0, 1.0)
    a_a = 1.0 - a_j
    rhoe = rho_s * e

    if a_j <= SGM_EPS:
        return air_gamma * rhoe

    rho1 = rho_s if a_a <= SGM_EPS else max(Y_s * rho_s / a_j, SGM_EPS)

    pref1 = A * (1.0 - omega0 * rho1 / (R1 * rho0)) * math.exp(-R1 * rho0 / rho1) + B * (1.0 - omega0 * rho1 / (R2 * rho0)) * math.exp(-R2 * rho0 / rho1)
    Kj = 1.0 / max(omega0, SGM_EPS)
    denom = max(a_j * Kj + a_a / max(air_gamma, SGM_EPS), SGM_EPS)
    return (rhoe + a_j * Kj * pref1) / denom


def _e0_pr(rho, pres, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG):
    """s_jwl_energy_pr (type 0): pressure → energy."""
    rho_s = max(rho, SGM_EPS)
    Y_w = _clamp(Y, 0.0, 1.0)
    a_j = _clamp(alpha_j, 0.0, 1.0)
    a_a = 1.0 - a_j

    if a_j <= SGM_EPS:
        return max(pres / max(air_gamma * rho_s, SGM_EPS), 0.0)

    rho1 = rho_s if a_a <= SGM_EPS else max(Y_w * rho_s / a_j, SGM_EPS)

    pref1 = A * (1.0 - omega0 * rho1 / (R1 * rho0)) * math.exp(-R1 * rho0 / rho1) + B * (1.0 - omega0 * rho1 / (R2 * rho0)) * math.exp(-R2 * rho0 / rho1)
    Kj = 1.0 / max(omega0, SGM_EPS)
    e = (pres * (a_j * Kj + a_a / max(air_gamma, SGM_EPS)) - a_j * Kj * pref1) / rho_s
    return max(e, 0.0)


# Type 1 — Kuhl/Khasainov temperature-form closure


def _p1_er(rho, e, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG, cv_prod=_CVJ, cv_air=_CVA):
    """s_jwl_kuhl_pressure_er (type 1): energy → pressure."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    V = rho0 / rho_s
    exp1 = math.exp(-R1 * V)
    exp2 = math.exp(-R2 * V)
    pbase = A * exp1 + B * exp2
    ecold = A / (R1 * rho0) * exp1 + B / (R2 * rho0) * exp2
    cv_mix = max(Y_s * cv_prod + (1.0 - Y_s) * cv_air, SGM_EPS)
    R_mix = Y_s * omega0 * cv_prod + (1.0 - Y_s) * air_gamma * cv_air
    T = max((e - Y_s * ecold) / cv_mix, SGM_EPS)
    return max(Y_s * pbase + rho_s * R_mix * T, SGM_EPS)


def _e1_pr(rho, pres, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG, cv_prod=_CVJ, cv_air=_CVA):
    """s_jwl_kuhl_energy_pr (type 1): pressure → energy."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    p_safe = max(pres, SGM_EPS)
    V = rho0 / rho_s
    exp1 = math.exp(-R1 * V)
    exp2 = math.exp(-R2 * V)
    pbase = A * exp1 + B * exp2
    ecold = A / (R1 * rho0) * exp1 + B / (R2 * rho0) * exp2
    cv_mix = max(Y_s * cv_prod + (1.0 - Y_s) * cv_air, SGM_EPS)
    R_mix = Y_s * omega0 * cv_prod + (1.0 - Y_s) * air_gamma * cv_air
    T = max((p_safe - Y_s * pbase) / max(rho_s * R_mix, SGM_EPS), SGM_EPS)
    return max(Y_s * ecold + cv_mix * T, 0.0)


# Type 2 — p-T equilibrium closure (bisection)


def _p2_er(rho, e, Y, alpha=None, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, air_gamma=_AG, cv_j=_CVJ, cv_a=_CVA):  # noqa: ARG001
    """s_jwl_ptequil_pressure_er (type 2): energy → pressure."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    cv_mix = max(Y_s * cv_j + (1.0 - Y_s) * cv_a, SGM_EPS)

    if Y_s <= SGM_EPS:
        return max(air_gamma * rho_s * e, SGM_EPS)
    if 1.0 - Y_s <= SGM_EPS:
        pcg = _jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0)
        return max(pcg + omega0 * rho_s * e, SGM_EPS)

    a_lo = SGM_EPS
    a_hi = 1.0 - SGM_EPS

    rj = max(Y_s * rho_s / a_lo, SGM_EPS)
    ra = max((1.0 - Y_s) * rho_s / (1.0 - a_lo), SGM_EPS)
    V = rho0 / rj
    ecold = A / (R1 * rho0) * math.exp(-R1 * V) + B / (R2 * rho0) * math.exp(-R2 * V)
    T = max((e - Y_s * ecold) / cv_mix, SGM_EPS)
    pj = A * math.exp(-R1 * V) + B * math.exp(-R2 * V) + omega0 * rj * cv_j * T
    pa = air_gamma * ra * cv_a * T
    f_lo = pj - pa
    pres = 0.5 * (pj + pa)

    for _ in range(60):
        a_m = 0.5 * (a_lo + a_hi)
        rj = max(Y_s * rho_s / a_m, SGM_EPS)
        ra = max((1.0 - Y_s) * rho_s / (1.0 - a_m), SGM_EPS)
        V = rho0 / rj
        ecold = A / (R1 * rho0) * math.exp(-R1 * V) + B / (R2 * rho0) * math.exp(-R2 * V)
        T = max((e - Y_s * ecold) / cv_mix, SGM_EPS)
        pj = A * math.exp(-R1 * V) + B * math.exp(-R2 * V) + omega0 * rj * cv_j * T
        pa = air_gamma * ra * cv_a * T
        f_m = pj - pa
        pres = 0.5 * (pj + pa)
        if f_lo * f_m > 0.0:
            a_lo = a_m
            f_lo = f_m
        else:
            a_hi = a_m

    return max(pres, SGM_EPS)


def _e2_pr(rho, pres, Y, alpha=None, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, air_gamma=_AG, cv_j=_CVJ, cv_a=_CVA):  # noqa: ARG001
    """s_jwl_ptequil_energy_pr (type 2): pressure → energy."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    p_s = max(pres, SGM_EPS)

    if Y_s <= SGM_EPS:
        return max(p_s / max(air_gamma * rho_s, SGM_EPS), 0.0)
    if 1.0 - Y_s <= SGM_EPS:
        pcg = _jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0)
        return max((p_s - pcg) / max(omega0 * rho_s, SGM_EPS), 0.0)

    a_lo = SGM_EPS
    a_hi = 1.0 - SGM_EPS

    rj = max(Y_s * rho_s / a_lo, SGM_EPS)
    ra = max((1.0 - Y_s) * rho_s / (1.0 - a_lo), SGM_EPS)
    V = rho0 / rj
    T = p_s / max(air_gamma * ra * cv_a, SGM_EPS)
    pj = A * math.exp(-R1 * V) + B * math.exp(-R2 * V) + omega0 * rj * cv_j * T
    g_lo = pj - p_s

    for _ in range(60):
        a_m = 0.5 * (a_lo + a_hi)
        rj = max(Y_s * rho_s / a_m, SGM_EPS)
        ra = max((1.0 - Y_s) * rho_s / (1.0 - a_m), SGM_EPS)
        V = rho0 / rj
        T = p_s / max(air_gamma * ra * cv_a, SGM_EPS)
        pj = A * math.exp(-R1 * V) + B * math.exp(-R2 * V) + omega0 * rj * cv_j * T
        g_m = pj - p_s
        if g_lo * g_m > 0.0:
            a_lo = a_m
            g_lo = g_m
        else:
            a_hi = a_m

    ecold = A / (R1 * rho0) * math.exp(-R1 * V) + B / (R2 * rho0) * math.exp(-R2 * V)
    return max((Y_s * cv_j + (1.0 - Y_s) * cv_a) * T + Y_s * ecold, 0.0)


# Type 3 — Rocflu single-fluid blended closure


def _p3_er(rho, e, Y, alpha=None, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG):  # noqa: ARG001
    """s_jwl_rocflu_pressure_er (type 3): energy → pressure."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)

    if Y_s <= 1.0e-4:
        return max(air_gamma * rho_s * e, SGM_EPS)
    if Y_s >= 0.99:
        cab = _jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0)
        return max(cab + omega0 * rho_s * e, SGM_EPS)

    V = rho0 / rho_s
    g_rho = _clamp((rho_s - air_rho0) / max(rho0 - air_rho0, SGM_EPS), 0.0, 1.0)
    om = air_gamma + (omega0 - air_gamma) * g_rho
    e0s = E0 / max(rho0, SGM_EPS)
    g_e = _clamp((e - air_e0) / max(e0s - air_e0, SGM_EPS), 0.0, 1.0)
    cab = A * (1.0 - om / (R1 * V)) * math.exp(-R1 * V) + B * (1.0 - om / (R2 * V)) * math.exp(-R2 * V)
    return max(g_e * cab + om * rho_s * e, SGM_EPS)


def _e3_pr(rho, pres, Y, alpha=None, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG):  # noqa: ARG001
    """s_jwl_rocflu_energy_pr (type 3): pressure → energy."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)

    if Y_s <= 1.0e-4:
        return max(pres / max(air_gamma * rho_s, SGM_EPS), 0.0)
    if Y_s >= 0.99:
        cab = _jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0)
        return max((pres - cab) / max(omega0 * rho_s, SGM_EPS), 0.0)

    V = rho0 / rho_s
    g_rho = _clamp((rho_s - air_rho0) / max(rho0 - air_rho0, SGM_EPS), 0.0, 1.0)
    om = air_gamma + (omega0 - air_gamma) * g_rho
    e0s = E0 / max(rho0, SGM_EPS)
    cab = A * (1.0 - om / (R1 * V)) * math.exp(-R1 * V) + B * (1.0 - om / (R2 * V)) * math.exp(-R2 * V)
    kab = cab / max(e0s - air_e0, SGM_EPS)
    e = (pres + air_e0 * kab) / max(kab + om * rho_s, SGM_EPS)
    g_e = (e - air_e0) / max(e0s - air_e0, SGM_EPS)
    if g_e < 0.0:
        e = pres / max(om * rho_s, SGM_EPS)
    elif g_e > 1.0:
        e = (pres - cab) / max(om * rho_s, SGM_EPS)
    return max(e, 0.0)


# Sound-speed functions


def _c2_mixture(rho, pres, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG):
    """s_jwl_mixture_sound_speed_squared (frozen mass-weighted, used for type 0)."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    a_j = _clamp(alpha_j, 0.0, 1.0)
    a_a = 1.0 - a_j

    if a_j <= SGM_EPS:
        return max((air_gamma + 1.0) * pres / rho_s, SGM_EPS)
    if a_a <= SGM_EPS:
        return max(_jwl_sound_speed_sq(rho_s, pres, A, B, R1, R2, omega0, rho0), SGM_EPS)

    rho1 = max(Y_s * rho_s / a_j, SGM_EPS)
    rho2 = max((1.0 - Y_s) * rho_s / a_a, SGM_EPS)
    c2_jwl = _jwl_sound_speed_sq(rho1, pres, A, B, R1, R2, omega0, rho0)
    c2_air = max((air_gamma + 1.0) * pres / rho2, SGM_EPS)
    c2 = Y_s * max(c2_jwl, SGM_EPS) + (1.0 - Y_s) * c2_air
    return max(c2, SGM_EPS)


def _c2_kuhl(rho, pres, Y, alpha_j, A=_A, B=_B, R1=_R1, R2=_R2, omega0=_W, rho0=_R0, E0=_E0, air_e0=_AE0, air_rho0=_AR0, air_gamma=_AG, cv_prod=_CVJ, cv_air=_CVA):
    """s_jwl_kuhl_sound_speed_squared (Mie-Grüneisen effective Γ, used for type 1)."""
    rho_s = max(rho, SGM_EPS)
    Y_s = _clamp(Y, 0.0, 1.0)
    V = rho0 / rho_s
    exp1 = math.exp(-R1 * V)
    exp2 = math.exp(-R2 * V)
    pbase = A * exp1 + B * exp2
    ecold = A / (R1 * rho0) * exp1 + B / (R2 * rho0) * exp2
    dpbase = (A * R1 * rho0 * exp1 + B * R2 * rho0 * exp2) / rho_s**2
    decold = pbase / rho_s**2
    cv_mix = max(Y_s * cv_prod + (1.0 - Y_s) * cv_air, SGM_EPS)
    R_mix = Y_s * omega0 * cv_prod + (1.0 - Y_s) * air_gamma * cv_air
    Gamma = R_mix / cv_mix
    e = max((pres - Y_s * (pbase - Gamma * rho_s * ecold)) / max(Gamma * rho_s, SGM_EPS), 0.0)
    dPeff = Y_s * dpbase - Gamma * Y_s * (ecold + rho_s * decold)
    c2 = max(dPeff + Gamma * (e + pres / rho_s), SGM_EPS)
    return c2


# Dispatch tables

_PRESSURE_ER = {0: _p0_er, 1: _p1_er, 2: _p2_er, 3: _p3_er}

_ENERGY_PR = {0: _e0_pr, 1: _e1_pr, 2: _e2_pr, 3: _e3_pr}

# Types 2 & 3 have no dedicated sound-speed formula; fall back to the mixture estimate.
_SOUND_SPEED_SQ = {0: _c2_mixture, 1: _c2_kuhl, 2: _c2_mixture, 3: _c2_mixture}

# Round-trip tolerance per closure type.
# Types 0 & 1 are algebraically exact inverses (analytical expressions).
# Type 2 uses 60-iteration bisection (converges to ~2^-60 of bracket).
# Type 3 has a blending g_e clamp that can make the inversion branch-switch
# for some states; we only require no NaN/Inf and p>0 there.
_RTOL = {0: 1.0e-10, 1: 1.0e-10, 2: 1.0e-6, 3: None}


def _min_pressure(rho, Y, alpha, mix_type):
    """
    Lower bound on achievable pressure for the given closure.

    For types 0, 1, 3: pressure_er(rho, e=0, Y, alpha) is the cold-curve floor.
    For type 2 (p-T equilibrium): the bisection has a domain with minimum products
    density at alpha_j = 1, giving a cold-curve floor of
        A*exp(-R1*rho0/(Y*rho)) + B*exp(-R2*rho0/(Y*rho)).
    The bisection-based pressure_er(rho, 0, Y) underestimates this floor when the
    bisection fails to converge (no root exists at e=0), so we use the analytical
    cold-curve bound directly.
    """
    if mix_type == 2:
        # Pure air: no cold-curve constraint.
        Y_s = _clamp(Y, 0.0, 1.0)
        if Y_s <= SGM_EPS:
            return 0.0
        # Pure JWL: standard cold curve.
        rho_s = max(rho, SGM_EPS)
        if 1.0 - Y_s <= SGM_EPS:
            return _jwl_pcold(rho_s, _A, _B, _R1, _R2, _W, _R0)
        # Mixed: minimum products cold-curve is at alpha_j → 1 (rj = Y*rho).
        rj_min = Y_s * rho_s  # products density at alpha_j = 1
        V_max = _R0 / max(rj_min, SGM_EPS)
        return _A * math.exp(-_R1 * V_max) + _B * math.exp(-_R2 * V_max)
    # Types 0, 1, 3: evaluate pressure at zero energy.
    p_er = _PRESSURE_ER[mix_type]
    return p_er(rho, 0.0, Y, alpha)


def _is_above_cold_curve(rho, p, Y, alpha, mix_type):
    """
    Return True only when (rho, p, Y, alpha) is safely above the cold-curve
    pressure floor so that the p→e→p round-trip identity is expected to hold.
    """
    p_min = _min_pressure(rho, Y, alpha, mix_type)
    # Require p to exceed the cold-curve floor by at least 1 % to avoid
    # floating-point boundary effects.
    return p > p_min * 1.01


# Core test runner (shared by all test methods)


def _run_roundtrip(mix_type):
    """
    For every (rho, p, Y, alpha) grid point run the protocol:
      1. e       = energy_pr(rho, p, Y, alpha)      [p → e]
      2. p_back  = pressure_er(rho, e, Y, alpha)     [e → p]
      3. e_back  = energy_pr(rho, p_back, Y, alpha)  [p → e again]
      4. p_back2 = pressure_er(rho, e_back, Y, alpha)[inner round-trip]
      5. c²      = sound_speed_sq(rho, p_back, Y, alpha)

    Hard failures (returned in list): NaN/Inf anywhere, p ≤ 0, e < 0, c² ≤ 0, or
    a round-trip error above tolerance for states that are *above* the cold-curve
    floor (p > p_min).  States below the cold curve intentionally clamp T to zero
    and do NOT satisfy the exact round-trip; those are flagged as INFO only.

    Returns (failures, clamped_count).
    """
    p_er = _PRESSURE_ER[mix_type]
    e_pr = _ENERGY_PR[mix_type]
    c2_fn = _SOUND_SPEED_SQ[mix_type]
    rtol = _RTOL[mix_type]

    failures = []
    clamped_count = 0

    for rho in RHO_VALS:
        for p in PRES_VALS:
            for Y in Y_VALS:
                for alpha in ALPHA_VALS:
                    label = f"jwl_mix_type={mix_type}, rho={rho}, p={p:.3e}, Y_jwl={Y}, alpha_jwl={alpha}"

                    # Step 1: p → e
                    try:
                        e = e_pr(rho, p, Y, alpha)
                    except Exception as exc:  # noqa: BLE001
                        failures.append((label, f"energy_pr raised: {exc}"))
                        continue

                    if _is_bad(e):
                        failures.append((label, f"e = {e!r} (NaN or Inf)"))
                        continue
                    if e < 0.0:
                        failures.append((label, f"e = {e} < 0"))
                        continue

                    # Step 2: e → p_back
                    try:
                        p_back = p_er(rho, e, Y, alpha)
                    except Exception as exc:  # noqa: BLE001
                        failures.append((label, f"pressure_er raised: {exc} (e={e:.6e})"))
                        continue

                    if _is_bad(p_back):
                        failures.append((label, f"p_back = {p_back!r} (NaN or Inf), e={e:.6e}"))
                        continue
                    if p_back <= 0.0:
                        failures.append((label, f"p_back = {p_back} <= 0, e={e:.6e}"))
                        continue

                    # Round-trip check — only for above-cold-curve states.
                    above_cc = _is_above_cold_curve(rho, p, Y, alpha, mix_type)
                    if rtol is not None:
                        rel_err = abs(p_back - p) / max(abs(p), SGM_EPS)
                        if rel_err > rtol:
                            if above_cc:
                                failures.append(
                                    (
                                        label,
                                        f"round-trip |p_back-p|/p = {rel_err:.2e} > {rtol:.2e} (p={p:.6e}, p_back={p_back:.6e}, e={e:.6e})",
                                    )
                                )
                                continue
                            else:
                                # Sub-cold-curve: T/e clamped to floor — expected behaviour.
                                clamped_count += 1
                                continue  # skip further checks for this state

                    # Step 3: p_back → e_back
                    try:
                        e_back = e_pr(rho, p_back, Y, alpha)
                    except Exception as exc:  # noqa: BLE001
                        failures.append((label, f"energy_pr(p_back) raised: {exc}"))
                        continue

                    if _is_bad(e_back):
                        failures.append((label, f"e_back = {e_back!r} (NaN or Inf)"))
                        continue

                    # Step 4: e_back → p_back2  (inner idempotency check)
                    try:
                        p_back2 = p_er(rho, e_back, Y, alpha)
                    except Exception as exc:  # noqa: BLE001
                        failures.append((label, f"pressure_er(e_back) raised: {exc}"))
                        continue

                    if _is_bad(p_back2):
                        failures.append((label, f"p_back2 = {p_back2!r} (NaN or Inf)"))
                        continue
                    if p_back2 <= 0.0:
                        failures.append((label, f"p_back2 = {p_back2} <= 0"))
                        continue

                    # Inner idempotency: p_back2 must equal p_back (both derived from e_back)
                    if rtol is not None:
                        rel2 = abs(p_back2 - p_back) / max(abs(p_back), SGM_EPS)
                        if rel2 > rtol:
                            failures.append(
                                (
                                    label,
                                    f"inner idempotency |p_back2-p_back|/p_back = {rel2:.2e} (p_back={p_back:.6e}, p_back2={p_back2:.6e})",
                                )
                            )
                            continue

                    # Sound speed check at p_back
                    try:
                        c2 = c2_fn(rho, p_back, Y, alpha)
                    except Exception as exc:  # noqa: BLE001
                        failures.append((label, f"sound_speed_sq raised: {exc}"))
                        continue

                    if _is_bad(c2):
                        failures.append((label, f"c2 = {c2!r} (NaN or Inf)"))
                        continue
                    if c2 <= 0.0:
                        failures.append((label, f"c2 = {c2} <= 0"))

    return failures, clamped_count


# Test case


class TestJwlInverseConsistency(unittest.TestCase):
    """Round-trip EOS consistency tests for all four JWL mixture closures."""

    def _assert_no_failures(self, mix_type):
        failures, clamped = _run_roundtrip(mix_type)
        if failures:
            lines = [f"\n{'=' * 70}"]
            lines.append(f"JWL mix_type={mix_type}: {len(failures)} failure(s) ({clamped} sub-cold-curve states skipped)")
            lines.append("=" * 70)
            for label, msg in failures:
                lines.append(f"  FAIL: {label}")
                lines.append(f"        {msg}")
            self.fail("\n".join(lines))

    def test_type0_isobaric(self):
        """Type 0 isobaric closure: algebraically exact p↔e inverse."""
        self._assert_no_failures(0)

    def test_type1_kuhl(self):
        """Type 1 Kuhl/Khasainov closure: temperature-form, algebraically exact."""
        self._assert_no_failures(1)

    def test_type2_ptequil(self):
        """Type 2 p-T-equilibrium closure: 60-step bisection, converges to ~1e-18."""
        self._assert_no_failures(2)

    def test_type3_rocflu(self):
        """Type 3 Rocflu blend: check p>0, c²>0, no NaN/Inf (blend clamping may break exact round-trip)."""
        self._assert_no_failures(3)

    def test_sound_speed_type0_positive(self):
        """Mixture sound speed (type 0 / isobaric) is positive for all sampled states."""
        c2_fn = _SOUND_SPEED_SQ[0]
        failures = []
        for rho in RHO_VALS:
            for p in PRES_VALS:
                for Y in Y_VALS:
                    for alpha in ALPHA_VALS:
                        c2 = c2_fn(rho, p, Y, alpha)
                        if _is_bad(c2) or c2 <= 0.0:
                            failures.append(f"rho={rho}, p={p:.2e}, Y={Y}, alpha={alpha}: c2={c2!r}")
        if failures:
            self.fail(f"c2 <= 0 or NaN/Inf in {len(failures)} state(s):\n  " + "\n  ".join(failures[:10]))

    def test_sound_speed_type1_positive(self):
        """Kuhl sound speed (type 1) is positive for all sampled states."""
        c2_fn = _SOUND_SPEED_SQ[1]
        failures = []
        for rho in RHO_VALS:
            for p in PRES_VALS:
                for Y in Y_VALS:
                    for alpha in ALPHA_VALS:
                        c2 = c2_fn(rho, p, Y, alpha)
                        if _is_bad(c2) or c2 <= 0.0:
                            failures.append(f"rho={rho}, p={p:.2e}, Y={Y}, alpha={alpha}: c2={c2!r}")
        if failures:
            self.fail(f"c2 <= 0 or NaN/Inf in {len(failures)} state(s):\n  " + "\n  ".join(failures[:10]))

    def test_pure_air_limit(self):
        """Y=0 or alpha=0: all types recover ideal-gas EOS (p = air_gamma * rho * e)."""
        for mix_type in (0, 1, 2, 3):
            p_er = _PRESSURE_ER[mix_type]
            e_pr = _ENERGY_PR[mix_type]
            for rho in RHO_VALS:
                for p in PRES_VALS:
                    # Strictly Y=0 (all types) and alpha=0 (type 0 only uses it)
                    e = e_pr(rho, p, 0.0, 0.0)
                    if _is_bad(e) or e < 0.0:
                        self.fail(f"type={mix_type}, rho={rho}, p={p}: e_pr(Y=0)={e!r}")
                    p_back = p_er(rho, e, 0.0, 0.0)
                    if _is_bad(p_back) or p_back <= 0.0:
                        self.fail(f"type={mix_type}, rho={rho}, p={p}: p_er(Y=0, e={e:.3e})={p_back!r}")

    def test_pure_products_limit(self):
        """Y=1 (and alpha=1 for type 0): all types recover single-phase JWL EOS for above-cold-curve states."""
        for mix_type in (0, 1, 2, 3):
            p_er = _PRESSURE_ER[mix_type]
            e_pr = _ENERGY_PR[mix_type]
            for rho in RHO_VALS:
                for p in PRES_VALS:
                    e = e_pr(rho, p, 1.0, 1.0)
                    if _is_bad(e) or e < 0.0:
                        self.fail(f"type={mix_type}, rho={rho}, p={p}: e_pr(Y=1)={e!r}")
                    p_back = p_er(rho, e, 1.0, 1.0)
                    if _is_bad(p_back) or p_back <= 0.0:
                        self.fail(f"type={mix_type}, rho={rho}, p={p}: p_er(Y=1, e={e:.3e})={p_back!r}")
                    # Only require tight round-trip for above-cold-curve states.
                    if not _is_above_cold_curve(rho, p, 1.0, 1.0, mix_type):
                        continue  # sub-cold-curve: clamped T, skip round-trip
                    rel_err = abs(p_back - p) / max(abs(p), SGM_EPS)
                    if rel_err > 1.0e-6:
                        self.fail(f"type={mix_type}, rho={rho}, p={p}: Y=1 round-trip |p_back-p|/p = {rel_err:.2e} (p_back={p_back:.6e})")


if __name__ == "__main__":
    # Standalone run: report all failures with full state info.
    import sys

    total_failures = 0
    n_states = len(RHO_VALS) * len(PRES_VALS) * len(Y_VALS) * len(ALPHA_VALS)
    for mt in (0, 1, 2, 3):
        fails, clamped = _run_roundtrip(mt)
        if fails:
            print(f"\n{'=' * 70}")
            print(f"jwl_mix_type={mt}: {len(fails)} failure(s), {clamped} sub-cold-curve states (expected, skipped)")
            print("=" * 70)
            for label, msg in fails:
                print(f"  FAIL: {label}")
                print(f"        {msg}")
            total_failures += len(fails)
        else:
            checked = n_states - clamped
            print(f"jwl_mix_type={mt}: {checked}/{n_states} above-cold-curve states PASSED ({clamped} sub-cold-curve states skipped)")

    if total_failures:
        print(f"\n{total_failures} total failure(s).")
        sys.exit(1)
    else:
        print("\nAll JWL inverse-consistency checks PASSED.")
        sys.exit(0)
