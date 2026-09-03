"""The Mie-Gruneisen reference curves as m_variables_conversion computes them, for tests and validation cases."""

import math


def mg_reference(rho, rho0, c0, s):
    """p_ref, e_ref and their d/drho for the linear-Hugoniot curve u_s = c0 + s u_p."""
    mu = rho / rho0 - 1.0
    if mu >= 0.0:
        d = 1.0 - (s - 1.0) * mu
        p = rho0 * c0**2 * mu * (1.0 + mu) / d**2
        dp_dmu = rho0 * c0**2 * ((1.0 + 2.0 * mu) * d + 2.0 * (s - 1.0) * mu * (1.0 + mu)) / d**3
    else:
        p = rho0 * c0**2 * mu
        dp_dmu = rho0 * c0**2
    e = p * mu / (2.0 * rho0 * (1.0 + mu))
    de_dmu = (dp_dmu * mu * (1.0 + mu) + p) / (2.0 * rho0 * (1.0 + mu) ** 2)
    return p, e, dp_dmu / rho0, de_dmu / rho0


def jwl_reference(rho, rho0, a, b, r1, r2):
    """p_ref, e_ref and their d/drho for the JWL curve p = A exp(-R1 V) + B exp(-R2 V), V = rho0/rho."""
    V = rho0 / rho
    ea, eb = a * math.exp(-r1 * V), b * math.exp(-r2 * V)
    p = ea + eb
    return p, (ea / r1 + eb / r2) / rho0, (rho0 / rho**2) * (r1 * ea + r2 * eb), p / rho**2


def coefficients_from_curve(rho, curve, gruneisen):
    """Gamma, Pi, dPi/drho from a reference curve (p, e, dp/drho, de/drho): what s_eos_coefficients returns."""
    p, e, dp, de = curve
    return 1.0 / gruneisen, rho * e - p / gruneisen, e + rho * de - dp / gruneisen


def eos_coefficients(rho, rho0, c0, s, gruneisen):
    return coefficients_from_curve(rho, mg_reference(rho, rho0, c0, s), gruneisen)


def jwl_coefficients(rho, rho0, a, b, r1, r2, omega):
    return coefficients_from_curve(rho, jwl_reference(rho, rho0, a, b, r1, r2), omega)


def sound_speed(rho, pres, gamma, pi, dpi):
    """c^2 = [((Gamma + 1) p + Pi)/rho - dPi/drho]/Gamma, the frozen single-phase speed the solver uses."""
    return ((((gamma + 1.0) * pres + pi) / rho - dpi) / gamma) ** 0.5


def jwl_isentrope(rho, rho0, a, b, r1, r2, omega, rho_start, p_start):
    """Pressure on the isentrope through (rho_start, p_start): p_ref(V) + C V^-(omega + 1), exact for JWL."""
    V, V0 = rho0 / rho, rho0 / rho_start
    C = (p_start - jwl_reference(rho_start, rho0, a, b, r1, r2)[0]) * V0 ** (omega + 1.0)
    return jwl_reference(rho, rho0, a, b, r1, r2)[0] + C * V ** (-(omega + 1.0))


def hugoniot_state(u_p, rho0, c0, s):
    """Shock speed and density behind a shock of particle velocity u_p driven into the reference state."""
    u_s = c0 + s * u_p
    return u_s, rho0 * u_s / (u_s - u_p)
