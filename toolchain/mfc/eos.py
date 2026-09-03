"""The Mie-Gruneisen reference curve as m_variables_conversion computes it, for tests and validation cases."""


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


def eos_coefficients(rho, rho0, c0, s, gruneisen):
    """Gamma, Pi, dPi/drho: the three numbers s_eos_coefficients returns."""
    p, e, dp, de = mg_reference(rho, rho0, c0, s)
    return 1.0 / gruneisen, rho * e - p / gruneisen, e + rho * de - dp / gruneisen


def sound_speed(rho, pres, rho0, c0, s, gruneisen):
    """c^2 = [((Gamma + 1) p + Pi)/rho - dPi/drho]/Gamma, the frozen single-phase speed the solver uses."""
    gamma, pi, dpi = eos_coefficients(rho, rho0, c0, s, gruneisen)
    return ((((gamma + 1.0) * pres + pi) / rho - dpi) / gamma) ** 0.5


def hugoniot_state(u_p, rho0, c0, s):
    """Shock speed and density behind a shock of particle velocity u_p driven into the reference state."""
    u_s = c0 + s * u_p
    return u_s, rho0 * u_s / (u_s - u_p)
