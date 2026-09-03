"""The Mie-Gruneisen reference curves as m_variables_conversion computes them, for tests and validation cases."""

import math


def mg_reference(rho, rho0, c0, s, s2=0.0, s3=0.0):
    """p_ref, e_ref and their d/drho for the Hugoniot u_s = c0 + s u_p + s2 u_p^2 + s3 u_p^3.

    The linear curve is closed form; with s2 or s3 the particle velocity behind the shock solves
    u_s(u_p) mu = u_p (1 + mu) by Newton from the linear solution, and dp/dmu follows by implicit
    differentiation. Release (mu < 0) is linear for every fit.
    """
    mu = rho / rho0 - 1.0
    if mu < 0.0:
        p, dp_dmu = rho0 * c0**2 * mu, rho0 * c0**2
    elif s2 == 0.0 and s3 == 0.0:
        d = 1.0 - (s - 1.0) * mu
        p = rho0 * c0**2 * mu * (1.0 + mu) / d**2
        dp_dmu = rho0 * c0**2 * ((1.0 + 2.0 * mu) * d + 2.0 * (s - 1.0) * mu * (1.0 + mu)) / d**3
    else:
        up = c0 * mu / (1.0 - (s - 1.0) * mu)
        for _ in range(8):
            us, dus = c0 + s * up + s2 * up**2 + s3 * up**3, s + 2.0 * s2 * up + 3.0 * s3 * up**2
            up -= (us * mu - up * (1.0 + mu)) / (dus * mu - (1.0 + mu))
        us, dus = c0 + s * up + s2 * up**2 + s3 * up**3, s + 2.0 * s2 * up + 3.0 * s3 * up**2
        dup_dmu = (up - us) / (dus * mu - (1.0 + mu))
        p = rho0 * us * up
        dp_dmu = rho0 * ((dus * up + us) * dup_dmu)
    e = p * mu / (2.0 * rho0 * (1.0 + mu))
    de_dmu = (dp_dmu * mu * (1.0 + mu) + p) / (2.0 * rho0 * (1.0 + mu) ** 2)
    return p, e, dp_dmu / rho0, de_dmu / rho0


def jwl_reference(rho, rho0, a, b, r1, r2):
    """p_ref, e_ref and their d/drho for the JWL curve p = A exp(-R1 V) + B exp(-R2 V), V = rho0/rho."""
    V = rho0 / rho
    ea, eb = a * math.exp(-r1 * V), b * math.exp(-r2 * V)
    p = ea + eb
    return p, (ea / r1 + eb / r2) / rho0, (rho0 / rho**2) * (r1 * ea + r2 * eb), p / rho**2


def coefficients_from_curve(rho, curve, gruneisen, dgruneisen=0.0):
    """Gamma, Pi, dPi/drho, dGamma/drho from a reference curve (p, e, dp/drho, de/drho) and Gamma_G(rho)."""
    p, e, dp, de = curve
    return 1.0 / gruneisen, rho * e - p / gruneisen, e + rho * de - dp / gruneisen + p * dgruneisen / gruneisen**2, -dgruneisen / gruneisen**2


def eos_coefficients(rho, rho0, c0, s, gruneisen, gruneisen_a=0.0, s2=0.0, s3=0.0):
    """Hugoniot-referenced Mie-Gruneisen with Gamma_G = Gamma_0 + a mu: what s_eos_coefficients returns."""
    mu = rho / rho0 - 1.0
    return coefficients_from_curve(rho, mg_reference(rho, rho0, c0, s, s2, s3), gruneisen + gruneisen_a * mu, gruneisen_a / rho0)


def jwl_coefficients(rho, rho0, a, b, r1, r2, omega):
    return coefficients_from_curve(rho, jwl_reference(rho, rho0, a, b, r1, r2), omega)


def sound_speed(rho, pres, gamma, pi, dpi, dgamma=0.0):
    """c^2 = [((Gamma + 1) p + Pi)/rho - dPi/drho - p dGamma/drho]/Gamma, the frozen single-phase speed the solver uses."""
    return ((((gamma + 1.0) * pres + pi) / rho - dpi - pres * dgamma) / gamma) ** 0.5


def reference_isentrope(reference, gruneisen, rho, rho_start, p_start):
    """Pressure on the isentrope through (rho_start, p_start) when the reference curve is itself an isentrope
    (JWL, Vinet) and Gamma_G is constant: p_ref(rho) + C (rho/rho_start)^(1 + Gamma_G), exact."""
    C = p_start - reference(rho_start)[0]
    return reference(rho)[0] + C * (rho / rho_start) ** (1.0 + gruneisen)


def hugoniot_state(u_p, rho0, c0, s, s2=0.0, s3=0.0):
    """Shock speed and density behind a shock of particle velocity u_p driven into the reference state."""
    u_s = c0 + s * u_p + s2 * u_p**2 + s3 * u_p**3
    return u_s, rho0 * u_s / (u_s - u_p)


def rk4(slope, x0, y0, x1, steps):
    """Classical fixed-step RK4 for y' = slope(x, y) from x0 to x1: the one stepper the device code mirrors."""
    x, y, h = x0, y0, (x1 - x0) / steps
    for _ in range(steps):
        k1 = slope(x, y)
        k2 = slope(x + 0.5 * h, y + 0.5 * h * k1)
        k3 = slope(x + 0.5 * h, y + 0.5 * h * k2)
        k4 = slope(x + h, y + h * k3)
        y += h * (k1 + 2.0 * (k2 + k3) + k4) / 6.0
        x += h
    return y


def isentrope_rk4(coefficients, rho_from, p_from, rho_to, steps=8):
    """p after an isentropic density change, dp/drho = c^2: what f_phase_pressure_on_isentrope does."""
    return rk4(lambda r, p: sound_speed(r, p, *coefficients(r)) ** 2, rho_from, p_from, rho_to, steps)


def reference_temperature(reference, gruneisen, rho0, t0, cv, rho, steps=16):
    """T_ref(rho) along a reference curve: dT/dV = (de_ref/dV + p_ref)/cv - Gamma_G T/V, integrated in V from rho0.

    `reference(rho)` returns (p, e, dp/drho, de/drho) and `gruneisen(rho)` returns Gamma_G. This is the
    Maxwell relation (de/dV)_T = T (dp/dT)_V - p applied to e = e_ref + cv (T - T_ref); an isentrope
    reference (JWL) makes the first term vanish and gives the closed form T0 (rho/rho0)^Gamma_G.
    """

    def slope(V, T):
        r = 1.0 / V
        p, _, _, de_drho = reference(r)
        return (-(r**2) * de_drho + p) / cv - gruneisen(r) * T / V

    return rk4(slope, 1.0 / rho0, t0, 1.0 / rho, steps)


def temperature(rho, pres, reference, gruneisen, rho0, t0, cv):
    """T = T_ref(rho) + (p - p_ref)/(rho Gamma_G cv), the Gamma/Pi form solved for e - e_ref."""
    return reference_temperature(reference, gruneisen, rho0, t0, cv, rho) + (pres - reference(rho)[0]) / (rho * gruneisen(rho) * cv)


def vinet_reference(rho, rho0, k0, k0p):
    """p_ref, e_ref and their d/drho for the Vinet cold curve, x = (rho0/rho)^(1/3), eta = 1.5 (K0' - 1).

    p_c = 3 K0 (1 - x)/x^2 exp(eta (1 - x)); e_c = -int p_c dv from rho0, which integrates in closed form.
    The curve is an isentrope (an isotherm at 0 K), so de_c/drho = p_c/rho^2.
    """
    eta = 1.5 * (k0p - 1.0)
    x = (rho0 / rho) ** (1.0 / 3.0)
    ex = math.exp(eta * (1.0 - x))
    p = 3.0 * k0 * (1.0 - x) / x**2 * ex
    e = 9.0 * k0 / (rho0 * eta**2) * (1.0 - (1.0 - eta * (1.0 - x)) * ex)
    dx_drho = -x / (3.0 * rho)
    dp_dx = 3.0 * k0 * ex * ((-1.0) / x**2 - 2.0 * (1.0 - x) / x**3 - eta * (1.0 - x) / x**2)
    return p, e, dp_dx * dx_drho, p / rho**2
