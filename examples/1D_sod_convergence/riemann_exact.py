"""Exact Riemann solver for the 1D Euler equations of an ideal gas (Toro, Riemann Solvers and Numerical Methods, ch. 4)."""

import numpy as np


def _f(p, rho, pk, ck, g):
    """Pressure function f_K and its derivative across the left or right wave."""
    if p > pk:  # shock
        a, b = 2.0 / ((g + 1.0) * rho), (g - 1.0) / (g + 1.0) * pk
        q = np.sqrt(a / (p + b))
        return (p - pk) * q, q * (1.0 - 0.5 * (p - pk) / (p + b))
    r = (p / pk) ** ((g - 1.0) / (2.0 * g))  # rarefaction
    return 2.0 * ck / (g - 1.0) * (r - 1.0), (p / pk) ** (-(g + 1.0) / (2.0 * g)) / (rho * ck)


def star_state(left, right, g):
    """Star-region pressure and velocity; left/right are (rho, u, p)."""
    (rl, ul, pl), (rr, ur, pr) = left, right
    cl, cr = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
    if 2.0 * (cl + cr) / (g - 1.0) <= ur - ul:
        raise ValueError("initial states generate vacuum")
    p = max(1e-12, 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr))  # PVRS guess
    for _ in range(100):
        fl, dl = _f(p, rl, pl, cl, g)
        fr, dr = _f(p, rr, pr, cr, g)
        p_new = max(1e-14, p - (fl + fr + ur - ul) / (dl + dr))
        if abs(p_new - p) < 1e-15 * (p_new + p):
            p = p_new
            break
        p = p_new
    fl, _ = _f(p, rl, pl, cl, g)
    fr, _ = _f(p, rr, pr, cr, g)
    return p, 0.5 * (ul + ur) + 0.5 * (fr - fl)


def sample(left, right, g, xi):
    """(rho, u, p) at the similarity coordinates xi = (x - x0) / t."""
    p_s, u_s = star_state(left, right, g)
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    out = np.empty((3, xi.size))
    for side, (rk, uk, pk), sgn in (("L", left, -1.0), ("R", right, 1.0)):
        ck = np.sqrt(g * pk / rk)
        sel = xi < u_s if side == "L" else xi >= u_s
        x = xi[sel]
        res = np.empty((3, x.size))
        if p_s > pk:  # shock
            r_s = rk * (p_s / pk + (g - 1.0) / (g + 1.0)) / ((g - 1.0) / (g + 1.0) * p_s / pk + 1.0)
            s = uk + sgn * ck * np.sqrt((g + 1.0) / (2.0 * g) * p_s / pk + (g - 1.0) / (2.0 * g))
            outside = sgn * (x - s) > 0.0
            res[:, outside] = np.array([[rk], [uk], [pk]])
            res[:, ~outside] = np.array([[r_s], [u_s], [p_s]])
        else:  # rarefaction
            r_s = rk * (p_s / pk) ** (1.0 / g)
            c_s = ck * (p_s / pk) ** ((g - 1.0) / (2.0 * g))
            head, tail = uk + sgn * ck, u_s + sgn * c_s
            outside = sgn * (x - head) > 0.0
            inside = sgn * (x - tail) < 0.0
            fan = ~outside & ~inside
            res[:, outside] = np.array([[rk], [uk], [pk]])
            res[:, inside] = np.array([[r_s], [u_s], [p_s]])
            xf = x[fan]
            u = 2.0 / (g + 1.0) * (-sgn * ck + 0.5 * (g - 1.0) * uk + xf)
            c = sgn * (xf - u)
            res[0, fan] = rk * (c / ck) ** (2.0 / (g - 1.0))
            res[1, fan] = u
            res[2, fan] = pk * (c / ck) ** (2.0 * g / (g - 1.0))
        out[:, sel] = res
    return out


def cell_averages(left, right, g, x0, t, edges, sub=64):
    """Exact cell averages of (rho, u, p) on cells with the given edges (midpoint rule on `sub` sub-cells)."""
    n = edges.size - 1
    frac = (np.arange(sub) + 0.5) / sub
    xs = edges[:-1, None] + np.diff(edges)[:, None] * frac[None, :]
    w = sample(left, right, g, (xs.ravel() - x0) / t).reshape(3, n, sub)
    return w.mean(axis=2)


def max_wave_speed(left, right, g):
    """Largest |signal speed| of the exact solution."""
    p_s, u_s = star_state(left, right, g)
    speeds = [abs(u_s)]
    for (rk, uk, pk), sgn in ((left, -1.0), (right, 1.0)):
        ck = np.sqrt(g * pk / rk)
        if p_s > pk:
            speeds.append(abs(uk + sgn * ck * np.sqrt((g + 1.0) / (2.0 * g) * p_s / pk + (g - 1.0) / (2.0 * g))))
        else:
            speeds += [abs(uk + sgn * ck), abs(u_s + sgn * ck * (p_s / pk) ** ((g - 1.0) / (2.0 * g)))]
    return max(speeds)
