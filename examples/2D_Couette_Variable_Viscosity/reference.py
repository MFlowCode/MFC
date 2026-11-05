#!/usr/bin/env python3
"""Exact steady solution of variable-viscosity planar Couette flow.

The incompressible steady state satisfies two coupled ODEs across the gap
0 <= y <= H:

    momentum:   d/dy [ mu(T) du/dy ] = 0          (uniform shear stress)
    energy:     d/dy [ k  dT/dy ] + mu(T) (du/dy)^2 = 0   (conduction + heating)

with mu(T) = exp(C + D/T), and boundary conditions

    u(0) = 0,  u(H) = U,   T(0) = T0,  T(H) = T1.

Writing the momentum flux tau = mu du/dy and heat flux q = k dT/dy as states
turns this into a first-order BVP that scipy.solve_bvp resolves to ~1e-10:

    u'   = tau / mu(T)
    tau' = 0
    T'   = q / k
    q'   = -tau^2 / mu(T)        (= -mu (du/dy)^2, the dissipation source)

The Brinkman-number -> 0 limit (negligible viscous heating) has the closed form
T(y) linear and u(y) = U * I(y)/I(H) with I(y) = integral_0^y dy'/mu(T(y')); we
use it as an independent check on the BVP solver. All physical constants come
from couette_config.py so this never drifts from the MFC case.
"""

import couette_config as cfg
import numpy as np
from scipy.integrate import solve_bvp


def solve_couette(npts=2001):
    """Solve the coupled BVP. Returns (y, u, T, tau, q_flux) on a fine mesh."""
    H, U, T0, T1, k = cfg.H, cfg.U, cfg.T0, cfg.T1, cfg.k_therm

    def rhs(y, s):
        u, tau, T, q = s
        mu = cfg.mu_of_T(T)
        return np.vstack([tau / mu, np.zeros_like(y), q / k, -(tau**2) / mu])

    def bc(sa, sb):
        # sa = state at y=0, sb = state at y=H
        return np.array([sa[0] - 0.0, sb[0] - U, sa[2] - T0, sb[2] - T1])

    y = np.linspace(0.0, H, 51)
    # Initial guess: linear u and T; stress and heat flux from those slopes.
    mu_mid = cfg.mu_of_T(cfg.T_ref)
    s0 = np.vstack(
        [
            U * y / H,
            np.full_like(y, mu_mid * U / H),
            T0 + (T1 - T0) * y / H,
            np.full_like(y, k * (T1 - T0) / H),
        ]
    )
    sol = solve_bvp(rhs, bc, y, s0, tol=1e-10, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(f"solve_bvp failed: {sol.message}")

    yf = np.linspace(0.0, H, npts)
    u, tau, T, q = sol.sol(yf)
    return yf, u, T, float(tau[0]), q


def u_decoupled(y):
    """Brinkman -> 0 closed form: linear T, u = U * I(y)/I(H)."""
    yg = np.linspace(0.0, cfg.H, 20001)
    integrand = 1.0 / cfg.mu_of_T(cfg.T_linear(yg))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(yg))])
    return cfg.U * np.interp(y, yg, cum) / cum[-1]


def u_at(y):
    """Exact (coupled-BVP) velocity sampled at arbitrary y."""
    yf, u, _T, _tau, _q = solve_couette()
    return np.interp(y, yf, u)


def T_at(y):
    """Exact (coupled-BVP) temperature sampled at arbitrary y."""
    yf, _u, T, _tau, _q = solve_couette()
    return np.interp(y, yf, T)


if __name__ == "__main__":
    yf, u, T, tau, q = solve_couette()
    # Cross-check the coupled solver against the Br->0 closed form. With the tiny
    # Brinkman number of this case the two velocity profiles agree closely; the
    # gap is exactly the viscous-heating contribution.
    ud = u_decoupled(yf)
    rel_u_gap = np.max(np.abs(u - ud)) / cfg.U
    Tlin = cfg.T_linear(yf)
    print(f"shear stress tau          = {tau:.6g}")
    print(f"max |u_coupled - u_decoup| = {rel_u_gap * cfg.U:.3e}  ({100 * rel_u_gap:.3f}% of U; ~ Brinkman effect, Br={cfg.Br:.3g})")
    print(f"max viscous-heating bump   = {np.max(T - Tlin):.4g} K above the linear profile")
    print(f"u(H/2) curved value        = {np.interp(cfg.H / 2, yf, u):.6g} (straight-line Couette would give {cfg.U / 2:.6g})")
