#!/usr/bin/env python3
"""Exact Riemann solver for the single-material JWL shock tube.

The EOS implemented here is byte-for-byte the same closure MFC uses for a pure
JWL material (alpha_j = 1) in src/common/m_variables_conversion.fpp:

    pref(rho) = A(1 - w*rho/(R1*rho0)) exp(-R1*rho0/rho)
              + B(1 - w*rho/(R2*rho0)) exp(-R2*rho0/rho)
    p(rho,e)  = pref(rho) + (w/rho0) * rho * e
    e(rho,p)  = (p - pref(rho)) * rho0 / (w * rho)

and the sound speed is the exact s_jwl_sound_speed_squared form.

Because the analytic solver and the solver share the SAME EOS, any discrepancy
in the comparison is attributable to the numerical scheme (reconstruction,
Riemann flux, time integration) -- this is a genuine verification, not a
physics-model check.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


class JWL:
    def __init__(self, A, B, R1, R2, omega, rho0):
        self.A, self.B, self.R1, self.R2 = A, B, R1, R2
        self.w, self.rho0 = omega, rho0

    def pref(self, rho):
        A, B, R1, R2, w, r0 = self.A, self.B, self.R1, self.R2, self.w, self.rho0
        return (A * (1 - w * rho / (R1 * r0)) * np.exp(-R1 * r0 / rho)
                + B * (1 - w * rho / (R2 * r0)) * np.exp(-R2 * r0 / rho))

    def p_of_rho_e(self, rho, e):
        return self.pref(rho) + (self.w / self.rho0) * rho * e

    def e_of_rho_p(self, rho, p):
        return (p - self.pref(rho)) * self.rho0 / (self.w * rho)

    def c2(self, rho, p):
        A, B, R1, R2, w, r0 = self.A, self.B, self.R1, self.R2, self.w, self.rho0
        v = r0 / rho
        e1, e2 = np.exp(-R1 * v), np.exp(-R2 * v)
        e = (v / w) * (p - A * (1 - w / (R1 * v)) * e1 - B * (1 - w / (R2 * v)) * e2)
        c2 = (A * e1 * (w / (R1 * v**2) - R1 * (1 - w / (R1 * v)))
              + B * e2 * (w / (R2 * v**2) - R2 * (1 - w / (R2 * v))))
        c2 = c2 * (-r0 / rho**2) + w * e / r0 + w * p / (rho * r0)
        return max(c2, 1e-12)

    def c(self, rho, p):
        return np.sqrt(self.c2(rho, p))

    # density behind a shock with post-shock pressure ps, pre-shock (rho_k,p_k)
    def rho_hugoniot(self, ps, rho_k, p_k):
        e_k = self.e_of_rho_p(rho_k, p_k)

        def hug(rho_s):
            e_s = self.e_of_rho_p(rho_s, ps)
            return (e_s - e_k) - 0.5 * (ps + p_k) * (1.0 / rho_k - 1.0 / rho_s)

        # compressed branch rho_s in (rho_k, ~8*rho0); scan for a sign change.
        grid = np.linspace(rho_k * (1 + 1e-9), 8.0 * self.rho0, 4000)
        vals = np.array([hug(r) for r in grid])
        sign = np.where(np.diff(np.sign(vals)) != 0)[0]
        if len(sign) == 0:
            raise RuntimeError("no Hugoniot root found")
        i = sign[0]
        return brentq(hug, grid[i], grid[i + 1], xtol=1e-10, rtol=1e-12)

    # isentrope ODE: integrate rho(p) and G(p)=int dp/(rho c) from p_k to p
    def _isentrope(self, p_k, rho_k, p_target):
        def rhs(p, y):
            rho = y[0]
            c = self.c(rho, p)
            return [1.0 / self.c2(rho, p), 1.0 / (rho * c)]
        sol = solve_ivp(rhs, [p_k, p_target], [rho_k, 0.0],
                        rtol=1e-9, atol=1e-12, dense_output=True, max_step=abs(p_k - p_target) / 50 + 1)
        return sol


def solve(jwlL, jwlR, WL, WR, x0, t, x):
    """WL=(rhoL,uL,pL), WR=(rhoR,uR,pR). Returns dict of exact profiles on x."""
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR
    J = jwlL  # same material both sides

    def F_and_state(ps, rho_k, p_k):
        """Return F_K(ps) where u* = u_k -/+ F_K, and rho behind wave."""
        if ps > p_k:  # shock
            rho_s = J.rho_hugoniot(ps, rho_k, p_k)
            F = np.sqrt((ps - p_k) * (1.0 / rho_k - 1.0 / rho_s))
            return F, rho_s
        else:  # rarefaction
            sol = J._isentrope(p_k, rho_k, ps)
            rho_s = sol.y[0, -1]
            G = sol.y[1, -1]          # int_{p_k}^{ps} dp/(rho c)  (negative since ps<p_k)
            return G, rho_s

    def velocity_mismatch(ps):
        FL, _ = F_and_state(ps, rhoL, pL)
        FR, _ = F_and_state(ps, rhoR, pR)
        # u*_L = uL - FL ; u*_R = uR + FR  ; match:
        return (uL - FL) - (uR + FR)

    # bracket p* by scanning for a sign change over a sensible pressure window.
    grid = np.linspace(0.2 * min(pL, pR), 3.0 * max(pL, pR), 2000)
    vals = np.array([velocity_mismatch(p) for p in grid])
    idx = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(idx) == 0:
        raise RuntimeError("no p* root found in scan window")
    k = idx[0]
    pstar = brentq(velocity_mismatch, grid[k], grid[k + 1], xtol=1.0, rtol=1e-12)

    FL, rhoLs = F_and_state(pstar, rhoL, pL)
    FR, rhoRs = F_and_state(pstar, rhoR, pR)
    ustar = 0.5 * ((uL - FL) + (uR + FR))

    cL, cR = J.c(rhoL, pL), J.c(rhoR, pR)
    cLs, cRs = J.c(rhoLs, pstar), J.c(rhoRs, pstar)

    # Build sampling helper for each fan via pressure-parametrised tables.
    def left_fan_tables():
        sol = J._isentrope(pL, rhoL, pstar)
        ps = np.linspace(pL, pstar, 400)
        rho = sol.sol(ps)[0]
        G = sol.sol(ps)[1]
        u = uL + (-(G))               # u(p)=uL - int_{p}^{pL} = uL + G(p) ; G negative -> wait
        # u* = uL - FL, FL=G(pstar)<0 => u*=uL-G(pstar). For p between: u=uL - G(p)
        u = uL - G
        c = np.array([J.c(r, p) for r, p in zip(rho, ps)])
        speed = u - c
        return ps, rho, u, c, speed

    def right_fan_tables():
        sol = J._isentrope(pR, rhoR, pstar)
        ps = np.linspace(pR, pstar, 400)
        rho = sol.sol(ps)[0]
        G = sol.sol(ps)[1]
        u = uR + G                    # u* = uR + FR, FR=G(pstar)<0
        c = np.array([J.c(r, p) for r, p in zip(rho, ps)])
        speed = u + c
        return ps, rho, u, c, speed

    rho_out = np.empty_like(x)
    u_out = np.empty_like(x)
    p_out = np.empty_like(x)

    left_is_shock = pstar > pL
    right_is_shock = pstar > pR
    if not left_is_shock:
        lp, lrho, lu, lc, lspeed = left_fan_tables()
    if not right_is_shock:
        rp, rrho, ru, rc, rspeed = right_fan_tables()

    SL_head = uL - cL
    SL_tail = ustar - cLs
    SR_tail = ustar + cRs
    SR_head = uR + cR
    if left_is_shock:
        jL = np.sqrt((pstar - pL) / (1.0 / rhoL - 1.0 / rhoLs))
        SL = uL - jL / rhoL
    if right_is_shock:
        jR = np.sqrt((pstar - pR) / (1.0 / rhoR - 1.0 / rhoRs))
        SR = uR + jR / rhoR

    for i, xi_x in enumerate(x):
        xi = (xi_x - x0) / t
        if xi < ustar:  # left of contact
            if left_is_shock:
                if xi < SL:
                    rho_out[i], u_out[i], p_out[i] = rhoL, uL, pL
                else:
                    rho_out[i], u_out[i], p_out[i] = rhoLs, ustar, pstar
            else:
                if xi <= SL_head:
                    rho_out[i], u_out[i], p_out[i] = rhoL, uL, pL
                elif xi >= SL_tail:
                    rho_out[i], u_out[i], p_out[i] = rhoLs, ustar, pstar
                else:
                    rho_out[i] = np.interp(xi, lspeed, lrho)
                    u_out[i] = np.interp(xi, lspeed, lu)
                    p_out[i] = np.interp(xi, lspeed, lp)
        else:  # right of contact
            if right_is_shock:
                if xi > SR:
                    rho_out[i], u_out[i], p_out[i] = rhoR, uR, pR
                else:
                    rho_out[i], u_out[i], p_out[i] = rhoRs, ustar, pstar
            else:
                if xi >= SR_head:
                    rho_out[i], u_out[i], p_out[i] = rhoR, uR, pR
                elif xi <= SR_tail:
                    rho_out[i], u_out[i], p_out[i] = rhoRs, ustar, pstar
                else:
                    # right fan tables go pR->pstar; speed monotonic
                    order = np.argsort(rspeed)
                    rho_out[i] = np.interp(xi, rspeed[order], rrho[order])
                    u_out[i] = np.interp(xi, rspeed[order], ru[order])
                    p_out[i] = np.interp(xi, rspeed[order], rp[order])

    return dict(rho=rho_out, u=u_out, p=p_out, pstar=pstar, ustar=ustar,
                rhoLs=rhoLs, rhoRs=rhoRs,
                left_shock=left_is_shock, right_shock=right_is_shock)
