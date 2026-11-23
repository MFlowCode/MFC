"""
Body-fitted 1D reference for a cylinder at Tw in quiescent gas at uniform pressure (low Mach): conduction only, with the
thermal conductivity lambda(T) and heat capacity c_p(T) of the real gas.

In the mass coordinate psi (d psi = rho r dr, the wall at psi = 0 since no gas crosses it) the expansion flow drops out:

    c_p dT/dt = d/dpsi (rho r^2 lambda dT/dpsi),    r^2 = R^2 + 2 int_0^psi dpsi'/rho,

a nonlinear diffusion solved by finite volumes (cells geometric from the wall) and scipy BDF. The wall heat into the body is
q = lambda dT/dr = rho r lambda dT/dpsi at psi = 0 and Q(t) = int 2 pi R q dt per unit depth is carried as one more unknown.
The initial profile is MFC's thermal layer, T = Tinf + (Tw - Tinf) (R/r)^(1/2) erfc((r - R)/delta0).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.special import erfc

RU = 8.314462618


def solve_props(Tw, Tinf, lam, cp, W, p, t_out, delta0, R, L, n=800, stretch=1.006):
    """lam(T) [W/m/K], cp(T) [J/kg/K], molar mass W [kg/mol], pressure p [Pa]. Returns rows (t, wall heat flux q, Q)."""
    rho = lambda T: p * W / (RU * T)  # noqa: E731
    rr = R + L * np.linspace(0, 1, 400001) ** 2
    T0r = Tinf + (Tw - Tinf) * np.sqrt(R / rr) * erfc((rr - R) / delta0)
    psi_r = np.concatenate([[0.0], np.cumsum(0.5 * (rho(T0r[1:]) * rr[1:] + rho(T0r[:-1]) * rr[:-1]) * np.diff(rr))])
    h = stretch ** np.arange(n)
    h *= psi_r[-1] / h.sum()
    pf = np.concatenate([[0.0], np.cumsum(h)])
    pc = 0.5 * (pf[1:] + pf[:-1])
    T = np.interp(pc, psi_r, T0r)
    dpc = np.diff(np.concatenate([[0.0], pc, [pf[-1]]]))

    def flux(T):
        r2f = R**2 + 2 * np.concatenate([[0.0], np.cumsum(h / rho(T))])
        Tf = np.concatenate([[Tw], 0.5 * (T[1:] + T[:-1]), [Tinf]])
        Tl = np.concatenate([[Tw], T])
        Tr = np.concatenate([T, [Tinf]])
        return rho(Tf) * r2f * lam(Tf) * (Tr - Tl) / dpc

    def rhs(t, y):
        F = flux(y[:-1])
        return np.concatenate([(F[1:] - F[:-1]) / (cp(y[:-1]) * h), [2 * np.pi * F[0]]])  # q = F[0] / R

    sp = diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(n + 1, n + 1)).tolil()
    sp[n, 0] = 1.0
    sp[n, 1] = 1.0
    sol = solve_ivp(rhs, (0, max(t_out)), np.concatenate([T, [0.0]]), method="BDF", t_eval=sorted(t_out), rtol=1e-8, atol=1e-14, jac_sparsity=sp.tocsr())
    return np.array([(t, flux(sol.y[:-1, j])[0] / R, sol.y[-1, j]) for j, t in enumerate(sol.t)])
