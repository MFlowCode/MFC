#!/usr/bin/env python3
"""
Independent reference for the mixing layer of case.py --mode mixing: 1D isobaric multicomponent diffusion without reactions, in
Lagrangian mass coordinates (a wall at x = 0).

    ref_mixing_layer.py M OUT.npz       # M uniform initial cells (a multiple of 3), output at 10 and 20 us

With psi = int rho dx and p fixed,

    dY_k/dt = -d j_k/d psi,    c_p dT/dt = -d q/d psi + sum_k h_k d j_k/d psi,
    j_k = -rho D_k (W_k/W) dX_k/dx + Y_k sum_j rho D_j (W_j/W) dX_j/dx,    q = -lambda dT/dx + sum_k h_k j_k,

with the mixture-averaged D_k, lambda, h_k and c_p of Cantera (Soret off). Finite volumes of fixed mass, two-point face gradients over
the current cell widths, arithmetic face averages (central face Y in the correction flux), zero flux at both ends; SSP-RK3 with
dt = 0.25 min dx^2 / max(D_k, lambda / (rho c_v)). The initial condition is the step of case.py in cell averages.
"""

import sys

import cantera as ct
import numpy as np

P0, LX, X0 = 101325.0, 3.0e-3, 1.0e-3
LEFT, RIGHT = (300.0, "H2:0.3,N2:0.7"), (1000.0, "O2:0.2,H2O:0.2,N2:0.6")
T_OUT = (1.0e-5, 2.0e-5)


class Lagrangian:
    def __init__(self, M, mech="h2o2.yaml"):
        self.g = ct.Solution(mech)
        self.g.transport_model = "mixture-averaged"
        self.K, self.M, self.W = self.g.n_species, M, self.g.molecular_weights
        xf = np.linspace(0.0, LX, M + 1)
        left = 0.5 * (xf[1:] + xf[:-1]) < X0
        self.a = ct.SolutionArray(self.g, M)
        T, Y = np.empty(M), np.empty((M, self.K))
        for sel, (Tk, Xk) in ((left, LEFT), (~left, RIGHT)):
            self.g.TPX = Tk, P0, Xk
            T[sel], Y[sel] = Tk, self.g.Y
        self.a.TPY = T, P0, Y
        self.m = self.a.density * np.diff(xf)
        self.y = np.concatenate([Y.T.ravel(), T / 1e3])
        self.dtmax = 0.0

    def rhs(self, y):
        K, M, a = self.K, self.M, self.a
        Y = y[: K * M].reshape(K, M).T
        a.TPY = y[K * M :] * 1e3, P0, Y
        rho, T, X, Wm = a.density, a.T, a.X, a.mean_molecular_weight
        D, lam = a.mix_diff_coeffs, a.thermal_conductivity
        hk = a.partial_molar_enthalpies / self.W[None, :]
        dx = self.m / rho
        self.dtmax = (dx**2 / np.maximum(D.max(1), lam / (rho * a.cv_mass))).min()
        dl = 0.5 * (dx[1:] + dx[:-1])
        av = lambda f: 0.5 * (f[1:] + f[:-1])  # noqa: E731
        jp = -(av(rho)[:, None] * av(D) * self.W[None, :] / av(Wm)[:, None]) * (np.diff(X, axis=0) / dl[:, None])
        j = jp - av(Y) * jp.sum(1)[:, None]
        q = -av(lam) * np.diff(T) / dl + (av(hk) * j).sum(1)
        Jf, Qf = np.zeros((M + 1, K)), np.zeros(M + 1)
        Jf[1:-1], Qf[1:-1] = j, q
        dY = -(Jf[1:] - Jf[:-1]) / self.m[:, None]
        dT = (-(Qf[1:] - Qf[:-1]) / self.m - (hk * dY).sum(1)) / a.cp_mass
        return np.concatenate([dY.T.ravel(), dT / 1e3])

    def solve(self, t_out=T_OUT, fourier=0.25):
        y, t, out = self.y.copy(), 0.0, []
        for to in t_out:
            while t < to * (1 - 1e-12):
                k1 = self.rhs(y)
                dt = min(fourier * self.dtmax, to - t)
                y1 = y + dt * k1
                y2 = 0.75 * y + 0.25 * (y1 + dt * self.rhs(y1))
                y = y / 3 + 2 / 3 * (y2 + dt * self.rhs(y2))
                t += dt
            Y = y[: self.K * self.M].reshape(self.K, self.M).T
            self.a.TPY = y[self.K * self.M :] * 1e3, P0, Y
            out.append((Y.copy(), self.a.T.copy()))
        psi = np.cumsum(self.m) - 0.5 * self.m
        return psi, out


if __name__ == "__main__":
    M = int(sys.argv[1])
    psi, out = Lagrangian(M).solve()
    np.savez(sys.argv[2], psi=psi, Y=[o[0] for o in out], T=[o[1] for o in out], t=T_OUT, species=Lagrangian(3).g.species_names)
