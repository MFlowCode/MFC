#!/usr/bin/env python3
"""
1D spherically (--geom 2) or cylindrically (--geom 1) symmetric low-Mach reference for nD_ibm_carbon_burn: a carbon body at fixed T_w
gasifying in quiescent hot moist air, with the gas-phase CO flame, from the same mechanism files as the MFC case.

    python3 examples/nD_ibm_carbon_burn/ref1d.py --geom 1 --tend 6e-5 --N 400 --out ref_cylinder.txt

Finite volumes on r in [R, R_out] (geometric stretching, first cell 0.002 R), cell unknowns T and Y_k; Cantera thermodynamics,
mixture-averaged transport (correction-velocity form), gas kinetics and the surface mechanism.
  rho DY_k/Dt = -div j_k + W_k w_k,   rho cp DT/Dt = -div q - sum_k j_k cp_k dT/dr - sum_k h_k W_k w_k
  div u = (1/T) DT/Dt + Wbar sum_k (1/W_k) DY_k/Dt      (ideal gas at constant p)  ->  r^g u = R^g u_w + int S r^g dr
Wall r = R: fixed T_w; the gas receives the surface-solve species fluxes W_k s_k (as ib_flux_wall); the wall gas state is the first cell;
Stefan velocity u_w = sum W_k s_k / rho_w. Outer face: Dirichlet T_inf, Y_inf. Time integration: scipy BDF with banded sparsity.
Output: t [us], gasification rate [kg/s, kg/m/s for the cylinder], T_max [K] and its radius over R.
"""

import argparse
import os

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

MECH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2D_ibm_reacting_surface")
GAS = os.path.join(MECH, "carbon_gasphase_reduced_gri11.yaml")
SURF = os.path.join(MECH, "carbon_surface_bradley_11species.yaml")


class Radial:
    def __init__(self, R=50e-6, Tw=2000.0, Tinf=2000.0, p=ct.one_atm, Yinf=None, Rout_over_R=60.0, N=160, dr0_over_R=2e-3, frozen=False, surface=True, gas_file=GAS, surf_file=SURF, geom=2):
        """geom 2: sphere (areas r^2 per steradian); geom 1: cylinder (areas r per radian, unit length)."""
        self.gas = ct.Solution(gas_file)
        if surface:
            graphite = ct.Solution(surf_file, "graphite")
            self.surf = ct.Interface(surf_file, "carbon_surface", [self.gas, graphite])
        self.surface, self.frozen = surface, frozen
        self.R, self.Tw, self.Tinf, self.p = R, Tw, Tinf, p
        self.K = self.gas.n_species
        self.W = self.gas.molecular_weights
        self.gas.TPY = Tinf, p, Yinf
        self.Yinf = self.gas.Y.copy()
        self.N = N
        # geometric faces: first cell dr0, last face at Rout
        Ro = Rout_over_R * R
        q = self._ratio(dr0_over_R * R, Ro - R, N)
        dr = dr0_over_R * R * q ** np.arange(N)
        self.rf = R + np.concatenate([[0.0], np.cumsum(dr)])
        self.rf[-1] = Ro
        self.rc = 0.5 * (self.rf[1:] + self.rf[:-1])
        self.geom = geom
        self.Af = self.rf**geom
        self.V = (self.rf[1:] ** (geom + 1) - self.rf[:-1] ** (geom + 1)) / (geom + 1)
        self.arr = ct.SolutionArray(self.gas, N)
        self.nrhs = 0

    @staticmethod
    def _ratio(dr0, L, N):
        lo, hi = 1.0 + 1e-12, 2.0
        for _ in range(200):
            m = 0.5 * (lo + hi)
            s = dr0 * (m**N - 1) / (m - 1)
            lo, hi = (m, hi) if s < L else (lo, m)
        return 0.5 * (lo + hi)

    def unpack(self, y):
        z = y.reshape(self.N, 1 + self.K)
        return z[:, 0], np.clip(z[:, 1:], 0.0, None)

    def wall_fluxes(self, T0, Y0):
        """Species mass fluxes into the gas [kg/m^2/s] from the surface solve at (T_w, wall gas = first cell)."""
        if not self.surface:
            return np.zeros(self.K)
        self.gas.TPY = self.Tw, self.p, Y0
        self.surf.TP = self.Tw, self.p
        return self.surf.get_net_production_rates(self.gas) * self.W

    def rhs(self, t, y):
        self.nrhs += 1
        N, K, W = self.N, self.K, self.W
        T, Y = self.unpack(y)
        Y = Y / Y.sum(axis=1, keepdims=True)
        a = self.arr
        a.TPY = T, self.p, Y
        rho, cp, lam = a.density, a.cp_mass, a.thermal_conductivity
        X, Wb, D = a.X, a.mean_molecular_weight, a.mix_diff_coeffs
        hk = a.partial_molar_enthalpies / W  # J/kg
        cpk = a.partial_molar_cp / W
        wdot = np.zeros((N, K)) if self.frozen else a.net_production_rates  # kmol/m^3/s

        # ghost at the outer boundary (Dirichlet, at the face): use the far state at r = Rout
        self.gas.TPY = self.Tinf, self.p, self.Yinf
        Xo, To = self.gas.X, self.Tinf
        rf, rc = self.rf, self.rc
        # interior faces 1..N-1 and outer face N
        dX = np.vstack([X[1:] - X[:-1], Xo - X[-1]])
        dT = np.concatenate([T[1:] - T[:-1], [To - T[-1]]])
        dn = np.concatenate([rc[1:] - rc[:-1], [rf[-1] - rc[-1]]])
        rhof = np.concatenate([0.5 * (rho[1:] + rho[:-1]), [rho[-1]]])
        Df = np.vstack([0.5 * (D[1:] + D[:-1]), D[-1:]])
        Yf = np.vstack([0.5 * (Y[1:] + Y[:-1]), Y[-1:]])
        Xf = np.vstack([0.5 * (X[1:] + X[:-1]), X[-1:]])
        Wbf = np.concatenate([0.5 * (Wb[1:] + Wb[:-1]), [Wb[-1]]])
        lamf = np.concatenate([0.5 * (lam[1:] + lam[:-1]), [lam[-1]]])
        j = -rhof[:, None] * Df * (W[None, :] / Wbf[:, None]) * dX / dn[:, None]
        j -= Yf * j.sum(axis=1, keepdims=True)  # correction velocity
        qf = -lamf * dT / dn
        hkf = np.vstack([0.5 * (hk[1:] + hk[:-1]), hk[-1:]])

        # wall face
        mw = self.wall_fluxes(T[0], Y[0])
        mdot = mw.sum()
        self.gas.TPY = self.Tw, self.p, Y[0]
        rhow = self.gas.density
        jw = mw - mdot * Y[0]  # diffusive part (wall gas = first cell)
        qw = -lam[0] * (T[0] - self.Tw) / (rc[0] - rf[0])

        J = np.vstack([jw[None, :], j])  # faces 0..N
        Q = np.concatenate([[qw], qf])
        A = self.Af
        divj = (A[1:, None] * J[1:] - A[:-1, None] * J[:-1]) / self.V[:, None]
        divq = (A[1:] * Q[1:] - A[:-1] * Q[:-1]) / self.V
        # sum_k j_k cp_k dT/dr at cell centres (face average of j times central gradient)
        Jc = 0.5 * (J[1:] + J[:-1])
        Tg = np.concatenate([[self.Tw], T, [To]])
        rg = np.concatenate([[rf[0]], rc, [rf[-1]]])
        dTdr = (Tg[2:] - Tg[:-2]) / (rg[2:] - rg[:-2])
        aY = (-divj + W[None, :] * wdot) / rho[:, None]
        aT = (-divq - (Jc * cpk).sum(axis=1) * dTdr - (hk * W[None, :] * wdot).sum(axis=1)) / (rho * cp)
        # wall enthalpy: products enter at T_w (energy flux sum m_k h_k(T_w)), the cell's enthalpy form needs the
        # difference to the cell's own h_k: rho cp DT/Dt gets -(1/V) A_0 sum m_k (h_k(T0) - h_k(Tw)) ... added below
        self.gas.TPY = self.Tw, self.p, Y[0]
        hkw = self.gas.partial_molar_enthalpies / W
        aT[0] += A[0] * (mw * (hkw - hk[0])).sum() / self.V[0] / (rho[0] * cp[0])

        S = aT / T + Wb * (aY / W[None, :]).sum(axis=1)
        uw = mdot / rhow
        r2u = A[0] * uw + np.concatenate([[0.0], np.cumsum(S * self.V)])  # r^2 u at faces 0..N
        uc = 0.5 * (r2u[1:] / A[1:] + r2u[:-1] / A[:-1])
        # upwind advection (u >= 0 expected); upstream neighbour: wall value for cell 0 is the cell itself
        Tup = np.concatenate([[T[0]], T[:-1]])
        Yup = np.vstack([Y[:1], Y[:-1]])
        dup = np.concatenate([[1.0], rc[1:] - rc[:-1]])
        Tdn = np.concatenate([T[1:], [To]])
        Ydn = np.vstack([Y[1:], self.Yinf[None, :]])
        ddn = dn
        pos = uc >= 0
        dTadv = np.where(pos, (T - Tup) / dup, (Tdn - T) / ddn)
        dYadv = np.where(pos[:, None], (Y - Yup) / dup[:, None], (Ydn - Y) / ddn[:, None])
        dT_dt = aT - uc * dTadv
        dY_dt = aY - uc[:, None] * dYadv
        self.last = dict(mdot=mdot, mw=mw, uw=uw, T=T, Y=Y, rho=rho)
        return np.column_stack([dT_dt, dY_dt]).ravel()

    def sparsity(self):
        n, B = self.N, 1 + self.K
        S = lil_matrix((n * B, n * B), dtype=np.int8)
        for i in range(n):
            for k in (i - 1, i, i + 1):
                if 0 <= k < n:
                    S[i * B : (i + 1) * B, k * B : (k + 1) * B] = 1
        return S.tocsr()

    def initial(self, T0=None, Y0=None):
        T = np.full(self.N, self.Tinf if T0 is None else T0)
        Y = np.tile(self.Yinf if Y0 is None else Y0, (self.N, 1))
        return np.column_stack([T, Y]).ravel()

    def run(self, t_out, y0=None, rtol=1e-6, atol=1e-10, max_step=np.inf, first_step=1e-12):
        y0 = self.initial() if y0 is None else y0
        out = solve_ivp(self.rhs, (0.0, t_out[-1]), y0, method="BDF", t_eval=t_out, rtol=rtol, atol=atol, jac_sparsity=self.sparsity(), max_step=max_step, first_step=first_step)
        return out

    def mdot_total(self, y):
        """Gasification rate [kg/s; kg/m/s for the cylinder] (4 pi R^2 or 2 pi R times sum W_k s_k) for state y."""
        T, Y = self.unpack(y)
        return (4 * np.pi * self.R**2 if self.geom == 2 else 2 * np.pi * self.R) * self.wall_fluxes(T[0], Y[0] / Y[0].sum()).sum()


AIR = {"N2": 0.760557, "O2": 0.233, "H2O": 0.0058, "CO2": 0.000643}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geom", type=int, default=1, choices=[1, 2], help="1 cylinder, 2 sphere")
    ap.add_argument("--R", type=float, default=50e-6)
    ap.add_argument("--p", type=float, default=10.0, help="atm")
    ap.add_argument("--Tw", type=float, default=2000.0)
    ap.add_argument("--Tinf", type=float, default=2000.0)
    ap.add_argument("--N", type=int, default=400)
    ap.add_argument("--Rout", type=float, default=200.0, help="outer radius over R")
    ap.add_argument("--tend", type=float, default=6e-5)
    ap.add_argument("--nout", type=int, default=120)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    s = Radial(R=a.R, Tw=a.Tw, Tinf=a.Tinf, p=a.p * ct.one_atm, Yinf=AIR, N=a.N, Rout_over_R=a.Rout, geom=a.geom)
    t = a.tend * np.linspace(0, 1, a.nout + 1)[1:]
    out = s.run(t, rtol=1e-6, atol=1e-12)
    rows = []
    for k in range(out.y.shape[1]):
        T, _ = s.unpack(out.y[:, k])
        i = int(T.argmax())
        rows.append([out.t[k] * 1e6, abs(s.mdot_total(out.y[:, k])), T[i], s.rc[i] / a.R])
    unit = "kg/m/s" if a.geom == 1 else "kg/s"
    np.savetxt(
        a.out,
        np.array(rows),
        fmt="%.6e",
        header=f"1D transient reference, geom {a.geom}, R = {a.R * 1e6:g} um, {a.p:g} atm, Tw = {a.Tw:g} K, Tinf = {a.Tinf:g} K, N = {a.N}, Rout = {a.Rout:g} R\nt_us mdot[{unit}] T_max[K] r_Tmax/R",
    )
    print(f"mdot(t_end) = {rows[-1][1]:.6e} {unit}")


if __name__ == "__main__":
    main()
