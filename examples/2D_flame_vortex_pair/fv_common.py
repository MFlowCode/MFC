"""Shared setup of the premixed H2/air flame + vortex cases (2D pair here, 3D ring in 3D_flame_vortex_ring): geometry, the Cantera
flame profile (flame_profile.txt, written by make_profile.py) and the numerics common to both case files."""

import json
import os

import cantera as ct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MECH = "h2o2.yaml"
P0, TU, PHI = 101325.0, 300.0, 1.0
L = 3.2e-3  # domain [0, L]^d, periodic in the transverse directions
XF = 1.8e-3  # initial flame position (T = 1200 K crossing)
DT = 2.632e-9  # fixed time step [s] at the 64^2 grid; the same on every grid so the first-order split-chemistry error is grid independent


def gas_unburnt():
    g = ct.Solution(MECH)
    g.TP = TU, P0
    g.set_equivalence_ratio(PHI, "H2", "O2:1, N2:3.76")
    return g


def profile():
    """Cantera FreeFlame profile (lab frame, reactants at rest) as a dict of arrays: x, T, u, Y (species x points), plus S_L and rho_u."""
    path = os.path.join(HERE, "flame_profile.txt")
    d = np.loadtxt(path, delimiter=",", skiprows=2)
    SL = float(open(path).readline().split("SL=")[1].split()[0])
    return {"x": d[:, 0], "T": d[:, 1], "u": d[:, 2], "Y": d[:, 3:].T, "SL": SL}


def flame_on(xc, f):
    """T, u (lab frame) and Y of the flame (T = 1200 K at XF) at the cell centres xc."""
    x1 = f["x"] - np.interp(1200.0, f["T"], f["x"]) + XF
    return np.interp(xc, x1, f["T"]), np.interp(xc, x1, f["u"]) - f["SL"], np.array([np.interp(xc, x1, y) for y in f["Y"]])


def numerics(g, mfc):
    """Case entries common to the 2D and 3D cases: WENO5-M, HLLC, SSP-RK3, viscous, mixture-averaged chemistry with split alpha-QSS reactions."""
    return {
        "run_time_info": "T",
        "cyl_coord": "F",
        "model_eqns": 2,
        "num_fluids": 1,
        "alt_soundspeed": "F",
        "mpp_lim": "F",
        "mixture_err": "T",
        "time_stepper": 3,
        "weno_order": 5,
        "weno_eps": 1e-16,
        "weno_Re_flux": "T",
        "weno_avg": "T",
        "avg_state": 2,
        "mapped_weno": "T",
        "null_weights": "F",
        "mp_weno": "F",
        "riemann_solver": 2,
        "wave_speeds": 1,
        "viscous": "T",
        "fd_order": 4,
        "chemistry": "T",
        "chem_params%diffusion": "T",
        "chem_params%reactions": "T",
        "chem_params%transport_model": 1,
        "chem_params%reaction_substeps": 2,
        "chem_params%adap_substeps": "T",
        "chem_params%reaction_substeps_max": 16,
        "cantera_file": MECH,
        "chem_wrt_T": "T",
        "fluid_pp(1)%gamma": 1.0 / 0.4,
        "fluid_pp(1)%eos": "stiffened_gas",
        "fluid_pp(1)%pi_inf": 0.0,
        "fluid_pp(1)%Re(1)": 1.0 / g.viscosity,
        "parallel_io": "T" if mfc.get("mpi", True) else "F",
        "format": "silo",
        "precision": "double",
        "prim_vars_wrt": "T",
        "bc_x%beg": -2,
        "bc_x%end": -8,
        "bc_x%grcbc_out": "T",
        "bc_x%pres_out": P0,
    }


def apply_set(case, items):
    """--set KEY=VAL overrides (numbers stay numbers: the toolchain checks types)."""
    for kv in items:
        k, v = kv.split("=", 1)
        try:
            case[k] = json.loads(v)
        except json.JSONDecodeError:
            case[k] = v
