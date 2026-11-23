#!/usr/bin/env python3
"""Writes flame_profile.txt: Cantera FreeFlame (h2o2.yaml, stoichiometric H2/air, 300 K, 1 atm, mixture-averaged transport, 2 cm wide,
refine ratio 2, slope 0.02, curve 0.04): x [m], T [K], u [m/s], Y_k. The header carries the flame speed S_L [m/s].

    python3 examples/2D_flame_vortex_pair/make_profile.py
"""

import os

import cantera as ct
import numpy as np
from fv_common import HERE, MECH, gas_unburnt

g = gas_unburnt()
f = ct.FreeFlame(g, width=0.02)
f.transport_model = "mixture-averaged"
f.set_refine_criteria(ratio=2, slope=0.02, curve=0.04)
f.solve(loglevel=0, auto=True)
# a 4 mm window around the flame carries the whole profile; the rest is uniform
T = f.T
x0 = np.interp(1200.0, T, f.grid)
keep = np.abs(f.grid - x0) < 2.0e-3
cols = np.column_stack([f.grid[keep], T[keep], f.velocity[keep], f.Y[:, keep].T])
head = f"SL={f.velocity[0]:.8f} thermal_thickness={(T[-1] - T[0]) / np.gradient(T, f.grid).max():.6e}\nx,T,u," + ",".join(g.species_names)
np.savetxt(os.path.join(HERE, "flame_profile.txt"), cols, delimiter=",", header=head, comments="", fmt="%.7e")
print(f"S_L = {f.velocity[0]:.4f} m/s, {keep.sum()} points")
