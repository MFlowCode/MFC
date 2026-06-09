#!/usr/bin/env python3
"""Lid-driven cavity with Herschel-Bulkley non-Newtonian viscosity.

A shear-thinning fluid (power-law index nn=0.5, consistency K=1e-2) fills a
unit square cavity with a moving lid at y_end (bc_y%ve1 = 0.5). The
Papanastasiou-regularized Herschel-Bulkley model (tau0=0, so pure power-law)
is used with viscosity clamped to [mu_min, mu_max].

This case is a qualitative demonstration of HB non-Newtonian viscosity in MFC
and is based on the shear-thinning lid-driven cavity of:

  Li, Z.-Y., Fang, M., Zhang, X., & Wu, Y.-L. (2015). Coupling lattice
  Boltzmann model for simulation of thermal flows on standard lattices.
  Phys. Rev. E, 85, 016710.

The grid (m=n=99) is intentionally coarse for a quick smoke-run.  For
production validation use m=n=499 or finer with a longer t_step_stop.

dt is set small (5e-6) to satisfy the viscous CFL limit imposed by mu_max=1.0
on the coarse 1/99 mesh.  Production runs at finer resolution should reduce dt
proportionally or use a larger pressure to increase the acoustic time scale.

Re_eff = 1/K = 100 (effective Reynolds number at unit shear rate).
"""

import json

K = 1e-2  # consistency index; Re_eff = 1/K = 100

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": 99,
            "n": 99,
            "p": 0,
            "dt": 5e-6,
            "t_step_start": 0,
            "t_step_stop": 4000,
            "t_step_save": 1000,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%ve1": 0.5,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(3)": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 5,
            "patch_icpp(1)%alpha_rho(1)": 0.5,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": 0.5,
            "patch_icpp(1)%alpha(2)": 0.5,
            # Fluids Physical Parameters
            # Re(1) = 1/K registers fluid as viscous; the HB model overrides
            # mu_eff at each cell via the Papanastasiou-regularized formula.
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / K,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%tau0": 0.0,
            "fluid_pp(1)%K": K,
            "fluid_pp(1)%nn": 0.5,
            "fluid_pp(1)%mu_max": 1.0,
            "fluid_pp(1)%mu_min": 1e-6,
            "fluid_pp(1)%hb_m": 1000.0,
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": 1.0 / K,
            "fluid_pp(2)%non_newtonian": "T",
            "fluid_pp(2)%tau0": 0.0,
            "fluid_pp(2)%K": K,
            "fluid_pp(2)%nn": 0.5,
            "fluid_pp(2)%mu_max": 1.0,
            "fluid_pp(2)%mu_min": 1e-6,
            "fluid_pp(2)%hb_m": 1000.0,
        }
    )
)
