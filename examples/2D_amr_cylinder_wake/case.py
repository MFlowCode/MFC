#!/usr/bin/env python3
# 2D viscous cylinder wake with block-structured AMR (see docs/documentation/amr.md).
#
# A low-Mach stream (U = 0.2, c = 1, Re = 100 on the diameter) sheds a von Karman street
# behind a no-slip immersed cylinder in a 36 D x 32 D domain (3 % blockage). One static 2:1
# block covers the body and the near wake, where the shedding is set; the far wake stays
# coarse. A small transverse-velocity patch behind the body starts the shedding early.
# Drag and lift are in D/ib_forces.dat; check_wake.py prints Cd and St (README.md).
#
# Variants (the "When AMR pays" table): set amr = False for a uniform grid at N cells per
# diameter (N = 32 is the grid the block resolves); subcycle = True for amr_subcycle.
# dt follows the finest cells; every variant ends at the same time t_end.
import json

N = 16  # level-0 cells per diameter
amr = True
subcycle = False

D, U, c, Re, gam = 1.0, 0.2, 1.0, 100.0, 1.4
x0, x1, y0, y1 = -12.0, 24.0, -16.0, 16.0
dx = D / N
t_end = 650.0

# CFL 0.384 on the cells each step advances at dt (the level-0 cells when subcycled)
h = dx / 2 if amr and not subcycle else dx
dt = 0.32 * h
steps = int(round(t_end / dt))
stride = int(round(0.2 / dt))  # force sample every 0.2 time units


def cell(x, lo):
    return int(round((x - lo) / dx))


case = {
    # Logistics
    "run_time_info": "T",
    # Computational domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "m": cell(x1, x0) - 1,
    "n": cell(y1, y0) - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": steps,
    # Simulation algorithm
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "viscous": "T",
    "fd_order": 2,
    # Dirichlet inflow holds the stream; extrapolation elsewhere
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "ib_force_wrt": "T",
    "ib_force_stride": stride,
    # Patch 1: uniform stream in +x
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (x0 + x1),
    "patch_icpp(1)%y_centroid": 0.5 * (y0 + y1),
    "patch_icpp(1)%length_x": x1 - x0,
    "patch_icpp(1)%length_y": y1 - y0,
    "patch_icpp(1)%vel(1)": U,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": c**2 / gam,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Patch 2: transverse kick above the near wake, to break the symmetry
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 1.5,
    "patch_icpp(2)%y_centroid": 0.5,
    "patch_icpp(2)%length_x": 2.0,
    "patch_icpp(2)%length_y": 1.0,
    "patch_icpp(2)%vel(1)": U,
    "patch_icpp(2)%vel(2)": 0.2 * U,
    "patch_icpp(2)%pres": c**2 / gam,
    "patch_icpp(2)%alpha_rho(1)": 1.0,
    "patch_icpp(2)%alpha(1)": 1.0,
    # Immersed cylinder
    "ib": "T",
    "num_ibs": 1,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.0,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": D / 2,
    "patch_ib(1)%slip": "F",
    # Fluid: ideal gas, mu = rho U D / Re
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": Re / (U * D),
}

if amr:
    case.update(
        {
            # One static block over the body and near wake: -1.5 <= x <= 6, |y| <= 2.5
            "amr": "T",
            "amr_block_beg(1)": cell(-1.5, x0),
            "amr_block_end(1)": cell(6.0, x0) - 1,
            "amr_block_beg(2)": cell(-2.5, y0),
            "amr_block_end(2)": cell(2.5, y0) - 1,
            "amr_regrid_int": 0,
            "amr_max_grid_size": 256,
            "amr_subcycle": "T" if subcycle else "F",
        }
    )

print(json.dumps(case))
