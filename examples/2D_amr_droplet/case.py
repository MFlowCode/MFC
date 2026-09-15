#!/usr/bin/env python3
# 2D advected density droplet with block-structured AMR (see docs/documentation/amr.md).
#
# A dense circular droplet in pressure equilibrium is carried across the domain by a
# uniform flow. Its sharp interface is the only feature worth resolving, so a 2:1 refined
# block tracks the moving droplet by dynamic regridding (amr_regrid_int > 0), subcycled at
# dt/2 (amr_subcycle), while the smooth surroundings stay coarse. This is the intended
# use of AMR: a compact, interior feature refined locally at bounded cost. The four
# required AMR settings are grouped at the bottom; delete that block for a uniform run.
import json

N = 128
dx = 1.0 / N
u = 0.8  # uniform advection velocity (+x)

# Coarse-grid CFL step: max wave speed ~ u + sqrt(1.4 * 1 / 0.125) ~ 0.8 + 3.7. The fine
# block subcycles at dt/2, so this coarse dt already satisfies the finest-cell CFL.
dt = 0.1 * dx / 4.5

# Initial refined block over the droplet, in level-0 cell indices (0-based, inclusive).
# It sits well inside the domain (a block must stay >= buff_size cells from every
# boundary) and spans at most half the grid per dimension.
bx0, bx1 = int(0.20 * N), int(0.44 * N)
by0, by1 = int(0.38 * N), int(0.62 * N)

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational domain
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": N - 1,
            "n": N - 1,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 600,
            "t_step_save": 30,
            # Simulation algorithm
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: ambient gas in uniform +x flow
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": u,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 0.125,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: dense droplet (same pressure and velocity - a passive contact)
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%x_centroid": 0.3,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%radius": 0.12,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": u,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%alpha_rho(1)": 1.0,
            "patch_icpp(2)%pres": 1.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluid (ideal gas, gamma = 1.4)
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            # --- AMR: the four required settings ---
            # amr + an initial block are mandatory; amr_regrid_int > 0 turns on dynamic
            # regridding (density-gradient tagging), amr_subcycle advances the fine block at
            # dt/2. Delete this block for a uniform-grid run.
            "amr": "T",
            "amr_block_beg(1)": bx0,
            "amr_block_end(1)": bx1,
            "amr_block_beg(2)": by0,
            "amr_block_end(2)": by1,
            "amr_regrid_int": 10,
            "amr_tag_eps": 0.1,
            "amr_buf": 4,
            "amr_max_blocks": 16,
            # A closed interface tags a thin ring; a looser clustering efficiency lets the
            # Berger-Rigoutsos clusterer cover it with a few boxes instead of many small tiles.
            "amr_cluster_eff": 0.35,
            "amr_subcycle": "T",
        }
    )
)
