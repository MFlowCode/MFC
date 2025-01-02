#!/usr/bin/env python3
#!/usr/bin/env python3
# Dependencies and Logistics
# Command to navigate between directories

import math
import json

myeps = 1.4 / 150.0

p_l, p_g = 1e06, 1e06
rho_l, rho_g = 1000.0, 1.0
v_l, v_g = 500.0, -500.0
c_l, c_g = math.sqrt((p_l + 4.0e8) / rho_l), math.sqrt(1.4 * p_g / rho_g)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.5,
            "x_domain%end": 0.5,
            "y_domain%beg": -0.5,
            "y_domain%end": 0.5,
            "cyl_coord": "F",
            "m": 99,
            "n": 99,
            "p": 0,
            "dt": 5.0e-10,
            "t_step_start": 0,
            "t_step_stop": 2000,
            "t_step_save": 20,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Top fluid, water
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0,
            "patch_icpp(1)%length_x": 1.0e00,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": v_l,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": p_l,
            "patch_icpp(1)%alpha_rho(1)": rho_l,
            "patch_icpp(1)%alpha_rho(2)": rho_l,
            "patch_icpp(1)%alpha(1)": 0.5e00,
            "patch_icpp(1)%alpha(2)": 0.5e00,
            # Patch 2: Main bottom fluid, water
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.25,
            "patch_icpp(2)%length_x": 1.0e00,
            "patch_icpp(2)%length_y": 0.5,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": v_g,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": p_g,
            "patch_icpp(2)%alpha_rho(1)": 0.0,
            "patch_icpp(2)%alpha_rho(2)": rho_g,
            "patch_icpp(2)%alpha(1)": 0.0,
            "patch_icpp(2)%alpha(2)": 1.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%Re(1)": 0.0001,
            "fluid_pp(1)%Re(2)": 0.0001,
            "fluid_pp(2)%Re(1)": 0.0001,
            "fluid_pp(2)%Re(2)": 0.0001,
        }
    )
)
