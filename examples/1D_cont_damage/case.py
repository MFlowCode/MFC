#!/usr/bin/python
import json

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": 100,
            "n": 0,
            "p": 0,
            "dt": 5e-10,
            "t_step_start": 0,
            "t_step_stop": 50000,
            "t_step_save": 50000,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 1,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Turning on Hypoelasticity
            "hypoelasticity": "T",
            "fd_order": 4,
            "cont_damage": "T",
            "tau_star": 0.0,
            "cont_damage_s": 2.0,
            "alpha_bar": 1e-4,
            # Acoustic Source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 1,
            "acoustic(1)%loc(1)": 0.1,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%npulse": 999,
            "acoustic(1)%dir": 1.0,
            "acoustic(1)%mag": 1000.0,
            "acoustic(1)%frequency": 1e4,
            "acoustic(1)%delay": 0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1 L
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 1e5,
            "patch_icpp(1)%alpha_rho(1)": 1000 * (1.0 - 1e-6),
            "patch_icpp(1)%alpha_rho(2)": 1000 * 1e-6,
            "patch_icpp(1)%alpha(1)": 1.0 - 1e-6,
            "patch_icpp(1)%alpha(2)": 1e-6,
            "patch_icpp(1)%tau_e(1)": 0.0,
            # Patch 2 R
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 1e5,
            "patch_icpp(2)%alpha_rho(1)": 1000 * 1e-6,
            "patch_icpp(2)%alpha_rho(2)": 1000 * (1.0 - 1e-6),
            "patch_icpp(2)%alpha(1)": 1e-6,
            "patch_icpp(2)%alpha(2)": 1.0 - 1e-6,
            "patch_icpp(2)%tau_e(1)": 0.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%G": 0e0,
            "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%G": 1.0e09,
        }
    )
)
