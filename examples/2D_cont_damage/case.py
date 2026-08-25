#!/usr/bin/env python3
# 2D continuum-damage demonstration: solid disk in water hit by a planar pulse.
#
# A planar high-pressure strip in the water launches a pulse that diffracts around
# and transmits through a damageable solid disk. Tension concentrations (poles of
# the disk during passage, interior wave focusing afterwards) exceed tau_star and
# accumulate damage; the surrounding water stays undamaged (damage is carried by
# the solid partial mass, U_D = m_s D).
#
# Parameters are deliberately sub-critical: with the tangent model, stress is
# retained as D grows, and driving D all the way to 1 under sustained tension
# leaves a zero-shear-stiffness cell with residual stress (a documented model
# limitation). tau_star and alpha_bar here keep max D around 0.2.
import json

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.01,
            "y_domain%beg": 0.0,
            "y_domain%end": 0.01,
            "m": 199,
            "n": 199,
            "p": 0,
            "dt": 5.0e-9,
            "t_step_start": 0,
            "t_step_stop": 1200,
            "t_step_save": 100,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hll",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            # Hypoelasticity + continuum damage
            "hypoelasticity": "T",
            "fd_order": 4,
            "cont_damage": "T",
            "tau_star": 2.0e7,
            "cont_damage_s": 2.0,
            "alpha_bar": 2.0e-5,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Background water (ambient)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.005,
            "patch_icpp(1)%y_centroid": 0.005,
            "patch_icpp(1)%length_x": 0.01,
            "patch_icpp(1)%length_y": 0.01,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0e5,
            "patch_icpp(1)%alpha_rho(1)": 1000.0 * (1.0 - 1e-6),
            "patch_icpp(1)%alpha(1)": 1.0 - 1e-6,
            "patch_icpp(1)%alpha_rho(2)": 1000.0 * 1e-6,
            "patch_icpp(1)%alpha(2)": 1e-6,
            "patch_icpp(1)%tau_e(1)": 0.0,
            "patch_icpp(1)%tau_e(2)": 0.0,
            "patch_icpp(1)%tau_e(3)": 0.0,
            # Driver strip (interior so the open boundaries never act as a reservoir)
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.00175,
            "patch_icpp(2)%y_centroid": 0.005,
            "patch_icpp(2)%length_x": 0.0005,
            "patch_icpp(2)%length_y": 0.01,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 4.0e8,
            "patch_icpp(2)%alpha_rho(1)": 1000.0 * (1.0 - 1e-6),
            "patch_icpp(2)%alpha(1)": 1.0 - 1e-6,
            "patch_icpp(2)%alpha_rho(2)": 1000.0 * 1e-6,
            "patch_icpp(2)%alpha(2)": 1e-6,
            "patch_icpp(2)%tau_e(1)": 0.0,
            "patch_icpp(2)%tau_e(2)": 0.0,
            "patch_icpp(2)%tau_e(3)": 0.0,
            # Damageable solid disk
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": 0.005,
            "patch_icpp(3)%y_centroid": 0.005,
            "patch_icpp(3)%radius": 0.0015,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": 1.0e5,
            "patch_icpp(3)%alpha_rho(1)": 1000.0 * 1e-6,
            "patch_icpp(3)%alpha(1)": 1e-6,
            "patch_icpp(3)%alpha_rho(2)": 1000.0 * (1.0 - 1e-6),
            "patch_icpp(3)%alpha(2)": 1.0 - 1e-6,
            "patch_icpp(3)%tau_e(1)": 0.0,
            "patch_icpp(3)%tau_e(2)": 0.0,
            "patch_icpp(3)%tau_e(3)": 0.0,
            # Fluids: 1 = water (no shear stiffness), 2 = damageable solid
            "fluid_pp(1)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(1)%G": 0.0,
            "fluid_pp(2)%gamma": 1.0e00 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 4.4e00 * 6.0e08 / (4.4e00 - 1.0e00),
            "fluid_pp(2)%G": 1.0e9,
        }
    )
)
