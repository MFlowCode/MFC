import json
import math

# Free stream
p_inf = 16281
rho_inf = 0.416
c_inf = math.sqrt(1.4 * p_inf / rho_inf)
vel_inf = 2.5 * c_inf

# Jet inlet conditions
p_J = 1.0 * p_inf
vel_J = 41.0
rho_J = 1000.0

# Viscosities
mu_inf = 1.85e-5
mu_J = 2.67e-3

d_jet = 1.2e-3

# Simulation parameters
Ny = 3999
Nx = 11999
time_end = 5e-4
cfl = 0.6
eps = 1e-8

# Case dictionary

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational domain parameters
            "x_domain%beg": -40 * d_jet,
            "x_domain%end": 60 * d_jet,
            "y_domain%beg": 0.0,
            "y_domain%end": 35 * d_jet,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "cfl_adap_dt": "T",
            "t_stop": time_end,
            "t_save": time_end / 50,
            "n_start": 0,
            "cfl_target": cfl,
            # Simulation algorithm parameters
            "num_patches": 2,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "F",
            "time_stepper": 3,
            "recon_type": 2,
            "muscl_order": 2,
            "muscl_lim": 1,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "elliptic_smoothing": "T",
            "elliptic_smoothing_iters": 10,  # 50,
            "bc_x%beg": -17,
            "bc_x%end": -12,
            "bc_y%beg": -17,
            "bc_y%end": -3,
            "num_bc_patches": 0,
            # Formatted Database File Structures
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Free stream
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 50 * d_jet,
            "patch_icpp(1)%length_x": 1000 * d_jet,
            "patch_icpp(1)%length_y": 1000 * d_jet,
            "patch_icpp(1)%vel(1)": f"{vel_inf} * tanh(y / {d_jet / 4})",
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_inf,
            "patch_icpp(1)%alpha_rho(1)": rho_inf * (1.0 - eps),
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha_rho(2)": eps,
            "patch_icpp(1)%alpha(2)": eps,
            # Patch 2: Jet
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": d_jet,
            "patch_icpp(2)%length_y": d_jet,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": vel_J,
            "patch_icpp(2)%pres": p_J,
            "patch_icpp(2)%alpha_rho(1)": eps,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_J,
            "patch_icpp(2)%alpha(2)": 1.0 - eps,
            # Fluid properties
            "fluid_pp(1)%gamma": 1.00 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(2)%gamma": 1.0 / (6.12 - 1.0),
            "fluid_pp(2)%pi_inf": 6.12 * 3.43e8 / (6.12 - 1),
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1 / mu_inf,
            "fluid_pp(2)%Re(1)": 1 / mu_J,
            "surface_tension": "T",
            "cf_wrt": "T",
            "sigma": 0.072,
        }
    )
)
