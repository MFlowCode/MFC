#!/usr/bin/env python3
# Benchmark viscosity_weno_Re_flux_T_weno_order_5_bubbles_T_bubble_mode_3_acoustic_T
# Additional Benchmarked Features
# - viscosity enabled
# - weno_Re_flux : T
# - weno_order : 5
# - bubbles_euler : T
# - bubble_model : 3
# - acoustic_source : T

import json, math, argparse

parser = argparse.ArgumentParser(prog="Benchmarking Case 2", description="This MFC case was created for the purposes of benchmarking MFC.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory per GPU.")

ARGS = vars(parser.parse_args())
DICT = ARGS["mfc"]

size = 1 if DICT["gpu"] else 0

ppg = 8000000 / 16.0
procs = DICT["nodes"] * DICT["tasks_per_node"]
ncells = math.floor(ppg * procs * ARGS["gbpp"])
s = math.floor((ncells / 2.0) ** (1 / 3))
Nx, Ny, Nz = 2 * s, s, s

x0 = 10.0e-04
y0 = 10.0e-04
z0 = 10.0e-04
p0 = 101325.0
p0ref = p0
rho0 = 1.0e03
rho0ref = rho0
u0 = math.sqrt(p0 / rho0)
p0ext = 101325.0

# water props
gam_l = 7.1
pi_inf_l = 306.0e06
mu_l = 1.002e-03  # viscosity
ss = 0.07275  # surface tension
pv = 2.3388e03  # vapor pressure

gamma_v = 1.33
M_v = 18.02
mu_v = 0.8816e-05
k_v = 0.019426

# air props
gam_g = 1.4

# reference bubble size
R0ref = 10.0e-06

pa = 0.1 * 1.0e06

# IC setup
vf0 = 0.00004
n0 = vf0 / (math.pi * 4.0e00 / 3.0e00)

cact = 1475.0
t0 = x0 / u0

cfl = 0.01
Ldomain = 20.0e-03
L = Ldomain / x0
dx = L / float(Nx)
dt = cfl * dx * u0 / cact
Lpulse = 0.3 * Ldomain
Tpulse = Lpulse / cact

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": -10.0e-03 / x0,
            "x_domain%end": 10.0e-03 / x0,
            "y_domain%beg": -5.0e-03 / y0,
            "y_domain%end": 5.0e-03 / y0,
            "z_domain%beg": -5.0e-03 / z0,
            "z_domain%end": 5.0e-03 / z0,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": int(15 * (5 * size + 5)),
            "t_step_save": int(15 * (5 * size + 5)),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 20.0e-03 / x0,
            "patch_icpp(1)%length_y": 10.0e-03 / y0,
            "patch_icpp(1)%length_z": 10.0e-03 / z0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": p0ext / p0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - 1.0e-12) * 1.0e03 / rho0,
            "patch_icpp(1)%alpha(1)": 1.0e-12,
            "patch_icpp(1)%r0": R0ref / x0,
            "patch_icpp(1)%v0": 0.0e00,
            # Patch 2 Screen
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%length_x": 5.0e-03 / x0,
            "patch_icpp(2)%length_y": 10.0e-03 / y0,
            "patch_icpp(2)%length_z": 10.0e-03 / z0,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": p0ext / p0,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - vf0) * 1.0e03 / rho0,
            "patch_icpp(2)%alpha(1)": vf0,
            "patch_icpp(2)%r0": R0ref / x0,
            "patch_icpp(2)%v0": 0.0e00,
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0e00 / (gam_l - 1.0e00),
            "fluid_pp(1)%pi_inf": gam_l * (pi_inf_l / p0) / (gam_l - 1.0),
            "fluid_pp(1)%Re(1)": 1e3,
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 3,
            "polytropic": "T",
            "polydisperse": "F",
            "thermal": 3,
            "nb": 1,
            "bub_pp%R0ref": R0ref / x0,
            "bub_pp%p0ref": p0ref / p0,
            "bub_pp%rho0ref": rho0ref / rho0,
            "bub_pp%ss": ss / (rho0 * x0 * u0 * u0),
            "bub_pp%pv": pv / p0,
            "bub_pp%mu_l": mu_l / (rho0 * x0 * u0),
            "bub_pp%gam_g": gam_g,
            # Acoustic source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 3,
            "acoustic(1)%loc(1)": -5.0e-03 / x0,
            "acoustic(1)%loc(2)": 0.0,
            "acoustic(1)%length": 1.0 / x0,
            "acoustic(1)%height": 1.0 / x0,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%dir": 1.0,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%mag": pa / p0,
            "acoustic(1)%wavelength": (1.0 / (300000.0)) * cact / x0,
        }
    )
)
