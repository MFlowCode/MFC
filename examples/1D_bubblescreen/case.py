#!/usr/bin/env python3
import math, json

x0 = 10.0e-06
p0 = 101325.0
rho0 = 1.0e03
c0 = math.sqrt(p0 / rho0)
patm = 1.0

# water props
n_tait = 7.1
B_tait = 306.0e06 / p0
mul0 = 1.002e-03  # viscosity
ss = 0.07275  # surface tension
pv = 2.3388e03  # vapor pressure

gamma_v = 1.33
M_v = 18.02
mu_v = 0.8816e-05
k_v = 0.019426

# air props
gamma_n = 1.4
M_n = 28.97
mu_n = 1.8e-05
k_n = 0.02556

# air props
gamma_gas = gamma_n

# reference bubble size
R0ref = 10.0e-06

pa = 0.1 * 1.0e06 / 101325.0

# Characteristic velocity
uu = math.sqrt(p0 / rho0)
# Cavitation number
Ca = (p0 - pv) / (rho0 * (uu**2.0))
# Weber number
We = rho0 * (uu**2.0) * R0ref / ss
# Inv. bubble Reynolds number
Re_inv = mul0 / (rho0 * uu * R0ref)

# IC setup
vf0 = 0.00004
n0 = vf0 / (math.pi * 4.0e00 / 3.0e00)

cact = 1475.0
t0 = x0 / c0

nbubbles = 1
myr0 = R0ref

cfl = 0.1
Nx = 100
Ldomain = 20.0e-03
L = Ldomain / x0
dx = L / float(Nx)
dt = cfl * dx * c0 / cact
Lpulse = 0.3 * Ldomain
Tpulse = Lpulse / cact
Tfinal = 0.25 * 10.0 * Tpulse * c0 / x0
Nt = int(Tfinal / dt)

Nfiles = 20.0
Nout = int(math.ceil(Nt / Nfiles))
Nt = int(Nout * Nfiles)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -10.0e-03 / x0,
            "x_domain%end": 10.0e-03 / x0,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nout,
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
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            "fd_order": 1,
            # 'schlieren_wrt'                :'T',
            "probe_wrt": "T",
            "num_probes": 1,
            "probe(1)%x": 0.0,
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%length_x": 20.0e-03 / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": patm,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - 1.0e-12) * 1.0e03 / rho0,
            "patch_icpp(1)%alpha(1)": 1.0e-12,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0e00,
            # Patch 2 Screen
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%length_x": 5.0e-03 / x0,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": patm,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - vf0) * 1.0e03 / rho0,
            "patch_icpp(2)%alpha(1)": vf0,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0e00,
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0e00 / (n_tait - 1.0e00),
            "fluid_pp(1)%pi_inf": n_tait * B_tait / (n_tait - 1.0),
            "fluid_pp(1)%mul0": mul0,
            "fluid_pp(1)%ss": ss,
            "fluid_pp(1)%pv": pv,
            "fluid_pp(1)%gamma_v": gamma_v,
            "fluid_pp(1)%M_v": M_v,
            "fluid_pp(1)%mu_v": mu_v,
            "fluid_pp(1)%k_v": k_v,
            # Last fluid_pp is always reserved for bubble gas state
            # if applicable
            "fluid_pp(2)%gamma": 1.0 / (gamma_gas - 1.0),
            "fluid_pp(2)%pi_inf": 0.0e00,
            "fluid_pp(2)%gamma_v": gamma_n,
            "fluid_pp(2)%M_v": M_n,
            "fluid_pp(2)%mu_v": mu_n,
            "fluid_pp(2)%k_v": k_n,
            # Non-polytropic gas compression model AND/OR Tait EOS
            "pref": p0,
            "rhoref": rho0,
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 3,
            "polytropic": "T",
            "polydisperse": "F",
            # 'poly_sigma'                   : 0.3,
            "thermal": 3,
            "R0ref": myr0,
            "nb": 1,
            "Ca": Ca,
            "Web": We,
            "Re_inv": Re_inv,
            # Acoustic source
            "acoustic_source": "T",
            "num_source": 1,
            "acoustic(1)%support": 1,
            "acoustic(1)%loc(1)": -5.0e-03 / x0,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%dir": 1.0,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%mag": pa,
            "acoustic(1)%wavelength": (1.0 / (300000.0)) * cact / x0,
        }
    )
)
