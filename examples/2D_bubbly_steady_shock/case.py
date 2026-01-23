import json
import math

Mu = 1.84e-05
gam_a = 1.4
gam_b = 1.1


x0 = 10e-06
p0 = 101325.0
rho0 = 1.0e03
c0 = math.sqrt(p0 / rho0)
patm = 1.0

# water props
n_tait = 7.1
B_tait = 3.43e05 / p0
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
# gamma_gas = gamma_n
gamma_gas = 1.4

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

vft = 1e-12
vf0 = 1e-03

cact = math.sqrt(n_tait * (p0 + p0 * B_tait) / ((1 - vf0) * rho0))
cfl = 0.3
Nx = 400
Ny = 200
dx = 6.0e-03 / (x0 * float(Nx))
dt = cfl * dx * c0 / cact

Min = 1.2
beta = (n_tait + 1) * Min**2 / (Min**2 * (n_tait - 1 + 2 * vf0) + 2 * (1 - vf0))
vel = Min * cact * (beta - 1) / beta
delta = (1 - vf0) + n_tait * Min**2 * (beta - 1) * (1 + B_tait) / beta

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 6.0e-03 / x0,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 3.0e-03 / x0,
            "cyl_coord": "F",
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 1000,  # 3000
            "t_step_save": 10,  # 10
            # Simulation Algorithm Parameters
            "num_patches": 2,
            # Use the 5 equation model
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            # Advect both volume fractions
            # 'adv_alphan'                   : 'T',
            # No need to ensure the volume fractions sum to unity at the end of each
            # time step
            "mpp_lim": "F",
            # Correct errors when computing speed of sound
            "mixture_err": "T",
            # Use TVD RK3 for time marching
            "time_stepper": 3,
            # Reconstruct the primitive variables to minimize spurious
            # Use WENO5
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "avg_state": 2,
            # Use the mapped WENO weights to maintain monotinicity
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            # Use the HLLC  Riemann solver
            "riemann_solver": 2,
            "wave_speeds": 1,
            "bc_x%beg": -6,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Set IB to True and add 1 patch
            "ib": "T",
            "num_ibs": 1,
            # Formatted Database Files Structure Parameters
            # Export primitive variables in double precision with parallel
            # I/O to minimize I/O computational time during large simulations
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            # 'ib_wrt'                       :'T',
            "fd_order": 1,
            "omega_wrt(3)": "T",
            "parallel_io": "T",
            # Ambient State
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 3.0e-03 / x0,
            "patch_icpp(1)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(1)%length_x": 6.0e-03 / x0,
            "patch_icpp(1)%length_y": 3.0e-03 / x0,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0),
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%vel(1)": 1.5,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0e00,
            # Shocked State
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5e-03 / x0,
            "patch_icpp(2)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(2)%length_x": 1.0e-03 / x0,
            "patch_icpp(2)%length_y": 3.0e-03 / x0,
            "patch_icpp(2)%alpha_rho(1)": beta,
            "patch_icpp(2)%alpha(1)": beta * vf0,
            "patch_icpp(2)%vel(1)": vel / c0,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": delta,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0e00,
            "patch_icpp(2)%alter_patch(1)": "T",
            # Patch: Cylinder Immersed Boundary
            "patch_ib(1)%geometry": 4,
            "patch_ib(1)%x_centroid": 1.5e-03 / x0,
            "patch_ib(1)%y_centroid": 1.5e-03 / x0,
            "patch_ib(1)%c": 1.0e-03 / x0,
            "patch_ib(1)%t": 0.15,
            "patch_ib(1)%p": 0.4,
            "patch_ib(1)%m": 0.02,
            "patch_ib(1)%slip": "F",
            "patch_ib(1)%theta": 15,
            "pref": p0,
            "rhoref": rho0,
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
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 2,
            "polytropic": "T",
            "polydisperse": "F",
            "poly_sigma": 0.3,
            "thermal": 3,
            "R0ref": x0,
            "nb": 1,
            "Ca": Ca,
            "Web": We,
            "Re_inv": Re_inv,
            "qbmm": "F",
            "dist_type": 1,
            "sigR": 0.1,
            "sigV": 0.3,
            "rhoRV": 0.0,
        }
    )
)
