import json
import math

Mu = 1.84e-05
gam_a = 1.4
gam_b = 1.1


x0 = 10e-06
p0 = 101325
rho0 = 1.0
c0 = math.sqrt(p0 / rho0)
patm = 1.0
rhoatm = 1.0

# air props
n_tait = 1.4
B_tait = 0.0


vf0 = 0.0

cact = math.sqrt(n_tait * (p0 + p0 * B_tait) / ((1 - vf0) * rho0))
cfl = 0.3
Nx = 400
Ny = 200
dx = 6.0e-03 / (x0 * float(Nx))
dt = cfl * dx * c0 / cact


vel = 1.5
vel_ac = vel * c0

Min = (n_tait + 1) * vel_ac + math.sqrt((n_tait + 1) ** 2 * vel_ac**2 + 16 * cact**2)
Min = Min / (4 * cact)
beta = (n_tait + 1) * Min**2 / (Min**2 * (n_tait - 1 + 2 * vf0) + 2 * (1 - vf0))
delta = (1 - vf0) + n_tait * Min**2 * (beta - 1) * (1 + B_tait) / beta

vel = Min * cact * (beta - 1) / (beta * c0)

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
            "bc_x%beg": -7,
            "bc_x%end": -8,
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
            "fd_order": 1,
            "omega_wrt(3)": "T",
            "parallel_io": "T",
            # Ambient State
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 3.0e-03 / x0,
            "patch_icpp(1)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(1)%length_x": 6.0e-03 / x0,
            "patch_icpp(1)%length_y": 3.0e-03 / x0,
            "patch_icpp(1)%alpha_rho(1)": rhoatm,
            "patch_icpp(1)%alpha(1)": 1,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": patm,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0e00,
            # Shocked State
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5e-03 / x0,
            "patch_icpp(2)%y_centroid": 1.50e-03 / x0,
            "patch_icpp(2)%length_x": 1.0e-03 / x0,
            "patch_icpp(2)%length_y": 3.0e-03 / x0,
            "patch_icpp(2)%alpha_rho(1)": beta * rhoatm,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": delta * patm,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0e00,
            "patch_icpp(2)%alter_patch(1)": "T",
            # CBC Inflow / Outflow
            "bc_x%grcbc_in": "T",
            "bc_x%grcbc_out": "F",
            "bc_x%grcbc_vel_out": "F",
            "bc_x%vel_in(1)": vel,
            "bc_x%vel_in(2)": 0,
            "bc_x%vel_in(3)": 0,
            "bc_x%pres_in": delta * patm,
            "bc_x%alpha_rho_in(1)": beta * rhoatm,
            "bc_x%alpha_in(1)": 1,
            "bc_x%vel_out(1)": vel,
            "bc_x%vel_out(2)": 0,
            "bc_x%vel_out(3)": 0,
            "bc_x%pres_out": 1.0,
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
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0e00 / (n_tait - 1.0e00),
            "fluid_pp(1)%pi_inf": n_tait * B_tait / (n_tait - 1.0),
            #'fluid_pp(1)%Re(1)' 	    : 67567,
        }
    )
)
