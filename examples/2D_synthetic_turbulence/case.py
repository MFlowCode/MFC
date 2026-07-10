import json
import math

Mu = 1.84e-05
gam_a = 1.4
gam_b = 1.1

free_stream_vel = 0.5

n_shells = 3
k_shells = [
    2 * math.pi / 0.6,  # shell 1 — large eddies
    2 * math.pi / 0.30,  # shell 2 — medium eddies
    2 * math.pi / 0.15,
]  # shell 3 — small eddies
amp_shells = [
    0.1 * free_stream_vel,  # m/s  ~5 % of U_inf
    0.07 * free_stream_vel,  # m/s  ~3 %
    0.05 * free_stream_vel,
]  # m/s  ~1 %
waves_per_shell = [4, 8, 12]  # random directions per shell

# Gaussian source zone
src_x, src_y = 1.0, 3.0  # center
src_Lx, src_Ly = 1.0, 4.0  # full extents fed to synth_L

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            # axial direction
            "x_domain%beg": 0.0e00,
            "x_domain%end": 12.0,
            # r direction
            "y_domain%beg": 0.0e00,
            "y_domain%end": 6.0,
            "cyl_coord": "F",
            "m": 500,
            "n": 250,
            "p": 0,
            "dt": 2.0e-3,
            "t_step_start": 0,
            "t_step_stop": 20000,  # 3000
            "t_step_save": 80,  # 10
            # Simulation Algorithm Parameters
            # Only one patches are necessary, the air tube
            "num_patches": 1,
            # Use the 5 equation model
            "model_eqns": 2,
            # 6 equations model does not need the K \div(u) term
            "alt_soundspeed": "F",
            # One fluids: air
            "num_fluids": 2,
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
            "mp_weno": "T",
            # Use the HLLC  Riemann solver
            "riemann_solver": 2,
            "wave_speeds": 1,
            # We use reflective boundary conditions at octant edges and
            # non-reflective boundary conditions at the domain edges
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            # Set IB to True and add 1 patch
            "ib": "T",
            "num_ibs": 1,
            "fd_order": 2,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "parallel_io": "T",
            # Patch: Constant Tube filled with air
            # Specify the cylindrical air tube grid geometry
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 6.0,
            # Uniform medium density, centroid is at the center of the domain
            "patch_icpp(1)%y_centroid": 3.0,
            "patch_icpp(1)%length_x": 12.0,
            "patch_icpp(1)%length_y": 6.0,
            # Specify the patch primitive variables
            "patch_icpp(1)%vel(1)": free_stream_vel,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 1.0e00,
            "patch_icpp(1)%alpha_rho(1)": 0.8e00,
            "patch_icpp(1)%alpha(1)": 0.8e00,
            "patch_icpp(1)%alpha_rho(2)": 0.2e00,
            "patch_icpp(1)%alpha(2)": 0.2e00,
            # Patch: Airfoil Immersed Boundary
            "patch_ib(1)%geometry": 4,
            "patch_ib(1)%x_centroid": 6.0,
            "patch_ib(1)%y_centroid": 3.0,
            "patch_ib(1)%airfoil_id": 1,
            "ib_airfoil(1)%c": 2.0,
            "ib_airfoil(1)%t": 0.15,
            "ib_airfoil(1)%p": 0.4,
            "ib_airfoil(1)%m": 0.02,
            "patch_ib(1)%angles(3)": -0.5235987756,  # 30 degrees clockwise rotation, in radians
            "patch_ib(1)%angular_vel(3)": "15.0 * 0.1 * pi * sin(0.1 * pi * t) * pi / 180.",
            "patch_ib(1)%moving_ibm": 1,
            # Fluids Physical Parameters
            # Use the same stiffness as the air bubble
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50 (Not 1.40)
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(2)%gamma": 1.0e00 / (gam_b - 1.0e00),  # 2.50 (Not 1.40)
            "fluid_pp(2)%pi_inf": 0,
            # -- Synthetic turbulence --
            "synthetic_turbulence": "T",
            "synth_seed": 42,
            "synth_n_shells": n_shells,
            "synth_U_inf": free_stream_vel,
            "num_turbulent_sources": 1,
            # Energy shells (wave-number magnitudes, amplitudes, wave counts)
            "synth_k_shell(1)": k_shells[0],
            "synth_k_shell(2)": k_shells[1],
            "synth_k_shell(3)": k_shells[2],
            "synth_amp_shell(1)": amp_shells[0],
            "synth_amp_shell(2)": amp_shells[1],
            "synth_amp_shell(3)": amp_shells[2],
            "synth_n_waves_per_shell(1)": waves_per_shell[0],
            "synth_n_waves_per_shell(2)": waves_per_shell[1],
            "synth_n_waves_per_shell(3)": waves_per_shell[2],
            # Gaussian forcing zone (source 1)
            "turb_pos(1,1)": src_x,
            "turb_pos(1,2)": src_y,
            "turb_pos(1,3)": 0.0,  # z not used in 2-D
            "synth_L(1,1)": src_Lx,
            "synth_L(1,2)": src_Ly,
            "synth_L(1,3)": 1.0,  # z extent irrelevant in 2-D; set nonzero to avoid div-by-zero
        }
    )
)
