"""
MFC Parameter Definitions (Compact).

Single file containing all ~3,300 parameter definitions using loops.
This replaces the definitions/ directory.
"""

from .schema import ParamDef, ParamType, Stage
from .registry import REGISTRY

# Index limits
NP, NF, NI, NA, NPR, NB = 10, 10, 10, 4, 10, 10  # patches, fluids, ibs, acoustic, probes, bc_patches

# Parameters that can be hard-coded for GPU case optimization
CASE_OPT_PARAMS = {
    "mapped_weno", "wenoz", "teno", "wenoz_q", "nb", "weno_order",
    "num_fluids", "mhd", "relativity", "igr_order", "viscous",
    "igr_iter_solver", "igr", "igr_pres_lim", "recon_type",
    "muscl_order", "muscl_lim"
}

# Parameter constraints (choices, min, max)
CONSTRAINTS = {
    # Reconstruction
    "weno_order": {"choices": [0, 1, 3, 5, 7]},  # 0 for MUSCL mode
    "recon_type": {"choices": [1, 2]},  # 1=WENO, 2=MUSCL
    "muscl_order": {"choices": [1, 2]},
    "muscl_lim": {"choices": [1, 2, 3, 4, 5]},  # minmod, MC, Van Albada, Van Leer, SUPERBEE

    # Time stepping
    "time_stepper": {"choices": [1, 2, 3]},  # 1=Euler, 2=TVD-RK2, 3=TVD-RK3

    # Riemann solver
    "riemann_solver": {"choices": [1, 2, 3, 4, 5]},  # HLL, HLLC, Exact, HLLD, LF
    "wave_speeds": {"choices": [1, 2]},  # direct, pressure
    "avg_state": {"choices": [1, 2]},  # Roe, arithmetic

    # Model equations
    "model_eqns": {"choices": [1, 2, 3, 4]},  # gamma-law, 5-eq, 6-eq, 4-eq

    # Bubbles
    "bubble_model": {"choices": [1, 2, 3]},  # Gilmore, Keller-Miksis, RP

    # Output
    "format": {"choices": [1, 2]},  # Silo, binary
    "precision": {"choices": [1, 2]},  # single, double

    # Counts (must be positive)
    "num_fluids": {"min": 1, "max": 10},
    "num_patches": {"min": 0, "max": 10},
    "num_ibs": {"min": 0, "max": 10},
    "m": {"min": 0},
    "n": {"min": 0},
    "p": {"min": 0},
}

# Parameter dependencies (requires, recommends)
DEPENDENCIES = {
    "bubbles_euler": {
        "when_true": {
            "recommends": ["nb", "R0ref", "polytropic"],
        }
    },
    "viscous": {
        "when_true": {
            "recommends": ["fluid_pp(1)%Re(1)"],
        }
    },
    "polydisperse": {
        "when_true": {
            "requires": ["nb"],
        }
    },
    "chemistry": {
        "when_true": {
            "requires": ["cantera_file"],
        }
    },
    "qbmm": {
        "when_true": {
            "recommends": ["bubbles_euler"],
        }
    },
    "ib": {
        "when_true": {
            "requires": ["num_ibs"],
        }
    },
    "acoustic_source": {
        "when_true": {
            "requires": ["num_source"],
        }
    },
    "probe_wrt": {
        "when_true": {
            "requires": ["num_probes"],
        }
    },
}

def _r(name, ptype, stages):
    """Register a parameter."""
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        stages=stages,
        case_optimization=(name in CASE_OPT_PARAMS),
        constraints=CONSTRAINTS.get(name),
        dependencies=DEPENDENCIES.get(name),
    ))

def _load():  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Load all parameter definitions."""
    C, P, S, POST = Stage.COMMON, Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS
    INT, REAL, LOG, STR = ParamType.INT, ParamType.REAL, ParamType.LOG, ParamType.STR
    A_REAL = ParamType.ANALYTIC_REAL

    # === COMMON ===
    # Truly common params (used in multiple stages via inheritance)
    for n in ["m", "n", "p", "model_eqns", "num_fluids", "weno_order",
              "thermal", "precision", "relax_model", "igr_order",
              "recon_type", "muscl_order", "num_bc_patches"]:
        _r(n, INT, {C})
    for n in ["pref", "poly_sigma", "R0ref", "nb", "rhoref", "sigma", "Bx0",
              "Web", "Re_inv", "Ca"]:
        _r(n, REAL, {C})
    for n in ["cyl_coord", "mhd", "hypoelasticity", "hyperelasticity", "parallel_io",
              "polytropic", "mpp_lim", "bubbles_euler", "polydisperse", "file_per_process",
              "relax", "adv_n", "cfl_adap_dt", "cfl_const_dt", "chemistry", "relativity",
              "cont_damage", "igr", "down_sample"]:
        _r(n, LOG, {C})
    for n in ["case_dir", "cantera_file"]:
        _r(n, STR, {C})

    # === PRE_PROCESS ===
    for n in ["t_step_old", "t_step_start", "perturb_flow_fluid", "perturb_sph_fluid",
              "dist_type", "mixlayer_perturb_nk", "elliptic_smoothing_iters", "num_patches"]:
        _r(n, INT, {P})
    for n in ["mixlayer_vel_coef", "mixlayer_domain", "mixlayer_perturb_k0",
              "perturb_flow_mag", "fluid_rho", "sigR", "sigV", "rhoRV"]:
        _r(n, REAL, {P})
    for n in ["old_grid", "old_ic", "mixlayer_vel_profile", "mixlayer_perturb",
              "perturb_flow", "perturb_sph", "stretch_x", "stretch_y", "stretch_z",
              "cfl_dt", "pre_stress", "elliptic_smoothing", "simplex_perturb"]:
        _r(n, LOG, {P})
    for d in ["x", "y", "z"]:
        _r(f"{d}_a", REAL, {P})
        _r(f"{d}_b", REAL, {P})
        _r(f"a_{d}", REAL, {P})
        _r(f"loops_{d}", INT, {P})
        # Domain and BC params (PRE+SIM for both)
        _r(f"{d}_domain%beg", REAL, {P, S})
        _r(f"{d}_domain%end", REAL, {P, S})
        _r(f"bc_{d}%beg", INT, {P, S})
        _r(f"bc_{d}%end", INT, {P, S})

    # PRE + SIM
    for n in ["n_start_old"]:
        _r(n, INT, {P})
    for n in ["palpha_eps", "ptgalpha_eps", "pi_fac"]:
        _r(n, REAL, {P, S})
    for n in ["viscous"]:
        _r(n, LOG, {P, S})

    # PRE + SIM + POST
    for n in ["n_start", "num_ibs"]:
        _r(n, INT, {P, S, POST})
    for n in ["qbmm", "ib", "surface_tension", "bubbles_lagrange", "fft_wrt"]:
        _r(n, LOG, {P, S, POST})

    # === SIMULATION ===
    for n in ["time_stepper", "t_step_start", "t_step_stop", "t_step_save", "t_step_print",
              "riemann_solver", "wave_speeds", "avg_state", "fd_order", "num_probes",
              "bubble_model", "num_source", "num_integrals", "adap_dt_max_iters", "low_Mach",
              "num_igr_iters", "num_igr_warm_start_iters", "igr_iter_solver", "muscl_lim",
              "nv_uvm_igr_temps_on_gpu", "t_step_old"]:
        _r(n, INT, {S})
    for n in ["dt", "cfl_target", "cfl_max", "t_tol", "weno_eps", "teno_CT", "wenoz_q",
              "adap_dt_tol", "t_stop", "t_save", "tau_star", "cont_damage_s", "alpha_bar",
              "alf_factor", "ic_eps", "ic_beta"]:
        _r(n, REAL, {S})
    for n in ["mapped_weno", "wenoz", "teno", "weno_Re_flux", "weno_avg", "mp_weno",
              "run_time_info", "prim_vars_wrt", "alt_soundspeed", "null_weights",
              "mixture_err", "probe_wrt", "acoustic_source", "integral_wrt", "rdma_mpi",
              "adap_dt", "powell", "igr_pres_lim", "int_comp", "nv_uvm_out_of_core",
              "nv_uvm_pref_gpu"]:
        _r(n, LOG, {S})
    # cfl_target also POST_PROCESS
    _r("cfl_target", REAL, {POST})

    # Body force (SIM)
    for d in ["x", "y", "z"]:
        for v in ["k", "w", "p", "g"]:
            _r(f"{v}_{d}", REAL, {S})
        _r(f"bf_{d}", LOG, {S})

    # === POST_PROCESS ===
    for n in ["format", "flux_lim", "t_step_start", "t_step_stop", "t_step_save"]:
        _r(n, INT, {POST})
    for d in ["x", "y", "z"]:
        _r(f"bc_{d}%beg", INT, {POST})
        _r(f"bc_{d}%end", INT, {POST})
    for n in ["schlieren_alpha"]:
        _r(n, REAL, {POST})
    for n in ["schlieren_wrt", "alpha_rho_wrt", "rho_wrt", "mom_wrt", "vel_wrt",
              "flux_wrt", "E_wrt", "pres_wrt", "alpha_wrt", "kappa_wrt", "gamma_wrt",
              "heat_ratio_wrt", "pi_inf_wrt", "pres_inf_wrt", "cons_vars_wrt", "c_wrt",
              "omega_wrt", "qm_wrt", "liutex_wrt", "cf_wrt", "sim_data", "output_partial_domain"]:
        _r(n, LOG, {POST})
    for d in ["x", "y", "z"]:
        _r(f"{d}_output%beg", REAL, {POST})
        _r(f"{d}_output%end", REAL, {POST})
    for n in ["prim_vars_wrt", "alt_soundspeed", "mixture_err", "t_stop", "t_save", "fd_order"]:
        _r(n, LOG if "wrt" in n or n in ["alt_soundspeed", "mixture_err"] else (INT if n == "fd_order" else REAL), {POST})
    for v in ["lag_header", "lag_txt_wrt", "lag_db_wrt", "lag_id_wrt", "lag_pos_wrt",
              "lag_pos_prev_wrt", "lag_vel_wrt", "lag_rad_wrt", "lag_rvel_wrt",
              "lag_r0_wrt", "lag_rmax_wrt", "lag_rmin_wrt", "lag_dphidt_wrt",
              "lag_pres_wrt", "lag_mv_wrt", "lag_mg_wrt", "lag_betaT_wrt", "lag_betaC_wrt"]:
        _r(v, LOG, {POST})

    # === INDEXED PARAMETERS ===

    # patch_icpp (PRE_PROCESS) - 10 patches
    for i in range(1, NP + 1):
        px = f"patch_icpp({i})%"
        for a in ["geometry", "smooth_patch_id", "hcid", "model_spc"]:
            _r(f"{px}{a}", INT, {P})
        for a in ["smoothen", "alter_patch"] if i >= 2 else ["smoothen"]:
            _r(f"{px}{a}", LOG, {P})
        for a in ["radius", "radii", "epsilon", "beta", "normal", "alpha_rho",
                  "non_axis_sym", "smooth_coeff", "rho", "vel", "alpha", "gamma",
                  "pi_inf", "r0", "v0", "p0", "m0", "cv", "qv", "qvp", "model_threshold"]:
            _r(f"{px}{a}", REAL, {P})
        for j in range(2, 10):
            _r(f"{px}a({j})", REAL, {P})
        for a in ["pres", "Bx", "By", "Bz", "cf_val"]:
            _r(f"{px}{a}", A_REAL, {P})
        for j in range(1, 101):
            _r(f"{px}Y({j})", A_REAL, {P})
        _r(f"{px}model_filepath", STR, {P})
        for t in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _r(f"{px}model_{t}({j})", REAL, {P})
        for d in ["x", "y", "z"]:
            _r(f"{px}{d}_centroid", REAL, {P})
            _r(f"{px}length_{d}", REAL, {P})
        for j in range(1, 4):
            _r(f"{px}radii({j})", REAL, {P})
            _r(f"{px}normal({j})", REAL, {P})
            _r(f"{px}vel({j})", A_REAL, {P})
        for f in range(1, NF + 1):
            _r(f"{px}alpha({f})", A_REAL, {P})
            _r(f"{px}alpha_rho({f})", A_REAL, {P})
        for j in range(1, 7):
            _r(f"{px}tau_e({j})", A_REAL, {P})
        if i >= 2:
            for j in range(1, i):
                _r(f"{px}alter_patch({j})", LOG, {P})

    # fluid_pp (PRE + SIM + POST)
    for f in range(1, NF + 1):
        px = f"fluid_pp({f})%"
        for a in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                  "mu_v", "k_v", "cp_v", "G", "cv", "qv", "qvp", "D_v"]:
            _r(f"{px}{a}", REAL, {P, S, POST})
        for j in [1, 2]:
            _r(f"{px}Re({j})", REAL, {S})

    # bub_pp (PRE + SIM + POST)
    for a in ["R0ref", "p0ref", "rho0ref", "T0ref", "ss", "pv", "vd",
              "mu_l", "mu_v", "mu_g", "gam_v", "gam_g",
              "M_v", "M_g", "k_v", "k_g", "cp_v", "cp_g", "R_v", "R_g"]:
        _r(f"bub_pp%{a}", REAL, {P, S, POST})

    # patch_ib (PRE + SIM)
    for i in range(1, NI + 1):
        px = f"patch_ib({i})%"
        for a in ["geometry", "moving_ibm"]:
            _r(f"{px}{a}", INT, {P, S})
        for a, pt in [("radius", REAL), ("theta", REAL), ("slip", LOG), ("c", REAL),
                      ("p", REAL), ("t", REAL), ("m", REAL), ("mass", REAL)]:
            _r(f"{px}{a}", pt, {P, S})
        for j in range(1, 4):
            _r(f"{px}angles({j})", REAL, {P, S})
        for d in ["x", "y", "z"]:
            _r(f"{px}{d}_centroid", REAL, {P, S})
            _r(f"{px}length_{d}", REAL, {P, S})
        for a, pt in [("model_filepath", STR), ("model_spc", INT), ("model_threshold", REAL)]:
            _r(f"{px}{a}", pt, {P})
        for t in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _r(f"{px}model_{t}({j})", REAL, {P})
        for j in range(1, 4):
            _r(f"{px}vel({j})", REAL, {P, S})  # Both PRE_PROCESS and SIMULATION
            _r(f"{px}angular_vel({j})", REAL, {P, S})  # Both PRE_PROCESS and SIMULATION

    # acoustic sources (SIM)
    for i in range(1, NA + 1):
        px = f"acoustic({i})%"
        for a in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            _r(f"{px}{a}", INT, {S})
        _r(f"{px}dipole", LOG, {S})
        for a in ["mag", "length", "height", "wavelength", "frequency",
                  "gauss_sigma_dist", "gauss_sigma_time", "npulse",
                  "dir", "delay", "foc_length", "aperture",
                  "element_spacing_angle", "element_polygon_ratio",
                  "rotate_angle", "bb_bandwidth", "bb_lowest_freq"]:
            _r(f"{px}{a}", REAL, {S})
        for j in range(1, 4):
            _r(f"{px}loc({j})", REAL, {S})

    # probes (SIM)
    for i in range(1, NPR + 1):
        for d in ["x", "y", "z"]:
            _r(f"probe({i})%{d}", REAL, {S})

    # integrals (SIM)
    for i in range(1, 6):
        for d in ["x", "y", "z"]:
            _r(f"integral({i})%{d}min", REAL, {S})
            _r(f"integral({i})%{d}max", REAL, {S})

    # Extended BC (PRE + SIM)
    for d in ["x", "y", "z"]:
        px = f"bc_{d}%"
        for a in ["vb1", "vb2", "vb3", "ve1", "ve2", "ve3", "pres_in", "pres_out"]:
            _r(f"{px}{a}", REAL, {P, S})
        for a in ["grcbc_in", "grcbc_out", "grcbc_vel_out"]:
            _r(f"{px}{a}", LOG, {P, S})
        for f in range(1, NF + 1):
            _r(f"{px}alpha_rho_in({f})", REAL, {P, S})
            _r(f"{px}alpha_in({f})", REAL, {P, S})
        for j in range(1, 4):
            _r(f"{px}vel_in({j})", REAL, {P, S})
            _r(f"{px}vel_out({j})", REAL, {P, S})

    # patch_bc (PRE)
    for i in range(1, NB + 1):
        px = f"patch_bc({i})%"
        for a in ["geometry", "type", "dir", "loc"]:
            _r(f"{px}{a}", INT, {P})
        for j in range(1, 4):
            _r(f"{px}centroid({j})", REAL, {P})
            _r(f"{px}length({j})", REAL, {P})
        _r(f"{px}radius", REAL, {P})

    # simplex_params (PRE)
    for f in range(1, NF + 1):
        _r(f"simplex_params%perturb_dens({f})", LOG, {P})
        _r(f"simplex_params%perturb_dens_freq({f})", REAL, {P})
        _r(f"simplex_params%perturb_dens_scale({f})", REAL, {P})
        for j in range(1, 4):
            _r(f"simplex_params%perturb_dens_offset({f}, {j})", REAL, {P})
    for d in range(1, 4):
        _r(f"simplex_params%perturb_vel({d})", LOG, {P})
        _r(f"simplex_params%perturb_vel_freq({d})", REAL, {P})
        _r(f"simplex_params%perturb_vel_scale({d})", REAL, {P})
        for j in range(1, 4):
            _r(f"simplex_params%perturb_vel_offset({d},{j})", REAL, {P})

    # lag_params (SIM)
    for a in ["heatTransfer_model", "massTransfer_model", "pressure_corrector",
              "write_bubbles", "write_bubbles_stats"]:
        _r(f"lag_params%{a}", LOG, {S})
    for a in ["solver_approach", "cluster_type", "smooth_type", "nBubs_glb"]:
        _r(f"lag_params%{a}", INT, {S})
    for a in ["epsilonb", "valmaxvoid", "charwidth", "c0", "rho0", "T0", "x0", "Thost"]:
        _r(f"lag_params%{a}", REAL, {S})

    # chem_params (SIM)
    for a in ["diffusion", "reactions"]:
        _r(f"chem_params%{a}", LOG, {S})
    for a in ["gamma_method", "transport_model"]:
        _r(f"chem_params%{a}", INT, {S})

    # Per-fluid post-process output (POST)
    for f in range(1, NF + 1):
        _r(f"schlieren_alpha({f})", REAL, {POST})
        for a in ["alpha_rho_wrt", "alpha_wrt", "kappa_wrt", "alpha_rho_e_wrt"]:
            _r(f"{a}({f})", LOG, {POST})
    for j in range(1, 4):
        for a in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
            _r(f"{a}({j})", LOG, {POST})

    # chem_wrt (POST)
    for j in range(1, 101):
        _r(f"chem_wrt_Y({j})", LOG, {POST})
    _r("chem_wrt_T", LOG, {POST})

    # fluid_rho (PRE)
    for f in range(1, NF + 1):
        _r(f"fluid_rho({f})", REAL, {P})


# Load definitions when module imported
_load()
