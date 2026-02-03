"""
MFC Parameter Definitions (Compact).

Single file containing all ~3,300 parameter definitions using loops.
This replaces the definitions/ directory.
"""

from typing import Dict, Any
from .schema import ParamDef, ParamType
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


# ============================================================================
# Schema Validation for Constraints and Dependencies
# ============================================================================
# Uses rapidfuzz for "did you mean?" suggestions when typos are detected

_VALID_CONSTRAINT_KEYS = {"choices", "min", "max"}
_VALID_DEPENDENCY_KEYS = {"when_true", "when_set"}
_VALID_CONDITION_KEYS = {"requires", "recommends"}


def _validate_constraint(param_name: str, constraint: Dict[str, Any]) -> None:
    """Validate a constraint dict has valid keys with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error  # pylint: disable=import-outside-toplevel

    invalid_keys = set(constraint.keys()) - _VALID_CONSTRAINT_KEYS
    if invalid_keys:
        # Get suggestion for the first invalid key
        first_invalid = next(iter(invalid_keys))
        raise ValueError(
            invalid_key_error(
                f"constraint for '{param_name}'",
                first_invalid,
                _VALID_CONSTRAINT_KEYS
            )
        )

    # Validate types
    if "choices" in constraint and not isinstance(constraint["choices"], list):
        raise ValueError(f"Constraint 'choices' for '{param_name}' must be a list")
    if "min" in constraint and not isinstance(constraint["min"], (int, float)):
        raise ValueError(f"Constraint 'min' for '{param_name}' must be a number")
    if "max" in constraint and not isinstance(constraint["max"], (int, float)):
        raise ValueError(f"Constraint 'max' for '{param_name}' must be a number")


def _validate_dependency(param_name: str, dependency: Dict[str, Any]) -> None:
    """Validate a dependency dict has valid structure with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error  # pylint: disable=import-outside-toplevel

    invalid_keys = set(dependency.keys()) - _VALID_DEPENDENCY_KEYS
    if invalid_keys:
        first_invalid = next(iter(invalid_keys))
        raise ValueError(
            invalid_key_error(
                f"dependency for '{param_name}'",
                first_invalid,
                _VALID_DEPENDENCY_KEYS
            )
        )

    for condition_key in ["when_true", "when_set"]:
        if condition_key in dependency:
            condition = dependency[condition_key]
            if not isinstance(condition, dict):
                raise ValueError(
                    f"Dependency '{condition_key}' for '{param_name}' must be a dict"
                )
            invalid_cond_keys = set(condition.keys()) - _VALID_CONDITION_KEYS
            if invalid_cond_keys:
                first_invalid = next(iter(invalid_cond_keys))
                raise ValueError(
                    invalid_key_error(
                        f"condition in '{condition_key}' for '{param_name}'",
                        first_invalid,
                        _VALID_CONDITION_KEYS
                    )
                )
            for req_key in ["requires", "recommends"]:
                if req_key in condition and not isinstance(condition[req_key], list):
                    raise ValueError(
                        f"Dependency '{condition_key}/{req_key}' for '{param_name}' "
                        "must be a list"
                    )


def _validate_all_constraints(constraints: Dict[str, Dict]) -> None:
    """Validate all constraint definitions."""
    for param_name, constraint in constraints.items():
        _validate_constraint(param_name, constraint)


def _validate_all_dependencies(dependencies: Dict[str, Dict]) -> None:
    """Validate all dependency definitions."""
    for param_name, dependency in dependencies.items():
        _validate_dependency(param_name, dependency)


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
            "recommends": ["nb", "polytropic"],
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

def _r(name, ptype, tags=None):
    """Register a parameter with optional feature tags."""
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        case_optimization=(name in CASE_OPT_PARAMS),
        constraints=CONSTRAINTS.get(name),
        dependencies=DEPENDENCIES.get(name),
        tags=tags if tags else set(),
    ))


def _load():  # pylint: disable=too-many-locals,too-many-statements
    """Load all parameter definitions."""
    INT, REAL, LOG, STR = ParamType.INT, ParamType.REAL, ParamType.LOG, ParamType.STR
    A_REAL = ParamType.ANALYTIC_REAL

    # ==========================================================================
    # SIMPLE PARAMETERS (non-indexed)
    # ==========================================================================

    # --- Grid ---
    for n in ["m", "n", "p"]:
        _r(n, INT, {"grid"})
    _r("cyl_coord", LOG, {"grid"})
    for n in ["stretch_x", "stretch_y", "stretch_z"]:
        _r(n, LOG, {"grid"})
    for d in ["x", "y", "z"]:
        _r(f"{d}_a", REAL, {"grid"})
        _r(f"{d}_b", REAL, {"grid"})
        _r(f"a_{d}", REAL, {"grid"})
        _r(f"loops_{d}", INT, {"grid"})
        _r(f"{d}_domain%beg", REAL, {"grid"})
        _r(f"{d}_domain%end", REAL, {"grid"})

    # --- Time stepping ---
    for n in ["time_stepper", "t_step_old", "t_step_start", "t_step_stop",
              "t_step_save", "t_step_print", "adap_dt_max_iters"]:
        _r(n, INT, {"time"})
    for n in ["dt", "cfl_target", "cfl_max", "t_tol", "adap_dt_tol", "t_stop", "t_save"]:
        _r(n, REAL, {"time"})
    for n in ["cfl_adap_dt", "cfl_const_dt", "cfl_dt", "adap_dt"]:
        _r(n, LOG, {"time"})

    # --- WENO/reconstruction ---
    _r("weno_order", INT, {"weno"})
    _r("recon_type", INT)
    _r("muscl_order", INT)
    _r("muscl_lim", INT)
    for n in ["weno_eps", "teno_CT", "wenoz_q"]:
        _r(n, REAL, {"weno"})
    for n in ["mapped_weno", "wenoz", "teno", "weno_avg", "mp_weno", "null_weights"]:
        _r(n, LOG, {"weno"})
    _r("weno_Re_flux", LOG, {"weno", "viscosity"})

    # --- Riemann solver ---
    for n in ["riemann_solver", "wave_speeds", "avg_state", "low_Mach"]:
        _r(n, INT, {"riemann"})

    # --- MHD ---
    _r("Bx0", REAL, {"mhd"})
    for n in ["hyper_cleaning_speed", "hyper_cleaning_tau"]:
        _r(n, REAL, {"mhd"})
    for n in ["mhd", "hyper_cleaning", "powell"]:
        _r(n, LOG, {"mhd"})

    # --- Bubbles ---
    for n in ["R0ref", "nb", "Web", "Ca"]:
        _r(n, REAL, {"bubbles"})
    _r("Re_inv", REAL, {"bubbles", "viscosity"})
    _r("bubble_model", INT, {"bubbles"})
    for n in ["polytropic", "bubbles_euler", "polydisperse", "qbmm", "bubbles_lagrange"]:
        _r(n, LOG, {"bubbles"})

    # --- Viscosity ---
    _r("viscous", LOG, {"viscosity"})

    # --- Elasticity ---
    for n in ["hypoelasticity", "hyperelasticity"]:
        _r(n, LOG, {"elasticity"})

    # --- Surface tension ---
    _r("sigma", REAL, {"surface_tension"})
    _r("surface_tension", LOG, {"surface_tension"})

    # --- Chemistry ---
    _r("cantera_file", STR, {"chemistry"})
    _r("chemistry", LOG, {"chemistry"})

    # --- Acoustic ---
    _r("num_source", INT, {"acoustic"})
    _r("acoustic_source", LOG, {"acoustic"})

    # --- Immersed boundary ---
    _r("num_ibs", INT, {"ib"})
    _r("ib", LOG, {"ib"})

    # --- Probes ---
    for n in ["num_probes", "num_integrals"]:
        _r(n, INT, {"probes"})
    _r("probe_wrt", LOG, {"output", "probes"})
    _r("integral_wrt", LOG, {"output", "probes"})

    # --- Output ---
    _r("precision", INT, {"output"})
    _r("format", INT, {"output"})
    _r("schlieren_alpha", REAL, {"output"})
    for n in ["parallel_io", "file_per_process", "run_time_info", "prim_vars_wrt",
              "cons_vars_wrt", "fft_wrt"]:
        _r(n, LOG, {"output"})
    for n in ["schlieren_wrt", "alpha_rho_wrt", "rho_wrt", "mom_wrt", "vel_wrt",
              "flux_wrt", "E_wrt", "pres_wrt", "alpha_wrt", "kappa_wrt", "gamma_wrt",
              "heat_ratio_wrt", "pi_inf_wrt", "pres_inf_wrt", "c_wrt",
              "omega_wrt", "qm_wrt", "liutex_wrt", "cf_wrt", "sim_data", "output_partial_domain"]:
        _r(n, LOG, {"output"})
    for d in ["x", "y", "z"]:
        _r(f"{d}_output%beg", REAL, {"output"})
        _r(f"{d}_output%end", REAL, {"output"})
    # Lagrangian output
    for v in ["lag_header", "lag_txt_wrt", "lag_db_wrt", "lag_id_wrt", "lag_pos_wrt",
              "lag_pos_prev_wrt", "lag_vel_wrt", "lag_rad_wrt", "lag_rvel_wrt",
              "lag_r0_wrt", "lag_rmax_wrt", "lag_rmin_wrt", "lag_dphidt_wrt",
              "lag_pres_wrt", "lag_mv_wrt", "lag_mg_wrt", "lag_betaT_wrt", "lag_betaC_wrt"]:
        _r(v, LOG, {"bubbles", "output"})

    # --- Boundary conditions ---
    for d in ["x", "y", "z"]:
        _r(f"bc_{d}%beg", INT, {"bc"})
        _r(f"bc_{d}%end", INT, {"bc"})

    # --- Relativity ---
    _r("relativity", LOG, {"relativity"})

    # --- Other (no specific feature tag) ---
    for n in ["model_eqns", "num_fluids", "thermal", "relax_model", "igr_order",
              "num_bc_patches", "num_patches", "perturb_flow_fluid", "perturb_sph_fluid",
              "dist_type", "mixlayer_perturb_nk", "elliptic_smoothing_iters",
              "n_start_old", "n_start", "fd_order", "num_igr_iters",
              "num_igr_warm_start_iters", "igr_iter_solver", "nv_uvm_igr_temps_on_gpu",
              "flux_lim"]:
        _r(n, INT)
    for n in ["pref", "poly_sigma", "rhoref", "mixlayer_vel_coef", "mixlayer_domain",
              "mixlayer_perturb_k0", "perturb_flow_mag", "fluid_rho", "sigR", "sigV",
              "rhoRV", "palpha_eps", "ptgalpha_eps", "pi_fac", "tau_star",
              "cont_damage_s", "alpha_bar", "alf_factor", "ic_eps", "ic_beta"]:
        _r(n, REAL)
    for n in ["mpp_lim", "relax", "adv_n", "cont_damage", "igr", "down_sample",
              "old_grid", "old_ic", "mixlayer_vel_profile", "mixlayer_perturb",
              "perturb_flow", "perturb_sph", "pre_stress", "elliptic_smoothing",
              "simplex_perturb", "alt_soundspeed", "mixture_err", "rdma_mpi",
              "igr_pres_lim", "int_comp", "nv_uvm_out_of_core", "nv_uvm_pref_gpu"]:
        _r(n, LOG)
    _r("case_dir", STR)

    # Body force
    for d in ["x", "y", "z"]:
        for v in ["k", "w", "p", "g"]:
            _r(f"{v}_{d}", REAL)
        _r(f"bf_{d}", LOG)

    # ==========================================================================
    # INDEXED PARAMETERS
    # ==========================================================================

    # --- patch_icpp (10 patches) ---
    for i in range(1, NP + 1):
        px = f"patch_icpp({i})%"
        for a in ["geometry", "smooth_patch_id", "hcid", "model_spc"]:
            _r(f"{px}{a}", INT)
        for a in ["smoothen", "alter_patch"] if i >= 2 else ["smoothen"]:
            _r(f"{px}{a}", LOG)
        for a in ["radius", "radii", "epsilon", "beta", "normal", "alpha_rho",
                  "non_axis_sym", "smooth_coeff", "rho", "vel", "alpha", "gamma",
                  "pi_inf", "cv", "qv", "qvp", "model_threshold"]:
            _r(f"{px}{a}", REAL)
        # Bubble fields
        for a in ["r0", "v0", "p0", "m0"]:
            _r(f"{px}{a}", REAL, {"bubbles"})
        for j in range(2, 10):
            _r(f"{px}a({j})", REAL)
        _r(f"{px}pres", A_REAL)
        _r(f"{px}cf_val", A_REAL)
        # MHD fields
        for a in ["Bx", "By", "Bz"]:
            _r(f"{px}{a}", A_REAL, {"mhd"})
        # Chemistry species
        for j in range(1, 101):
            _r(f"{px}Y({j})", A_REAL, {"chemistry"})
        _r(f"{px}model_filepath", STR)
        for t in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _r(f"{px}model_{t}({j})", REAL)
        for d in ["x", "y", "z"]:
            _r(f"{px}{d}_centroid", REAL)
            _r(f"{px}length_{d}", REAL)
        for j in range(1, 4):
            _r(f"{px}radii({j})", REAL)
            _r(f"{px}normal({j})", REAL)
            _r(f"{px}vel({j})", A_REAL)
        for f in range(1, NF + 1):
            _r(f"{px}alpha({f})", A_REAL)
            _r(f"{px}alpha_rho({f})", A_REAL)
        # Elasticity stress tensor
        for j in range(1, 7):
            _r(f"{px}tau_e({j})", A_REAL, {"elasticity"})
        if i >= 2:
            for j in range(1, i):
                _r(f"{px}alter_patch({j})", LOG)

    # --- fluid_pp (10 fluids) ---
    for f in range(1, NF + 1):
        px = f"fluid_pp({f})%"
        for a in ["gamma", "pi_inf", "cv", "qv", "qvp"]:
            _r(f"{px}{a}", REAL)
        _r(f"{px}mul0", REAL, {"viscosity"})
        _r(f"{px}ss", REAL, {"surface_tension"})
        for a in ["pv", "gamma_v", "M_v", "mu_v", "k_v", "cp_v", "D_v"]:
            _r(f"{px}{a}", REAL, {"bubbles"})
        _r(f"{px}G", REAL, {"elasticity"})
        for j in [1, 2]:
            _r(f"{px}Re({j})", REAL, {"viscosity"})

    # --- bub_pp (bubble properties) ---
    for a in ["R0ref", "p0ref", "rho0ref", "T0ref", "ss", "pv", "vd",
              "mu_l", "mu_v", "mu_g", "gam_v", "gam_g",
              "M_v", "M_g", "k_v", "k_g", "cp_v", "cp_g", "R_v", "R_g"]:
        _r(f"bub_pp%{a}", REAL, {"bubbles"})

    # --- patch_ib (10 immersed boundaries) ---
    for i in range(1, NI + 1):
        px = f"patch_ib({i})%"
        for a in ["geometry", "moving_ibm"]:
            _r(f"{px}{a}", INT, {"ib"})
        for a, pt in [("radius", REAL), ("theta", REAL), ("slip", LOG), ("c", REAL),
                      ("p", REAL), ("t", REAL), ("m", REAL), ("mass", REAL)]:
            _r(f"{px}{a}", pt, {"ib"})
        for j in range(1, 4):
            _r(f"{px}angles({j})", REAL, {"ib"})
        for d in ["x", "y", "z"]:
            _r(f"{px}{d}_centroid", REAL, {"ib"})
            _r(f"{px}length_{d}", REAL, {"ib"})
        for a, pt in [("model_filepath", STR), ("model_spc", INT), ("model_threshold", REAL)]:
            _r(f"{px}{a}", pt, {"ib"})
        for t in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _r(f"{px}model_{t}({j})", REAL, {"ib"})
        for j in range(1, 4):
            _r(f"{px}vel({j})", REAL, {"ib"})
            _r(f"{px}angular_vel({j})", REAL, {"ib"})

    # --- acoustic sources (4 sources) ---
    for i in range(1, NA + 1):
        px = f"acoustic({i})%"
        for a in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            _r(f"{px}{a}", INT, {"acoustic"})
        _r(f"{px}dipole", LOG, {"acoustic"})
        for a in ["mag", "length", "height", "wavelength", "frequency",
                  "gauss_sigma_dist", "gauss_sigma_time", "npulse",
                  "dir", "delay", "foc_length", "aperture",
                  "element_spacing_angle", "element_polygon_ratio",
                  "rotate_angle", "bb_bandwidth", "bb_lowest_freq"]:
            _r(f"{px}{a}", REAL, {"acoustic"})
        for j in range(1, 4):
            _r(f"{px}loc({j})", REAL, {"acoustic"})

    # --- probes (10 probes) ---
    for i in range(1, NPR + 1):
        for d in ["x", "y", "z"]:
            _r(f"probe({i})%{d}", REAL, {"probes"})

    # --- integrals (5 integral regions) ---
    for i in range(1, 6):
        for d in ["x", "y", "z"]:
            _r(f"integral({i})%{d}min", REAL, {"probes"})
            _r(f"integral({i})%{d}max", REAL, {"probes"})

    # --- Extended BC ---
    for d in ["x", "y", "z"]:
        px = f"bc_{d}%"
        for a in ["vb1", "vb2", "vb3", "ve1", "ve2", "ve3", "pres_in", "pres_out"]:
            _r(f"{px}{a}", REAL, {"bc"})
        for a in ["grcbc_in", "grcbc_out", "grcbc_vel_out"]:
            _r(f"{px}{a}", LOG, {"bc"})
        for f in range(1, NF + 1):
            _r(f"{px}alpha_rho_in({f})", REAL, {"bc"})
            _r(f"{px}alpha_in({f})", REAL, {"bc"})
        for j in range(1, 4):
            _r(f"{px}vel_in({j})", REAL, {"bc"})
            _r(f"{px}vel_out({j})", REAL, {"bc"})

    # --- patch_bc (10 BC patches) ---
    for i in range(1, NB + 1):
        px = f"patch_bc({i})%"
        for a in ["geometry", "type", "dir", "loc"]:
            _r(f"{px}{a}", INT, {"bc"})
        for j in range(1, 4):
            _r(f"{px}centroid({j})", REAL, {"bc"})
            _r(f"{px}length({j})", REAL, {"bc"})
        _r(f"{px}radius", REAL, {"bc"})

    # --- simplex_params ---
    for f in range(1, NF + 1):
        _r(f"simplex_params%perturb_dens({f})", LOG)
        _r(f"simplex_params%perturb_dens_freq({f})", REAL)
        _r(f"simplex_params%perturb_dens_scale({f})", REAL)
        for j in range(1, 4):
            _r(f"simplex_params%perturb_dens_offset({f}, {j})", REAL)
    for d in range(1, 4):
        _r(f"simplex_params%perturb_vel({d})", LOG)
        _r(f"simplex_params%perturb_vel_freq({d})", REAL)
        _r(f"simplex_params%perturb_vel_scale({d})", REAL)
        for j in range(1, 4):
            _r(f"simplex_params%perturb_vel_offset({d},{j})", REAL)

    # --- lag_params (Lagrangian bubbles) ---
    for a in ["heatTransfer_model", "massTransfer_model", "pressure_corrector",
              "write_bubbles", "write_bubbles_stats"]:
        _r(f"lag_params%{a}", LOG, {"bubbles"})
    for a in ["solver_approach", "cluster_type", "smooth_type", "nBubs_glb"]:
        _r(f"lag_params%{a}", INT, {"bubbles"})
    for a in ["epsilonb", "valmaxvoid", "charwidth", "c0", "rho0", "T0", "x0", "Thost"]:
        _r(f"lag_params%{a}", REAL, {"bubbles"})

    # --- chem_params ---
    for a in ["diffusion", "reactions"]:
        _r(f"chem_params%{a}", LOG, {"chemistry"})
    for a in ["gamma_method", "transport_model"]:
        _r(f"chem_params%{a}", INT, {"chemistry"})

    # --- Per-fluid output arrays ---
    for f in range(1, NF + 1):
        _r(f"schlieren_alpha({f})", REAL, {"output"})
        for a in ["alpha_rho_wrt", "alpha_wrt", "kappa_wrt", "alpha_rho_e_wrt"]:
            _r(f"{a}({f})", LOG, {"output"})
    for j in range(1, 4):
        for a in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
            _r(f"{a}({j})", LOG, {"output"})

    # --- chem_wrt (chemistry output) ---
    for j in range(1, 101):
        _r(f"chem_wrt_Y({j})", LOG, {"chemistry", "output"})
    _r("chem_wrt_T", LOG, {"chemistry", "output"})

    # --- fluid_rho ---
    for f in range(1, NF + 1):
        _r(f"fluid_rho({f})", REAL)


# Load definitions when module imported and freeze registry
def _init_registry():
    """Initialize and freeze the registry. Called once at module import."""
    try:
        # Validate constraint and dependency schemas first
        # This catches typos like "choises" instead of "choices"
        _validate_all_constraints(CONSTRAINTS)
        _validate_all_dependencies(DEPENDENCIES)

        # Load all parameter definitions
        _load()

        # Freeze registry to prevent further modifications
        REGISTRY.freeze()
    except Exception as e:
        # Re-raise with context to help debugging initialization failures
        raise RuntimeError(
            f"Failed to initialize parameter registry: {e}\n"
            "This is likely a bug in the parameter definitions."
        ) from e

_init_registry()
