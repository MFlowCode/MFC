"""
MFC Parameter Definitions (Compact).

Single file containing all ~3,300 parameter definitions using loops.
This replaces the definitions/ directory.
"""

import re
from typing import Any, Dict

from .namelist_parser import get_fortran_constants
from .registry import REGISTRY, IndexedFamily
from .schema import ParamDef, ParamType

# Index limits — sourced from Fortran compile-time constants (m_constants.fpp).
# Falls back to the inline default when src/ is unavailable (e.g. Homebrew).
# Default must match src/common/m_constants.fpp — enforced by co-location.
_FC = get_fortran_constants()


def _fc(name: str, default: int) -> int:
    """Get a Fortran constant, using the inline default when m_constants.fpp is unavailable."""
    if _FC:
        if name not in _FC:
            raise RuntimeError(f"Fortran constant '{name}' not found in m_constants.fpp. Toolchain is out of sync with Fortran source.")
        return _FC[name]
    return default


NF = _fc("num_fluids_max", 10)  # fluid_pp
NPR = _fc("num_probes_max", 10)  # probe, acoustic, integral
NB = _fc("num_bc_patches_max", 10)  # patch_bc
NUM_PATCHES_MAX = _fc("num_patches_max", 10)  # patch_icpp (Fortran array bound)
NIB = _fc("num_ib_patches_max", 50000)  # patch_ib (Fortran array bound)
# Enumeration limits for families not yet converted to IndexedFamily.
# These are smaller than the Fortran array bounds to keep the registry compact.
# The CONSTRAINTS dict below uses the Fortran constants for validation.
NP = 10  # patch_icpp: has per-index variations, can't easily be IndexedFamily
NA = 4  # acoustic sources: enumerated individually


# Parameters that can be hard-coded for GPU case optimization
CASE_OPT_PARAMS = {
    "mapped_weno",
    "wenoz",
    "teno",
    "wenoz_q",
    "nb",
    "weno_order",
    "num_fluids",
    "mhd",
    "relativity",
    "igr_order",
    "viscous",
    "igr_iter_solver",
    "igr",
    "igr_pres_lim",
    "recon_type",
    "muscl_order",
    "muscl_lim",
}


# Data-driven Annotations for Doc Generation
# These dicts are the single source of truth for parameter hints in the docs.
# To annotate a new param, add an entry here instead of editing docs_gen.py.

HINTS = {
    "bc": {
        "grcbc_in": "Enables GRCBC subsonic inflow (bc type -7)",
        "grcbc_out": "Enables GRCBC subsonic outflow (bc type -8)",
        "grcbc_vel_out": "GRCBC velocity outlet (requires `grcbc_out`)",
        "vel_in": "Inlet velocity component (used with `grcbc_in`)",
        "vel_out": "Outlet velocity component (used with `grcbc_vel_out`)",
        "pres_in": "Inlet pressure (used with `grcbc_in`)",
        "pres_out": "Outlet pressure (used with `grcbc_out`)",
        "alpha_rho_in": "Inlet partial density per fluid (used with `grcbc_in`)",
        "alpha_in": "Inlet volume fraction per fluid (used with `grcbc_in`)",
        "vb1": "Boundary velocity component 1 at domain begin",
        "vb2": "Boundary velocity component 2 at domain begin",
        "vb3": "Boundary velocity component 3 at domain begin",
        "ve1": "Boundary velocity component 1 at domain end",
        "ve2": "Boundary velocity component 2 at domain end",
        "ve3": "Boundary velocity component 3 at domain end",
        "Twall_in": "Temperature of the entrance-side isothermal wall.",
        "Twall_out": "Temperature of the exit-side isothermal wall.",
    },
    "patch_bc": {
        "geometry": "Patch shape: 1=line, 2=circle, 3=rectangle",
        "type": "BC type applied within patch region",
        "dir": "Patch normal direction (1=x, 2=y, 3=z)",
        "loc": "Domain boundary (-1=begin, 1=end)",
        "centroid": "Patch center coordinate",
        "length": "Patch dimension",
        "radius": "Patch radius (geometry=2)",
    },
    "simplex_params": {
        "perturb_dens": "Enable simplex density perturbation",
        "perturb_dens_freq": "Density perturbation frequency",
        "perturb_dens_scale": "Density perturbation amplitude",
        "perturb_dens_offset": "Density perturbation offset seed",
        "perturb_vel": "Enable simplex velocity perturbation",
        "perturb_vel_freq": "Velocity perturbation frequency",
        "perturb_vel_scale": "Velocity perturbation amplitude",
        "perturb_vel_offset": "Velocity perturbation offset seed",
    },
    "fluid_pp": {
        "gamma": "Specific heat ratio (EOS)",
        "pi_inf": "Stiffness pressure (EOS)",
        "cv": "Specific heat at constant volume",
        "qv": "Heat of formation",
        "qvp": "Heat of formation derivative",
    },
}

# Tag → display name for docs. Dict order = priority when a param has multiple tags.
TAG_DISPLAY_NAMES = {
    "bubbles": "Bubble model",
    "mhd": "MHD",
    "chemistry": "Chemistry",
    "time": "Time-stepping",
    "grid": "Grid",
    "weno": "WENO",
    "viscosity": "Viscosity",
    "elasticity": "Elasticity",
    "surface_tension": "Surface tension",
    "acoustic": "Acoustic",
    "ib": "Immersed boundary",
    "probes": "Probe/integral",
    "riemann": "Riemann solver",
    "relativity": "Relativity",
    "output": "Output",
    "bc": "Boundary condition",
}

# Prefix → hint for untagged simple params
PREFIX_HINTS = {
    "mixlayer_": "Mixing layer parameter",
    "nv_uvm_": "GPU memory management",
    "ic_": "Initial condition parameter",
}


def _lookup_hint(name):
    """Auto-derive constraint hint from HINTS dict using family+attribute matching."""
    if "%" not in name:
        # Check PREFIX_HINTS for simple params
        for prefix, label in PREFIX_HINTS.items():
            if name.startswith(prefix):
                return label
        return ""
    # Compound name: extract family and attribute
    prefix, attr_full = name.split("%", 1)
    # Normalize family: "bc_x" → "bc", "patch_bc(1)" → "patch"
    family = re.sub(r"[_(].*", "", prefix)
    if family not in HINTS:
        # Fallback: keep underscores — "patch_bc" → "patch_bc", "simplex_params" → "simplex_params"
        m = re.match(r"^[a-zA-Z_]+", prefix)
        family = m.group(0) if m else ""
    if family not in HINTS:
        return ""
    # Strip index from attr: "vel_in(1)" → "vel_in"
    m = re.match(r"^[a-zA-Z_0-9]+", attr_full)
    if not m:
        return ""
    attr = m.group(0)
    return HINTS[family].get(attr, "")


# Schema Validation for Constraints and Dependencies
# Uses rapidfuzz for "did you mean?" suggestions when typos are detected

_VALID_CONSTRAINT_KEYS = {"choices", "min", "max", "value_labels"}
_VALID_DEPENDENCY_KEYS = {"when_true", "when_set", "when_value"}
_VALID_CONDITION_KEYS = {"requires", "recommends", "requires_value"}


def _validate_constraint(param_name: str, constraint: Dict[str, Any]) -> None:
    """Validate a constraint dict has valid keys with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error

    invalid_keys = set(constraint.keys()) - _VALID_CONSTRAINT_KEYS
    if invalid_keys:
        # Get suggestion for the first invalid key
        first_invalid = next(iter(invalid_keys))
        raise ValueError(invalid_key_error(f"constraint for '{param_name}'", first_invalid, _VALID_CONSTRAINT_KEYS))

    # Validate types
    if "choices" in constraint and not isinstance(constraint["choices"], list):
        raise ValueError(f"Constraint 'choices' for '{param_name}' must be a list")
    if "min" in constraint and not isinstance(constraint["min"], (int, float)):
        raise ValueError(f"Constraint 'min' for '{param_name}' must be a number")
    if "max" in constraint and not isinstance(constraint["max"], (int, float)):
        raise ValueError(f"Constraint 'max' for '{param_name}' must be a number")
    if "value_labels" in constraint:
        if not isinstance(constraint["value_labels"], dict):
            raise ValueError(f"Constraint 'value_labels' for '{param_name}' must be a dict")
        if "choices" in constraint:
            for key in constraint["value_labels"]:
                if key not in constraint["choices"]:
                    raise ValueError(f"value_labels key {key!r} for '{param_name}' not in choices {constraint['choices']}")


def _validate_dependency(param_name: str, dependency: Dict[str, Any]) -> None:
    """Validate a dependency dict has valid structure with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error

    invalid_keys = set(dependency.keys()) - _VALID_DEPENDENCY_KEYS
    if invalid_keys:
        first_invalid = next(iter(invalid_keys))
        raise ValueError(invalid_key_error(f"dependency for '{param_name}'", first_invalid, _VALID_DEPENDENCY_KEYS))

    def _validate_condition(cond_label: str, condition: Any) -> None:
        """Validate a condition dict (shared by when_true, when_set, when_value entries)."""
        if not isinstance(condition, dict):
            raise ValueError(f"Dependency '{cond_label}' for '{param_name}' must be a dict")
        invalid_cond_keys = set(condition.keys()) - _VALID_CONDITION_KEYS
        if invalid_cond_keys:
            first_invalid = next(iter(invalid_cond_keys))
            raise ValueError(invalid_key_error(f"condition in '{cond_label}' for '{param_name}'", first_invalid, _VALID_CONDITION_KEYS))
        for req_key in ["requires", "recommends"]:
            if req_key in condition and not isinstance(condition[req_key], list):
                raise ValueError(f"Dependency '{cond_label}/{req_key}' for '{param_name}' must be a list")
        if "requires_value" in condition:
            rv = condition["requires_value"]
            if not isinstance(rv, dict):
                raise ValueError(f"Dependency '{cond_label}/requires_value' for '{param_name}' must be a dict")
            for rv_param, rv_vals in rv.items():
                if not isinstance(rv_vals, list):
                    raise ValueError(f"Dependency '{cond_label}/requires_value/{rv_param}' for '{param_name}' must be a list")

    for condition_key in ["when_true", "when_set"]:
        if condition_key in dependency:
            _validate_condition(condition_key, dependency[condition_key])

    if "when_value" in dependency:
        wv = dependency["when_value"]
        if not isinstance(wv, dict):
            raise ValueError(f"Dependency 'when_value' for '{param_name}' must be a dict")
        for val, condition in wv.items():
            _validate_condition(f"when_value/{val}", condition)


def _validate_all_constraints(constraints: Dict[str, Dict]) -> None:
    """Validate all constraint definitions."""
    for param_name, constraint in constraints.items():
        _validate_constraint(param_name, constraint)


def _validate_all_dependencies(dependencies: Dict[str, Dict]) -> None:
    """Validate all dependency definitions."""
    for param_name, dependency in dependencies.items():
        _validate_dependency(param_name, dependency)


def get_value_label(param_name: str, value: int) -> str:
    """Look up the human-readable label for a parameter's integer code.

    Returns the label string, or ``str(value)`` when no label is defined.
    This is the single source of truth for value ↔ label mappings.
    """
    constraint = CONSTRAINTS.get(param_name)
    if constraint is None:
        return str(value)
    labels = constraint.get("value_labels")
    if labels is None:
        return str(value)
    return labels.get(value, str(value))


# Parameter constraints (choices, min, max)
CONSTRAINTS = {
    # Reconstruction
    "weno_order": {
        "choices": [0, 1, 3, 5, 7],
        "value_labels": {0: "MUSCL mode", 1: "1st order", 3: "WENO3", 5: "WENO5", 7: "WENO7"},
    },
    "recon_type": {
        "choices": [1, 2],
        "value_labels": {1: "WENO", 2: "MUSCL"},
    },
    "muscl_order": {
        "choices": [1, 2],
        "value_labels": {1: "1st order", 2: "2nd order"},
    },
    "muscl_lim": {
        "choices": [0, 1, 2, 3, 4, 5],
        "value_labels": {0: "unlimited", 1: "minmod", 2: "MC", 3: "Van Albada", 4: "Van Leer", 5: "SUPERBEE"},
    },
    "int_comp": {
        "choices": [0, 1, 2],
        "value_labels": {0: "off", 1: "THINC", 2: "MTHINC"},
    },
    # Time stepping
    "time_stepper": {
        "choices": [1, 2, 3],
        "value_labels": {1: "RK1 (Forward Euler)", 2: "RK2", 3: "RK3 (SSP)"},
    },
    # Riemann solver
    "riemann_solver": {
        "choices": [1, 2, 4, 5],
        "value_labels": {1: "HLL", 2: "HLLC", 4: "HLLD", 5: "Lax-Friedrichs"},
    },
    "wave_speeds": {
        "choices": [1, 2],
        "value_labels": {1: "direct", 2: "pressure"},
    },
    "avg_state": {
        "choices": [1, 2],
        "value_labels": {1: "Roe", 2: "arithmetic"},
    },
    # Model equations
    "model_eqns": {
        "choices": [1, 2, 3, 4],
        "value_labels": {1: "Gamma-law", 2: "5-Equation", 3: "6-Equation", 4: "4-Equation"},
    },
    # Bubbles
    "bubble_model": {
        "choices": [1, 2, 3],
        "value_labels": {1: "Gilmore", 2: "Keller-Miksis", 3: "Rayleigh-Plesset"},
    },
    # Output
    "format": {
        "choices": [1, 2],
        "value_labels": {1: "Silo", 2: "binary"},
    },
    "precision": {
        "choices": [1, 2],
        "value_labels": {1: "single", 2: "double"},
    },
    # Time stepping (must be positive)
    "dt": {"min": 0},
    "t_stop": {"min": 0},
    "t_save": {"min": 0},
    "t_step_save": {"min": 1},
    "t_step_print": {"min": 1},
    "cfl_target": {"min": 0},
    # WENO
    "weno_eps": {"min": 0},
    # MUSCL
    "muscl_eps": {"min": 0},
    # Physics (must be non-negative)
    "R0ref": {"min": 0},
    "sigma": {"min": 0},
    # Counts (must be positive)
    "num_fluids": {"min": 1, "max": NF},
    "num_patches": {"min": 0, "max": NUM_PATCHES_MAX},
    "num_ibs": {"min": 0},
    "num_source": {"min": 1},
    "num_probes": {"min": 1},
    "num_integrals": {"min": 1},
    "nb": {"min": 1},
    "m": {"min": 0},
    "n": {"min": 0},
    "p": {"min": 0},
}

# Parameter dependencies (requires, recommends)
DEPENDENCIES = {
    "bubbles_euler": {
        "when_true": {
            "recommends": ["nb", "polytropic"],
            "requires_value": {
                "model_eqns": [2, 4],
                "riemann_solver": [2],
                "avg_state": [2],
            },
        }
    },
    "model_eqns": {
        "when_value": {
            2: {"requires": ["num_fluids"]},
            3: {"requires_value": {"riemann_solver": [2], "avg_state": [2], "wave_speeds": [1]}},
            4: {"requires": ["rhoref", "pref"], "requires_value": {"num_fluids": [1]}},
        }
    },
    "viscous": {
        "when_true": {
            "recommends": ["fluid_pp(1)%Re(1)"],
        }
    },
    "polydisperse": {
        "when_true": {
            "requires": ["nb", "poly_sigma"],
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
    "collision_model": {
        "when_set": {
            "requires": ["ib", "coefficient_of_restitution", "collision_time"],
        }
    },
    "acoustic_source": {
        "when_true": {
            "requires": ["num_source"],
        }
    },
    "probe_wrt": {
        "when_true": {
            "requires": ["num_probes", "fd_order"],
        }
    },
    "stretch_x": {
        "when_true": {
            "requires": ["a_x", "x_a", "x_b"],
        }
    },
    "stretch_y": {
        "when_true": {
            "requires": ["a_y", "y_a", "y_b"],
        }
    },
    "stretch_z": {
        "when_true": {
            "requires": ["a_z", "z_a", "z_b"],
        }
    },
    "bf_x": {
        "when_true": {
            "requires": ["k_x", "w_x", "p_x", "g_x"],
        }
    },
    "bf_y": {
        "when_true": {
            "requires": ["k_y", "w_y", "p_y", "g_y"],
        }
    },
    "bf_z": {
        "when_true": {
            "requires": ["k_z", "w_z", "p_z", "g_z"],
        }
    },
    "teno": {
        "when_true": {
            "requires": ["teno_CT"],
        }
    },
    "recon_type": {
        "when_value": {
            2: {"recommends": ["muscl_order", "muscl_lim"]},
        }
    },
    "surface_tension": {
        "when_true": {
            "requires": ["sigma"],
        }
    },
    "relativity": {
        "when_true": {
            "requires": ["mhd"],
        }
    },
    "schlieren_wrt": {
        "when_true": {
            "requires": ["fd_order"],
        }
    },
    "cfl_adap_dt": {
        "when_true": {
            "recommends": ["cfl_target"],
        }
    },
    "cfl_dt": {
        "when_true": {
            "recommends": ["cfl_target"],
        }
    },
    "integral_wrt": {
        "when_true": {
            "requires": ["fd_order"],
        }
    },
}


def _r(name, ptype, tags=None, desc=None, hint=None, math=None, str_len=None, storage_precision=False):
    """Register a parameter with optional feature tags and description."""
    if hint is None:
        hint = _lookup_hint(name)
    if desc is None:
        from .descriptions import get_description

        desc = get_description(name)
    constraint = CONSTRAINTS.get(name)
    if constraint and "value_labels" in constraint:
        labels = constraint["value_labels"]
        suffix = ", ".join(f"{v}={labels[v]}" for v in sorted(labels))
        desc = f"{desc} ({suffix})".strip()
    REGISTRY.register(
        ParamDef(
            name=name,
            param_type=ptype,
            description=desc,
            case_optimization=(name in CASE_OPT_PARAMS),
            constraints=constraint,
            dependencies=DEPENDENCIES.get(name),
            tags=tags if tags else set(),
            hint=hint,
            math_symbol=math or "",
            str_len=str_len if str_len is not None else "name_len",
            storage_precision=storage_precision,
        )
    )


def _load():
    """Load all parameter definitions."""
    INT, REAL, LOG, STR = ParamType.INT, ParamType.REAL, ParamType.LOG, ParamType.STR
    A_REAL = ParamType.ANALYTIC_REAL

    # SIMPLE PARAMETERS (non-indexed)

    # Grid
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

    # Time stepping
    for n in ["time_stepper", "t_step_old", "t_step_start", "t_step_stop", "t_step_save", "t_step_print", "adap_dt_max_iters"]:
        _r(n, INT, {"time"})
    _r("dt", REAL, {"time"}, math=r"\f$\Delta t\f$")
    _r("cfl_target", REAL, {"time"}, math=r"\f$\mathrm{CFL}\f$")
    for n in ["adap_dt_tol", "t_stop", "t_save"]:
        _r(n, REAL, {"time"})
    for n in ["cfl_adap_dt", "cfl_const_dt", "cfl_dt", "adap_dt"]:
        _r(n, LOG, {"time"})

    # WENO/reconstruction
    _r("weno_order", INT, {"weno"})
    _r("recon_type", INT)
    _r("muscl_order", INT)
    _r("muscl_lim", INT)
    _r("muscl_eps", REAL)
    _r("weno_eps", REAL, {"weno"}, math=r"\f$\varepsilon\f$")
    _r("teno_CT", REAL, {"weno"}, math=r"\f$C_T\f$")
    _r("wenoz_q", REAL, {"weno"})
    for n in ["mapped_weno", "wenoz", "teno", "weno_avg", "mp_weno", "null_weights"]:
        _r(n, LOG, {"weno"})
    _r("weno_Re_flux", LOG, {"weno", "viscosity"})

    # Riemann solver
    for n in ["riemann_solver", "wave_speeds", "avg_state", "low_Mach"]:
        _r(n, INT, {"riemann"})

    # MHD
    _r("Bx0", REAL, {"mhd"}, math=r"\f$B_{x,0}\f$")
    _r("hyper_cleaning_speed", REAL, {"mhd"}, math=r"\f$c_h\f$")
    _r("hyper_cleaning_tau", REAL, {"mhd"})
    for n in ["mhd", "hyper_cleaning"]:
        _r(n, LOG, {"mhd"})

    # Bubbles
    _r("R0ref", REAL, {"bubbles"}, math=r"\f$R_0\f$")
    _r("nb", INT, {"bubbles"}, math=r"\f$N_b\f$")
    _r("Web", REAL, {"bubbles"}, math=r"\f$\mathrm{We}\f$")
    _r("Ca", REAL, {"bubbles"}, math=r"\f$\mathrm{Ca}\f$")
    _r("Re_inv", REAL, {"bubbles", "viscosity"}, math=r"\f$\mathrm{Re}^{-1}\f$")
    _r("bubble_model", INT, {"bubbles"})
    for n in ["polytropic", "bubbles_euler", "polydisperse", "qbmm", "bubbles_lagrange"]:
        _r(n, LOG, {"bubbles"})

    # Viscosity
    _r("viscous", LOG, {"viscosity"})

    # Elasticity
    for n in ["hypoelasticity", "hyperelasticity"]:
        _r(n, LOG, {"elasticity"})

    # Surface tension
    _r("sigma", REAL, {"surface_tension"}, math=r"\f$\sigma\f$")
    _r("surface_tension", LOG, {"surface_tension"})

    # Chemistry
    _r("cantera_file", STR, {"chemistry"})
    _r("chemistry", LOG, {"chemistry"})

    # Acoustic
    _r("num_source", INT, {"acoustic"})
    _r("acoustic_source", LOG, {"acoustic"})

    # Immersed boundary
    _r("num_ibs", INT, {"ib"})
    _r("ib", LOG, {"ib"})
    _r("collision_model", INT, {"ib"})
    _r("coefficient_of_restitution", REAL, {"ib"})
    _r("collision_time", REAL, {"ib"})
    _r("ib_coefficient_of_friction", REAL, {"ib"})

    # Probes
    for n in ["num_probes", "num_integrals"]:
        _r(n, INT, {"probes"})
    _r("probe_wrt", LOG, {"output", "probes"})
    _r("integral_wrt", LOG, {"output", "probes"})

    # Output
    _r("precision", INT, {"output"})
    _r("format", INT, {"output"})
    for n in ["parallel_io", "file_per_process", "run_time_info", "prim_vars_wrt", "cons_vars_wrt", "fft_wrt", "ib_state_wrt"]:
        _r(n, LOG, {"output"})
    for n in [
        "schlieren_wrt",
        "alpha_wrt",
        "rho_wrt",
        "E_wrt",
        "pres_wrt",
        "gamma_wrt",
        "heat_ratio_wrt",
        "pi_inf_wrt",
        "pres_inf_wrt",
        "c_wrt",
        "qm_wrt",
        "liutex_wrt",
        "cf_wrt",
        "sim_data",
        "output_partial_domain",
    ]:
        _r(n, LOG, {"output"})
    for d in ["x", "y", "z"]:
        _r(f"{d}_output%beg", REAL, {"output"})
        _r(f"{d}_output%end", REAL, {"output"})
    # Lagrangian output
    for v in [
        "lag_header",
        "lag_txt_wrt",
        "lag_db_wrt",
        "lag_id_wrt",
        "lag_pos_wrt",
        "lag_pos_prev_wrt",
        "lag_vel_wrt",
        "lag_rad_wrt",
        "lag_rvel_wrt",
        "lag_r0_wrt",
        "lag_rmax_wrt",
        "lag_rmin_wrt",
        "lag_dphidt_wrt",
        "lag_pres_wrt",
        "lag_mv_wrt",
        "lag_mg_wrt",
        "lag_betaT_wrt",
        "lag_betaC_wrt",
    ]:
        _r(v, LOG, {"bubbles", "output"})

    # Boundary conditions
    for d in ["x", "y", "z"]:
        _r(f"bc_{d}%beg", INT, {"bc"})
        _r(f"bc_{d}%end", INT, {"bc"})

    # Relativity
    _r("relativity", LOG, {"relativity"})

    # Other (no specific feature tag)
    for n in [
        "model_eqns",
        "num_fluids",
        "thermal",
        "relax_model",
        "igr_order",
        "num_bc_patches",
        "num_patches",
        "perturb_flow_fluid",
        "perturb_sph_fluid",
        "dist_type",
        "mixlayer_perturb_nk",
        "elliptic_smoothing_iters",
        "n_start_old",
        "n_start",
        "fd_order",
        "num_igr_iters",
        "num_igr_warm_start_iters",
        "igr_iter_solver",
        "nv_uvm_igr_temps_on_gpu",
        "flux_lim",
    ]:
        _r(n, INT)
    _r("pref", REAL, math=r"\f$p_\text{ref}\f$")
    _r("poly_sigma", REAL, math=r"\f$\sigma_\text{poly}\f$")
    _r("rhoref", REAL, math=r"\f$\rho_\text{ref}\f$")
    _r("palpha_eps", REAL, math=r"\f$\varepsilon_\alpha\f$")
    _r("ptgalpha_eps", REAL, math=r"\f$\varepsilon_\alpha\f$")
    _r("pi_fac", REAL, math=r"\f$\pi\text{-factor}\f$")
    for n in [
        "mixlayer_vel_coef",
        "mixlayer_perturb_k0",
        "perturb_flow_mag",
        "sigR",
        "sigV",
        "rhoRV",
        "tau_star",
        "cont_damage_s",
        "alpha_bar",
        "alf_factor",
        "ic_eps",
        "ic_beta",
    ]:
        _r(n, REAL)
    for n in [
        "mpp_lim",
        "relax",
        "adv_n",
        "cont_damage",
        "igr",
        "down_sample",
        "old_grid",
        "old_ic",
        "mixlayer_vel_profile",
        "mixlayer_perturb",
        "perturb_flow",
        "perturb_sph",
        "pre_stress",
        "elliptic_smoothing",
        "simplex_perturb",
        "alt_soundspeed",
        "mixture_err",
        "rdma_mpi",
        "igr_pres_lim",
        "nv_uvm_out_of_core",
        "nv_uvm_pref_gpu",
    ]:
        _r(n, LOG)
    _r("int_comp", INT)
    _r("case_dir", STR, str_len="path_len")

    # Body force
    for d in ["x", "y", "z"]:
        _r(f"g_{d}", REAL, math=r"\f$g_" + d + r"\f$")
        _r(f"k_{d}", REAL, math=r"\f$k_" + d + r"\f$")
        _r(f"w_{d}", REAL, math=r"\f$\omega_" + d + r"\f$")
        _r(f"p_{d}", REAL, math=r"\f$\phi_" + d + r"\f$")
        _r(f"bf_{d}", LOG)

    # INDEXED PARAMETERS

    # patch_icpp (10 patches)
    for i in range(1, NP + 1):
        px = f"patch_icpp({i})%"
        for a in ["geometry", "smooth_patch_id", "hcid", "model_spc"]:
            _r(f"{px}{a}", INT)
        for a in ["smoothen", "alter_patch"] if i >= 2 else ["smoothen"]:
            _r(f"{px}{a}", LOG)
        for a, sym in [("rho", r"\f$\rho\f$"), ("gamma", r"\f$\gamma\f$"), ("pi_inf", r"\f$\pi_\infty\f$"), ("cv", r"\f$c_v\f$"), ("qv", r"\f$q_v\f$"), ("qvp", r"\f$q'_v\f$")]:
            _r(f"{px}{a}", REAL, math=sym)
        for a in ["radius", "radii", "epsilon", "beta", "normal", "alpha_rho", "non_axis_sym", "smooth_coeff", "vel", "alpha", "model_threshold"]:
            _r(f"{px}{a}", REAL)
        # Bubble fields
        for a in ["r0", "v0", "p0", "m0"]:
            _r(f"{px}{a}", REAL, {"bubbles"})
        for j in range(2, 10):
            _r(f"{px}a({j})", REAL)
        _r(f"{px}pres", A_REAL, math=r"\f$p\f$")
        _r(f"{px}cf_val", A_REAL)
        # MHD fields
        for a, sym in [("Bx", r"\f$B_x\f$"), ("By", r"\f$B_y\f$"), ("Bz", r"\f$B_z\f$")]:
            _r(f"{px}{a}", A_REAL, {"mhd"}, math=sym)
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
            _r(f"{px}vel({j})", A_REAL, math=r"\f$u_" + str(j) + r"\f$")
        for f in range(1, NF + 1):
            _r(f"{px}alpha({f})", A_REAL, math=r"\f$\alpha_" + str(f) + r"\f$")
            _r(f"{px}alpha_rho({f})", A_REAL, math=r"\f$\alpha \rho\f$")
        # Elasticity stress tensor
        for j in range(1, 7):
            _r(f"{px}tau_e({j})", A_REAL, {"elasticity"}, math=r"\f$\tau_e\f$")
        if i >= 2:
            for j in range(1, i):
                _r(f"{px}alter_patch({j})", LOG)
        # 2D modal (geometry 13): Fourier modes and options
        for j in range(1, 11):
            _r(f"{px}fourier_cos({j})", REAL)
            _r(f"{px}fourier_sin({j})", REAL)
        _r(f"{px}modal_clip_r_to_min", LOG)
        _r(f"{px}modal_r_min", REAL)
        _r(f"{px}modal_use_exp_form", LOG)
        # 3D spherical harmonic (geometry 14): coeffs (l, m), l=0..5, m=-l..l
        for ll in range(0, 6):
            for mm in range(-ll, ll + 1):
                _r(f"{px}sph_har_coeff({ll},{mm})", REAL)

    # fluid_pp (10 fluids)
    for f in range(1, NF + 1):
        px = f"fluid_pp({f})%"
        for a, sym in [("gamma", r"\f$\gamma_k\f$"), ("pi_inf", r"\f$\pi_{\infty,k}\f$"), ("cv", r"\f$c_{v,k}\f$"), ("qv", r"\f$q_{v,k}\f$"), ("qvp", r"\f$q'_{v,k}\f$")]:
            _r(f"{px}{a}", REAL, math=sym)
        _r(f"{px}mul0", REAL, {"viscosity"}, math=r"\f$\mu_{l,k}\f$")
        _r(f"{px}ss", REAL, {"surface_tension"}, math=r"\f$\sigma_k\f$")
        for a in ["pv", "gamma_v", "M_v", "mu_v", "k_v", "cp_v", "D_v"]:
            _r(f"{px}{a}", REAL, {"bubbles"})
        _r(f"{px}G", REAL, {"elasticity"}, math=r"\f$G_k\f$")
        _r(f"{px}Re(1)", REAL, {"viscosity"}, math=r"\f$\mathrm{Re}_k\f$ (shear)")
        _r(f"{px}Re(2)", REAL, {"viscosity"}, math=r"\f$\mathrm{Re}_k\f$ (bulk)")

    # bub_pp (bubble properties)
    for a, sym in [
        ("R0ref", r"\f$R_0\f$"),
        ("p0ref", r"\f$p_0\f$"),
        ("rho0ref", r"\f$\rho_l\f$"),
        ("T0ref", r"\f$T_0\f$"),
        ("ss", r"\f$\sigma\f$"),
        ("pv", r"\f$p_v\f$"),
        ("vd", r"\f$D\f$"),
        ("mu_l", r"\f$\mu_l\f$"),
        ("mu_v", r"\f$\mu_v\f$"),
        ("mu_g", r"\f$\mu_g\f$"),
        ("gam_v", r"\f$\gamma_v\f$"),
        ("gam_g", r"\f$\gamma_g\f$"),
        ("M_v", r"\f$M_v\f$"),
        ("M_g", r"\f$M_g\f$"),
        ("k_v", r"\f$k_v\f$"),
        ("k_g", r"\f$k_g\f$"),
        ("cp_v", r"\f$c_{p,v}\f$"),
        ("cp_g", r"\f$c_{p,g}\f$"),
        ("R_v", r"\f$R_v\f$"),
        ("R_g", r"\f$R_g\f$"),
    ]:
        _r(f"bub_pp%{a}", REAL, {"bubbles"}, math=sym)

    # patch_ib (immersed boundaries) — registered as indexed family for O(1) lookup.
    # max_index is None so the parameter registry stays compact (no enumeration).
    # The Fortran-side upper bound (num_ib_patches_max in m_constants.fpp) is parsed
    # and enforced by the case_validator, not by max_index here.
    _ib_tags = {"ib"}
    _ib_attrs: Dict[str, tuple] = {}
    for a in ["geometry", "moving_ibm"]:
        _ib_attrs[a] = (INT, _ib_tags)
    for a, pt in [("radius", REAL), ("theta", REAL), ("slip", LOG), ("c", REAL), ("p", REAL), ("t", REAL), ("m", REAL), ("mass", REAL)]:
        _ib_attrs[a] = (pt, _ib_tags)
    for j in range(1, 4):
        _ib_attrs[f"angles({j})"] = (REAL, _ib_tags)
    for d in ["x", "y", "z"]:
        _ib_attrs[f"{d}_centroid"] = (REAL, _ib_tags)
        _ib_attrs[f"length_{d}"] = (REAL, _ib_tags)
    for a, pt in [("model_filepath", STR), ("model_spc", INT), ("model_threshold", REAL)]:
        _ib_attrs[a] = (pt, _ib_tags)
    for t in ["translate", "scale", "rotate"]:
        for j in range(1, 4):
            _ib_attrs[f"model_{t}({j})"] = (REAL, _ib_tags)
    for j in range(1, 4):
        _ib_attrs[f"vel({j})"] = (A_REAL, _ib_tags)
        _ib_attrs[f"angular_vel({j})"] = (A_REAL, _ib_tags)
    REGISTRY.register_family(
        IndexedFamily(
            base_name="patch_ib",
            attrs=_ib_attrs,
            tags=_ib_tags,
            max_index=NIB,
        )
    )

    # acoustic sources (4 sources)
    for i in range(1, NA + 1):
        px = f"acoustic({i})%"
        for a in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            _r(f"{px}{a}", INT, {"acoustic"})
        _r(f"{px}dipole", LOG, {"acoustic"})
        for a in [
            "mag",
            "length",
            "height",
            "wavelength",
            "frequency",
            "gauss_sigma_dist",
            "gauss_sigma_time",
            "npulse",
            "dir",
            "delay",
            "foc_length",
            "aperture",
            "element_spacing_angle",
            "element_polygon_ratio",
            "rotate_angle",
            "bb_bandwidth",
            "bb_lowest_freq",
        ]:
            _r(f"{px}{a}", REAL, {"acoustic"})
        for j in range(1, 4):
            _r(f"{px}loc({j})", REAL, {"acoustic"})

    # probes (10 probes)
    for i in range(1, NPR + 1):
        for d in ["x", "y", "z"]:
            _r(f"probe({i})%{d}", REAL, {"probes"})

    # integrals (5 integral regions)
    for i in range(1, 6):
        for d in ["x", "y", "z"]:
            _r(f"integral({i})%{d}min", REAL, {"probes"})
            _r(f"integral({i})%{d}max", REAL, {"probes"})

    # Extended BC
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

        for a in ["Twall_in", "Twall_out"]:
            _r(f"{px}{a}", REAL, {"bc"})
        for a in ["isothermal_in", "isothermal_out"]:
            _r(f"{px}{a}", LOG, {"bc"})

    # patch_bc (10 BC patches)
    for i in range(1, NB + 1):
        px = f"patch_bc({i})%"
        for a in ["geometry", "type", "dir", "loc"]:
            _r(f"{px}{a}", INT, {"bc"})
        for j in range(1, 4):
            _r(f"{px}centroid({j})", REAL, {"bc"})
            _r(f"{px}length({j})", REAL, {"bc"})
        _r(f"{px}radius", REAL, {"bc"})

    # simplex_params
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

    # lag_params (Lagrangian bubbles)
    for a in ["heatTransfer_model", "massTransfer_model", "pressure_corrector", "write_bubbles", "write_bubbles_stats"]:
        _r(f"lag_params%{a}", LOG, {"bubbles"})
    for a in ["solver_approach", "cluster_type", "smooth_type", "nBubs_glb"]:
        _r(f"lag_params%{a}", INT, {"bubbles"})
    for a in ["epsilonb", "valmaxvoid", "charwidth", "c0", "rho0", "T0", "x0", "Thost"]:
        _r(f"lag_params%{a}", REAL, {"bubbles"})

    # chem_params
    for a in ["diffusion", "reactions"]:
        _r(f"chem_params%{a}", LOG, {"chemistry"})
    for a in ["gamma_method", "transport_model"]:
        _r(f"chem_params%{a}", INT, {"chemistry"})

    # Per-fluid output arrays
    for f in range(1, NF + 1):
        _r(f"schlieren_alpha({f})", REAL, {"output"})
        for a in ["alpha_rho_wrt", "alpha_wrt", "alpha_rho_e_wrt"]:
            _r(f"{a}({f})", LOG, {"output"})
    for j in range(1, 4):
        for a in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
            _r(f"{a}({j})", LOG, {"output"})

    # chem_wrt (chemistry output)
    for j in range(1, 101):
        _r(f"chem_wrt_Y({j})", LOG, {"chemistry", "output"})
    _r("chem_wrt_T", LOG, {"chemistry", "output"})

    # fluid_rho
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
        raise RuntimeError(f"Failed to initialize parameter registry: {e}\nThis is likely a bug in the parameter definitions.") from e


_init_registry()

# Namelist target mapping for Fortran codegen.
# Maps each Fortran namelist root variable to the set of MFC executables whose
# namelist it appears in. Used by fortran_gen.py to generate per-target .fpp files.
#
# When adding a new parameter:
#   1. Add to definitions.py (type, constraints, etc.) — you are here
#   2. Add the namelist root variable to NAMELIST_VARS with its target set
#   3. Re-run cmake to regenerate the .fpp files (cmake reconfigure)

NAMELIST_VARS: dict[str, set[str]] = {}

# Maps indexed-family base names to their Fortran dimension expression.
# The generator emits `{type}, dimension({dim}) :: {name}` for each entry.
# Add here whenever a new array param needs no manual Fortran declaration.
FORTRAN_ARRAY_DIMS: dict[str, str] = {
    "fluid_rho": "num_fluids_max",
    "alpha_rho_wrt": "num_fluids_max",
    "alpha_rho_e_wrt": "num_fluids_max",
    "alpha_wrt": "num_fluids_max",
    "schlieren_alpha": "num_fluids_max",
    "chem_wrt_Y": "num_species",
    "flux_wrt": "3",
    "mom_wrt": "3",
    "omega_wrt": "3",
    "vel_wrt": "3",
}


def _nv(targets: set, *names: str) -> None:
    for n in names:
        NAMELIST_VARS[n] = set(targets)


_ALL, _PRE_SIM, _SIM_POST = {"pre", "sim", "post"}, {"pre", "sim"}, {"sim", "post"}
_PRE_POST = {"pre", "post"}
_SIM = {"sim"}
_PRE = {"pre"}
_POST = {"post"}

_nv(
    _ALL,
    "m",
    "n",
    "p",
    "cyl_coord",
    "bc_x",
    "bc_y",
    "bc_z",
    "num_bc_patches",
    "case_dir",
    "t_step_start",
    "cfl_adap_dt",
    "cfl_const_dt",
    "n_start",
    "model_eqns",
    "mpp_lim",
    "relax",
    "relax_model",
    "fluid_pp",
    "bub_pp",
    "rhoref",
    "pref",
    "bubbles_euler",
    "bubbles_lagrange",
    "R0ref",
    "polytropic",
    "thermal",
    "Ca",
    "Web",
    "Re_inv",
    "polydisperse",
    "poly_sigma",
    "qbmm",
    "sigma",
    "adv_n",
    "hypoelasticity",
    "hyperelasticity",
    "surface_tension",
    "relativity",
    "ib",
    "num_ibs",
    "cont_damage",
    "hyper_cleaning",
    "Bx0",
    "precision",
    "parallel_io",
    "file_per_process",
    "fft_wrt",
    "down_sample",
)
_nv(
    _SIM_POST,
    "t_step_stop",
    "t_step_save",
    "t_stop",
    "t_save",
    "cfl_target",
    "avg_state",
    "prim_vars_wrt",
    "alt_soundspeed",
    "mixture_err",
    "fd_order",
    "ib_state_wrt",
)
_nv(
    _PRE_SIM,
    "x_domain",
    "y_domain",
    "z_domain",
    "x_a",
    "y_a",
    "z_a",
    "x_b",
    "y_b",
    "z_b",
    "palpha_eps",
    "ptgalpha_eps",
    "t_step_old",
    "patch_ib",
    "pi_fac",
)
_nv(_PRE_POST, "num_fluids", "weno_order", "recon_type", "muscl_order", "mhd", "nb", "sigR", "igr", "igr_order")
_nv(
    _SIM,
    "dt",
    "t_step_print",
    "time_stepper",
    "adap_dt",
    "adap_dt_tol",
    "adap_dt_max_iters",
    "weno_eps",
    "teno_CT",
    "wenoz_q",
    "mp_weno",
    "weno_avg",
    "weno_Re_flux",
    "null_weights",
    "muscl_eps",
    "int_comp",
    "ic_eps",
    "ic_beta",
    "riemann_solver",
    "wave_speeds",
    "low_Mach",
    "hyper_cleaning_speed",
    "hyper_cleaning_tau",
    "run_time_info",
    "bubble_model",
    "lag_params",
    "probe_wrt",
    "num_probes",
    "probe",
    "integral_wrt",
    "num_integrals",
    "integral",
    "acoustic_source",
    "num_source",
    "acoustic",
    "chem_params",
    "bf_x",
    "bf_y",
    "bf_z",
    "k_x",
    "k_y",
    "k_z",
    "w_x",
    "w_y",
    "w_z",
    "p_x",
    "p_y",
    "p_z",
    "g_x",
    "g_y",
    "g_z",
    "collision_model",
    "coefficient_of_restitution",
    "collision_time",
    "ib_coefficient_of_friction",
    "tau_star",
    "cont_damage_s",
    "alpha_bar",
    "rdma_mpi",
    "alf_factor",
    "num_igr_iters",
    "num_igr_warm_start_iters",
    "igr_iter_solver",
    "igr_pres_lim",
    "nv_uvm_out_of_core",
    "nv_uvm_igr_temps_on_gpu",
    "nv_uvm_pref_gpu",
)
_nv(
    _PRE,
    "stretch_x",
    "stretch_y",
    "stretch_z",
    "a_x",
    "a_y",
    "a_z",
    "loops_x",
    "loops_y",
    "loops_z",
    "n_start_old",
    "num_patches",
    "patch_icpp",
    "patch_bc",
    "sigV",
    "dist_type",
    "rhoRV",
    "viscous",
    "old_grid",
    "old_ic",
    "perturb_flow",
    "perturb_flow_fluid",
    "perturb_flow_mag",
    "perturb_sph",
    "perturb_sph_fluid",
    "fluid_rho",
    "mixlayer_vel_profile",
    "mixlayer_vel_coef",
    "mixlayer_perturb",
    "mixlayer_perturb_nk",
    "mixlayer_perturb_k0",
    "pre_stress",
    "elliptic_smoothing",
    "elliptic_smoothing_iters",
    "simplex_perturb",
    "simplex_params",
)
_nv(
    _POST,
    "x_output",
    "y_output",
    "z_output",
    "format",
    "output_partial_domain",
    "sim_data",
    "G",
    "flux_lim",
    "cons_vars_wrt",
    "rho_wrt",
    "E_wrt",
    "pres_wrt",
    "c_wrt",
    "gamma_wrt",
    "heat_ratio_wrt",
    "pi_inf_wrt",
    "pres_inf_wrt",
    "omega_wrt",
    "qm_wrt",
    "liutex_wrt",
    "schlieren_wrt",
    "schlieren_alpha",
    "alpha_rho_wrt",
    "mom_wrt",
    "vel_wrt",
    "flux_wrt",
    "alpha_wrt",
    "cf_wrt",
    "chem_wrt_T",
    "chem_wrt_Y",
    "alpha_rho_e_wrt",
    "lag_header",
    "lag_txt_wrt",
    "lag_db_wrt",
    "lag_id_wrt",
    "lag_pos_wrt",
    "lag_pos_prev_wrt",
    "lag_vel_wrt",
    "lag_rad_wrt",
    "lag_rvel_wrt",
    "lag_r0_wrt",
    "lag_rmax_wrt",
    "lag_rmin_wrt",
    "lag_dphidt_wrt",
    "lag_pres_wrt",
    "lag_mv_wrt",
    "lag_mg_wrt",
    "lag_betaT_wrt",
    "lag_betaC_wrt",
)

# Case-optimization params appear in the sim namelist under #:if not MFC_CASE_OPTIMIZATION.
for _v in CASE_OPT_PARAMS:
    NAMELIST_VARS.setdefault(_v, set()).add("sim")
