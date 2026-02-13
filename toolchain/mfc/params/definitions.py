"""
MFC Parameter Definitions (Compact).

Single file containing all ~3,300 parameter definitions using loops.
This replaces the definitions/ directory.
"""  # pylint: disable=too-many-lines

import re
from typing import Dict, Any
from .schema import ParamDef, ParamType
from .registry import REGISTRY

# Index limits
NP, NF, NI, NA, NPR, NB = 10, 10, 10, 4, 10, 10  # patches, fluids, ibs, acoustic, probes, bc_patches


# =============================================================================
# Auto-generated Descriptions
# =============================================================================
# Descriptions are auto-generated from parameter names using naming conventions.
# Override with explicit desc= parameter when auto-generation is inadequate.

# Prefix descriptions for indexed parameter families
_PREFIX_DESCS = {
    "patch_icpp": "initial condition patch",
    "patch_ib": "immersed boundary",
    "patch_bc": "boundary condition patch",
    "fluid_pp": "fluid",
    "acoustic": "acoustic source",
    "probe": "probe",
    "integral": "integral region",
}

# Attribute descriptions (suffix after %)
_ATTR_DESCS = {
    # Geometry/position
    "geometry": "Geometry type",
    "x_centroid": "X-coordinate of centroid",
    "y_centroid": "Y-coordinate of centroid",
    "z_centroid": "Z-coordinate of centroid",
    "length_x": "X-dimension length",
    "length_y": "Y-dimension length",
    "length_z": "Z-dimension length",
    "radius": "Radius",
    "radii": "Radii array",
    "normal": "Normal direction",
    "theta": "Theta angle",
    "angles": "Orientation angles",
    # Physics
    "vel": "Velocity",
    "pres": "Pressure",
    "rho": "Density",
    "alpha": "Volume fraction",
    "alpha_rho": "Partial density",
    "gamma": "Specific heat ratio",
    "pi_inf": "Stiffness pressure",
    "cv": "Specific heat (const. volume)",
    "qv": "Heat of formation",
    "qvp": "Heat of formation prime",
    "G": "Shear modulus",
    "Re": "Reynolds number",
    "mul0": "Reference viscosity",
    "ss": "Surface tension",
    "pv": "Vapor pressure",
    # MHD
    "Bx": "Magnetic field (x-component)",
    "By": "Magnetic field (y-component)",
    "Bz": "Magnetic field (z-component)",
    # Model/smoothing
    "smoothen": "Enable smoothing",
    "smooth_patch_id": "Patch ID to smooth against",
    "smooth_coeff": "Smoothing coefficient",
    "alter_patch": "Alter with another patch",
    "model_filepath": "STL model file path",
    "model_spc": "Model spacing",
    "model_threshold": "Model threshold",
    "model_translate": "Model translation",
    "model_scale": "Model scale",
    "model_rotate": "Model rotation",
    # Bubbles
    "r0": "Initial bubble radius",
    "v0": "Initial bubble velocity",
    "p0": "Initial bubble pressure",
    "m0": "Initial bubble mass",
    # IB specific
    "slip": "Enable slip condition",
    "moving_ibm": "Enable moving boundary",
    "angular_vel": "Angular velocity",
    "mass": "Mass",
    # BC specific
    "vel_in": "Inlet velocity",
    "vel_out": "Outlet velocity",
    "alpha_rho_in": "Inlet partial density",
    "alpha_in": "Inlet volume fraction",
    "pres_in": "Inlet pressure",
    "pres_out": "Outlet pressure",
    "grcbc_in": "Enable GRCBC inlet",
    "grcbc_out": "Enable GRCBC outlet",
    "grcbc_vel_out": "Enable GRCBC velocity outlet",
    # Acoustic
    "loc": "Location",
    "mag": "Magnitude",
    "pulse": "Pulse type",
    "support": "Support type",
    "frequency": "Frequency",
    "wavelength": "Wavelength",
    "length": "Length",
    "height": "Height",
    "delay": "Delay time",
    "dipole": "Enable dipole",
    "dir": "Direction",
    # Output
    "x": "X-coordinate",
    "y": "Y-coordinate",
    "z": "Z-coordinate",
    "xmin": "X minimum",
    "xmax": "X maximum",
    "ymin": "Y minimum",
    "ymax": "Y maximum",
    "zmin": "Z minimum",
    "zmax": "Z maximum",
    # Chemistry
    "Y": "Species mass fraction",
    # Shape coefficients
    "a": "Shape coefficient",
    # Elasticity
    "tau_e": "Elastic stress component",
    # Misc
    "cf_val": "Color function value",
    "hcid": "Hard-coded ID",
    "epsilon": "Interface thickness",
    "beta": "Shape parameter beta",
    "non_axis_sym": "Non-axisymmetric parameter",
}

# Simple parameter descriptions (non-indexed)
_SIMPLE_DESCS = {
    # Grid
    "m": "Grid cells in x-direction",
    "n": "Grid cells in y-direction",
    "p": "Grid cells in z-direction",
    "cyl_coord": "Enable cylindrical coordinates",
    "stretch_x": "Enable grid stretching in x",
    "stretch_y": "Enable grid stretching in y",
    "stretch_z": "Enable grid stretching in z",
    "a_x": "Grid stretching rate in x",
    "a_y": "Grid stretching rate in y",
    "a_z": "Grid stretching rate in z",
    "x_a": "Stretching start (negative x)",
    "x_b": "Stretching start (positive x)",
    "y_a": "Stretching start (negative y)",
    "y_b": "Stretching start (positive y)",
    "z_a": "Stretching start (negative z)",
    "z_b": "Stretching start (positive z)",
    "loops_x": "Stretching iterations in x",
    "loops_y": "Stretching iterations in y",
    "loops_z": "Stretching iterations in z",
    # Time
    "dt": "Time step size",
    "t_step_start": "Starting time step",
    "t_step_stop": "Ending time step",
    "t_step_save": "Save interval (steps)",
    "t_step_print": "Print interval (steps)",
    "t_stop": "Stop time",
    "t_save": "Save interval (time)",
    "time_stepper": "Time integration scheme",
    "cfl_target": "Target CFL number",
    "cfl_max": "Maximum CFL number",
    "cfl_adap_dt": "Enable adaptive CFL time stepping",
    "cfl_const_dt": "Use constant CFL time stepping",
    "cfl_dt": "Enable CFL-based time stepping",
    "adap_dt": "Enable adaptive time stepping",
    "adap_dt_tol": "Adaptive time stepping tolerance",
    "adap_dt_max_iters": "Max iterations for adaptive dt",
    "t_tol": "Time tolerance",
    # Model
    "model_eqns": "Model equations",
    "num_fluids": "Number of fluids",
    "num_patches": "Number of IC patches",
    "mpp_lim": "Mixture pressure positivity limiter",
    # WENO
    "weno_order": "WENO reconstruction order",
    "weno_eps": "WENO epsilon parameter",
    "mapped_weno": "Enable mapped WENO",
    "wenoz": "Enable WENO-Z",
    "teno": "Enable TENO",
    "mp_weno": "Enable monotonicity-preserving WENO",
    # Riemann
    "riemann_solver": "Riemann solver",
    "wave_speeds": "Wave speed estimate method",
    "avg_state": "Average state",
    # Physics toggles
    "viscous": "Enable viscous effects",
    "mhd": "Enable magnetohydrodynamics",
    "hyper_cleaning": "Enable hyperbolic divergence cleaning",
    "hyper_cleaning_speed": "Divergence cleaning wave speed",
    "hyper_cleaning_tau": "Divergence cleaning damping time",
    "powell": "Enable Powell source terms for MHD",
    "bubbles_euler": "Enable Euler bubble model",
    "bubbles_lagrange": "Enable Lagrangian bubbles",
    "polytropic": "Enable polytropic gas",
    "polydisperse": "Enable polydisperse bubbles",
    "qbmm": "Enable QBMM",
    "chemistry": "Enable chemistry",
    "surface_tension": "Enable surface tension",
    "hypoelasticity": "Enable hypoelastic model",
    "hyperelasticity": "Enable hyperelastic model",
    "relativity": "Enable special relativity",
    "ib": "Enable immersed boundaries",
    "acoustic_source": "Enable acoustic sources",
    # Output
    "parallel_io": "Enable parallel I/O",
    "probe_wrt": "Write probe data",
    "prim_vars_wrt": "Write primitive variables",
    "cons_vars_wrt": "Write conservative variables",
    "run_time_info": "Print runtime info",
    # Misc
    "case_dir": "Case directory path",
    "cantera_file": "Cantera mechanism file",
    "num_ibs": "Number of immersed boundaries",
    "num_source": "Number of acoustic sources",
    "num_probes": "Number of probes",
    "num_integrals": "Number of integral regions",
    "nb": "Number of bubble bins",
    "R0ref": "Reference bubble radius",
    "sigma": "Surface tension coefficient",
    "Bx0": "Background magnetic field (x)",
    "old_grid": "Load grid from previous simulation",
    "old_ic": "Load initial conditions from previous",
    "t_step_old": "Time step to restart from",
    "fd_order": "Finite difference order",
    "recon_type": "Reconstruction type",
    "muscl_order": "MUSCL reconstruction order",
    "muscl_lim": "MUSCL limiter type",
    "low_Mach": "Low Mach number correction",
    "bubble_model": "Bubble dynamics model",
    "Ca": "Cavitation number",
    "Web": "Weber number",
    "Re_inv": "Inverse Reynolds number",
    "format": "Output format",
    "precision": "Output precision",
    # Body forces
    "bf_x": "Enable body force in x",
    "bf_y": "Enable body force in y",
    "bf_z": "Enable body force in z",
    "k_x": "Body force wavenumber in x",
    "k_y": "Body force wavenumber in y",
    "k_z": "Body force wavenumber in z",
    "w_x": "Body force frequency in x",
    "w_y": "Body force frequency in y",
    "w_z": "Body force frequency in z",
    "p_x": "Body force phase in x",
    "p_y": "Body force phase in y",
    "p_z": "Body force phase in z",
    "g_x": "Gravitational acceleration in x",
    "g_y": "Gravitational acceleration in y",
    "g_z": "Gravitational acceleration in z",
    # More output
    "E_wrt": "Write energy field",
    "c_wrt": "Write sound speed field",
    "rho_wrt": "Write density field",
    "pres_wrt": "Write pressure field",
    "schlieren_wrt": "Write schlieren images",
    "cf_wrt": "Write color function",
    "omega_wrt": "Write vorticity",
    "qm_wrt": "Write Q-criterion",
    "liutex_wrt": "Write Liutex vortex field",
    "gamma_wrt": "Write gamma field",
    "heat_ratio_wrt": "Write heat capacity ratio",
    "pi_inf_wrt": "Write pi_inf field",
    "pres_inf_wrt": "Write reference pressure",
    "fft_wrt": "Write FFT output",
    "kappa_wrt": "Write curvature field",
    "chem_wrt_T": "Write temperature (chemistry)",
    # Misc physics
    "alt_soundspeed": "Alternative sound speed formulation",
    "mixture_err": "Enable mixture error checking",
    "cont_damage": "Enable continuum damage model",
}


def _auto_describe(name: str) -> str:
    """Auto-generate description from parameter name."""
    # Check simple params first
    if name in _SIMPLE_DESCS:
        return _SIMPLE_DESCS[name]

    # Handle indexed params: prefix(N)%attr or prefix(N)%attr(M)
    match = re.match(r"([a-z_]+)\((\d+)\)%(.+)", name)
    if match:
        prefix, idx, attr = match.group(1), match.group(2), match.group(3)
        prefix_desc = _PREFIX_DESCS.get(prefix, prefix.replace("_", " "))

        # Check for nested index: attr(M) or attr(M, K)
        attr_match = re.match(r"([a-z_]+)\((\d+)(?:,\s*(\d+))?\)", attr)
        if attr_match:
            attr_base = attr_match.group(1)
            idx2 = attr_match.group(2)
            attr_desc = _ATTR_DESCS.get(attr_base, attr_base.replace("_", " "))
            return f"{attr_desc} {idx2} for {prefix_desc} {idx}"

        attr_desc = _ATTR_DESCS.get(attr, attr.replace("_", " "))
        return f"{attr_desc} for {prefix_desc} {idx}"

    # Handle bc_x%attr style (no index in prefix)
    if "%" in name:
        prefix, attr = name.split("%", 1)
        # Check for indexed attr
        attr_match = re.match(r"([a-z_]+)\((\d+)\)", attr)
        if attr_match:
            attr_base, idx = attr_match.group(1), attr_match.group(2)
            attr_desc = _ATTR_DESCS.get(attr_base, attr_base.replace("_", " "))
            return f"{attr_desc} {idx} for {prefix.replace('_', ' ')}"

        attr_desc = _ATTR_DESCS.get(attr, "")
        if attr_desc:
            return f"{attr_desc} for {prefix.replace('_', ' ')}"
        # Fallback: just clean up the name
        return f"{attr.replace('_', ' ').title()} for {prefix.replace('_', ' ')}"

    # Handle suffix-indexed: name(N) or name(N, M)
    match = re.match(r"([a-z_]+)\((\d+)(?:,\s*(\d+))?\)", name)
    if match:
        base, idx = match.group(1), match.group(2)
        # Handle _wrt patterns
        if base.endswith("_wrt"):
            field = base[:-4].replace("_", " ")
            return f"Write {field} for component {idx}"
        return f"{base.replace('_', ' ').title()} {idx}"

    # Fallback patterns
    if name.endswith("_wrt"):
        return f"Write {name[:-4].replace('_', ' ')}"
    if name.startswith("num_"):
        return f"Number of {name[4:].replace('_', ' ')}"

    # Last resort: clean up the name
    return name.replace("_", " ").replace("%", " ")

# Parameters that can be hard-coded for GPU case optimization
CASE_OPT_PARAMS = {
    "mapped_weno", "wenoz", "teno", "wenoz_q", "nb", "weno_order",
    "num_fluids", "mhd", "relativity", "igr_order", "viscous",
    "igr_iter_solver", "igr", "igr_pres_lim", "recon_type",
    "muscl_order", "muscl_lim"
}


# =============================================================================
# Data-driven Annotations for Doc Generation
# =============================================================================
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
    if '%' not in name:
        # Check PREFIX_HINTS for simple params
        for prefix, label in PREFIX_HINTS.items():
            if name.startswith(prefix):
                return label
        return ""
    # Compound name: extract family and attribute
    prefix, attr_full = name.split('%', 1)
    # Normalize family: "bc_x" → "bc", "patch_bc(1)" → "patch"
    family = re.sub(r'[_(].*', '', prefix)
    if family not in HINTS:
        # Fallback: keep underscores — "patch_bc" → "patch_bc", "simplex_params" → "simplex_params"
        m = re.match(r'^[a-zA-Z_]+', prefix)
        family = m.group(0) if m else ""
    if family not in HINTS:
        return ""
    # Strip index from attr: "vel_in(1)" → "vel_in"
    m = re.match(r'^[a-zA-Z_0-9]+', attr_full)
    if not m:
        return ""
    attr = m.group(0)
    return HINTS[family].get(attr, "")


# ============================================================================
# Schema Validation for Constraints and Dependencies
# ============================================================================
# Uses rapidfuzz for "did you mean?" suggestions when typos are detected

_VALID_CONSTRAINT_KEYS = {"choices", "min", "max", "value_labels"}
_VALID_DEPENDENCY_KEYS = {"when_true", "when_set", "when_value"}
_VALID_CONDITION_KEYS = {"requires", "recommends", "requires_value"}


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
    if "value_labels" in constraint:
        if not isinstance(constraint["value_labels"], dict):
            raise ValueError(f"Constraint 'value_labels' for '{param_name}' must be a dict")
        if "choices" in constraint:
            for key in constraint["value_labels"]:
                if key not in constraint["choices"]:
                    raise ValueError(
                        f"value_labels key {key!r} for '{param_name}' "
                        f"not in choices {constraint['choices']}"
                    )


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

    def _validate_condition(cond_label: str, condition: Any) -> None:
        """Validate a condition dict (shared by when_true, when_set, when_value entries)."""
        if not isinstance(condition, dict):
            raise ValueError(
                f"Dependency '{cond_label}' for '{param_name}' must be a dict"
            )
        invalid_cond_keys = set(condition.keys()) - _VALID_CONDITION_KEYS
        if invalid_cond_keys:
            first_invalid = next(iter(invalid_cond_keys))
            raise ValueError(
                invalid_key_error(
                    f"condition in '{cond_label}' for '{param_name}'",
                    first_invalid,
                    _VALID_CONDITION_KEYS
                )
            )
        for req_key in ["requires", "recommends"]:
            if req_key in condition and not isinstance(condition[req_key], list):
                raise ValueError(
                    f"Dependency '{cond_label}/{req_key}' for '{param_name}' "
                    "must be a list"
                )
        if "requires_value" in condition:
            rv = condition["requires_value"]
            if not isinstance(rv, dict):
                raise ValueError(
                    f"Dependency '{cond_label}/requires_value' for '{param_name}' "
                    "must be a dict"
                )
            for rv_param, rv_vals in rv.items():
                if not isinstance(rv_vals, list):
                    raise ValueError(
                        f"Dependency '{cond_label}/requires_value/{rv_param}' "
                        f"for '{param_name}' must be a list"
                    )

    for condition_key in ["when_true", "when_set"]:
        if condition_key in dependency:
            _validate_condition(condition_key, dependency[condition_key])

    if "when_value" in dependency:
        wv = dependency["when_value"]
        if not isinstance(wv, dict):
            raise ValueError(
                f"Dependency 'when_value' for '{param_name}' must be a dict"
            )
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
        "choices": [1, 2, 3, 4, 5],
        "value_labels": {1: "minmod", 2: "MC", 3: "Van Albada", 4: "Van Leer", 5: "SUPERBEE"},
    },

    # Time stepping
    "time_stepper": {
        "choices": [1, 2, 3],
        "value_labels": {1: "RK1 (Forward Euler)", 2: "RK2", 3: "RK3 (SSP)"},
    },

    # Riemann solver
    "riemann_solver": {
        "choices": [1, 2, 3, 4, 5],
        "value_labels": {1: "HLL", 2: "HLLC", 3: "Exact", 4: "HLLD", 5: "Lax-Friedrichs"},
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
    "cfl_max": {"min": 0},

    # WENO
    "weno_eps": {"min": 0},

    # Physics (must be non-negative)
    "R0ref": {"min": 0},
    "sigma": {"min": 0},

    # Counts (must be positive)
    "num_fluids": {"min": 1, "max": 10},
    "num_patches": {"min": 0, "max": 10},
    "num_ibs": {"min": 0, "max": 10},
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

def _r(name, ptype, tags=None, desc=None, hint=None):
    """Register a parameter with optional feature tags and description."""
    if hint is None:
        hint = _lookup_hint(name)
    description = desc if desc else _auto_describe(name)
    constraint = CONSTRAINTS.get(name)
    if constraint and "value_labels" in constraint:
        labels = constraint["value_labels"]
        suffix = ", ".join(f"{v}={labels[v]}" for v in sorted(labels))
        description = f"{description} ({suffix})"
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        description=description,
        case_optimization=(name in CASE_OPT_PARAMS),
        constraints=constraint,
        dependencies=DEPENDENCIES.get(name),
        tags=tags if tags else set(),
        hint=hint,
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
