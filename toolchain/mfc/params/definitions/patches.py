"""
Patch Parameter Definitions.

Defines patch_icpp parameters - the largest family of indexed parameters.
Uses loop-based generation to match case_dicts.py structure.
"""

from ..schema import ParamDef, ParamType, Stage
from ..registry import REGISTRY

# Constants matching case_dicts.py
NUM_PATCHES = 10
NUM_FLUIDS = 10


def _reg(name: str, ptype: ParamType, stages: set, category: str = "patch_icpp"):
    """Helper to register a parameter."""
    REGISTRY.register(ParamDef(
        name=name,
        param_type=ptype,
        stages=stages,
        category=category,
    ))


def register_patch_icpp_params():
    """Register all patch_icpp parameters matching case_dicts.py PRE_PROCESS."""

    for p_id in range(1, NUM_PATCHES + 1):
        prefix = f"patch_icpp({p_id})%"

        # Integer/Log attributes
        for attr, ptype in [("geometry", ParamType.INT),
                            ("smoothen", ParamType.LOG),
                            ("smooth_patch_id", ParamType.INT),
                            ("hcid", ParamType.INT)]:
            _reg(f"{prefix}{attr}", ptype, {Stage.PRE_PROCESS})

        # Real attributes (simple)
        for attr in ["radius", "radii", "epsilon", "beta", "normal", "alpha_rho",
                     "non_axis_sym", "smooth_coeff", "rho", "vel",
                     "alpha", "gamma", "pi_inf", "r0", "v0", "p0", "m0", "cv",
                     "qv", "qvp"]:
            _reg(f"{prefix}{attr}", ParamType.REAL, {Stage.PRE_PROCESS})

        # a(2) through a(9)
        for a_id in range(2, 10):
            _reg(f"{prefix}a({a_id})", ParamType.REAL, {Stage.PRE_PROCESS})

        # Analytic parameters (can be expressions)
        _reg(f"{prefix}pres", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})
        _reg(f"{prefix}Bx", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})
        _reg(f"{prefix}By", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})
        _reg(f"{prefix}Bz", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # Y(0) through Y(99) - species mass fractions
        for y_id in range(100):
            _reg(f"{prefix}Y({y_id})", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # Model STL parameters
        _reg(f"{prefix}model_filepath", ParamType.STR, {Stage.PRE_PROCESS})
        _reg(f"{prefix}model_spc", ParamType.INT, {Stage.PRE_PROCESS})
        _reg(f"{prefix}model_threshold", ParamType.REAL, {Stage.PRE_PROCESS})

        for transform in ["translate", "scale", "rotate"]:
            for j in range(1, 4):
                _reg(f"{prefix}model_{transform}({j})", ParamType.REAL, {Stage.PRE_PROCESS})

        # Centroid and length for x, y, z
        for cmp in ["x", "y", "z"]:
            _reg(f"{prefix}{cmp}_centroid", ParamType.REAL, {Stage.PRE_PROCESS})
            _reg(f"{prefix}length_{cmp}", ParamType.REAL, {Stage.PRE_PROCESS})

        # Indexed real attributes: radii, normal, vel (1-3)
        for cmp_id in range(1, 4):
            _reg(f"{prefix}radii({cmp_id})", ParamType.REAL, {Stage.PRE_PROCESS})
            _reg(f"{prefix}normal({cmp_id})", ParamType.REAL, {Stage.PRE_PROCESS})
            _reg(f"{prefix}vel({cmp_id})", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # Multi-fluid: alpha and alpha_rho for each fluid
        for f_id in range(1, NUM_FLUIDS + 1):
            _reg(f"{prefix}alpha({f_id})", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})
            _reg(f"{prefix}alpha_rho({f_id})", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # tau_e (elastic stress tensor components) - 6 components
        for tau_id in range(1, 7):
            _reg(f"{prefix}tau_e({tau_id})", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # cf_val (continuous damage field)
        _reg(f"{prefix}cf_val", ParamType.ANALYTIC_REAL, {Stage.PRE_PROCESS})

        # alter_patch - only for patches 2 and above
        if p_id >= 2:
            _reg(f"{prefix}alter_patch", ParamType.LOG, {Stage.PRE_PROCESS})
            # Can alter patches 1 through p_id-1
            for alter_id in range(1, p_id):
                _reg(f"{prefix}alter_patch({alter_id})", ParamType.LOG, {Stage.PRE_PROCESS})


def register_fluid_pp_params():
    """Register fluid_pp parameters for all stages."""

    for f_id in range(1, NUM_FLUIDS + 1):
        prefix = f"fluid_pp({f_id})%"

        # Common real attributes (PRE_PROCESS, SIMULATION, POST_PROCESS all have these)
        for attr in ["gamma", "pi_inf", "G", "cv", "qv", "qvp"]:
            # Register once with all stages - registry will merge
            REGISTRY.register(ParamDef(
                name=f"{prefix}{attr}",
                param_type=ParamType.REAL,
                stages={Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS},
                category="fluid_pp",
            ))

        # Re (Reynolds number) - SIMULATION only
        for re_id in [1, 2]:
            REGISTRY.register(ParamDef(
                name=f"{prefix}Re({re_id})",
                param_type=ParamType.REAL,
                stages={Stage.SIMULATION},
                category="fluid_pp",
            ))


def register_bub_pp_params():
    """Register bub_pp (bubble properties) parameters."""

    bub_vars = ["R0ref", "p0ref", "rho0ref", "T0ref", "ss", "pv", "vd",
                "mu_l", "mu_v", "mu_g", "gam_v", "gam_g",
                "M_v", "M_g", "k_v", "k_g", "cp_v", "cp_g", "R_v", "R_g"]

    for var in bub_vars:
        # Present in PRE_PROCESS, SIMULATION, POST_PROCESS
        REGISTRY.register(ParamDef(
            name=f"bub_pp%{var}",
            param_type=ParamType.REAL,
            stages={Stage.PRE_PROCESS, Stage.SIMULATION, Stage.POST_PROCESS},
            category="bub_pp",
        ))


# Auto-register when module is imported
register_patch_icpp_params()
register_fluid_pp_params()
register_bub_pp_params()
