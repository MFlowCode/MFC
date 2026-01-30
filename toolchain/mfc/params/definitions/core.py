"""
Core Parameter Definitions.

Defines fundamental MFC parameters:
- Grid dimensions: m, n, p
- Model configuration: model_eqns, num_fluids, num_patches
"""

from ..schema import ParamDef, ConstraintRule, ParamType, Stage, ConstraintType
from ..registry import REGISTRY


# =============================================================================
# Grid Dimension Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="m",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Number of grid cells in x-direction (must be positive)",
    category="grid",
))

REGISTRY.register(ParamDef(
    name="n",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    default=0,
    description="Number of grid cells in y-direction (0 for 1D)",
    category="grid",
))

REGISTRY.register(ParamDef(
    name="p",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    default=0,
    description="Number of grid cells in z-direction (0 for 1D/2D)",
    category="grid",
))


# =============================================================================
# Model Configuration Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="model_eqns",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Model equations: 1=pi-gamma, 2=5-equation, 3=6-equation, 4=4-equation",
    category="model",
))

REGISTRY.register(ParamDef(
    name="num_fluids",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Number of fluid components",
    category="model",
))

REGISTRY.register(ParamDef(
    name="num_patches",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Number of initial condition patches",
    category="model",
))


# =============================================================================
# Grid Constraints
# =============================================================================

REGISTRY.register_constraint(ConstraintRule(
    rule_id="M_REQUIRED",
    params=["m"],
    constraint_type=ConstraintType.REQUIRED,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="m must be set",
    check_method="check_simulation_domain",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="M_POSITIVE",
    params=["m"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="m must be positive",
    check_method="check_simulation_domain",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="N_NON_NEGATIVE",
    params=["n"],
    constraint_type=ConstraintType.NON_NEGATIVE,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="n must be non-negative",
    check_method="check_simulation_domain",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="P_NON_NEGATIVE",
    params=["p"],
    constraint_type=ConstraintType.NON_NEGATIVE,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="p must be non-negative",
    check_method="check_simulation_domain",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="P_REQUIRES_N",
    params=["n", "p"],
    constraint_type=ConstraintType.DEPENDENCY,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    if_param="p",
    if_condition="> 0",
    then_param="n",
    then_condition="> 0",
    message="3D (p > 0) requires 2D (n > 0)",
    check_method="check_simulation_domain",
    priority=30,
))


# =============================================================================
# Model Constraints
# =============================================================================

REGISTRY.register_constraint(ConstraintRule(
    rule_id="MODEL_EQNS_CHOICES",
    params=["model_eqns"],
    constraint_type=ConstraintType.CHOICES,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    valid_values=[1, 2, 3, 4],
    message="model_eqns must be 1, 2, 3, or 4",
    check_method="check_model_eqns_and_num_fluids",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="NUM_FLUIDS_POSITIVE",
    params=["num_fluids"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="num_fluids must be positive",
    check_method="check_model_eqns_and_num_fluids",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="NUM_PATCHES_POSITIVE",
    params=["num_patches"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.PRE_PROCESS},
    message="num_patches must be positive",
    check_method="check_patches",
    priority=20,
))
