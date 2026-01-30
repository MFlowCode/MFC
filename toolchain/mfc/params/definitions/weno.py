"""
WENO Parameter Definitions.

Defines WENO (Weighted Essentially Non-Oscillatory) scheme parameters:
- weno_order: Order of WENO reconstruction
- mapped_weno, wenoz, teno: WENO variants (mutually exclusive)
- weno_eps: Small regularization parameter
"""

from ..schema import ParamDef, ConstraintRule, ParamType, Stage, ConstraintType
from ..registry import REGISTRY


# =============================================================================
# WENO Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="weno_order",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="WENO reconstruction order: 1, 3, 5, or 7",
    category="weno",
    case_optimization=True,
))

REGISTRY.register(ParamDef(
    name="mapped_weno",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Use mapped WENO (improves accuracy near discontinuities)",
    category="weno",
))

REGISTRY.register(ParamDef(
    name="wenoz",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Use WENO-Z (improved smoothness indicators)",
    category="weno",
))

REGISTRY.register(ParamDef(
    name="teno",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Use TENO (Targeted Essentially Non-Oscillatory)",
    category="weno",
))

REGISTRY.register(ParamDef(
    name="weno_eps",
    param_type=ParamType.REAL,
    stages={Stage.SIMULATION},
    default=1e-6,
    description="WENO regularization parameter (prevents division by zero)",
    category="weno",
))

REGISTRY.register(ParamDef(
    name="weno_Re_flux",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Apply WENO to viscous flux reconstruction",
    category="weno",
))

REGISTRY.register(ParamDef(
    name="weno_avg",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Use averaged WENO reconstruction",
    category="weno",
))


# =============================================================================
# WENO Constraints
# =============================================================================

REGISTRY.register_constraint(ConstraintRule(
    rule_id="WENO_ORDER_CHOICES",
    params=["weno_order"],
    constraint_type=ConstraintType.CHOICES,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    valid_values=[1, 3, 5, 7],
    message="weno_order must be 1, 3, 5, or 7",
    check_method="check_weno",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="WENO_VARIANTS_MUTEX",
    params=["mapped_weno", "wenoz", "teno"],
    constraint_type=ConstraintType.MUTEX,
    stages={Stage.SIMULATION},
    mutex_params=["mapped_weno", "wenoz", "teno"],
    message="Only one of mapped_weno, wenoz, or teno can be enabled",
    check_method="check_weno",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="WENO_EPS_POSITIVE",
    params=["weno_eps"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.SIMULATION},
    message="weno_eps must be positive",
    check_method="check_weno",
    priority=30,
))
