"""
Domain Parameter Definitions.

Defines domain geometry and boundary condition parameters:
- x_domain%beg, x_domain%end: Domain bounds in x
- y_domain%beg, y_domain%end: Domain bounds in y
- z_domain%beg, z_domain%end: Domain bounds in z
- bc_x%beg, bc_x%end, etc.: Boundary conditions
"""

from ..schema import ParamDef, ConstraintRule, ParamType, Stage, ConstraintType
from ..registry import REGISTRY


# =============================================================================
# Domain Bounds Parameters
# =============================================================================

# X domain
REGISTRY.register(ParamDef(
    name="x_domain%beg",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain start in x-direction",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="x_domain%end",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain end in x-direction",
    category="domain",
))

# Y domain
REGISTRY.register(ParamDef(
    name="y_domain%beg",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain start in y-direction",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="y_domain%end",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain end in y-direction",
    category="domain",
))

# Z domain
REGISTRY.register(ParamDef(
    name="z_domain%beg",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain start in z-direction",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="z_domain%end",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Domain end in z-direction",
    category="domain",
))


# =============================================================================
# Boundary Condition Parameters
# =============================================================================

# X boundary conditions
REGISTRY.register(ParamDef(
    name="bc_x%beg",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at x_domain%beg",
    category="boundary_conditions",
))

REGISTRY.register(ParamDef(
    name="bc_x%end",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at x_domain%end",
    category="boundary_conditions",
))

# Y boundary conditions
REGISTRY.register(ParamDef(
    name="bc_y%beg",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at y_domain%beg",
    category="boundary_conditions",
))

REGISTRY.register(ParamDef(
    name="bc_y%end",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at y_domain%end",
    category="boundary_conditions",
))

# Z boundary conditions
REGISTRY.register(ParamDef(
    name="bc_z%beg",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at z_domain%beg",
    category="boundary_conditions",
))

REGISTRY.register(ParamDef(
    name="bc_z%end",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Boundary condition at z_domain%end",
    category="boundary_conditions",
))


# =============================================================================
# Additional Domain Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="cyl_coord",
    param_type=ParamType.LOG,
    stages={Stage.COMMON},
    default='F',
    description="Use cylindrical/axisymmetric coordinates",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="stretch_x",
    param_type=ParamType.LOG,
    stages={Stage.PRE_PROCESS},
    default='F',
    description="Enable grid stretching in x-direction",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="stretch_y",
    param_type=ParamType.LOG,
    stages={Stage.PRE_PROCESS},
    default='F',
    description="Enable grid stretching in y-direction",
    category="domain",
))

REGISTRY.register(ParamDef(
    name="stretch_z",
    param_type=ParamType.LOG,
    stages={Stage.PRE_PROCESS},
    default='F',
    description="Enable grid stretching in z-direction",
    category="domain",
))


# =============================================================================
# Domain Constraints
# =============================================================================

REGISTRY.register_constraint(ConstraintRule(
    rule_id="X_DOMAIN_BEG_REQUIRED",
    params=["x_domain%beg"],
    constraint_type=ConstraintType.REQUIRED,
    stages={Stage.PRE_PROCESS},
    message="x_domain%beg must be set",
    check_method="check_simulation_domain",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="X_DOMAIN_END_REQUIRED",
    params=["x_domain%end"],
    constraint_type=ConstraintType.REQUIRED,
    stages={Stage.PRE_PROCESS},
    message="x_domain%end must be set",
    check_method="check_simulation_domain",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="BC_X_BEG_REQUIRED",
    params=["bc_x%beg"],
    constraint_type=ConstraintType.REQUIRED,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="bc_x%beg must be set",
    check_method="check_boundary_conditions",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="BC_X_END_REQUIRED",
    params=["bc_x%end"],
    constraint_type=ConstraintType.REQUIRED,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    message="bc_x%end must be set",
    check_method="check_boundary_conditions",
    priority=10,
))

# Domain consistency: end > beg
REGISTRY.register_constraint(ConstraintRule(
    rule_id="X_DOMAIN_ORDER",
    params=["x_domain%beg", "x_domain%end"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.PRE_PROCESS},
    predicate=lambda p: p.get("x_domain%beg") is None or p.get("x_domain%end") is None or p.get("x_domain%end") > p.get("x_domain%beg"),
    message="x_domain%end must be greater than x_domain%beg",
    check_method="check_simulation_domain",
    priority=20,
))

# Y domain required if n > 0
REGISTRY.register_constraint(ConstraintRule(
    rule_id="Y_DOMAIN_REQUIRED_2D",
    params=["n", "y_domain%beg", "y_domain%end"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.PRE_PROCESS},
    predicate=lambda p: p.get("n", 0) == 0 or (p.get("y_domain%beg") is not None and p.get("y_domain%end") is not None),
    message="y_domain%beg and y_domain%end must be set for 2D/3D (n > 0)",
    check_method="check_simulation_domain",
    priority=30,
))

# Y boundary conditions required if n > 0
REGISTRY.register_constraint(ConstraintRule(
    rule_id="BC_Y_REQUIRED_2D",
    params=["n", "bc_y%beg", "bc_y%end"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    predicate=lambda p: p.get("n", 0) == 0 or (p.get("bc_y%beg") is not None and p.get("bc_y%end") is not None),
    message="bc_y%beg and bc_y%end must be set for 2D/3D (n > 0)",
    check_method="check_boundary_conditions",
    priority=30,
))

# Z domain required if p > 0
REGISTRY.register_constraint(ConstraintRule(
    rule_id="Z_DOMAIN_REQUIRED_3D",
    params=["p", "z_domain%beg", "z_domain%end"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.PRE_PROCESS},
    predicate=lambda p: p.get("p", 0) == 0 or (p.get("z_domain%beg") is not None and p.get("z_domain%end") is not None),
    message="z_domain%beg and z_domain%end must be set for 3D (p > 0)",
    check_method="check_simulation_domain",
    priority=30,
))

# Cylindrical coordinates constraint: p must be odd
REGISTRY.register_constraint(ConstraintRule(
    rule_id="CYL_COORD_P_ODD",
    params=["cyl_coord", "p"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.PRE_PROCESS, Stage.SIMULATION},
    predicate=lambda par: par.get("cyl_coord") not in ('T', 't', True) or par.get("p", 0) == 0 or par.get("p", 0) % 2 == 1,
    message="With cyl_coord=T, p must be 0 or odd",
    check_method="check_simulation_domain",
    priority=40,
))
