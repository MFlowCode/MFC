"""
Time Stepping Parameter Definitions.

Defines time integration parameters:
- dt: Time step size
- time_stepper: Time integration method
- t_step_start, t_step_stop, t_step_save: Time step control
- cfl_target, cfl_adap_dt: Adaptive time stepping
"""

from ..schema import ParamDef, ConstraintRule, ParamType, Stage, ConstraintType
from ..registry import REGISTRY


# =============================================================================
# Time Step Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="dt",
    param_type=ParamType.REAL,
    stages={Stage.COMMON},
    description="Time step size (if None, uses adaptive time stepping)",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="time_stepper",
    param_type=ParamType.INT,
    stages={Stage.SIMULATION},
    default=3,
    description="Time integration method: 1=Euler, 2=RK2, 3=RK3, 4=RK4, 5=RK5",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="t_step_start",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    default=0,
    description="Starting time step index",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="t_step_stop",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Ending time step index",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="t_step_save",
    param_type=ParamType.INT,
    stages={Stage.COMMON},
    description="Save data every N time steps",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="t_step_print",
    param_type=ParamType.INT,
    stages={Stage.SIMULATION},
    description="Print status every N time steps",
    category="time_stepping",
))


# =============================================================================
# Adaptive Time Stepping Parameters
# =============================================================================

REGISTRY.register(ParamDef(
    name="cfl_adap_dt",
    param_type=ParamType.LOG,
    stages={Stage.SIMULATION},
    default='F',
    description="Enable adaptive time stepping based on CFL",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="cfl_target",
    param_type=ParamType.REAL,
    stages={Stage.SIMULATION, Stage.POST_PROCESS},
    default=0.5,
    description="Target CFL number for adaptive time stepping (0 < cfl <= 1)",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="cfl_max",
    param_type=ParamType.REAL,
    stages={Stage.SIMULATION},
    description="Maximum allowed CFL number",
    category="time_stepping",
))

REGISTRY.register(ParamDef(
    name="t_tol",
    param_type=ParamType.REAL,
    stages={Stage.SIMULATION},
    description="Time tolerance for adaptive stepping",
    category="time_stepping",
))


# =============================================================================
# Time Stepping Constraints
# =============================================================================

REGISTRY.register_constraint(ConstraintRule(
    rule_id="DT_POSITIVE",
    params=["dt"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.SIMULATION},
    message="dt must be positive (or None for adaptive time stepping)",
    check_method="check_time_stepping",
    priority=10,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="TIME_STEPPER_CHOICES",
    params=["time_stepper"],
    constraint_type=ConstraintType.CHOICES,
    stages={Stage.SIMULATION},
    valid_values=[1, 2, 3, 4, 5],
    message="time_stepper must be 1 (Euler), 2 (RK2), 3 (RK3), 4 (RK4), or 5 (RK5)",
    check_method="check_time_stepping",
    priority=20,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="T_STEP_START_NON_NEGATIVE",
    params=["t_step_start"],
    constraint_type=ConstraintType.NON_NEGATIVE,
    stages={Stage.COMMON},
    message="t_step_start must be non-negative",
    check_method="check_time_stepping",
    priority=30,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="T_STEP_STOP_GE_START",
    params=["t_step_start", "t_step_stop"],
    constraint_type=ConstraintType.CUSTOM,
    stages={Stage.COMMON},
    predicate=lambda p: p.get("t_step_stop") is None or p.get("t_step_start") is None or p.get("t_step_stop") >= p.get("t_step_start"),
    message="t_step_stop must be >= t_step_start",
    check_method="check_time_stepping",
    priority=40,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="T_STEP_SAVE_POSITIVE",
    params=["t_step_save"],
    constraint_type=ConstraintType.POSITIVE,
    stages={Stage.COMMON},
    message="t_step_save must be positive",
    check_method="check_time_stepping",
    priority=50,
))

REGISTRY.register_constraint(ConstraintRule(
    rule_id="CFL_TARGET_RANGE",
    params=["cfl_target"],
    constraint_type=ConstraintType.RANGE,
    stages={Stage.SIMULATION},
    min_value=0,
    min_exclusive=True,
    max_value=1,
    message="cfl_target must be in (0, 1]",
    check_method="check_time_stepping",
    priority=60,
))
