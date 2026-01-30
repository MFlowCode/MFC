"""
Parameter Schema Definitions.

This module defines the core dataclasses for MFC's parameter schema:
- ParamDef: Defines a single parameter (name, type, stages, default, etc.)
- ConstraintRule: Defines a validation constraint (required, range, dependency, etc.)

These dataclasses serve as the single source of truth for parameter metadata,
replacing the scattered definitions across case_dicts.py and case_validator.py.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Optional, Any, Set, Callable, Dict, Union


class ParamType(Enum):
    """
    Parameter types matching MFC's Fortran types.

    These correspond to the type tags used in case_dicts.py:
    - INT: Integer values (Fortran INTEGER)
    - REAL: Floating point values (Fortran REAL)
    - LOG: Logical/boolean values (Fortran LOGICAL, uses 'T'/'F' strings)
    - STR: String values (Fortran CHARACTER)
    - ANALYTIC_INT: Integer expressions evaluated at runtime (e.g., "2*m")
    - ANALYTIC_REAL: Real expressions evaluated at runtime (e.g., "0.5*sin(x)")
    """
    INT = "int"
    REAL = "real"
    LOG = "log"
    STR = "str"
    ANALYTIC_INT = "analytic:int"
    ANALYTIC_REAL = "analytic:real"


class Stage(Enum):
    """
    MFC execution stages where parameters apply.

    Parameters may be valid in one or more stages:
    - COMMON: Used by all stages (pre_process, simulation, post_process)
    - PRE_PROCESS: Only used during mesh/IC generation
    - SIMULATION: Only used during the actual simulation
    - POST_PROCESS: Only used during post-processing/visualization
    """
    COMMON = auto()
    PRE_PROCESS = auto()
    SIMULATION = auto()
    POST_PROCESS = auto()


class ConstraintType(Enum):
    """
    Types of validation constraints.

    Declarative constraints (~80%) can be auto-generated:
    - REQUIRED: Parameter must be set (not None)
    - RANGE: Value must be within min/max bounds
    - CHOICES: Value must be one of a discrete set
    - POSITIVE: Value must be > 0
    - NON_NEGATIVE: Value must be >= 0
    - DEPENDENCY: If condition A, then B must be true
    - MUTEX: Only one of a set of parameters can be true

    Custom constraints (~20%) need hand-coded predicates:
    - CUSTOM: Complex multi-parameter logic with a callable predicate
    """
    REQUIRED = auto()
    RANGE = auto()
    CHOICES = auto()
    POSITIVE = auto()
    NON_NEGATIVE = auto()
    DEPENDENCY = auto()
    MUTEX = auto()
    CUSTOM = auto()


@dataclass
class ParamDef:
    """
    Definition of a single MFC parameter.

    This captures all metadata needed to:
    1. Generate type schemas for case_dicts.py
    2. Generate documentation
    3. Enable IDE autocompletion

    Attributes:
        name: Parameter name as used in case files (e.g., "m", "weno_order")
        param_type: Data type (INT, REAL, LOG, STR, ANALYTIC_*)
        stages: Set of stages where this parameter is valid
        default: Default value (None means required)
        description: Human-readable description for docs
        category: Grouping category (e.g., "grid", "weno", "physics")
        case_optimization: If True, can be hard-coded for GPU optimization
        fortran_name: Actual Fortran variable name if different from name
    """
    name: str
    param_type: ParamType
    stages: Set[Stage] = field(default_factory=lambda: {Stage.COMMON})
    default: Any = None
    description: str = ""
    category: Optional[str] = None
    case_optimization: bool = False
    fortran_name: Optional[str] = None

    def __post_init__(self):
        # Ensure stages is a set
        if not isinstance(self.stages, set):
            self.stages = set(self.stages)

    @property
    def is_required(self) -> bool:
        """Parameter is required if it has no default."""
        return self.default is None

    @property
    def type_tag(self) -> str:
        """Return the type tag used in case_dicts.py schema."""
        return self.param_type.value


@dataclass
class ConstraintRule:
    """
    Definition of a validation constraint.

    Constraints enforce rules on parameter values and combinations.
    They can be declarative (auto-generated) or custom (hand-coded).

    Attributes:
        rule_id: Unique identifier for this constraint (e.g., "M_POSITIVE")
        params: List of parameter names involved in this constraint
        constraint_type: Type of constraint (REQUIRED, RANGE, etc.)
        stages: Stages where this constraint applies
        message: Error message to show when constraint is violated

        # For RANGE constraints:
        min_value: Minimum allowed value (inclusive)
        max_value: Maximum allowed value (inclusive)
        min_exclusive: If True, min_value is exclusive (> instead of >=)
        max_exclusive: If True, max_value is exclusive (< instead of <=)

        # For CHOICES constraints:
        valid_values: List of allowed values

        # For DEPENDENCY constraints:
        if_param: Parameter to check condition on
        if_condition: Condition on if_param (e.g., "> 0", "== 'T'", "is not None")
        then_param: Parameter that must satisfy then_condition
        then_condition: Required condition on then_param

        # For MUTEX constraints:
        mutex_params: List of parameters where only one can be True/'T'

        # For CUSTOM constraints:
        predicate: Callable[[dict], bool] that returns True if valid
        check_method: Name of the check method this belongs to (for organization)

        priority: Lower values run first (default 100)
    """
    rule_id: str
    params: List[str]
    constraint_type: ConstraintType
    stages: Set[Stage] = field(default_factory=lambda: {Stage.COMMON})
    message: str = ""

    # RANGE constraint fields
    min_value: Optional[Union[int, float]] = None
    max_value: Optional[Union[int, float]] = None
    min_exclusive: bool = False
    max_exclusive: bool = False

    # CHOICES constraint fields
    valid_values: Optional[List[Any]] = None

    # DEPENDENCY constraint fields
    if_param: Optional[str] = None
    if_condition: Optional[str] = None
    then_param: Optional[str] = None
    then_condition: Optional[str] = None

    # MUTEX constraint fields
    mutex_params: Optional[List[str]] = None

    # CUSTOM constraint fields
    predicate: Optional[Callable[[Dict[str, Any]], bool]] = None
    check_method: Optional[str] = None

    priority: int = 100

    def __post_init__(self):
        # Ensure stages is a set
        if not isinstance(self.stages, set):
            self.stages = set(self.stages)

        # Validate constraint-specific fields
        if self.constraint_type == ConstraintType.RANGE:
            if self.min_value is None and self.max_value is None:
                raise ValueError(f"RANGE constraint {self.rule_id} needs min_value or max_value")

        elif self.constraint_type == ConstraintType.CHOICES:
            if not self.valid_values:
                raise ValueError(f"CHOICES constraint {self.rule_id} needs valid_values")

        elif self.constraint_type == ConstraintType.DEPENDENCY:
            if not all([self.if_param, self.if_condition, self.then_param, self.then_condition]):
                raise ValueError(f"DEPENDENCY constraint {self.rule_id} needs if_param, "
                                 "if_condition, then_param, then_condition")

        elif self.constraint_type == ConstraintType.MUTEX:
            if not self.mutex_params or len(self.mutex_params) < 2:
                raise ValueError(f"MUTEX constraint {self.rule_id} needs at least 2 mutex_params")

        elif self.constraint_type == ConstraintType.CUSTOM:
            if self.predicate is None:
                raise ValueError(f"CUSTOM constraint {self.rule_id} needs a predicate function")

    def applies_to_stage(self, stage: Stage) -> bool:
        """Check if this constraint applies to the given stage."""
        return Stage.COMMON in self.stages or stage in self.stages


@dataclass
class PatternDef:
    """
    Definition of a parameterized parameter pattern.

    MFC uses many indexed parameters like patch_icpp(1)%pres, fluid_pp(2)%gamma.
    This class defines the pattern for generating these parameters.

    Attributes:
        pattern: Pattern string with {i}, {j}, etc. placeholders
        param_type: Data type for generated parameters
        stages: Stages where generated parameters are valid
        indices: Dict mapping placeholder to (start, end) range or list of values
        description: Description template (can use {i}, {j} placeholders)
        category: Category for generated parameters
    """
    pattern: str
    param_type: ParamType
    stages: Set[Stage] = field(default_factory=lambda: {Stage.COMMON})
    indices: Dict[str, Union[tuple, list]] = field(default_factory=dict)
    description: str = ""
    category: Optional[str] = None
    default: Any = None

    def expand(self, index_limits: Dict[str, int]) -> List[ParamDef]:
        """
        Expand the pattern into concrete ParamDef instances.

        Args:
            index_limits: Dict mapping index name to max value
                          e.g., {"i": 10, "j": 3} for patch_icpp(i)%alpha(j)

        Returns:
            List of ParamDef instances for each combination
        """
        params = []

        # Build ranges for each index
        ranges = {}
        for idx_name, spec in self.indices.items():
            if isinstance(spec, list):
                ranges[idx_name] = spec
            else:
                start, end_key = spec
                end = index_limits.get(end_key, 1)
                ranges[idx_name] = list(range(start, end + 1))

        # Generate all combinations
        from itertools import product
        idx_names = list(ranges.keys())
        if idx_names:
            for combo in product(*[ranges[n] for n in idx_names]):
                substitutions = dict(zip(idx_names, combo))
                name = self.pattern
                desc = self.description
                for idx_name, val in substitutions.items():
                    name = name.replace(f"{{{idx_name}}}", str(val))
                    desc = desc.replace(f"{{{idx_name}}}", str(val))

                params.append(ParamDef(
                    name=name,
                    param_type=self.param_type,
                    stages=self.stages,
                    default=self.default,
                    description=desc,
                    category=self.category,
                ))
        else:
            # No indices - just a single parameter
            params.append(ParamDef(
                name=self.pattern,
                param_type=self.param_type,
                stages=self.stages,
                default=self.default,
                description=self.description,
                category=self.category,
            ))

        return params
