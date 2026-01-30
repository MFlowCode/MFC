"""
Parameter Schema Definitions.

Provides ParamDef for parameter metadata, constraints, and dependencies.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Set, Any, Optional, Dict, List


class ParamType(Enum):
    """Parameter types matching MFC's Fortran types."""
    INT = "int"
    REAL = "real"
    LOG = "log"
    STR = "str"
    ANALYTIC_INT = "analytic:int"
    ANALYTIC_REAL = "analytic:real"


class Stage(Enum):
    """MFC execution stages."""
    COMMON = auto()
    PRE_PROCESS = auto()
    SIMULATION = auto()
    POST_PROCESS = auto()


@dataclass
class ParamDef:
    """
    Definition of a single MFC parameter.

    Attributes:
        name: Parameter name
        param_type: Type (INT, REAL, LOG, STR, ANALYTIC_*)
        stages: Which stages this param applies to
        description: Human-readable description
        case_optimization: Can be hard-coded for GPU builds
        constraints: Validation constraints (choices, min, max)
        dependencies: Related params (requires, recommends)
    """
    name: str
    param_type: ParamType
    stages: Set[Stage] = field(default_factory=lambda: {Stage.COMMON})
    description: str = ""
    case_optimization: bool = False
    constraints: Optional[Dict[str, Any]] = None  # {"choices": [...], "min": N, "max": N}
    dependencies: Optional[Dict[str, Any]] = None  # {"requires": [...], "recommends": [...]}

    def __post_init__(self):
        if not isinstance(self.stages, set):
            self.stages = set(self.stages)

    @property
    def type_tag(self) -> str:
        return self.param_type.value

    def validate_value(self, value: Any) -> List[str]:
        """
        Validate a value against this parameter's constraints.

        Returns list of error messages (empty if valid).
        """
        errors = []
        if self.constraints is None:
            return errors

        # Check choices constraint
        if "choices" in self.constraints:
            choices = self.constraints["choices"]
            if value not in choices:
                errors.append(
                    f"{self.name} must be one of {choices}, got {value}"
                )

        # Check numeric range constraints
        if "min" in self.constraints and value < self.constraints["min"]:
            errors.append(
                f"{self.name} must be >= {self.constraints['min']}, got {value}"
            )
        if "max" in self.constraints and value > self.constraints["max"]:
            errors.append(
                f"{self.name} must be <= {self.constraints['max']}, got {value}"
            )

        return errors
