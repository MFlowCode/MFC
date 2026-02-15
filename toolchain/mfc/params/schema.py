"""
Parameter Schema Definitions.

Provides ParamDef for parameter metadata, constraints, and dependencies.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Set, Any, Optional, Dict, List

from .errors import constraint_error


class ParamType(Enum):
    """Parameter types matching MFC's Fortran types with JSON schema support."""
    INT = "int"
    REAL = "real"
    LOG = "log"
    STR = "str"
    ANALYTIC_INT = "analytic:int"
    ANALYTIC_REAL = "analytic:real"

    @property
    def json_schema(self) -> Dict[str, Any]:
        """
        Return JSON schema fragment for this parameter type.

        Used by fastjsonschema for case file validation.
        """
        schemas = {
            ParamType.INT: {"type": "integer"},
            ParamType.REAL: {"type": "number"},
            ParamType.LOG: {"enum": ["T", "F"]},
            ParamType.STR: {"type": "string"},
            # Analytic types accept either the base type or a string expression
            ParamType.ANALYTIC_INT: {"type": ["integer", "string"]},
            ParamType.ANALYTIC_REAL: {"type": ["number", "string"]},
        }
        return schemas[self]


@dataclass
class ParamDef:  # pylint: disable=too-many-instance-attributes
    """
    Definition of a single MFC parameter.

    Attributes:
        name: Parameter name
        param_type: Type (INT, REAL, LOG, STR, ANALYTIC_*)
        description: Human-readable description
        case_optimization: Can be hard-coded for GPU builds
        constraints: Validation constraints (choices, min, max)
        dependencies: Related params (requires, recommends)
        tags: Feature tags for grouping (e.g., "mhd", "bubbles", "weno")
    """
    name: str
    param_type: ParamType
    description: str = ""
    case_optimization: bool = False
    constraints: Optional[Dict[str, Any]] = None  # {"choices": [...], "min": N, "max": N}
    dependencies: Optional[Dict[str, Any]] = None  # {"requires": [...], "recommends": [...]}
    tags: Set[str] = field(default_factory=set)  # Feature tags: "mhd", "bubbles", etc.
    hint: str = ""  # Constraint/usage hint for docs (e.g. "Used with grcbc_in")
    math_symbol: str = ""  # LaTeX math symbol (Doxygen format, e.g. "\\f$\\gamma_k\\f$")

    def __post_init__(self):
        # Validate name
        if not self.name or not isinstance(self.name, str):
            raise ValueError("ParamDef name must be a non-empty string")

    @property
    def type_tag(self) -> str:
        return self.param_type.value

    def validate_value(self, value: Any) -> List[str]:
        """
        Validate a value against this parameter's constraints.

        Returns list of error messages (empty if valid).
        Uses consistent error formatting from errors.py.
        """
        errors = []
        if self.constraints is None:
            return errors

        # Check choices constraint
        if "choices" in self.constraints:
            choices = self.constraints["choices"]
            if value not in choices:
                errors.append(constraint_error(self.name, "choices", choices, value))

        # Check numeric range constraints (only for numeric values, not analytic strings)
        if isinstance(value, (int, float)):
            if "min" in self.constraints and value < self.constraints["min"]:
                errors.append(
                    constraint_error(self.name, "min", self.constraints["min"], value)
                )
            if "max" in self.constraints and value > self.constraints["max"]:
                errors.append(
                    constraint_error(self.name, "max", self.constraints["max"], value)
                )

        return errors
