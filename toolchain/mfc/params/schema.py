"""
Parameter Schema Definitions (Minimal).

Provides ParamDef for parameter metadata and ParamType enum.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Set, Any, Optional


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
    """Definition of a single MFC parameter."""
    name: str
    param_type: ParamType
    stages: Set[Stage] = field(default_factory=lambda: {Stage.COMMON})
    description: str = ""

    def __post_init__(self):
        if not isinstance(self.stages, set):
            self.stages = set(self.stages)

    @property
    def type_tag(self) -> str:
        return self.param_type.value
