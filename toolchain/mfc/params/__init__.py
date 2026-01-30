"""
MFC Parameter Schema Package.

This package provides a unified, single-source-of-truth architecture for MFC's
~3,300 case parameters. It eliminates manual synchronization between case_dicts.py,
case_validator.py, and documentation.

Architecture:
- schema.py: ParamDef, ConstraintRule dataclasses
- registry.py: Central ParamRegistry coordinator
- definitions/: Parameter definitions by domain (core, weno, time_stepping, etc.)
- generators/: Generate case_dicts.py and validator code from registry

Usage:
    from mfc.params import REGISTRY
    from mfc.params.schema import ParamDef, ConstraintRule, ParamType, Stage

    # Query parameters
    params = REGISTRY.get_params_by_stage(Stage.SIMULATION)

    # Query constraints
    constraints = REGISTRY.get_constraints_for_stage(Stage.SIMULATION)
"""

from .registry import REGISTRY
from .schema import ParamDef, ConstraintRule, ParamType, Stage, ConstraintType

__all__ = [
    'REGISTRY',
    'ParamDef',
    'ConstraintRule',
    'ParamType',
    'Stage',
    'ConstraintType',
]
