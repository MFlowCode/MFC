"""
MFC Parameter Schema Package (Minimal).

Single source of truth for MFC's ~3,300 case parameters.
"""

from .registry import REGISTRY
from .schema import ParamDef, ParamType, Stage

__all__ = ['REGISTRY', 'ParamDef', 'ParamType', 'Stage']
