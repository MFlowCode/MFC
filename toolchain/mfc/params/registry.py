"""
Parameter Registry (Minimal).

Central storage for MFC parameter definitions.
"""

from typing import Dict, Set
from collections import defaultdict

from .schema import ParamDef, Stage


class ParamRegistry:
    """Central registry for MFC parameters."""

    def __init__(self):
        self._params: Dict[str, ParamDef] = {}
        self._by_stage: Dict[Stage, Set[str]] = defaultdict(set)

    def register(self, param: ParamDef) -> None:
        """Register a parameter, merging stages if it already exists."""
        if param.name in self._params:
            existing = self._params[param.name]
            if existing.param_type != param.param_type:
                raise ValueError(f"Type mismatch for '{param.name}'")
            existing.stages.update(param.stages)
            for stage in param.stages:
                self._by_stage[stage].add(param.name)
            return

        self._params[param.name] = param
        for stage in param.stages:
            self._by_stage[stage].add(param.name)

    @property
    def all_params(self) -> Dict[str, ParamDef]:
        """Get all registered parameters."""
        return self._params

    def get_params_by_stage(self, stage: Stage) -> Dict[str, ParamDef]:
        """Get parameters for a stage (includes COMMON)."""
        result = {}
        for name in self._by_stage.get(Stage.COMMON, set()):
            result[name] = self._params[name]
        if stage != Stage.COMMON:
            for name in self._by_stage.get(stage, set()):
                result[name] = self._params[name]
        return result


# Global registry instance
REGISTRY = ParamRegistry()
