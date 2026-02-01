"""
Parameter Registry.

Central storage for MFC parameter definitions. This module provides the
ParamRegistry class which serves as the single source of truth for all
~3,300 MFC parameters.

Usage
-----
The global REGISTRY instance is populated by importing the definitions module.
Once populated, parameters can be queried by name or by stage:

    from mfc.params import REGISTRY

    # Get a specific parameter
    param = REGISTRY.all_params.get('m')

    # Get all parameters for simulation stage
    sim_params = REGISTRY.get_params_by_stage(Stage.SIMULATION)

Thread Safety
-------------
The registry is populated once at import time and frozen (made immutable).
After freezing, it is safe to read from multiple threads. Attempts to
register new parameters after freezing will raise RuntimeError.
"""

from typing import Dict, Set, Mapping
from types import MappingProxyType
from collections import defaultdict

from .schema import ParamDef, Stage


class RegistryFrozenError(RuntimeError):
    """Raised when attempting to modify a frozen registry."""


class ParamRegistry:
    """
    Central registry for MFC parameters.

    This class stores parameter definitions and provides lookup methods
    for retrieving parameters by name or by execution stage.

    The registry can be frozen after initialization to prevent further
    modifications, ensuring thread-safety for read operations.

    Attributes:
        _params: Dictionary mapping parameter names to ParamDef instances.
        _by_stage: Dictionary mapping stages to sets of parameter names.
        _frozen: Whether the registry has been frozen (immutable).
    """

    def __init__(self):
        """Initialize an empty registry."""
        self._params: Dict[str, ParamDef] = {}
        self._by_stage: Dict[Stage, Set[str]] = defaultdict(set)
        self._frozen: bool = False
        self._params_proxy: Mapping[str, ParamDef] = None

    def freeze(self) -> None:
        """
        Freeze the registry, preventing further modifications.

        After calling this method:
        - register() will raise RegistryFrozenError
        - all_params returns a read-only view (MappingProxyType)

        This method is idempotent (safe to call multiple times).
        """
        if not self._frozen:
            self._frozen = True
            self._params_proxy = MappingProxyType(self._params)

    @property
    def is_frozen(self) -> bool:
        """Return True if the registry has been frozen."""
        return self._frozen

    def register(self, param: ParamDef) -> None:
        """
        Register a parameter definition.

        If a parameter with the same name already exists, the stages are
        merged (i.e., the new stages are added to the existing ones).
        This allows parameters to be defined incrementally across multiple
        definition files.

        Args:
            param: The parameter definition to register.

        Raises:
            RegistryFrozenError: If the registry has been frozen.
            ValueError: If a parameter with the same name exists but has
                a different type (type mismatch is not allowed).
        """
        if self._frozen:
            raise RegistryFrozenError(
                f"Cannot register '{param.name}': registry is frozen. "
                "All parameters must be registered during module initialization."
            )

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
    def all_params(self) -> Mapping[str, ParamDef]:
        """
        Get all registered parameters.

        Returns:
            Mapping of parameter names to their definitions.
            If the registry is frozen, returns a read-only view.
            If not frozen, returns the internal dict (mutable).
        """
        if self._frozen and self._params_proxy is not None:
            return self._params_proxy
        return self._params

    def get_params_by_stage(self, stage: Stage) -> Dict[str, ParamDef]:
        """
        Get parameters applicable to a specific execution stage.

        Parameters in Stage.COMMON are included for all stages.

        Args:
            stage: The execution stage (PRE_PROCESS, SIMULATION, or POST_PROCESS).

        Returns:
            Dictionary mapping parameter names to their definitions,
            including both stage-specific and COMMON parameters.
        """
        result = {}
        for name in self._by_stage.get(Stage.COMMON, set()):
            result[name] = self._params[name]
        if stage != Stage.COMMON:
            for name in self._by_stage.get(stage, set()):
                result[name] = self._params[name]
        return result


# Global registry instance - populated when definitions module is imported
REGISTRY = ParamRegistry()
