"""
Parameter Registry.

Central storage for MFC parameter definitions. This module provides the
ParamRegistry class which serves as the single source of truth for all
~3,300 MFC parameters.

Usage
-----
The global REGISTRY instance is populated by importing the definitions module.
Once populated, parameters can be queried by name or by tag:

    from mfc.params import REGISTRY

    # Get a specific parameter
    param = REGISTRY.all_params.get('m')

    # Get parameters by feature tag
    mhd_params = REGISTRY.get_params_by_tag('mhd')

Thread Safety
-------------
The registry is populated once at import time and frozen (made immutable).
After freezing, it is safe to read from multiple threads. Attempts to
register new parameters after freezing will raise RuntimeError.
"""

from typing import Dict, Set, Mapping, Any
from types import MappingProxyType
from collections import defaultdict
from functools import lru_cache

from .schema import ParamDef


class RegistryFrozenError(RuntimeError):
    """Raised when attempting to modify a frozen registry."""


class ParamRegistry:
    """
    Central registry for MFC parameters.

    This class stores parameter definitions and provides lookup methods
    for retrieving parameters by name or by feature tag.

    The registry can be frozen after initialization to prevent further
    modifications, ensuring thread-safety for read operations.

    Attributes:
        _params: Dictionary mapping parameter names to ParamDef instances.
        _by_tag: Dictionary mapping tags to sets of parameter names.
        _frozen: Whether the registry has been frozen (immutable).
    """

    def __init__(self):
        """Initialize an empty registry."""
        self._params: Dict[str, ParamDef] = {}
        self._by_tag: Dict[str, Set[str]] = defaultdict(set)
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

        If a parameter with the same name already exists, the tags are
        merged. This allows parameters to be defined incrementally across
        multiple definition files.

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
            existing.tags.update(param.tags)
            for tag in param.tags:
                self._by_tag[tag].add(param.name)
            return

        self._params[param.name] = param
        for tag in param.tags:
            self._by_tag[tag].add(param.name)

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

    def get_params_by_tag(self, tag: str) -> Dict[str, ParamDef]:
        """
        Get parameters with a specific feature tag.

        Args:
            tag: The feature tag (e.g., "mhd", "bubbles", "weno").

        Returns:
            Dictionary mapping parameter names to their definitions.
        """
        return {name: self._params[name] for name in self._by_tag.get(tag, set())}

    def get_all_tags(self) -> Set[str]:
        """
        Get all feature tags used in the registry.

        Returns:
            Set of all tag names.
        """
        return set(self._by_tag.keys())

    def get_json_schema(self) -> Dict[str, Any]:
        """
        Generate JSON schema for case file validation.

        Returns:
            JSON schema dict compatible with fastjsonschema.
        """
        properties = {
            name: param.param_type.json_schema
            for name, param in self.all_params.items()
        }

        return {
            "type": "object",
            "properties": properties,
            "additionalProperties": False
        }

    def get_validator(self):
        """
        Get a cached JSON schema validator for all parameters.

        Returns:
            Compiled fastjsonschema validator function.
        """
        # Use module-level cache since registry is frozen
        return _get_cached_validator(id(self))


@lru_cache(maxsize=1)
def _get_cached_validator(registry_id: int):  # pylint: disable=unused-argument
    """Cache the validator at module level (registry is immutable after freeze).

    Note: registry_id is used as cache key to invalidate when registry changes.
    """
    import fastjsonschema  # pylint: disable=import-outside-toplevel
    return fastjsonschema.compile(REGISTRY.get_json_schema())


# Global registry instance - populated when definitions module is imported
REGISTRY = ParamRegistry()
