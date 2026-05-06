"""
Parameter Registry.

Central storage for MFC parameter definitions. This module provides the
ParamRegistry class which serves as the single source of truth for all
MFC parameters.

Usage
-----
The global REGISTRY instance is populated by importing the definitions module.
Once populated, parameters can be queried by name or by tag:

    from mfc.params import REGISTRY

    # Get a specific parameter
    param = REGISTRY.get_param_def('m')

    # Check if a parameter name is valid (including indexed families)
    REGISTRY.is_known_param('patch_ib(500)%geometry')  # True

    # Get parameters by feature tag
    mhd_params = REGISTRY.get_params_by_tag('mhd')

Thread Safety
-------------
The registry is populated once at import time and frozen (made immutable).
After freezing, it is safe to read from multiple threads. Attempts to
register new parameters after freezing will raise RuntimeError.
"""

import re
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from types import MappingProxyType
from typing import Any, Dict, Iterator, Optional, Set, Tuple

from .schema import ParamDef, ParamType


class RegistryFrozenError(RuntimeError):
    """Raised when attempting to modify a frozen registry."""


# Regex for parsing indexed family parameter names:
#   patch_ib(123)%vel(1)  -> base="patch_ib", index=123, attr="vel(1)"
#   patch_ib(1)%geometry  -> base="patch_ib", index=1,   attr="geometry"
_INDEXED_RE = re.compile(r"^([a-zA-Z_]\w*)\((\d+)\)%(.+)$")


def _resolve_family(name: str, families: Dict[str, "IndexedFamily"]) -> Optional[Tuple[ParamType, Set[str]]]:
    """
    Resolve a parameter name against indexed families.

    Returns (ParamType, tags) if the name matches a registered family
    attribute, or None otherwise.  This is the single implementation of
    family pattern-matching used by both _FamilyAwareMapping and
    ParamRegistry.
    """
    m = _INDEXED_RE.match(name)
    if m is None:
        return None
    base, idx_str, attr = m.groups()
    fam = families.get(base)
    if fam is None:
        return None
    idx = int(idx_str)
    if idx < 1:
        return None
    if fam.max_index is not None and idx > fam.max_index:
        return None
    entry = fam.attrs.get(attr)
    return entry if entry is not None else None


@dataclass(frozen=True)
class IndexedFamily:
    """
    Template for an indexed parameter family like patch_ib(N)%attr.

    Instead of registering every index individually (e.g., patch_ib(1)%geometry
    through patch_ib(N)%geometry for all attributes), we store one template and
    validate parameter names via pattern matching.

    Attributes:
        base_name: Family prefix (e.g., "patch_ib")
        attrs: Mapping of attribute name to (ParamType, tags) — attribute names
               may include sub-indices like "vel(1)", "angles(3)".
        tags: Metadata-only tags for this family (not used in resolution;
              per-attribute tags in ``attrs`` are what get returned).
        max_index: Upper bound on the index (1-based). None = unlimited.
    """

    base_name: str
    attrs: Dict[str, Tuple[ParamType, Set[str]]] = field(default_factory=dict)
    tags: Set[str] = field(default_factory=set)
    max_index: Optional[int] = None

    def __post_init__(self):
        if not self.base_name or not re.match(r"^[a-zA-Z_]\w*$", self.base_name):
            raise ValueError(f"Invalid base_name: {self.base_name!r}")
        if self.max_index is not None and self.max_index < 1:
            raise ValueError(f"max_index must be >= 1 or None, got {self.max_index}")


class _FamilyAwareMapping(Mapping):
    """
    Read-only mapping that combines scalar params with indexed family lookups.

    For containment checks and item access, indexed family params like
    ``patch_ib(500)%geometry`` are resolved via pattern matching against
    registered families — no enumeration needed.

    For iteration (items/keys/values/len), only scalar params and one
    representative example per family attribute (index=1) are yielded.
    This keeps iteration bounded regardless of max_index.

    keys(), items(), values() are inherited from collections.abc.Mapping
    and return proper KeysView/ItemsView/ValuesView objects.
    """

    __slots__ = ("_scalars", "_families", "_examples")

    def __init__(
        self,
        scalars: Dict[str, ParamDef],
        families: Dict[str, IndexedFamily],
    ):
        self._scalars = scalars
        self._families = families
        # Pre-build one example per family attr for iteration/docs
        self._examples: Dict[str, ParamDef] = {}
        for fam in families.values():
            for attr_name, (ptype, tags) in fam.attrs.items():
                key = f"{fam.base_name}(1)%{attr_name}"
                self._examples[key] = ParamDef(name=key, param_type=ptype, tags=set(tags))

    def _make_param_def(self, name: str) -> Optional[ParamDef]:
        """Build a ParamDef from a family match, or return None."""
        result = _resolve_family(name, self._families)
        if result is None:
            return None
        ptype, tags = result
        return ParamDef(name=name, param_type=ptype, tags=set(tags))

    def __getitem__(self, key: str) -> ParamDef:
        try:
            return self._scalars[key]
        except KeyError:
            pass
        result = self._make_param_def(key)
        if result is not None:
            return result
        raise KeyError(key)

    def __contains__(self, key: object) -> bool:
        # Note: this intentionally deviates from the standard Mapping contract.
        # `key in self` may be True for family params (e.g., patch_ib(500)%geometry)
        # that do NOT appear in iter(self) (which only yields index=1 examples).
        if key in self._scalars:
            return True
        if isinstance(key, str):
            return _resolve_family(key, self._families) is not None
        return False

    def __iter__(self) -> Iterator[str]:
        yield from self._scalars
        yield from self._examples

    def __len__(self) -> int:
        return len(self._scalars) + len(self._examples)


class ParamRegistry:
    """
    Central registry for MFC parameters.

    Supports two kinds of parameters:
    1. Scalar/small-indexed params — stored individually in _params.
    2. Indexed families — stored as templates in _families, matched by pattern.

    The registry can be frozen after initialization to prevent further
    modifications, ensuring thread-safety for read operations.
    """

    def __init__(self):
        """Initialize an empty registry."""
        self._params: Dict[str, ParamDef] = {}
        self._families: Dict[str, IndexedFamily] = {}
        self._by_tag: Dict[str, Set[str]] = defaultdict(set)
        self._frozen: bool = False
        self._all_params_view: Optional[_FamilyAwareMapping] = None

    def freeze(self) -> None:
        """
        Freeze the registry, preventing further modifications.

        After calling this method:
        - register() and register_family() will raise RegistryFrozenError
        - all_params returns a family-aware read-only mapping

        This method is idempotent (safe to call multiple times).
        """
        if not self._frozen:
            self._frozen = True
            self._all_params_view = _FamilyAwareMapping(self._params, self._families)

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
            raise RegistryFrozenError(f"Cannot register '{param.name}': registry is frozen. All parameters must be registered during module initialization.")

        if param.name in self._params:
            existing = self._params[param.name]
            if existing.param_type != param.param_type:
                raise ValueError(
                    f"Type mismatch for '{param.name}': "
                    f"existing type is {existing.param_type!r}, "
                    f"new type is {param.param_type!r}"
                )
            existing.tags.update(param.tags)
            for tag in param.tags:
                self._by_tag[tag].add(param.name)
            return

        self._params[param.name] = param
        for tag in param.tags:
            self._by_tag[tag].add(param.name)

    def register_family(self, family: IndexedFamily) -> None:
        """
        Register an indexed parameter family.

        Instead of registering N*attrs individual params, this stores a
        single template that is matched by pattern. This makes validation
        O(1) per parameter regardless of max_index.

        Args:
            family: The indexed family definition.

        Raises:
            RegistryFrozenError: If the registry has been frozen.
        """
        if self._frozen:
            raise RegistryFrozenError(f"Cannot register family '{family.base_name}': registry is frozen.")
        self._families[family.base_name] = family
        # Register tags for the family (using example names)
        for attr_name, (_, tags) in family.attrs.items():
            example = f"{family.base_name}(1)%{attr_name}"
            for tag in tags:
                self._by_tag[tag].add(example)

    @property
    def families(self) -> Mapping[str, IndexedFamily]:
        """Get all registered indexed families (read-only view)."""
        return MappingProxyType(self._families)

    @property
    def all_params(self) -> Mapping[str, ParamDef]:
        """
        Get all registered parameters as a mapping.

        Returns a family-aware mapping that supports:
        - Containment: ``'patch_ib(500)%geometry' in registry.all_params``
        - Lookup: ``registry.all_params.get('patch_ib(500)%geometry')``
        - Iteration: yields scalar params + one example per family attr

        If the registry is frozen, returns an immutable view.
        If not frozen and no families are registered, returns the internal dict.
        If not frozen but families exist, raises RuntimeError (the plain dict
        cannot resolve family params — call freeze() first).
        """
        if self._frozen and self._all_params_view is not None:
            return self._all_params_view
        if self._families:
            raise RuntimeError("Cannot access all_params before freeze() when indexed families are registered. Call freeze() first.")
        return self._params

    def is_known_param(self, name: str) -> bool:
        """
        Check if a parameter name is valid (scalar or indexed family).

        This is the fast path for validating user-provided parameter names.
        O(1) for both scalar params and indexed family params.
        """
        if name in self._params:
            return True
        return _resolve_family(name, self._families) is not None

    def get_param_def(self, name: str) -> Optional[ParamDef]:
        """
        Get the ParamDef for a parameter name, resolving indexed families.

        Returns None if the name is not recognized.
        """
        if name in self._params:
            return self._params[name]
        result = _resolve_family(name, self._families)
        if result is None:
            return None
        ptype, tags = result
        return ParamDef(name=name, param_type=ptype, tags=set(tags))

    def get_params_by_tag(self, tag: str) -> Dict[str, ParamDef]:
        """
        Get parameters with a specific feature tag.

        Args:
            tag: The feature tag (e.g., "mhd", "bubbles", "weno").

        Returns:
            Dictionary mapping parameter names to their definitions.
        """
        result = {}
        for name in self._by_tag.get(tag, set()):
            if name in self._params:
                result[name] = self._params[name]
            else:
                # Family example — resolve it
                param_def = self.get_param_def(name)
                if param_def is not None:
                    result[name] = param_def
        return result

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

        Indexed parameter families (e.g., patch_ib(N)%radius) are
        represented as patternProperties regexes.

        Returns:
            JSON schema dict compatible with fastjsonschema.
        """
        properties = {}
        pattern_props = {}

        # Scalar and small-indexed params
        for name, param in self._params.items():
            if "(" not in name:
                properties[name] = param.param_type.json_schema
            else:
                # Small indexed param — collapse into pattern
                pattern = re.sub(r"\(\d+\)", "__IDX__", name)
                pattern = re.escape(pattern).replace("__IDX__", r"\(\d+\)")
                pattern = f"^{pattern}$"
                if pattern not in pattern_props:
                    pattern_props[pattern] = param.param_type.json_schema

        # Indexed families — generate one pattern per attribute
        for fam in self._families.values():
            base_esc = re.escape(fam.base_name)
            for attr_name, (ptype, _tags) in fam.attrs.items():
                # Escape the attr name but replace sub-indices with \(\d+\)
                attr_pattern = re.sub(r"\(\d+\)", "__IDX__", attr_name)
                attr_pattern = re.escape(attr_pattern).replace("__IDX__", r"\(\d+\)")
                pattern = f"^{base_esc}\\([1-9]\\d*\\)%{attr_pattern}$"
                if pattern not in pattern_props:
                    pattern_props[pattern] = ptype.json_schema

        return {
            "type": "object",
            "properties": properties,
            "patternProperties": pattern_props,
            "additionalProperties": False,
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
def _get_cached_validator(registry_id: int):
    """Cache the validator at module level (registry is immutable after freeze).

    Note: registry_id is used as cache key to invalidate when registry changes.
    """
    import fastjsonschema

    return fastjsonschema.compile(REGISTRY.get_json_schema())


# Global registry instance - populated when definitions module is imported
REGISTRY = ParamRegistry()
