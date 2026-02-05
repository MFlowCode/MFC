"""
MFC Parameter Schema Package (Minimal).

Single source of truth for MFC's ~3,300 case parameters.

Import Order
------------
The imports below follow a specific order:
1. REGISTRY is imported first (empty at this point)
2. Schema classes (ParamDef, ParamType) for type definitions
3. definitions module is imported LAST to populate and freeze REGISTRY

The definitions import is a side-effect import that registers all ~3,300
parameters with REGISTRY and then freezes it. This must happen at package
import time so that any code importing from this package gets a fully
populated, immutable registry.

After initialization, REGISTRY.is_frozen is True and any attempt to
register new parameters will raise RegistryFrozenError.
"""

from .registry import REGISTRY, RegistryFrozenError
from .schema import ParamDef, ParamType

# IMPORTANT: This import populates REGISTRY with all parameter definitions
# and freezes it. It must come after REGISTRY is imported and must not be removed.
from . import definitions  # noqa: F401  pylint: disable=unused-import

__all__ = ['REGISTRY', 'RegistryFrozenError', 'ParamDef', 'ParamType']
