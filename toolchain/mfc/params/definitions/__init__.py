"""
Parameter Definitions Package.

This package contains parameter definitions organized by domain:
- core: Grid dimensions, model equations, num_fluids
- weno: WENO scheme parameters
- time_stepping: dt, time_stepper, t_step_*, cfl_*
- domain: Domain bounds and boundary conditions

Each module registers its parameters and constraints with the global REGISTRY
when imported.
"""

from typing import List


def load_all_definitions() -> List[str]:
    """
    Load all parameter definition modules.

    This imports each definition module, causing their parameters
    and constraints to be registered with the global REGISTRY.

    Returns:
        List of loaded module names
    """
    loaded = []

    # Import in dependency order (core first, then domain-specific)
    from . import core
    loaded.append('core')

    from . import weno
    loaded.append('weno')

    from . import time_stepping
    loaded.append('time_stepping')

    from . import domain
    loaded.append('domain')

    from . import patches
    loaded.append('patches')

    from . import simulation
    loaded.append('simulation')

    from . import post_process
    loaded.append('post_process')

    return loaded


__all__ = ['load_all_definitions']
