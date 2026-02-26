"""Bounded FIFO step cache shared by tui.py and interactive.py.

Keeps up to CACHE_MAX assembled timesteps in memory, evicting the oldest
entry when the cap is reached.  The module-level state is intentional:
both the TUI and the interactive server are single-instance; a per-session
cache avoids redundant disk reads while bounding peak memory usage.
"""

from typing import Callable

CACHE_MAX: int = 50
_cache: dict = {}
_cache_order: list = []


def load(step: int, read_func: Callable) -> object:
    """Return cached data for *step*, calling *read_func* on a miss."""
    if step not in _cache:
        if len(_cache) >= CACHE_MAX:
            evict = _cache_order.pop(0)
            _cache.pop(evict, None)
        _cache[step] = read_func(step)
        _cache_order.append(step)
    return _cache[step]


def seed(step: int, data: object) -> None:
    """Clear the cache and pre-populate it with already-loaded data."""
    _cache.clear()
    _cache_order.clear()
    _cache[step] = data
    _cache_order.append(step)
