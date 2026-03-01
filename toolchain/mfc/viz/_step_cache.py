"""Bounded FIFO step cache shared by tui.py and interactive.py.

Keeps up to CACHE_MAX assembled timesteps in memory, evicting the oldest
entry when the cap is reached.  The module-level state is intentional:
both the TUI and the interactive server are single-instance; a per-session
cache avoids redundant disk reads while bounding peak memory usage.

Thread safety: Dash runs Flask with threading enabled.  All mutations are
protected by _lock so concurrent callbacks from multiple browser tabs cannot
corrupt _cache_order or exceed CACHE_MAX.
"""

import threading
from typing import Callable

CACHE_MAX: int = 50
_cache: dict = {}
_cache_order: list = []
_lock = threading.Lock()


def load(step: int, read_func: Callable) -> object:
    """Return cached data for *step*, calling *read_func* on a miss.

    read_func is called *before* eviction so that a failed read (e.g. a
    missing or corrupt file) does not discard a valid cache entry.
    """
    with _lock:
        if step in _cache:
            return _cache[step]
        # Read outside-the-lock would allow concurrent loads of the same
        # step; keeping it inside is simpler and safe since read_func is
        # plain file I/O that never calls load() recursively.
        data = read_func(step)
        # Evict only after a successful read.
        if len(_cache) >= CACHE_MAX:
            evict = _cache_order.pop(0)
            _cache.pop(evict, None)
        _cache[step] = data
        _cache_order.append(step)
        return data


def seed(step: int, data: object) -> None:
    """Clear the cache and pre-populate it with already-loaded data."""
    with _lock:
        _cache.clear()
        _cache_order.clear()
        _cache[step] = data
        _cache_order.append(step)


def clear() -> None:
    """Reset the cache to empty (useful for test teardown)."""
    with _lock:
        _cache.clear()
        _cache_order.clear()
