"""Bounded FIFO step cache shared by tui.py and interactive.py.

Keeps up to CACHE_MAX assembled timesteps in memory, evicting the oldest
entry when the cap is reached.  The module-level state is intentional:
both the TUI and the interactive server are single-instance; a per-session
cache avoids redundant disk reads while bounding peak memory usage.

Thread safety: Dash runs Flask with threading enabled.  All mutations are
protected by _lock.

The disk read in load() is performed *outside* the lock so concurrent
callbacks loading different steps proceed in parallel.  Two threads may
both read the same step on a cold-cache miss; the second writer simply
re-inserts identical data (safe but slightly wasteful).

prefetch() submits background loads for adjacent steps into a small thread
pool.  It should only be called when loads are fast — for large grids where
each read takes multiple seconds, concurrent prefetch workers compete for
I/O bandwidth and make direct reads slower.  Callers are responsible for
gating prefetch on measured load time (see interactive.py).
"""

import atexit
import logging
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List, Optional

logger = logging.getLogger(__name__)

CACHE_MAX: int = 40
_cache: dict = {}
_cache_order: list = []
_in_flight: set = set()  # steps currently being prefetched
_lock = threading.Lock()

_prefetch_pool: Optional[ThreadPoolExecutor] = None
_prefetch_pool_lock = threading.Lock()


def _get_prefetch_pool() -> ThreadPoolExecutor:
    """Return the prefetch pool, creating it lazily on first use."""
    global _prefetch_pool  # noqa: PLW0603
    with _prefetch_pool_lock:
        if _prefetch_pool is None:
            _prefetch_pool = ThreadPoolExecutor(max_workers=3, thread_name_prefix="mfc_prefetch")
            atexit.register(_prefetch_pool.shutdown, wait=False)
        return _prefetch_pool


def load(step: int, read_func: Callable) -> object:
    """Return cached data for *step*, calling *read_func* on a miss.

    The disk read is performed outside the lock so concurrent callbacks
    for different steps do not serialize behind a single I/O operation.
    """
    with _lock:
        if step in _cache:
            return _cache[step]

    # Read outside lock — concurrent reads of the same step are safe.
    data = read_func(step)

    with _lock:
        # Re-check: another thread may have loaded this step while we read.
        if step in _cache:
            return _cache[step]
        if len(_cache) >= CACHE_MAX:
            evict = _cache_order.pop(0)
            _cache.pop(evict, None)
        _cache[step] = data
        _cache_order.append(step)

    return data


def prefetch(current_step: int, all_steps: List[int], read_func: Callable) -> None:
    """Eagerly load neighbours of *current_step* into the cache.

    Submits background tasks for the next two and previous one steps.
    No-ops if the step is already cached or a load is already in flight.

    WARNING: only call this when direct loads are fast.  For large grids,
    concurrent prefetch I/O competes with direct reads and causes slowdowns.
    Gate this call on measured load time in the caller.
    """
    try:
        idx = all_steps.index(current_step)
    except ValueError:
        return

    candidates = []
    if idx + 1 < len(all_steps):
        candidates.append(all_steps[idx + 1])
    if idx + 2 < len(all_steps):
        candidates.append(all_steps[idx + 2])
    if idx > 0:
        candidates.append(all_steps[idx - 1])

    for s in candidates:
        with _lock:
            if s in _cache or s in _in_flight:
                continue
            _in_flight.add(s)
        _get_prefetch_pool().submit(_bg_load, s, read_func)


def _bg_load(key: object, read_func: Callable) -> None:
    """Background worker: load *key* and insert into cache."""
    try:
        data = read_func(key)
        with _lock:
            if key not in _cache:
                if len(_cache) >= CACHE_MAX:
                    evict = _cache_order.pop(0)
                    _cache.pop(evict, None)
                _cache[key] = data
                _cache_order.append(key)
    except Exception:
        logger.debug("Prefetch failed for key %s", key, exc_info=True)
    finally:
        with _lock:
            _in_flight.discard(key)


def prefetch_one(key: object, read_func: Callable) -> None:
    """Submit a single background load for *key* if not cached or in-flight.

    Works with any hashable key (int step, (step, var) tuple, etc.).
    No-ops immediately if the key is already cached or loading.
    """
    with _lock:
        if key in _cache or key in _in_flight:
            return
        _in_flight.add(key)
    _get_prefetch_pool().submit(_bg_load, key, read_func)


def seed(step: int, data: object) -> None:
    """Clear the cache and pre-populate it with already-loaded data."""
    with _lock:
        _cache.clear()
        _cache_order.clear()
        _in_flight.clear()
        _cache[step] = data
        _cache_order.append(step)


def clear() -> None:
    """Reset the cache to empty (useful for test teardown)."""
    with _lock:
        _cache.clear()
        _cache_order.clear()
        _in_flight.clear()
