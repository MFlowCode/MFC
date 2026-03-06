"""
Silo-HDF5 reader for MFC post-processed output.

Silo files produced by MFC are valid HDF5 underneath.  Each Silo object
is stored as an HDF5 Named Datatype whose ``silo`` compound attribute
carries the metadata (mesh name, data-array path, dimensions, etc.).
Actual data lives in numbered datasets under the ``.silo/`` group.

This reader uses h5py to navigate that structure.

Performance note
----------------
The HDF5 file structure (which internal paths correspond to which
variables) is **identical for every timestep of the same case**.  On
the first read for each processor rank directory, ``_get_structure``
parses the metadata and caches a ``_SiloStructure``.  All subsequent
reads for that rank skip attribute iteration entirely and go directly
to the HDF5 dataset paths, reducing per-step overhead from ~0.5–3 s to
near-raw-I/O speed (~0.15 s).
"""

import atexit
import os
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

from .reader import AssembledData, ProcessorData, assemble_from_proc_data

# Silo type constants (from silo.h)
_DB_QUADMESH = 130
_DB_QUADVAR = 501


# ---------------------------------------------------------------------------
# File-structure cache
# ---------------------------------------------------------------------------

@dataclass
class _SiloStructure:
    """Cached HDF5 layout for one processor's silo files.

    All timestep files in the same rank directory have the same internal
    path assignments, so we only need to parse this once.
    """
    ndims: int
    coord_paths: List[bytes]          # HDF5 paths for x, y[, z] node coords
    var_paths: Dict[str, bytes] = field(default_factory=dict)  # varname → data path


_struct_cache: Dict[str, _SiloStructure] = {}   # key = rank directory path
_struct_lock = threading.Lock()


def _parse_structure(f) -> _SiloStructure:
    """Walk the HDF5 file once and extract the silo layout.

    This is the expensive operation (~0.5–3 s per rank).  The result is
    cached in ``_struct_cache`` so it is only paid once per rank.
    """
    ndims = 0
    coord_paths: List[bytes] = []
    var_paths: Dict[str, bytes] = {}

    for key, obj in f.items():
        if key in ("..", ".silo"):
            continue
        if not isinstance(obj, h5py.Datatype):
            continue
        silo_type = obj.attrs.get("silo_type")
        if silo_type is None:
            continue
        silo_type_int = int(silo_type)

        if silo_type_int == _DB_QUADMESH and "silo" in obj.attrs:
            mesh_attr = obj.attrs["silo"]
            ndims = int(mesh_attr["ndims"])
            for i in range(ndims):
                raw = mesh_attr[f"coord{i}"]
                coord_paths.append(bytes(raw))

        elif silo_type_int == _DB_QUADVAR and "silo" in obj.attrs:
            attr = obj.attrs["silo"]
            try:
                raw = attr["value0"]
                var_paths[key] = bytes(raw)
            except (KeyError, ValueError):
                warnings.warn(
                    f"Variable '{key}' missing 'value0' in silo attr, skipping",
                    stacklevel=4,
                )

    if not coord_paths:
        raise ValueError("No rectilinear mesh found in file")

    return _SiloStructure(ndims=ndims, coord_paths=coord_paths, var_paths=var_paths)


def _get_structure(rank_dir: str, example_path: str) -> _SiloStructure:
    """Return the cached structure for *rank_dir*, parsing *example_path* on a miss."""
    with _struct_lock:
        if rank_dir in _struct_cache:
            return _struct_cache[rank_dir]

    # Parse outside lock — two threads may both parse on a cold miss, which is
    # safe (identical result) and far cheaper than holding the lock during I/O.
    with h5py.File(example_path, "r") as f:
        struct = _parse_structure(f)

    with _struct_lock:
        # Re-check in case another thread beat us here.
        if rank_dir not in _struct_cache:
            _struct_cache[rank_dir] = struct
        return _struct_cache[rank_dir]


def clear_structure_cache() -> None:
    """Invalidate the structure cache (useful when case directory changes)."""
    with _struct_lock:
        _struct_cache.clear()


# ---------------------------------------------------------------------------
# File reader
# ---------------------------------------------------------------------------

def _resolve_path(h5file, path_bytes):
    """Resolve a silo internal path (e.g. b'/.silo/#000003') to a numpy array."""
    path = path_bytes.decode() if isinstance(path_bytes, (bytes, np.bytes_)) else str(path_bytes)
    return np.array(h5file[path])


def read_silo_file(  # pylint: disable=too-many-locals
    path: str,
    var_filter: Optional[str] = None,
    rank_dir: Optional[str] = None,
) -> ProcessorData:
    """
    Read a single Silo-HDF5 file produced by MFC post_process.

    On the first call for a given *rank_dir* the file structure is parsed
    and cached.  Subsequent calls for the same rank skip all metadata work
    and read data directly.

    Args:
        path:       Path to the ``.silo`` file.
        var_filter: If given, only load this variable (case-sensitive).
        rank_dir:   Directory containing this file (used as structure cache
                    key).  Defaults to ``os.path.dirname(path)``.
    """
    if rank_dir is None:
        rank_dir = os.path.dirname(path)

    struct = _get_structure(rank_dir, path)

    with h5py.File(path, "r") as f:
        # Coordinates — read directly by cached path
        coords = [_resolve_path(f, cp) for cp in struct.coord_paths]

        x_cb = coords[0]
        y_cb = coords[1] if struct.ndims >= 2 else np.array([0.0])
        z_cb = coords[2] if struct.ndims >= 3 else np.array([0.0])

        m = len(x_cb) - 1
        n = (len(y_cb) - 1) if struct.ndims >= 2 else 0
        p = (len(z_cb) - 1) if struct.ndims >= 3 else 0

        # Variables — read directly by cached path, skip metadata iteration
        variables: Dict[str, np.ndarray] = {}
        for var_name, data_path in struct.var_paths.items():
            if var_filter is not None and var_name != var_filter:
                continue
            data = _resolve_path(f, data_path).astype(np.float64)

            # MFC's DBPUTQV1 passes the Fortran column-major array as a flat
            # buffer.  HDF5 stores it row-major.  Reinterpret in Fortran order
            # so data[i,j,k] = value at (x_i, y_j, z_k).
            if data.ndim >= 2:
                data = np.ascontiguousarray(data).ravel().reshape(data.shape, order='F')
            variables[var_name] = data

    return ProcessorData(
        m=m, n=n, p=p, x_cb=x_cb, y_cb=y_cb, z_cb=z_cb, variables=variables
    )


# ---------------------------------------------------------------------------
# Persistent thread pool for parallel rank reads
# ---------------------------------------------------------------------------

_READ_POOL: Optional[ThreadPoolExecutor] = None
_POOL_LOCK = threading.Lock()


def _get_pool() -> ThreadPoolExecutor:
    """Return a module-level thread pool, creating it on first use."""
    global _READ_POOL  # pylint: disable=global-statement
    with _POOL_LOCK:
        if _READ_POOL is None:
            _READ_POOL = ThreadPoolExecutor(
                max_workers=32, thread_name_prefix="mfc_silo"
            )
            atexit.register(_READ_POOL.shutdown, wait=False)
        return _READ_POOL


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble_silo(
    case_dir: str,
    step: int,
    var: Optional[str] = None,
) -> AssembledData:
    """
    Read and assemble multi-processor Silo-HDF5 data for a given timestep.
    """

    base = os.path.join(case_dir, "silo_hdf5")
    if not os.path.isdir(base):
        raise FileNotFoundError(f"Silo-HDF5 directory not found: {base}")

    ranks: List[int] = []
    for entry in os.listdir(base):
        if entry.startswith("p") and entry[1:].isdigit():
            ranks.append(int(entry[1:]))
    ranks.sort()

    if not ranks:
        raise FileNotFoundError(f"No processor directories in {base}")

    # Validate all paths exist synchronously so errors are immediate.
    rank_paths: List[tuple] = []
    for rank in ranks:
        rank_dir = os.path.join(base, f"p{rank}")
        silo_file = os.path.join(rank_dir, f"{step}.silo")
        if not os.path.isfile(silo_file):
            raise FileNotFoundError(
                f"Processor file not found: {silo_file}. "
                "Incomplete output (missing rank) would produce incorrect data."
            )
        rank_paths.append((rank, rank_dir, silo_file))

    def _read_one(args):
        rank, rank_dir, silo_file = args
        return rank, read_silo_file(silo_file, var_filter=var, rank_dir=rank_dir)

    pool = _get_pool()
    results = list(pool.map(_read_one, rank_paths))

    proc_data: List[Tuple[int, ProcessorData]] = list(results)

    if not proc_data:
        raise FileNotFoundError(f"No Silo data found for step {step}")

    return assemble_from_proc_data(proc_data)
