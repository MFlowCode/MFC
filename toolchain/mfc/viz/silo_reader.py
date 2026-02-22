"""
Silo-HDF5 reader for MFC post-processed output.

Silo files produced by MFC are valid HDF5 underneath.  Each Silo object
is stored as an HDF5 Named Datatype whose ``silo`` compound attribute
carries the metadata (mesh name, data-array path, dimensions, etc.).
Actual data lives in numbered datasets under the ``.silo/`` group.

This reader uses h5py to navigate that structure.

Requires: h5py (optional dependency).
"""

import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from .reader import AssembledData, ProcessorData

try:
    import h5py

    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

# Silo type constants (from silo.h)
_DB_QUADMESH = 130
_DB_QUADVAR = 501


def _check_h5py():
    if not HAS_H5PY:
        raise ImportError(
            "h5py is required to read Silo-HDF5 files.\n"
            "Install it with: pip install h5py\n"
            "Or re-run post_process with format=2 to produce binary output."
        )


def _read_silo_object(h5file, name):
    """Read a Silo Named-Datatype object and return its ``silo`` attribute."""
    obj = h5file[name]
    if not isinstance(obj, h5py.Datatype):
        return None
    if "silo" not in obj.attrs:
        return None
    return obj.attrs["silo"]


def _resolve_path(h5file, path_bytes):
    """Resolve a silo internal path (e.g. b'/.silo/#000003') to a dataset."""
    path = path_bytes.decode() if isinstance(path_bytes, bytes) else str(path_bytes)
    return np.array(h5file[path])


def read_silo_file(  # pylint: disable=too-many-locals
    path: str, var_filter: Optional[str] = None
) -> ProcessorData:
    """
    Read a single Silo-HDF5 file produced by MFC post_process.

    Args:
        path: Path to the ``.silo`` file.
        var_filter: If given, only load this variable (case-sensitive).

    Returns:
        ProcessorData with grid coordinates and variable arrays.
    """
    _check_h5py()

    with h5py.File(path, "r") as f:
        # --- locate the mesh ------------------------------------------------
        mesh_name = None
        mesh_attr = None
        for key, obj in f.items():
            if key in ("..", ".silo"):
                continue
            if not isinstance(obj, h5py.Datatype):
                continue
            silo_type = obj.attrs.get("silo_type")
            if silo_type is not None and int(silo_type) == _DB_QUADMESH:
                mesh_name = key
                mesh_attr = obj.attrs["silo"]
                break

        if mesh_attr is None:
            raise ValueError(f"No rectilinear mesh found in {path}")

        ndims = int(mesh_attr["ndims"])
        coord_paths = []
        for i in range(ndims):
            coord_paths.append(mesh_attr[f"coord{i}"])

        coords = [_resolve_path(f, cp) for cp in coord_paths]

        x_cb = coords[0]
        y_cb = coords[1] if ndims >= 2 else np.array([0.0])
        z_cb = coords[2] if ndims >= 3 else np.array([0.0])

        # Grid dimensions: node counts minus 1 give cell counts
        m = len(x_cb) - 1
        n = (len(y_cb) - 1) if ndims >= 2 else 0
        p = (len(z_cb) - 1) if ndims >= 3 else 0

        # --- locate variables ------------------------------------------------
        variables: Dict[str, np.ndarray] = {}
        for key, obj in f.items():
            if key in ("..", ".silo", mesh_name):
                continue
            if not isinstance(obj, h5py.Datatype):
                continue
            silo_type = obj.attrs.get("silo_type")
            if silo_type is None or int(silo_type) != _DB_QUADVAR:
                continue

            # Apply variable filter
            if var_filter is not None and key != var_filter:
                continue

            attr = obj.attrs["silo"]
            data_path = attr["value0"]
            data = _resolve_path(f, data_path).astype(np.float64)

            # MFC's DBPUTQV1 passes the Fortran column-major array as a
            # flat buffer.  HDF5 stores it row-major.  Reinterpret the
            # bytes in Fortran order so data[i,j,k] = value at (x_i,y_j,z_k),
            # matching the binary reader convention.
            if data.ndim >= 2:
                data = np.ascontiguousarray(data).ravel().reshape(data.shape, order='F')
            variables[key] = data

    return ProcessorData(
        m=m, n=n, p=p, x_cb=x_cb, y_cb=y_cb, z_cb=z_cb, variables=variables
    )


def discover_timesteps_silo(case_dir: str) -> List[int]:
    """Return sorted list of available timesteps from ``silo_hdf5/`` directory."""
    p0_dir = os.path.join(case_dir, "silo_hdf5", "p0")
    if not os.path.isdir(p0_dir):
        return []
    steps = set()
    for fname in os.listdir(p0_dir):
        if fname.endswith(".silo") and not fname.startswith("collection"):
            try:
                steps.add(int(fname[:-5]))
            except ValueError:
                pass
    return sorted(steps)


def assemble_silo(  # pylint: disable=too-many-locals,too-many-statements,too-many-branches
    case_dir: str,
    step: int,
    var: Optional[str] = None,
) -> AssembledData:
    """
    Read and assemble multi-processor Silo-HDF5 data for a given timestep.
    """
    _check_h5py()

    base = os.path.join(case_dir, "silo_hdf5")
    ranks: List[int] = []
    for entry in os.listdir(base):
        if entry.startswith("p") and entry[1:].isdigit():
            ranks.append(int(entry[1:]))
    ranks.sort()

    if not ranks:
        raise FileNotFoundError(f"No processor directories in {base}")

    proc_data: List[Tuple[int, ProcessorData]] = []
    for rank in ranks:
        silo_file = os.path.join(base, f"p{rank}", f"{step}.silo")
        if not os.path.isfile(silo_file):
            continue
        pdata = read_silo_file(silo_file, var_filter=var)
        proc_data.append((rank, pdata))

    if not proc_data:
        raise FileNotFoundError(f"No Silo data found for step {step}")

    # --- single processor — fast path ------------------------------------
    if len(proc_data) == 1:
        _, pd = proc_data[0]
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (
            (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        )
        z_cc = (
            (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        )
        ndim = 1 + (pd.n > 0) + (pd.p > 0)
        return AssembledData(
            ndim=ndim,
            x_cc=x_cc,
            y_cc=y_cc,
            z_cc=z_cc,
            variables=pd.variables,
        )

    # --- multi-processor assembly ----------------------------------------
    sample = proc_data[0][1]
    ndim = 1 + (sample.n > 0) + (sample.p > 0)

    # Compute cell centers for each processor
    proc_centers: list = []
    for rank, pd in proc_data:
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (
            (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        )
        z_cc = (
            (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        )
        proc_centers.append((rank, pd, x_cc, y_cc, z_cc))

    # Build unique sorted global coordinate arrays (handles ghost overlap)
    all_x = np.concatenate([xc for _, _, xc, _, _ in proc_centers])
    global_x = np.unique(np.round(all_x, 12))
    if ndim >= 2:
        all_y = np.concatenate([yc for _, _, _, yc, _ in proc_centers])
        global_y = np.unique(np.round(all_y, 12))
    else:
        global_y = np.array([0.0])
    if ndim >= 3:
        all_z = np.concatenate([zc for _, _, _, _, zc in proc_centers])
        global_z = np.unique(np.round(all_z, 12))
    else:
        global_z = np.array([0.0])

    varnames = list(proc_data[0][1].variables.keys())
    nx, ny, nz = len(global_x), len(global_y), len(global_z)

    global_vars: Dict[str, np.ndarray] = {}
    for vn in varnames:
        if ndim == 3:
            global_vars[vn] = np.zeros((nx, ny, nz))
        elif ndim == 2:
            global_vars[vn] = np.zeros((nx, ny))
        else:
            global_vars[vn] = np.zeros(nx)

    # Place each processor's data using per-cell coordinate lookup
    # (handles ghost/buffer cell overlap between processors)
    for _rank, pd, x_cc, y_cc, z_cc in proc_centers:
        xi = np.searchsorted(global_x, np.round(x_cc, 12))
        yi = np.searchsorted(global_y, np.round(y_cc, 12)) if ndim >= 2 else np.array([0])
        zi = np.searchsorted(global_z, np.round(z_cc, 12)) if ndim >= 3 else np.array([0])

        for vn, data in pd.variables.items():
            if vn not in global_vars:
                continue
            if ndim == 3:
                global_vars[vn][np.ix_(xi, yi, zi)] = data
            elif ndim == 2:
                global_vars[vn][np.ix_(xi, yi)] = data
            else:
                global_vars[vn][xi] = data

    return AssembledData(
        ndim=ndim,
        x_cc=global_x,
        y_cc=global_y,
        z_cc=global_z,
        variables=global_vars,
    )
