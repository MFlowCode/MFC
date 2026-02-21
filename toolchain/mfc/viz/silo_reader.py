"""
Silo-HDF5 reader for MFC post-processed output.

Silo files produced by MFC are valid HDF5 underneath. This reader
uses h5py to navigate the HDF5 tree and extract mesh coordinates
and variable arrays.

Requires: h5py (optional dependency).
"""

import os
from typing import List, Optional, Tuple

import numpy as np

from .reader import AssembledData, ProcessorData

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


def _check_h5py():
    if not HAS_H5PY:
        raise ImportError(
            "h5py is required to read Silo-HDF5 files.\n"
            "Install it with: pip install h5py\n"
            "Or re-run post_process with format=2 to produce binary output."
        )


def _find_mesh_and_vars(h5file):
    """Navigate the HDF5 tree to find mesh coordinates and variables."""
    mesh_coords = {}
    variables = {}

    # Silo stores data in a nested structure. Common patterns:
    # /<meshname>/ contains coordinate arrays
    # Variables are stored at the top level or in subdirectories
    for key in h5file.keys():
        obj = h5file[key]
        if isinstance(obj, h5py.Dataset):
            variables[key] = np.array(obj)
        elif isinstance(obj, h5py.Group):
            # Check for mesh data
            for subkey in obj.keys():
                subobj = obj[subkey]
                if isinstance(subobj, h5py.Dataset):
                    arr = np.array(subobj)
                    if subkey in ('coord0', 'coord1', 'coord2'):
                        mesh_coords[subkey] = arr
                    elif subkey.startswith('_coord'):
                        mesh_coords[subkey] = arr
                    else:
                        variables[subkey] = arr

    return mesh_coords, variables


def read_silo_file(path: str, var_filter: Optional[str] = None) -> ProcessorData:
    """
    Read a single Silo-HDF5 file.

    Args:
        path: Path to the .silo file.
        var_filter: If given, only load this variable.

    Returns:
        ProcessorData with grid and variable data.
    """
    _check_h5py()

    with h5py.File(path, 'r') as f:
        mesh_coords, raw_vars = _find_mesh_and_vars(f)

    # Extract coordinates
    x_cb = mesh_coords.get('coord0', np.array([0.0, 1.0]))
    y_cb = mesh_coords.get('coord1', np.array([0.0]))
    z_cb = mesh_coords.get('coord2', np.array([0.0]))

    m = len(x_cb) - 2 if len(x_cb) > 1 else 0
    n = len(y_cb) - 2 if len(y_cb) > 1 else 0
    p = len(z_cb) - 2 if len(z_cb) > 1 else 0

    variables = {}
    for name, data in raw_vars.items():
        if var_filter is not None and name != var_filter:
            continue
        variables[name] = data.astype(np.float64)

    return ProcessorData(m=m, n=n, p=p, x_cb=x_cb, y_cb=y_cb, z_cb=z_cb, variables=variables)


def discover_timesteps_silo(case_dir: str) -> List[int]:
    """Return sorted list of available timesteps from silo_hdf5/ directory."""
    p0_dir = os.path.join(case_dir, 'silo_hdf5', 'p0')
    if not os.path.isdir(p0_dir):
        return []
    steps = set()
    for dname in os.listdir(p0_dir):
        if dname.startswith('t_step='):
            try:
                steps.add(int(dname.split('=')[1]))
            except (ValueError, IndexError):
                pass
    return sorted(steps)


def assemble_silo(case_dir: str, step: int,  # pylint: disable=too-many-locals,too-many-statements
                  var: Optional[str] = None) -> AssembledData:
    """
    Read and assemble multi-processor Silo-HDF5 data for a given timestep.
    """
    _check_h5py()

    base = os.path.join(case_dir, 'silo_hdf5')
    ranks = []
    for entry in os.listdir(base):
        if entry.startswith('p') and entry[1:].isdigit():
            ranks.append(int(entry[1:]))
    ranks.sort()

    if not ranks:
        raise FileNotFoundError(f"No processor directories in {base}")

    proc_data: List[Tuple[int, ProcessorData]] = []
    for rank in ranks:
        silo_dir = os.path.join(base, f'p{rank}', f't_step={step}')
        if not os.path.isdir(silo_dir):
            continue
        silo_file = os.path.join(silo_dir, f'{step}.silo')
        if not os.path.isfile(silo_file):
            # Try finding any .silo file in the directory
            for f in os.listdir(silo_dir):
                if f.endswith('.silo'):
                    silo_file = os.path.join(silo_dir, f)
                    break
        if not os.path.isfile(silo_file):
            continue
        pdata = read_silo_file(silo_file, var_filter=var)
        proc_data.append((rank, pdata))

    if not proc_data:
        raise FileNotFoundError(f"No Silo data found for step {step}")

    # For single processor, return directly
    if len(proc_data) == 1:
        _, pd = proc_data[0]
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        z_cc = (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        ndim = 1
        if pd.n > 0:
            ndim = 2
        if pd.p > 0:
            ndim = 3
        return AssembledData(ndim=ndim, x_cc=x_cc, y_cc=y_cc, z_cc=z_cc,
                             variables=pd.variables)

    # Multi-processor assembly — simplified version of binary reader's logic

    sample = proc_data[0][1]
    ndim = 1
    if sample.n > 0:
        ndim = 2
    if sample.p > 0:
        ndim = 3

    proc_centers = []
    for rank, pd in proc_data:
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        z_cc = (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        proc_centers.append((rank, pd, x_cc, y_cc, z_cc))

    # Build global coordinates by sorting and concatenating unique chunks
    x_chunks = {}
    y_chunks = {}
    z_chunks = {}

    for rank, pd, x_cc, y_cc, z_cc in proc_centers:
        x_key = round(x_cc[0], 12)
        y_key = round(y_cc[0], 12) if ndim >= 2 else 0.0
        z_key = round(z_cc[0], 12) if ndim >= 3 else 0.0
        if x_key not in x_chunks:
            x_chunks[x_key] = (len(x_cc), x_cc)
        if y_key not in y_chunks:
            y_chunks[y_key] = (len(y_cc), y_cc)
        if z_key not in z_chunks:
            z_chunks[z_key] = (len(z_cc), z_cc)

    sorted_x_keys = sorted(x_chunks.keys())
    sorted_y_keys = sorted(y_chunks.keys())
    sorted_z_keys = sorted(z_chunks.keys())

    global_x = np.concatenate([x_chunks[k][1] for k in sorted_x_keys])
    global_y = np.concatenate([y_chunks[k][1] for k in sorted_y_keys]) if ndim >= 2 else np.array([0.0])
    global_z = np.concatenate([z_chunks[k][1] for k in sorted_z_keys]) if ndim >= 3 else np.array([0.0])

    x_offsets = {}
    off = 0
    for k in sorted_x_keys:
        x_offsets[k] = off
        off += x_chunks[k][0]

    y_offsets = {}
    off = 0
    for k in sorted_y_keys:
        y_offsets[k] = off
        off += y_chunks[k][0]

    z_offsets = {}
    off = 0
    for k in sorted_z_keys:
        z_offsets[k] = off
        off += z_chunks[k][0]

    varnames = list(proc_data[0][1].variables.keys())
    nx, ny, nz = len(global_x), len(global_y), len(global_z)

    global_vars = {}
    for vn in varnames:
        if ndim == 3:
            global_vars[vn] = np.zeros((nx, ny, nz))
        elif ndim == 2:
            global_vars[vn] = np.zeros((nx, ny))
        else:
            global_vars[vn] = np.zeros(nx)

    for rank, pd, x_cc, y_cc, z_cc in proc_centers:
        x_key = round(x_cc[0], 12)
        y_key = round(y_cc[0], 12) if ndim >= 2 else 0.0
        z_key = round(z_cc[0], 12) if ndim >= 3 else 0.0

        xi = x_offsets[x_key]
        yi = y_offsets[y_key] if ndim >= 2 else 0
        zi = z_offsets[z_key] if ndim >= 3 else 0

        for vn, data in pd.variables.items():
            if vn not in global_vars:
                continue
            if ndim == 3:
                global_vars[vn][xi:xi + pd.m + 1, yi:yi + pd.n + 1, zi:zi + pd.p + 1] = data
            elif ndim == 2:
                global_vars[vn][xi:xi + pd.m + 1, yi:yi + pd.n + 1] = data
            else:
                global_vars[vn][xi:xi + pd.m + 1] = data

    return AssembledData(ndim=ndim, x_cc=global_x, y_cc=global_y, z_cc=global_z,
                         variables=global_vars)
