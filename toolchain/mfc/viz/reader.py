"""
Binary format reader for MFC post-processed output.

Reads Fortran unformatted binary files produced by MFC's post_process
with format=2. Each Fortran `write` produces one record:
  [4-byte record-marker][payload][4-byte record-marker]

File layout per processor:
  Record 1 (header):   m(int32), n(int32), p(int32), dbvars(int32)
  Record 2 (grid):     x_cb [, y_cb [, z_cb]]  (float32 or float64)
  Records 3..N (vars): varname(50-char) + data((m+1)*(n+1)*(p+1) floats)
"""

import os
import struct
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


NAME_LEN = 50  # Fortran character length for variable names


@dataclass
class ProcessorData:
    """Data from a single processor file."""
    m: int
    n: int
    p: int
    x_cb: np.ndarray
    y_cb: np.ndarray
    z_cb: np.ndarray
    variables: Dict[str, np.ndarray] = field(default_factory=dict)


@dataclass
class AssembledData:
    """Assembled multi-processor data on a global grid."""
    ndim: int
    x_cc: np.ndarray
    y_cc: np.ndarray
    z_cc: np.ndarray
    variables: Dict[str, np.ndarray] = field(default_factory=dict)


def read_record(f) -> bytes:
    """Read one Fortran unformatted record, returning the payload bytes."""
    raw = f.read(4)
    if len(raw) < 4:
        raise EOFError("Unexpected end of file reading record marker")
    rec_len = struct.unpack('<i', raw)[0]
    if rec_len < 0:
        # Try big-endian
        rec_len = struct.unpack('>i', raw)[0]
        if rec_len < 0:
            raise ValueError(f"Invalid record length: {rec_len}")
    payload = f.read(rec_len)
    if len(payload) < rec_len:
        raise EOFError("Unexpected end of file reading record payload")
    f.read(4)  # trailing marker
    return payload


def _detect_endianness(path: str) -> str:
    """Detect endianness from the first record marker (should be 16 for header)."""
    with open(path, 'rb') as f:
        raw = f.read(4)
    le = struct.unpack('<i', raw)[0]
    if le == 16:
        return '<'
    be = struct.unpack('>i', raw)[0]
    if be == 16:
        return '>'
    raise ValueError(
        f"Cannot detect endianness: first record marker is {le} (LE) / {be} (BE), expected 16"
    )


def _read_record_endian(f, endian: str) -> bytes:
    """Read one Fortran unformatted record with known endianness."""
    raw = f.read(4)
    if len(raw) < 4:
        raise EOFError("Unexpected end of file reading record marker")
    rec_len = struct.unpack(f'{endian}i', raw)[0]
    payload = f.read(rec_len)
    if len(payload) < rec_len:
        raise EOFError("Unexpected end of file reading record payload")
    f.read(4)  # trailing marker
    return payload


def read_binary_file(path: str, var_filter: Optional[str] = None) -> ProcessorData:  # pylint: disable=too-many-locals,too-many-statements
    """
    Read a single MFC binary post-process file.

    Args:
        path: Path to the .dat file.
        var_filter: If given, only load this variable (skip others).

    Returns:
        ProcessorData with grid and variable data.
    """
    endian = _detect_endianness(path)

    with open(path, 'rb') as f:
        # Record 1: header [m, n, p, dbvars] — 4 int32
        hdr = _read_record_endian(f, endian)
        m, n, p, dbvars = struct.unpack(f'{endian}4i', hdr)

        # Record 2: grid coordinates — all in one record
        grid_raw = _read_record_endian(f, endian)
        grid_bytes = len(grid_raw)

        # Determine number of grid values
        if p > 0:
            n_vals = (m + 2) + (n + 2) + (p + 2)
        elif n > 0:
            n_vals = (m + 2) + (n + 2)
        else:
            n_vals = m + 2

        # Auto-detect grid precision from record size
        bytes_per_val = grid_bytes / n_vals
        if abs(bytes_per_val - 8.0) < 0.5:
            grid_dtype = np.dtype(f'{endian}f8')
        elif abs(bytes_per_val - 4.0) < 0.5:
            grid_dtype = np.dtype(f'{endian}f4')
        else:
            raise ValueError(
                f"Cannot determine grid precision: {grid_bytes} bytes for {n_vals} values "
                f"({bytes_per_val:.1f} bytes/value)"
            )

        grid_arr = np.frombuffer(grid_raw, dtype=grid_dtype)

        # Split into x_cb, y_cb, z_cb
        offset = 0
        x_cb = grid_arr[offset:offset + m + 2].astype(np.float64)
        offset += m + 2
        if n > 0:
            y_cb = grid_arr[offset:offset + n + 2].astype(np.float64)
            offset += n + 2
        else:
            y_cb = np.array([0.0])
        if p > 0:
            z_cb = grid_arr[offset:offset + p + 2].astype(np.float64)
        else:
            z_cb = np.array([0.0])

        # Records 3..N: variables
        variables: Dict[str, np.ndarray] = {}
        data_size = (m + 1) * max(n + 1, 1) * max(p + 1, 1)

        for _ in range(dbvars):
            var_raw = _read_record_endian(f, endian)
            varname = var_raw[:NAME_LEN].decode('ascii', errors='replace').strip()

            if var_filter is not None and varname != var_filter:
                continue

            # Auto-detect variable data precision from record size
            data_bytes = len(var_raw) - NAME_LEN
            var_bpv = data_bytes / data_size
            if abs(var_bpv - 8.0) < 0.5:
                var_dtype = np.dtype(f'{endian}f8')
            elif abs(var_bpv - 4.0) < 0.5:
                var_dtype = np.dtype(f'{endian}f4')
            else:
                raise ValueError(
                    f"Cannot determine variable precision for '{varname}': "
                    f"{data_bytes} bytes for {data_size} values ({var_bpv:.1f} bytes/value)"
                )

            data = np.frombuffer(var_raw[NAME_LEN:], dtype=var_dtype).astype(np.float64)

            # Reshape for multi-dimensional data (Fortran column-major order)
            if p > 0:
                data = data.reshape((m + 1, n + 1, p + 1), order='F')
            elif n > 0:
                data = data.reshape((m + 1, n + 1), order='F')

            variables[varname] = data

    return ProcessorData(m=m, n=n, p=p, x_cb=x_cb, y_cb=y_cb, z_cb=z_cb, variables=variables)


def discover_format(case_dir: str) -> str:
    """Detect whether case has binary or silo_hdf5 output."""
    if os.path.isdir(os.path.join(case_dir, 'binary')):
        return 'binary'
    if os.path.isdir(os.path.join(case_dir, 'silo_hdf5')):
        return 'silo'
    raise FileNotFoundError(
        f"No 'binary/' or 'silo_hdf5/' directory found in {case_dir}. "
        "Run post_process with format=1 (Silo) or format=2 (binary) first."
    )


def discover_timesteps(case_dir: str, fmt: str) -> List[int]:
    """Return sorted list of available timesteps."""
    if fmt == 'binary':
        # Check root/ first (1D), then p0/
        root_dir = os.path.join(case_dir, 'binary', 'root')
        if os.path.isdir(root_dir):
            steps = set()
            for fname in os.listdir(root_dir):
                if fname.endswith('.dat'):
                    try:
                        steps.add(int(fname[:-4]))
                    except ValueError:
                        pass
            if steps:
                return sorted(steps)

        # Multi-dimensional: look in p0/
        p0_dir = os.path.join(case_dir, 'binary', 'p0')
        if os.path.isdir(p0_dir):
            steps = set()
            for fname in os.listdir(p0_dir):
                if fname.endswith('.dat'):
                    try:
                        steps.add(int(fname[:-4]))
                    except ValueError:
                        pass
            return sorted(steps)

    elif fmt == 'silo':
        p0_dir = os.path.join(case_dir, 'silo_hdf5', 'p0')
        if os.path.isdir(p0_dir):
            steps = set()
            for fname in os.listdir(p0_dir):
                if fname.endswith('.silo') and not fname.startswith('collection'):
                    try:
                        steps.add(int(fname[:-5]))
                    except ValueError:
                        pass
            return sorted(steps)

    return []


def _discover_processors(case_dir: str, fmt: str) -> List[int]:
    """Return sorted list of processor ranks."""
    if fmt == 'binary':
        base = os.path.join(case_dir, 'binary')
    else:
        base = os.path.join(case_dir, 'silo_hdf5')

    ranks = []
    if not os.path.isdir(base):
        return ranks
    for entry in os.listdir(base):
        if entry.startswith('p') and entry[1:].isdigit():
            ranks.append(int(entry[1:]))
    return sorted(ranks)


def _is_1d(case_dir: str) -> bool:
    """Check if the output is 1D (has binary/root/ directory)."""
    return os.path.isdir(os.path.join(case_dir, 'binary', 'root'))


def assemble(case_dir: str, step: int, fmt: str = 'binary',  # pylint: disable=too-many-locals,too-many-statements
             var: Optional[str] = None) -> AssembledData:
    """
    Read and assemble multi-processor data for a given timestep.

    For 1D, reads the root file directly.
    For 2D/3D, reads all processor files and assembles into global arrays.
    """
    if fmt != 'binary':
        raise ValueError(f"Format '{fmt}' not supported by binary reader. Use silo_reader.")

    # 1D case: read root file directly
    if _is_1d(case_dir):
        root_path = os.path.join(case_dir, 'binary', 'root', f'{step}.dat')
        if not os.path.isfile(root_path):
            raise FileNotFoundError(f"Root file not found: {root_path}")
        pdata = read_binary_file(root_path, var_filter=var)
        x_cc = (pdata.x_cb[:-1] + pdata.x_cb[1:]) / 2.0
        return AssembledData(
            ndim=1, x_cc=x_cc,
            y_cc=np.array([0.0]), z_cc=np.array([0.0]),
            variables=pdata.variables,
        )

    # Multi-dimensional: read all processor files
    ranks = _discover_processors(case_dir, fmt)
    if not ranks:
        raise FileNotFoundError(f"No processor directories found in {case_dir}/binary/")

    # Read all processor data
    proc_data: List[Tuple[int, ProcessorData]] = []
    for rank in ranks:
        fpath = os.path.join(case_dir, 'binary', f'p{rank}', f'{step}.dat')
        if not os.path.isfile(fpath):
            continue
        pdata = read_binary_file(fpath, var_filter=var)
        if pdata.m == 0 and pdata.n == 0 and pdata.p == 0:
            continue
        proc_data.append((rank, pdata))

    if not proc_data:
        raise FileNotFoundError(f"No valid processor data found for step {step}")

    ndim = 1
    sample = proc_data[0][1]
    if sample.n > 0:
        ndim = 2
    if sample.p > 0:
        ndim = 3

    # Compute cell centers for each processor
    proc_centers = []
    for rank, pd in proc_data:
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        z_cc = (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
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
        ndim=ndim, x_cc=global_x, y_cc=global_y, z_cc=global_z,
        variables=global_vars,
    )
