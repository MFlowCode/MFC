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

import glob
import itertools
import os
import struct
import warnings
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np


NAME_LEN = 50  # Fortran character length for variable names


@dataclass
class ProcessorData:
    """Data from a single processor file.

    m, n, p store dimension metadata but their exact semantics differ
    between readers:
      - Binary: m = Fortran header value, x_cb has m+2 elements,
        cell count is m+1.
      - Silo: m = cell count = len(x_cb) - 1.
    Assembly code intentionally avoids using m/n/p for array sizing —
    it derives everything from x_cb/y_cb/z_cb lengths.  If future code
    needs m directly, this discrepancy must be resolved.
    """
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



def _detect_endianness(path: str) -> str:
    """Detect endianness from the first record marker.

    The header record contains 4 int32s (m, n, p, dbvars) = 16 bytes,
    so the leading Fortran record marker must be 16.
    """
    with open(path, 'rb') as f:
        raw = f.read(4)
    if len(raw) < 4:
        raise EOFError(f"File too short to detect endianness: {path}")
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
    if rec_len < 0:
        raise ValueError(f"Invalid Fortran record length: {rec_len}")
    payload = f.read(rec_len)
    if len(payload) < rec_len:
        raise EOFError("Unexpected end of file reading record payload")
    trail = f.read(4)
    if len(trail) < 4:
        raise EOFError("Unexpected end of file reading trailing record marker")
    trail_len = struct.unpack(f'{endian}i', trail)[0]
    if trail_len != rec_len:
        raise ValueError(
            f"Fortran record marker mismatch: leading={rec_len}, trailing={trail_len}. "
            "File may be corrupted."
        )
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
        if m < 0 or n < 0 or p < 0 or dbvars < 0:
            raise ValueError(
                f"Invalid header in {path}: m={m}, n={n}, p={p}, dbvars={dbvars}"
            )

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
        if grid_bytes == n_vals * 8:
            grid_dtype = np.dtype(f'{endian}f8')
        elif grid_bytes == n_vals * 4:
            grid_dtype = np.dtype(f'{endian}f4')
        else:
            bytes_per_val = grid_bytes / n_vals if n_vals else 0
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
            # Read leading record-length marker and variable name only.
            # If this variable is filtered out we can seek past the data
            # instead of reading the full payload into memory.
            raw_len = f.read(4)
            if len(raw_len) < 4:
                raise EOFError("Unexpected end of file reading variable record marker")
            rec_len = struct.unpack(f'{endian}i', raw_len)[0]
            if rec_len < NAME_LEN:
                raise ValueError(f"Variable record too short: {rec_len} bytes")

            name_raw = f.read(NAME_LEN)
            if len(name_raw) < NAME_LEN:
                raise EOFError("Unexpected end of file reading variable name")
            varname = name_raw.decode('ascii', errors='replace').strip()

            data_bytes = rec_len - NAME_LEN
            if var_filter is not None and varname != var_filter:
                # Seek past remaining data + trailing record-length marker.
                f.seek(data_bytes + 4, 1)
                continue

            data_raw = f.read(data_bytes)
            if len(data_raw) < data_bytes:
                raise EOFError("Unexpected end of file reading variable data")
            trail = f.read(4)
            if len(trail) < 4:
                raise EOFError("Unexpected end of file reading trailing variable record marker")
            trail_len = struct.unpack(f'{endian}i', trail)[0]
            if trail_len != rec_len:
                raise ValueError(
                    f"Fortran record marker mismatch for '{varname}': "
                    f"leading={rec_len}, trailing={trail_len}"
                )

            # Auto-detect variable data precision from record size
            if data_bytes == data_size * 8:
                var_dtype = np.dtype(f'{endian}f8')
            elif data_bytes == data_size * 4:
                var_dtype = np.dtype(f'{endian}f4')
            else:
                var_bpv = data_bytes / data_size if data_size else 0
                raise ValueError(
                    f"Cannot determine variable precision for '{varname}': "
                    f"{data_bytes} bytes for {data_size} values ({var_bpv:.1f} bytes/value)"
                )

            data = np.frombuffer(data_raw, dtype=var_dtype).astype(np.float64)

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
    if fmt not in ('binary', 'silo'):
        raise ValueError(f"Unknown format '{fmt}'. Supported: 'binary', 'silo'.")

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

    return []  # no timestep files found in expected directories


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
    """Check if the output is 1D (binary/root/ exists with .dat files, no p0/ present)."""
    root = os.path.join(case_dir, 'binary', 'root')
    p0 = os.path.join(case_dir, 'binary', 'p0')
    return (os.path.isdir(root)
            and any(f.endswith('.dat') for f in os.listdir(root))
            and not os.path.isdir(p0))


def assemble_from_proc_data(  # pylint: disable=too-many-locals,too-many-statements
    proc_data: List[Tuple[int, ProcessorData]],
) -> AssembledData:
    """
    Assemble multi-processor data into a single global grid.

    Shared helper used by both binary and silo assembly paths.
    Handles ghost/buffer cell overlap between processors by using
    per-cell coordinate lookup (np.unique + np.searchsorted + np.ix_).
    """
    if not proc_data:
        raise ValueError("No processor data to assemble")

    # Single processor — fast path
    if len(proc_data) == 1:
        _, pd = proc_data[0]
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        z_cc = (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        ndim = 1 + (pd.n > 0) + (pd.p > 0)
        return AssembledData(
            ndim=ndim, x_cc=x_cc, y_cc=y_cc, z_cc=z_cc,
            variables=pd.variables,
        )

    # Multi-processor assembly
    sample = proc_data[0][1]
    ndim = 1 + (sample.n > 0) + (sample.p > 0)

    # Compute cell centers for each processor
    proc_centers = []
    for rank, pd in proc_data:
        x_cc = (pd.x_cb[:-1] + pd.x_cb[1:]) / 2.0
        y_cc = (pd.y_cb[:-1] + pd.y_cb[1:]) / 2.0 if pd.n > 0 else np.array([0.0])
        z_cc = (pd.z_cb[:-1] + pd.z_cb[1:]) / 2.0 if pd.p > 0 else np.array([0.0])
        proc_centers.append((rank, pd, x_cc, y_cc, z_cc))

    # Build unique sorted global coordinate arrays (handles ghost overlap).
    # Normalize each axis by its extent before rounding so that precision is
    # always 12 significant digits *relative to the domain size*.  This is
    # correct for both micro-scale domains (extent ~ 1e-10) and large-scale
    # domains (extent > 1e12) where the old formula (decimals = -log10(extent)
    # + 12) could go negative, causing np.round to round to the nearest 10 or
    # 100 and incorrectly merging distinct cell centers.
    def _dedup(arr):
        extent = arr.max() - arr.min()
        if extent > 0:
            origin = arr.min()
            norm = np.round((arr - origin) / extent, 12)
            return origin + np.unique(norm) * extent, origin, extent
        return np.unique(arr), arr.min(), 0.0

    def _norm_round(arr, origin, extent):
        """Round *arr* with the same relative tolerance used by _dedup."""
        if extent > 0:
            return origin + np.round((arr - origin) / extent, 12) * extent
        return arr

    all_x = np.concatenate([xc for _, _, xc, _, _ in proc_centers])
    global_x, x_orig, x_ext = _dedup(all_x)
    if ndim >= 2:
        all_y = np.concatenate([yc for _, _, _, yc, _ in proc_centers])
        global_y, y_orig, y_ext = _dedup(all_y)
    else:
        global_y, y_orig, y_ext = np.array([0.0]), 0.0, 0.0
    if ndim >= 3:
        all_z = np.concatenate([zc for _, _, _, _, zc in proc_centers])
        global_z, z_orig, z_ext = _dedup(all_z)
    else:
        global_z, z_orig, z_ext = np.array([0.0]), 0.0, 0.0

    varnames = sorted({vn for _, pd in proc_data for vn in pd.variables})
    nx, ny, nz = len(global_x), len(global_y), len(global_z)

    global_vars: Dict[str, np.ndarray] = {}
    for vn in varnames:
        if ndim == 3:
            global_vars[vn] = np.zeros((nx, ny, nz))
        elif ndim == 2:
            global_vars[vn] = np.zeros((nx, ny))
        else:
            global_vars[vn] = np.zeros(nx)

    # Place each processor's data using per-cell coordinate lookup.
    # Apply the same normalized rounding used by _dedup so that lookup
    # coordinates match the global grid entries exactly.
    for _rank, pd, x_cc, y_cc, z_cc in proc_centers:
        xi = np.clip(np.searchsorted(global_x, _norm_round(x_cc, x_orig, x_ext)), 0, nx - 1)
        yi = np.clip(np.searchsorted(global_y, _norm_round(y_cc, y_orig, y_ext)), 0, ny - 1) if ndim >= 2 else np.array([0])
        zi = np.clip(np.searchsorted(global_z, _norm_round(z_cc, z_orig, z_ext)), 0, nz - 1) if ndim >= 3 else np.array([0])

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


def assemble(case_dir: str, step: int, fmt: str = 'binary',  # pylint: disable=too-many-locals
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

    proc_data: List[Tuple[int, ProcessorData]] = []
    for rank in ranks:
        fpath = os.path.join(case_dir, 'binary', f'p{rank}', f'{step}.dat')
        if not os.path.isfile(fpath):
            raise FileNotFoundError(
                f"Processor file not found: {fpath}. "
                "Incomplete output (missing rank) would produce incorrect data."
            )
        pdata = read_binary_file(fpath, var_filter=var)
        if pdata.m == 0 and pdata.n == 0 and pdata.p == 0:
            warnings.warn(f"Processor p{rank} has zero dimensions, skipping", stacklevel=2)
            continue
        proc_data.append((rank, pdata))

    if not proc_data:
        raise FileNotFoundError(f"No valid processor data found for step {step}")

    return assemble_from_proc_data(proc_data)


# ---------------------------------------------------------------------------
# Lagrange bubble position reader
# ---------------------------------------------------------------------------

@lru_cache(maxsize=32)
def _nBubs_per_step(path: str) -> int:
    """Count how many bubble rows share the first time value in *path*.

    The result is cached so repeated calls for the same file (across
    different steps in an MP4 render) only scan the file once.
    """
    with open(path, encoding='ascii', errors='replace') as f:
        f.readline()  # skip header
        first = f.readline()
        if not first.strip():
            return 0
        t0 = first.split()[0]
        n = 1
        for line in f:
            parts = line.split()
            if parts and parts[0] == t0:
                n += 1
            else:
                break
    return n


def read_lag_bubbles_at_step(case_dir: str, step: int) -> Optional[np.ndarray]:
    """Read Lagrange bubble positions at a given save-step index.

    Reads ``D/lag_bubble_evol_<rank>.dat`` files written by MFC when
    ``lag_params%write_bubbles = T``.  Each file contains one row per
    bubble per simulation time-step (appended after every completed RK
    stage).  This function seeks efficiently to the block for *step* by
    counting the fixed number of bubbles per rank.

    Returns an ``(N, 4)`` float64 array of ``(x, y, z, r)`` in
    simulation-normalized units across all MPI ranks, or ``None`` when
    no bubble data is found.
    """
    d_dir = os.path.join(case_dir, 'D')
    if not os.path.isdir(d_dir):
        return None

    files = sorted(glob.glob(os.path.join(d_dir, 'lag_bubble_evol_*.dat')))
    if not files:
        return None

    chunks: List[np.ndarray] = []
    for fpath in files:
        try:
            nBubs = _nBubs_per_step(fpath)
            if nBubs == 0:
                continue
            # Line layout: 1 header + step * nBubs prior data rows.
            # MFC writes one bubble block per simulation timestep (at the last
            # RK stage), so block index == MFC timestep integer.  This is
            # correct for fresh runs (t_step_start=0).  For restarts where
            # t_step_start>0 the lag file starts at 0 but step numbers begin
            # at t_step_start — seeking would overshoot; restart support is
            # not yet implemented.
            skip = 1 + step * nBubs
            rows = []
            with open(fpath, encoding='ascii', errors='replace') as f:
                for _ in itertools.islice(f, skip):
                    pass
                for line in itertools.islice(f, nBubs):
                    parts = line.split()
                    if len(parts) >= 8:
                        # cols: time id x y z mv conc r [rdot p]
                        rows.append((float(parts[2]), float(parts[3]),
                                     float(parts[4]), float(parts[7])))
            if rows:
                chunks.append(np.array(rows, dtype=np.float64))
        except (OSError, ValueError):
            continue

    return np.concatenate(chunks, axis=0) if chunks else None


def has_lag_bubble_evol(case_dir: str) -> bool:
    """Return True if ``D/lag_bubble_evol_*.dat`` files exist in *case_dir*."""
    d_dir = os.path.join(case_dir, 'D')
    return bool(glob.glob(os.path.join(d_dir, 'lag_bubble_evol_*.dat')))
