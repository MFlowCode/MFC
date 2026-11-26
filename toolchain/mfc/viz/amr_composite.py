"""
Composite (finest available) cell data of a one-level AMR run written with parallel_io = T.

restart_data/lustre_<step>.dat holds the level-0 conservative fields (restricted fine data under the blocks) and
restart_data/lustre_amr_<step>.dat the fine blocks (layout: src/common/m_amr_restart_io.fpp). A uniform-grid run has no AMR file.
"""

import os
import re

import numpy as np


def read_inp(run_dir):
    """{name: value string} from pre_process.inp and simulation.inp."""
    out = {}
    for name in ("pre_process.inp", "simulation.inp"):
        with open(os.path.join(run_dir, name)) as f:
            out.update({k: v.strip().strip("'\"") for k, v in re.findall(r"^\s*([\w%()]+)\s*=\s*([^,\n]+)", f.read(), re.M)})
    return out


def parse_num(s):
    return float(s.replace("d", "e").replace("D", "e"))


def grid(inp):
    """Level-0 cell counts, lower domain corners and cell widths per active dimension."""
    cells, lo, dx = [], [], []
    for k, ax in (("m", "x"), ("n", "y"), ("p", "z")):
        c = int(inp.get(k, "0")) + 1
        if c == 1:
            break
        beg, end = parse_num(inp[f"{ax}_domain%beg"]), parse_num(inp[f"{ax}_domain%end"])
        cells.append(c)
        lo.append(beg)
        dx.append((end - beg) / c)
    return cells, np.array(lo), np.array(dx)


def leaves(run_dir, step, var=0):
    """Cells of the finest level at each point: a list of (centers, cell width, values) per patch, with patch 0 the uncovered
    level-0 cells (centers is a tuple of flat coordinate arrays). var is the 0-based conservative variable index."""
    inp = read_inp(run_dir)
    cells, lo, dx = grid(inp)
    ncell = int(np.prod(cells))
    raw = np.fromfile(os.path.join(run_dir, "restart_data", f"lustre_{step}.dat"))  # double precision only (checked below)
    if raw.size % ncell:
        raise SystemExit(f"{run_dir}: lustre_{step}.dat does not hold whole double-precision fields (single or mixed precision run?)")
    nvar = raw.size // ncell
    q0 = raw.reshape(nvar, ncell)[var].reshape(cells, order="F")
    centers0 = np.meshgrid(*[lo[d] + (np.arange(cells[d]) + 0.5) * dx[d] for d in range(len(cells))], indexing="ij")
    covered = np.zeros(cells, dtype=bool)
    patches = []
    for blo, bhi, lvl, ext, data in _blocks(os.path.join(run_dir, "restart_data", f"lustre_amr_{step}.dat"), len(cells)):
        if lvl != 1:
            raise SystemExit("amr_composite reads one refined level only")
        rr = (ext[0] + 1) // (bhi[0] - blo[0] + 1)
        fdx = dx / rr
        c = np.meshgrid(*[lo[d] + blo[d] * dx[d] + (np.arange(ext[d] + 1) + 0.5) * fdx[d] for d in range(len(cells))], indexing="ij")
        patches.append((tuple(a.ravel(order="F") for a in c), fdx, data[var].ravel(order="F"), (blo, bhi)))
        covered[tuple(slice(blo[d], bhi[d] + 1) for d in range(len(cells)))] = True
    keep = ~covered
    base = (tuple(a[keep] for a in centers0), dx, q0[keep], None)
    return [base] + patches, q0


def _blocks(path, ndim):
    """(lo, hi, level, extents, data[nvar, ...]) of each block in a parallel_io AMR restart file; none for a uniform run."""
    if not os.path.isfile(path):
        return
    with open(path, "rb") as f:
        hdr = np.fromfile(f, dtype=np.int32, count=6)
        nblk, nvar, wbytes = int(hdr[1]), int(hdr[2]), int(hdr[5])
        if wbytes != 8:
            raise SystemExit(f"{path}: word size {wbytes}; amr_composite reads double-precision runs only")
        for _ in range(nblk):
            rec = np.fromfile(f, dtype=np.int32, count=11)
            blo, bhi, lvl, own, ext = rec[0:3], rec[3:6], int(rec[6]), int(rec[7]), rec[8:11]
            if own <= 0:
                raise SystemExit(f"{path}: block without an owner")
            shape = tuple(int(e) + 1 for e in ext[:ndim])
            data = np.fromfile(f, dtype=np.float64, count=nvar * int(np.prod(shape)))
            yield blo[:ndim], bhi[:ndim], lvl, ext[:ndim], data.reshape((nvar,) + shape[::-1]).transpose([0] + list(range(ndim, 0, -1)))
