"""
MFC amr_check command - compare the interior of an AMR block, cell for cell, with a uniform fine-grid run.

Away from its edges a block advances the same cells as the uniform fine grid, so the two must agree to roundoff. The
block's initial data is interpolated from the coarse grid, so this holds only when that interpolation reproduces the
fine-grid initial condition (piecewise-constant patches on coarse cell faces do); the command checks that first.
"""

import glob
import itertools
import json
import math
import os
import subprocess
import sys

from .common import MFC_ROOT_DIR, MFCException, create_directory, file_write, system
from .printer import cons
from .run import input as run_input
from .state import ARG

DIMS = ("m", "n", "p")


def _active(params):
    return [d for d, k in enumerate(DIMS) if (params.get(k) or 0) > 0]


def block(params):
    """(lo, hi) level-0 cell indices of the block per active dimension: the case's own, else the central half."""
    out = []
    for d in _active(params):
        cells = params[DIMS[d]] + 1
        lo = params.get(f"amr_block_beg({d + 1})", cells // 4)
        out.append((lo, params.get(f"amr_block_end({d + 1})", lo + cells // 2 - 1)))
    return out


def variants(params, steps):
    """The case with one static block, and the same case on a uniform grid at the block's resolution."""
    if any(params.get(k) == "T" for k in ("cfl_dt", "cfl_adap_dt", "cfl_const_dt")):
        raise MFCException("amr_check needs a fixed dt: both runs must take the same time steps.")
    if params.get("t_step_start", 0) != 0:
        raise MFCException("amr_check starts both runs from the initial condition: set t_step_start = 0.")
    both = {**params, "t_step_stop": steps, "t_step_save": steps, "parallel_io": "F"}
    amr = {**both, "amr": "T", "amr_regrid_int": 0}
    for d, (lo, hi) in zip(_active(params), block(params)):
        amr.update({f"amr_block_beg({d + 1})": lo, f"amr_block_end({d + 1})": hi})
    fine = {k: v for k, v in both.items() if not k.startswith("amr")}
    fine.update({DIMS[d]: 2 * params[DIMS[d]] + 1 for d in _active(params)})
    if params.get("amr_subcycle") == "T":
        fine.update({"dt": params["dt"] / 2, "t_step_stop": 2 * steps, "t_step_save": 2 * steps})
    return amr, fine


def reach(params, steps):
    """Coarse cells from a block edge that the edge can influence in `steps` steps: one stencil per Runge-Kutta stage."""
    stencil = (2 if params.get("recon_type") in (2, "muscl") else params.get("weno_order", 5) // 2 + 1) + 1
    stages = params["time_stepper"] if isinstance(params.get("time_stepper"), int) else 3
    return math.ceil(steps * stages * (2 if params.get("amr_subcycle") == "T" else 1) * stencil / 2) + 1


def load_fields(dirpath, step):
    """{variable: {cell index tuple: value}} from the serial D/cons.* files of one save, all ranks."""
    rows = {}
    for path in glob.glob(os.path.join(dirpath, "D", f"cons.*.{step:06d}.dat")):
        with open(path) as f:
            rows.setdefault(int(os.path.basename(path).split(".")[1]), []).extend(tuple(map(float, line.split())) for line in f if line.strip())
    fields = {}
    for var, data in rows.items():
        index = [{x: i for i, x in enumerate(sorted({r[d] for r in data}))} for d in range(len(data[0]) - 1)]
        fields[var] = {tuple(index[d][r[d]] for d in range(len(index))): r[-1] for r in data}
    return fields


def _cells(blk):
    return itertools.product(*(range(lo, hi + 1) for lo, hi in blk))


def restrict(fine, blk):
    """Mean of the fine cells under each block cell."""
    corners = list(itertools.product((0, 1), repeat=len(blk)))
    return {var: {c: sum(f[tuple(2 * i + o for i, o in zip(c, corner))] for corner in corners) / len(corners) for c in _cells(blk)} for var, f in fine.items()}


def depth_errors(a, b, blk):
    """Largest |a - b| / max|b| over the variables, per depth (cells from the nearest block edge)."""
    worst = {}
    for var, ref in b.items():
        scale = max(map(abs, ref.values())) or 1.0
        for c in _cells(blk):
            depth = min(min(i - lo, hi - i) for i, (lo, hi) in zip(c, blk))
            worst[depth] = max(worst.get(depth, 0.0), abs(a[var][c] - ref[c]) / scale)
    return worst


def verdict(initial, final, active, reach, tol):  # pylint: disable=redefined-outer-name
    """(passed, reason). The core is the block cells at least `reach` cells from every block edge."""
    core = [d for d in final if d >= reach]
    if not core:
        return False, f"the block is too small: no cell is {reach} cells from its edges (fewer --steps or a larger block)"
    if max(initial.values()) > tol:
        return False, f"the initial data differs by {max(initial.values()):.1e}: interpolating the coarse initial condition does not reproduce the fine one"
    if max(active[d] for d in core) == 0.0:
        return False, "the core did not change during the run: put part of the flow inside the block"
    err = max(final[d] for d in core)
    return err <= tol, f"core differs from the fine-grid run by {err:.1e} (tolerance {tol:.1e})"


def absolute_paths(params, dirpath):
    """The parameters, with each file named relative to the case directory made absolute: the runs use other directories."""
    return {k: os.path.join(dirpath, v) if isinstance(v, str) and os.path.isfile(os.path.join(dirpath, v)) else v for k, v in params.items()}


def write_case(dirpath, params):
    create_directory(dirpath)
    case = os.path.join(dirpath, "case.json")
    file_write(case, json.dumps(params, indent=4) + "\n")
    return case


def _run(dirpath, params):
    case = write_case(dirpath, params)
    log = os.path.join(dirpath, "out.txt")
    with open(log, "w") as out:
        cmd = [os.path.join(MFC_ROOT_DIR, "mfc.sh"), "run", case, "-t", "pre_process", "simulation", *ARG("--")]
        if system(cmd, print_cmd=False, stdout=out, stderr=subprocess.STDOUT).returncode != 0:
            raise MFCException(f"amr_check: the run failed, see {log}.")


def amr_check():
    """Run the case with a block and on the uniform fine grid, then compare the block's interior."""
    case = run_input.load(ARG("input"), do_print=False)
    steps, tol = ARG("steps"), ARG("tol")
    blk, far = block(case.params), reach(case.params, steps)
    work = os.path.join(case.dirpath, "amr_check")
    saves = {}
    for name, params in zip(("amr", "fine"), variants(absolute_paths(case.params, case.dirpath), steps)):
        cons.print(f"Running the {name} case in {os.path.join(work, name)}...")
        _run(os.path.join(work, name), params)
        saves[name] = [load_fields(os.path.join(work, name), t) for t in (0, params["t_step_stop"])]
    fine = [restrict(f, blk) for f in saves["fine"]]
    initial, final = (depth_errors(a, f, blk) for a, f in zip(saves["amr"], fine))
    active = depth_errors(fine[1], fine[0], blk)
    passed, reason = verdict(initial, final, active, far, tol)
    cons.print(f"Block {blk}, {steps} step(s); core = cells at least {far} from every block edge.")
    cons.print("  depth  |amr - fine|")
    for depth in sorted(d for d in final if d < far):
        cons.print(f"  {depth:5d}  {final[depth]:.1e}")
    if any(d >= far for d in final):
        cons.print(f"  >={far:3d}  {max(v for d, v in final.items() if d >= far):.1e}  core")
    cons.print(f"[bold {'green]✓ PASS' if passed else 'red]✗ FAIL'}[/]: {reason}")
    sys.exit(0 if passed else 1)
