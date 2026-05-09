"""Shared helpers for MFC convergence/order test runners."""

import json
import math
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

MFC = "./mfc.sh"
CONS_TOL = 1e-10


def read_cons_var(run_dir, step, var_idx, num_ranks=1, expected_size=None):
    """Read q_cons_vf{var_idx} from all ranks; rank-order concatenation."""
    chunks = []
    for rank in range(num_ranks):
        path = os.path.join(run_dir, "p_all", f"p{rank}", str(step), f"q_cons_vf{var_idx}.dat")
        with open(path, "rb") as f:
            rec_len = struct.unpack("i", f.read(4))[0]
            data = np.frombuffer(f.read(rec_len), dtype=np.float64)
            f.read(4)
        chunks.append(data.copy())
    combined = np.concatenate(chunks)
    if expected_size is not None and combined.size != expected_size:
        raise ValueError(f"Expected {expected_size} values across {num_ranks} ranks, got {combined.size}")
    return combined


def conservation_errors(run_dir, Nt, cell_vol, var_list, num_ranks, expected_size=None):
    """|Σq(T) - Σq(0)| / |Σq(0)| for each named variable."""
    errs = {}
    for name, idx in var_list:
        q0 = read_cons_var(run_dir, 0, idx, num_ranks, expected_size)
        qT = read_cons_var(run_dir, Nt, idx, num_ranks, expected_size)
        s0 = float(np.sum(q0)) * cell_vol
        sT = float(np.sum(qT)) * cell_vol
        errs[name] = abs(sT - s0) / (abs(s0) + 1e-300)
    return errs


def l2_norm(diff, scale):
    """sqrt(sum(diff^2) * scale)."""
    return float(np.sqrt(np.sum(diff**2) * scale))


def fit_rate(errors, h_values):
    """Least-squares slope of log(error) vs log(h)."""
    log_h = np.log(np.array(h_values, dtype=float))
    log_err = np.log(np.array(errors, dtype=float))
    slope, _ = np.polyfit(log_h, log_err, 1)
    return float(slope)


def pairwise_rates(errors, h_values):
    """Pairwise rates aligned with errors (first entry is None)."""
    rates = [None]
    for i in range(1, len(errors)):
        log_h0 = math.log(h_values[i - 1])
        log_h1 = math.log(h_values[i])
        rates.append((math.log(errors[i]) - math.log(errors[i - 1])) / (log_h1 - log_h0))
    return rates


def run_mfc_case(case_path, tmpdir, run_tag, case_args, num_ranks=1):
    """Run case.py once and copy p_all to tmpdir/run_tag. Returns (cfg_dict, run_dir)."""
    result = subprocess.run(
        [sys.executable, case_path, "--mfc", "{}"] + case_args,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"case.py failed:\n{result.stderr}")
    cfg = json.loads(result.stdout)

    cmd = [MFC, "run", case_path, "-t", "pre_process", "simulation", "-n", str(num_ranks), "--"] + case_args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd(), check=False)
    if result.returncode != 0:
        print(result.stdout[-3000:])
        print(result.stderr)
        raise RuntimeError(f"./mfc.sh run failed for {run_tag}")

    case_dir = os.path.dirname(case_path)
    src = os.path.join(case_dir, "p_all")
    dst = os.path.join(tmpdir, run_tag, "p_all")
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    shutil.rmtree(src, ignore_errors=True)
    shutil.rmtree(os.path.join(case_dir, "D"), ignore_errors=True)

    return cfg, os.path.join(tmpdir, run_tag)


def print_conservation_check(all_cons_errs, var_list, tol=CONS_TOL):
    """Print per-variable max conservation error and return overall pass."""
    print(f"\n  Conservation (need rel. error < {tol:.0e}):")
    passed = True
    for name, _ in var_list:
        max_err = max(ce[name] for ce in all_cons_errs)
        ok = max_err < tol
        print(f"    {name:<14}: max = {max_err:.2e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            passed = False
    return passed


def print_summary(results, label_width=18):
    """Print PASS/FAIL summary table; return overall bool."""
    print(f"\n{'=' * 60}")
    print("  Summary")
    print(f"{'=' * 60}")
    all_pass = True
    for label, passed in results.items():
        print(f"  {label:<{label_width}} {'PASS' if passed else 'FAIL'}")
        if not passed:
            all_pass = False
    return all_pass


def run_with_traceback(label, fn, *args, **kwargs):
    """Run a test function, print traceback on failure, return pass/fail bool."""
    try:
        return fn(*args, **kwargs)
    except Exception as e:
        import traceback

        print(f"  ERROR ({label}): {e}")
        traceback.print_exc()
        return False
