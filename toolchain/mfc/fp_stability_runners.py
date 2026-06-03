"""Verrou subprocess runners for the FP-stability suite.

Each routine drives the verrou/valgrind binary and returns parsed results.  Pure
parsing / metric helpers live in fp_stability_metrics, which this module imports.
"""

import glob
import os
import shutil
import subprocess
import tempfile

from .common import MFC_ROOT_DIR, MFCException
from .fp_stability_metrics import (
    VPREC_MANTISSA_BITS,
    _max_diff_np,
    _parse_cancel_gen,
    _parse_vg_error_locs,
)
from .printer import cons


def _has_verrou_tool(valgrind_bin: str, env: dict = None) -> bool:
    """True if this valgrind actually provides the 'verrou' tool. A plain system
    valgrind does not — accepting one would only fail later at run time. Pass env
    (with VALGRIND_LIB) to verify a relocated prebuilt tree, which cannot load its
    tool without it."""
    try:
        return subprocess.run([valgrind_bin, "--tool=verrou", "--version"], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False).returncode == 0
    except OSError:
        return False


def _find_verrou() -> str:
    verrou_home = os.environ.get("VERROU_HOME", os.path.join(os.path.expanduser("~"), ".local", "verrou"))
    candidate = os.path.join(verrou_home, "bin", "valgrind")
    # Require the $VERROU_HOME tree to actually run the verrou tool (with VALGRIND_LIB
    # for a relocated prebuilt). A broken/stale/non-Verrou tree there must read as
    # "absent" so it gets reinstalled, not used until it fails on every run.
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK) and _has_verrou_tool(candidate, _verrou_env(candidate)):
        return candidate
    # Fall back to a valgrind on PATH only if it is Verrou-enabled; a bare system
    # valgrind must read as "Verrou absent" so it gets installed, not misused. Verify
    # with VALGRIND_LIB too, so a relocated prebuilt on PATH (env.sh not sourced) isn't
    # wrongly judged absent.
    path_vg = shutil.which("valgrind")
    if path_vg and _has_verrou_tool(path_vg, _verrou_env(path_vg)):
        return path_vg
    return ""


def _find_binary(name: str) -> str:
    install_dir = os.path.join(MFC_ROOT_DIR, "build", "install")
    candidates = glob.glob(os.path.join(install_dir, "*", "bin", name))
    return max(candidates, key=os.path.getmtime) if candidates else ""


def _verrou_env(verrou_bin: str) -> dict:
    """os.environ plus VALGRIND_LIB, so a relocated install tree (e.g. a prebuilt
    artifact extracted to a new prefix) can locate its tool — Valgrind bakes its
    build prefix into the binary otherwise. Harmless for a source-built tree, where
    VALGRIND_LIB just equals the compiled-in path. A VALGRIND_LIB already in the
    environment (user sourced env.sh) is left untouched."""
    env = os.environ.copy()
    libdir = os.path.join(os.path.dirname(os.path.dirname(verrou_bin)), "libexec", "valgrind")
    if "VALGRIND_LIB" not in env and os.path.isdir(libdir):
        env["VALGRIND_LIB"] = libdir
    return env


def _write_inp(params: dict, target_name: str, work_dir: str) -> None:
    """Write a Fortran namelist .inp file from a Python params dict."""
    from .run import case_dicts

    master_keys = case_dicts.get_input_dict_keys(target_name)
    lines = [f"{k} = {v}" for k, v in params.items() if k in master_keys]
    with open(os.path.join(work_dir, f"{target_name}.inp"), "w") as fh:
        fh.write("&user_inputs\n" + "\n".join(lines) + "\n&end/\n")


def _run_preprocess(pp_bin: str, pre_params: dict, work_dir: str):
    _write_inp(pre_params, "pre_process", work_dir)
    with open(os.path.join(work_dir, "pre.log"), "w") as f:
        result = subprocess.run([pp_bin], cwd=work_dir, stdout=f, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise MFCException(f"pre_process failed (rc={result.returncode}). See {work_dir}/pre.log")


def _run_simulation_verrou(
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    run_dir: str,
    rounding_mode: str = None,
    extra_flags: list = None,
):
    """Copy ICs into a fresh tmpdir, run simulation under verrou, collect D/ output.

    rounding_mode is passed as --rounding-mode=<mode> when not None.
    extra_flags are appended before the binary (e.g. --backend=vprec ...).
    """
    with tempfile.TemporaryDirectory(prefix="mfc-fps-") as tmpdir:
        for fname in ["simulation.inp", "indices.dat", "pre_time_data.dat", "io_time_data.dat"]:
            src = os.path.join(work_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, tmpdir)
        shutil.copytree(os.path.join(work_dir, "p_all"), os.path.join(tmpdir, "p_all"))
        os.makedirs(os.path.join(tmpdir, "D"))

        log_path = os.path.join(run_dir, "verrou.log")
        cmd = [verrou_bin, "--tool=verrou", "--error-limit=no", f"--log-file={log_path}"]
        if rounding_mode:
            cmd.append(f"--rounding-mode={rounding_mode}")
        cmd.extend(extra_flags or [])
        cmd.append(sim_bin)

        with open(os.path.join(run_dir, "sim.out"), "w") as f:
            result = subprocess.run(cmd, cwd=tmpdir, env=_verrou_env(verrou_bin), stdout=f, stderr=subprocess.STDOUT, check=False)

        if result.returncode != 0:
            tag = rounding_mode or "vprec"
            raise MFCException(f"simulation ({tag}) exited {result.returncode}. See {run_dir}/sim.out")

        os.makedirs(run_dir, exist_ok=True)
        for fn in os.listdir(os.path.join(tmpdir, "D")):
            shutil.copy2(os.path.join(tmpdir, "D", fn), run_dir)


def _run_cancellation_check(verrou_bin: str, sim_bin: str, work_dir: str, threshold: int = 10):
    """Run --check-cancellation at the given bit threshold; return [(fname, line)]
    of MFC cancellation sites (subtractions losing >= `threshold` significant bits),
    or None if the run itself failed (distinct from [] = ran and found none)."""
    tag = f"cancellation_{threshold}"
    run_dir = os.path.join(work_dir, tag)
    os.makedirs(run_dir, exist_ok=True)
    gen_path = os.path.join(run_dir, "cancel_gen.txt")
    flags = [
        "--check-cancellation=yes",
        f"--cc-threshold-double={threshold}",
        f"--cc-gen-file={gen_path}",
    ]
    try:
        _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="nearest", extra_flags=flags)
    except MFCException as exc:
        cons.print(f"  [yellow]cancellation run (threshold {threshold}) failed: {exc}[/yellow]")
        return None
    return _parse_cancel_gen(gen_path)


def _run_float_max_check(verrou_bin: str, sim_bin: str, work_dir: str):
    """Run with --check-max-float=yes; return [(fname, line)] of overflow sites,
    or None if the run failed (distinct from [] = ran and found none)."""
    run_dir = os.path.join(work_dir, "float_max")
    os.makedirs(run_dir, exist_ok=True)
    try:
        _run_simulation_verrou(
            verrou_bin,
            sim_bin,
            work_dir,
            run_dir,
            rounding_mode="nearest",
            extra_flags=["--check-max-float=yes"],
        )
    except MFCException as exc:
        cons.print(f"  [yellow]float-max run failed: {exc}[/yellow]")
        return None
    return _parse_vg_error_locs(os.path.join(run_dir, "verrou.log"), "Max float")


def _run_float_proxy(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, ref_dir: str) -> float:
    """One run with --rounding-mode=float; returns L∞ deviation from nearest-ref."""
    run_dir = os.path.join(work_dir, "float_proxy")
    os.makedirs(run_dir)
    _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, rounding_mode="float")
    return _max_diff_np(ref_dir, run_dir, case["compare"])


def _run_vprec_sweep(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, ref_dir: str) -> list:
    """Run at each mantissa-bit level. Returns [(bits, dev), ...]."""
    results = []
    for bits in VPREC_MANTISSA_BITS:
        run_dir = os.path.join(work_dir, f"vprec_{bits}")
        os.makedirs(run_dir)
        flags = [
            "--backend=vprec",
            "--vprec-mode=full",
            f"--vprec-precision-binary64={bits}",
            "--vprec-range-binary64=11",
        ]
        try:
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, extra_flags=flags)
            dev = _max_diff_np(ref_dir, run_dir, case["compare"])
        except MFCException:
            dev = float("inf")
        results.append((bits, dev))
    return results
