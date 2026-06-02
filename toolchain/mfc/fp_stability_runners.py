"""Verrou subprocess runners for the FP-stability suite.

Each routine drives the verrou/valgrind binary (or the verrou_dd_* delta-debug
tools) and returns parsed results.  Pure parsing / metric helpers live in
fp_stability_metrics, which this module imports.
"""

import glob
import math
import os
import shutil
import stat
import subprocess
import tempfile
import textwrap

from .common import MFC_ROOT_DIR, MFCException
from .fp_stability_metrics import (
    _DD_FALLBACK_THRESHOLD,
    VPREC_MANTISSA_BITS,
    _build_source_filter,
    _confirm_decision,
    _is_arithmetic_loc,
    _max_abs_np,
    _max_diff_np,
    _parse_cancel_gen,
    _parse_rddmin_locs,
    _parse_rddmin_syms,
    _parse_vg_error_locs,
    _rank_locs,
)
from .printer import cons


def _find_verrou() -> str:
    verrou_home = os.environ.get("VERROU_HOME", os.path.join(os.path.expanduser("~"), ".local", "verrou"))
    candidate = os.path.join(verrou_home, "bin", "valgrind")
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return shutil.which("valgrind") or ""


def _find_binary(name: str) -> str:
    install_dir = os.path.join(MFC_ROOT_DIR, "build", "install")
    candidates = glob.glob(os.path.join(install_dir, "*", "bin", name))
    return max(candidates, key=os.path.getmtime) if candidates else ""


def _find_dd_tool(verrou_bin: str, tool: str) -> str:
    """Path to a verrou_dd_* tool (e.g. 'verrou_dd_sym') next to the verrou binary,
    or '' if absent."""
    c = os.path.join(os.path.dirname(verrou_bin), tool)
    return c if os.path.isfile(c) else ""


def _verrou_pythonpath(verrou_bin: str) -> str:
    """Path that must be on PYTHONPATH for verrou_dd_* imports (valgrind/ subdir)."""
    verrou_home = os.path.dirname(os.path.dirname(verrou_bin))
    matches = glob.glob(os.path.join(verrou_home, "lib", "python*", "site-packages", "valgrind"))
    return matches[0] if matches else ""


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
            result = subprocess.run(cmd, cwd=tmpdir, stdout=f, stderr=subprocess.STDOUT, check=False)

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
    raw = _parse_cancel_gen(gen_path)
    filtered = [(f, ln) for f, ln in raw if _is_arithmetic_loc(f, ln, ln)]
    skipped = len(raw) - len(filtered)
    if skipped and threshold == 10:
        cons.print(f"  [dim]cancellation: filtered {skipped} control-flow boundary site(s)[/dim]")
    return filtered


def _run_mca_samples(
    case: dict,
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    ref_dir: str,
    n_mca: int,
) -> tuple:
    """Run N mcaquad samples; return (max_dev, sig_bits_lower_bound, n_ok) where
    n_ok is how many samples actually completed (0 => no usable measurement)."""
    compare = case["compare"]
    ref_scale = _max_abs_np(ref_dir, compare)
    max_dev = 0.0
    n_ok = 0
    flags = ["--backend=mcaquad", "--mca-mode=mca"]
    for i in range(n_mca):
        run_dir = os.path.join(work_dir, f"mca_{i:02d}")
        os.makedirs(run_dir, exist_ok=True)
        try:
            _run_simulation_verrou(verrou_bin, sim_bin, work_dir, run_dir, extra_flags=flags)
            max_dev = max(max_dev, _max_diff_np(ref_dir, run_dir, compare))
            n_ok += 1
        except MFCException as exc:
            cons.print(f"  [dim]MCA sample {i} failed: {exc}[/dim]")
    sig_bits = None
    if n_ok and max_dev > 0.0 and ref_scale > 0.0:
        sig_bits = max(0, int(math.floor(-math.log2(max_dev / ref_scale))))
    return max_dev, sig_bits, n_ok


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


def _write_dd_run_sh(path: str, verrou_bin: str, sim_bin: str, ic_dir: str):
    """Generate dd_run.sh for verrou_dd_sym / verrou_dd_line.

    verrou_dd_* calls: dd_run.sh RUNDIR and injects function/line exclusion via
    VERROU_EXCLUDE / VERROU_SOURCE environment variables.  For test runs, we use
    --rounding-mode=float (deterministic, same deviation every call, --nruns=1 suffices).
    For the reference run, verrou_dd_sym sets VERROU_ROUNDING_MODE=nearest in the
    environment — we honour that so the reference is a stable nearest-rounding baseline
    to compare against.  CLI --rounding-mode would override the env var and break the
    reference, so we pass the mode via ${VERROU_ROUNDING_MODE:-float} instead.
    """
    content = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        # Generated by mfc.sh fp-stability — do not edit by hand.
        VERROU_BIN={verrou_bin!r}
        SIM_BIN={sim_bin!r}
        IC_DIR={ic_dir!r}

        RUNDIR="$1"
        TMPDIR_RUN=$(mktemp -d)
        trap 'rm -rf "$TMPDIR_RUN"' EXIT

        cp -r "$IC_DIR/p_all" "$TMPDIR_RUN/p_all"
        cp "$IC_DIR/simulation.inp" "$TMPDIR_RUN/simulation.inp"
        for fname in indices.dat pre_time_data.dat io_time_data.dat; do
            [ -f "$IC_DIR/$fname" ] && cp "$IC_DIR/$fname" "$TMPDIR_RUN/"
        done
        mkdir -p "$TMPDIR_RUN/D"

        # verrou_dd_sym sets VERROU_ROUNDING_MODE=nearest for its reference run and
        # leaves it unset for test runs.  Defaulting to float gives deterministic
        # test steps while letting the reference use nearest-rounding.
        ROUND="${{VERROU_ROUNDING_MODE:-float}}"

        # verrou_dd_sym injects VERROU_EXCLUDE (symbols to exclude from perturbation).
        # verrou_dd_line injects VERROU_SOURCE (source lines to restrict perturbation to).
        # Forward them as valgrind flags when set.
        EXTRA=""
        [ -n "${{VERROU_EXCLUDE:-}}" ] && EXTRA="$EXTRA --exclude=$VERROU_EXCLUDE"
        [ -n "${{VERROU_SOURCE:-}}" ]  && EXTRA="$EXTRA --source=$VERROU_SOURCE"

        cd "$TMPDIR_RUN"
        "$VERROU_BIN" --tool=verrou --error-limit=no --rounding-mode="$ROUND" $EXTRA "$SIM_BIN"
        rc=$?

        [ -d "$TMPDIR_RUN/D" ] && cp -a "$TMPDIR_RUN/D/." "$RUNDIR/"
        exit $rc
    """)
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _write_dd_cmp_py(path: str, compare_files: list, threshold: float):
    """Generate dd_cmp.py for verrou_dd_sym / verrou_dd_line.

    verrou_dd_* calls: dd_cmp.py REF_DIR RUN_DIR
    Exits 0 (stable) or 1 (unstable) based on threshold.
    """
    content = textwrap.dedent(f"""\
        #!/usr/bin/env python3
        # Generated by mfc.sh fp-stability — do not edit by hand.
        import sys, os, numpy as np

        COMPARE_FILES = {compare_files!r}
        THRESHOLD = {threshold!r}

        ref_dir, run_dir = sys.argv[1], sys.argv[2]
        max_dev = 0.0
        for fname in COMPARE_FILES:
            ref_p = os.path.join(ref_dir, fname)
            run_p = os.path.join(run_dir, fname)
            if not os.path.exists(ref_p) or not os.path.exists(run_p):
                print(f"MISSING: {{fname}}")
                sys.exit(1)
            ref = np.atleast_2d(np.loadtxt(ref_p))[:, 1]
            run = np.atleast_2d(np.loadtxt(run_p))[:, 1]
            dev = float(np.max(np.abs(ref - run)))
            max_dev = max(max_dev, dev)

        print(f"max_dev={{max_dev:.3e}}  threshold={{THRESHOLD:.0e}}")
        sys.exit(0 if max_dev <= THRESHOLD else 1)
    """)
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _dd_env(verrou_bin: str) -> dict:
    """Environment with PYTHONPATH set for verrou_dd_* imports."""
    py_pkg = _verrou_pythonpath(verrou_bin)
    env = os.environ.copy()
    if py_pkg:
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = ":".join(filter(None, [py_pkg, existing]))
    return env


def _run_dd_tool(
    dd_bin: str,
    dd_dir: str,
    dd_run_sh: str,
    dd_cmp_py: str,
    env: dict,
    log_name: str,
    summary_subdir: str,
    label: str,
) -> list:
    """Generic runner for verrou_dd_sym / verrou_dd_line. Returns raw summary lines."""
    log_file = os.path.join(dd_dir, log_name)
    cmd = [dd_bin, "--nruns=1", "--rddmin=d", "--reference-rounding=nearest", dd_run_sh, dd_cmp_py]
    cons.print(f"  [dim]running {label} (--nruns=1 float-mode --rddmin=d)...[/dim]")
    with open(log_file, "w") as f:
        result = subprocess.run(cmd, cwd=dd_dir, env=env, stdout=f, stderr=subprocess.STDOUT, check=False)
    summary_path = os.path.join(dd_dir, summary_subdir, "rddmin_summary")
    summary_lines = []
    if result.returncode == 0:
        if os.path.isfile(summary_path):
            with open(summary_path) as f:
                summary_lines = f.readlines()
            cons.print(f"  [bold yellow]{label} result[/bold yellow]:")
            for line in summary_lines:
                cons.print(f"    {line.rstrip()}")
        else:
            cons.print(f"  [dim]{label} done; see {log_file}[/dim]")
    else:
        cons.print(f"  [bold yellow]{label} exited {result.returncode}[/bold yellow] (see {log_file})")
    return summary_lines


def _setup_dd_run(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, dd_dir: str, threshold: float):
    """Write dd_run.sh and dd_cmp.py for a verrou_dd_* run into dd_dir; return their
    paths.  The threshold falls back to _DD_FALLBACK_THRESHOLD when unset."""
    os.makedirs(dd_dir, exist_ok=True)
    dd_run_sh = os.path.join(dd_dir, "dd_run.sh")
    dd_cmp_py = os.path.join(dd_dir, "dd_cmp.py")
    _write_dd_run_sh(dd_run_sh, verrou_bin, sim_bin, work_dir)
    _write_dd_cmp_py(dd_cmp_py, case["compare"], threshold if threshold is not None else _DD_FALLBACK_THRESHOLD)
    return dd_run_sh, dd_cmp_py


def _run_dd_sym(case: dict, verrou_bin: str, sim_bin: str, work_dir: str, log_dir: str, threshold: float = None) -> list:
    """Run verrou_dd_sym; return list of responsible symbol names."""
    dd_bin = _find_dd_tool(verrou_bin, "verrou_dd_sym")
    if not dd_bin:
        cons.print("  [dim]verrou_dd_sym not found; skipping delta-debug[/dim]")
        return []

    dd_dir = os.path.join(log_dir, case["name"])
    dd_run_sh, dd_cmp_py = _setup_dd_run(case, verrou_bin, sim_bin, work_dir, dd_dir, threshold)
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_sym.log", "dd.sym", "verrou_dd_sym")
    cons.print(f"  [dim]dd_sym logs: {dd_dir}[/dim]")
    return _parse_rddmin_syms(os.path.join(dd_dir, "dd.sym", "rddmin_summary"))


def _run_dd_line(
    case: dict,
    verrou_bin: str,
    sim_bin: str,
    work_dir: str,
    log_dir: str,
    threshold: float = None,
) -> list:
    """Run verrou_dd_line; return [{path, start, end, macro}] location dicts."""
    dd_bin = _find_dd_tool(verrou_bin, "verrou_dd_line")
    if not dd_bin:
        cons.print("  [dim]verrou_dd_line not found; skipping line-level debug[/dim]")
        return []

    dd_dir = os.path.join(log_dir, case["name"])
    dd_run_sh, dd_cmp_py = _setup_dd_run(case, verrou_bin, sim_bin, work_dir, dd_dir, threshold)
    _run_dd_tool(dd_bin, dd_dir, dd_run_sh, dd_cmp_py, _dd_env(verrou_bin), "dd_line.log", "dd.line", "verrou_dd_line")
    return _parse_rddmin_locs(os.path.join(dd_dir, "dd.line", "rddmin_summary"))


def _source_perturb_dev(verrou_bin, sim_bin, work_dir, ref_dir, conf_dir, src_lines, compare, tag):
    """Perturb only the lines in src_lines (deterministic float mode) and return
    the L-inf deviation from the nearest-rounding reference, or None on failure."""
    src_path = os.path.join(conf_dir, f"source_{tag}.txt")
    with open(src_path, "w") as fh:
        fh.writelines(src_lines)
    run_dir = os.path.join(conf_dir, f"perturb_{tag}")
    os.makedirs(run_dir, exist_ok=True)
    try:
        _run_simulation_verrou(
            verrou_bin,
            sim_bin,
            work_dir,
            run_dir,
            rounding_mode="float",
            extra_flags=[f"--source={src_path}"],
        )
    except MFCException:
        return None
    return _max_diff_np(ref_dir, run_dir, compare)


def _capture_gen_source(verrou_bin, sim_bin, work_dir, run_dir, gen_path):
    """Run nearest-rounding with --gen-source to capture the symbol-correct
    executed source lines (FILE\\tLINE\\tSYMBOL); return them, or None on failure."""
    try:
        _run_simulation_verrou(
            verrou_bin,
            sim_bin,
            work_dir,
            run_dir,
            rounding_mode="nearest",
            extra_flags=[f"--gen-source={gen_path}"],
        )
    except MFCException:
        return None
    if not os.path.isfile(gen_path):
        return None
    with open(gen_path) as fh:
        return fh.readlines()


def _run_confirmation(case, verrou_bin, sim_bin, work_dir, ref_dir, dd_line_locs, dd_threshold, float_proxy):
    """Positive control for dd_line: perturb ONLY the suspect lines and confirm
    the instability reproduces, then rank each line by its individual share.

    Verrou's --source matches file+line+symbol (not file+line alone), so we first
    capture the symbol-correct executed source lines via --gen-source, filter them
    to the suspect set, then run deterministic float-mode restricted to just those
    lines.  If the suspect-only deviation reaches dd_threshold the attribution is
    confirmed; if it stays near zero the reported lines do not actually carry the
    instability (e.g. a #:for-expanded line blamed for the wrong instance).

    Each line is then perturbed alone so its 'share_dev' (and 'share' of
    float_proxy) shows which computation dominates.

    Returns (confirmed, suspect_dev, ranked_locs).
    """
    if not dd_line_locs:
        return None, None, dd_line_locs
    conf_dir = os.path.join(work_dir, "confirm")
    os.makedirs(conf_dir, exist_ok=True)
    gen_lines = _capture_gen_source(verrou_bin, sim_bin, work_dir, conf_dir, os.path.join(conf_dir, "gen_source.txt"))
    if gen_lines is None:
        return None, None, dd_line_locs
    compare = case["compare"]

    # whole-set positive control
    suspects = [(loc["path"], loc["start"], loc["end"]) for loc in dd_line_locs]
    set_src = _build_source_filter(gen_lines, suspects)
    if not set_src:
        # none of the reported lines performs an instrumented FP op -> not reproduced
        return False, 0.0, dd_line_locs
    set_dev = _source_perturb_dev(verrou_bin, sim_bin, work_dir, ref_dir, conf_dir, set_src, compare, "set")
    confirmed = _confirm_decision(set_dev, dd_threshold)

    # per-line ranking (a single line trivially owns the whole set deviation)
    if len(dd_line_locs) == 1:
        dd_line_locs[0]["share_dev"] = set_dev
    else:
        for i, loc in enumerate(dd_line_locs):
            one = _build_source_filter(gen_lines, [(loc["path"], loc["start"], loc["end"])])
            loc["share_dev"] = _source_perturb_dev(verrou_bin, sim_bin, work_dir, ref_dir, conf_dir, one, compare, f"line{i:02d}") if one else 0.0
    ranked = _rank_locs(dd_line_locs, total=(float_proxy or set_dev))
    return confirmed, set_dev, ranked
