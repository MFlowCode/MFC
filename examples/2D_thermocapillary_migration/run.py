#!/usr/bin/env python3
"""Run the thermocapillary validation cases and measure them.

Each grid variant runs in its own runs/<name>/ dir (MFC writes output next to the case
file), is measured by measure.py, and its result is collected into runs/<summary>.json.

Targets:
  fig5  case_Ma_0p001.py -- 2D zero-Ma rise, grid convergence (v/v_YGB -> ~0.80)
  fig7  case_Ma_20.py    -- 2D low-Ma migration (U*/U_r, peak ~0.13)
  all   fig5 + fig7

Usage (invokes mpirun -- run from a normal shell):
    python3 run.py <fig5|fig7|all> [run|remeasure]   (default: fig5 run)

  run        run each variant, then measure it
  remeasure  re-measure existing runs/ WITHOUT running any simulation

The canonical case files have no grid knob, so a variant's grid is set by copying the case
into its run dir and rewriting its `Nx = <n>` line; the committed case is never touched.
"""

import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")

# Per target: the case file, the measure.py mode, the summary file, and the grid variants to run.
# Each variant is (run-dir name, MPI ranks, Nx override or None for the case default).
TARGETS = {
    "fig5": dict(
        case="case_Ma_0p001.py",
        mode="fig5",
        summary="summary.json",
        variants=[
            ("tc1/w064", 6, 64),
            ("tc1/w128", 16, 128),
            ("tc1/w256", 64, 256),  # Ma=0.001 is conduction-dt-limited; w256 is very expensive
        ],
    ),
    "fig7": dict(
        case="case_Ma_20.py",
        mode="fig7",
        summary="fig7_summary.json",
        variants=[
            ("tc2/w064", 8, 64),
            ("tc2/w128", 16, 128),
        ],
    ),
}

# Headline numbers to show per mode in the end-of-run table.
TABLE_KEYS = {
    "fig5": ("cells_per_D", "ratio_plateau", "overshoot", "ratio_final"),
    "fig7": ("cells_per_D", "peak", "t_peak_tr", "terminal"),
}


def run_variant(case_file, name, ranks, nx):
    """Run one variant in a fresh runs/<name>/ dir. Returns True on success.

    The case is copied in; if nx is given, the copy's `Nx = <n>` line is rewritten. The
    committed case file is never touched.
    """
    wd = os.path.join(RUNS, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)

    dst = os.path.join(wd, os.path.basename(case_file))
    shutil.copy(os.path.join(HERE, case_file), dst)
    if nx is not None:
        text, hits = re.subn(r"(?m)^Nx = \d+", f"Nx = {nx}", open(dst).read(), count=1)
        if hits != 1:
            print(f"  WARNING: no `Nx = <int>` line in {name} -- using case default")
        open(dst, "w").write(text)

    print(f"\n>>> {name}: {os.path.basename(case_file)} ranks={ranks} Nx={nx or 'default'}", flush=True)
    rel = os.path.relpath(dst, REPO)
    p = subprocess.run(
        ["./mfc.sh", "run", rel, "-n", str(ranks), "--case-optimization"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}):")
        print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
        return False
    return True


def measure(wd, mode):
    """Run measure.py <wd> <mode> and return its RESULT_JSON dict (or None on failure)."""
    m = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), wd, mode], capture_output=True, text=True, check=False)
    print("\n".join(m.stdout.splitlines()[-4:]))
    if m.returncode != 0:
        print(f"  MEASURE FAILED: {m.stderr[-500:]}")
        return None
    for line in m.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def process_target(target, run_mode):
    """Run (or just re-measure) every variant of one target and update its summary file."""
    cfg = TARGETS[target]
    spath = os.path.join(RUNS, cfg["summary"])
    summary = json.load(open(spath)) if os.path.isfile(spath) else {}

    for name, ranks, nx in cfg["variants"]:
        wd = os.path.join(RUNS, name)
        if run_mode == "run":
            if not run_variant(cfg["case"], name, ranks, nx):
                continue
        elif not os.path.isfile(os.path.join(wd, "simulation.inp")):
            continue  # remeasure: skip variants that were never run
        res = measure(wd, cfg["mode"])
        if res is not None:
            summary[name] = res
            json.dump(summary, open(spath, "w"), indent=2)  # checkpoint after each variant

    keys = TABLE_KEYS[target]
    print(f"\n=== {target} summary ({len(summary)} runs) -> {spath} ===")
    for name, res in sorted(summary.items()):
        print(f"  {name:>16}: " + "  ".join(f"{k}={res[k]}" for k in keys if k in res))


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "fig5"
    run_mode = sys.argv[2] if len(sys.argv) > 2 else "run"
    targets = ["fig5", "fig7"] if target == "all" else [target]
    if any(t not in TARGETS for t in targets):
        sys.exit(f"unknown target {target!r}; choose from {', '.join(TARGETS)}, all")

    os.makedirs(RUNS, exist_ok=True)
    for t in targets:
        process_target(t, run_mode)


if __name__ == "__main__":
    main()
