#!/usr/bin/env bash
# Run the variable-viscosity Couette grid-refinement sweep.
#
# The first grid is a normal building run (the analytic density IC compiles into
# pre_process, so the build is case-specific); once it finishes, the remaining
# grids reuse that binary via --no-build and run concurrently, each in its own
# directory (runs/n<N>/) on a single rank so they cannot collide and never touch
# the MPI domain decomposition. Pass grid sizes as arguments; default is
# "32 64 96".
#
#   ./run_suite.sh                # 32 64 96
#   ./run_suite.sh 32 64 96 128   # add a finer grid (n=128 is viscous-CFL
#                                 # limited: ~221k steps, hours on one rank)
#
# After this finishes, analyze with:  python3 validate.py
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
GRIDS="${*:-32 64 96}"

pids=()
first=1
for N in $GRIDS; do
    d="$HERE/runs/n$N"
    rm -rf "$d"
    mkdir -p "$d"
    cp "$HERE/case.py" "$HERE/couette_config.py" "$d/"
    if [ "$first" -eq 1 ]; then
        first=0
        echo "running n=$N -> $d (builds MFC, then the other grids launch)"
        (
            cd "$ROOT"
            COUETTE_N="$N" ./mfc.sh run "$d/case.py" -n 1 >"$d/run.log" 2>&1
        )
        continue
    fi
    echo "launching n=$N -> $d"
    (
        cd "$ROOT"
        COUETTE_N="$N" ./mfc.sh run "$d/case.py" -n 1 --no-build >"$d/run.log" 2>&1
    ) &
    pids+=("$!")
    sleep 5  # stagger so the per-run input-file generation does not overlap
done

echo "waiting on ${#pids[@]} runs: ${pids[*]}"
rc=0
for p in "${pids[@]}"; do
    wait "$p" || rc=1
done
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one grid failed; check runs/n*/run.log" >&2
    exit 1
fi
echo "all grids done"
