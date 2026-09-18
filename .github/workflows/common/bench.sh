#!/bin/bash
# Unified benchmark script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Expects env vars: $job_device, $job_interface, $job_slug, $job_cluster

set -euo pipefail

source .github/scripts/bench-preamble.sh

# Cap parallel jobs at 64 to avoid overwhelming MPI daemons on large nodes
# (GNR nodes have 192 cores but nproc is too aggressive for build).
n_jobs=$(( $(nproc) > 64 ? 64 : $(nproc) ))

# --- Phoenix TMPDIR setup ---
if [ "$job_cluster" = "phoenix" ]; then
    tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build
    currentdir=$tmpbuild/run-$(( RANDOM % 9000 ))
    mkdir -p $tmpbuild
    mkdir -p $currentdir
    export TMPDIR=$currentdir
    trap 'rm -rf "$currentdir" || true' EXIT
fi

# --- Build ---
# Phoenix builds everything inside SLURM (no login-node build step).
# Frontier/Frontier AMD: deps already fetched on login node via --deps-only;
# source code is built here on the compute node.
# Phoenix: always nuke stale builds (heterogeneous compute nodes → ISA mismatch risk).
if [ "$job_cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
fi

source .github/scripts/retry-build.sh
retry_build ./mfc.sh build -j $n_jobs $build_opts || exit 1

# Probe this node before spending the allocation on it. Placed after the build,
# not in the sbatch template: these scripts nuke and rebuild build/ themselves
# (Phoenix does so precisely because its compute nodes are heterogeneous), so a
# probe running earlier would test a stale binary from a previous job -- likely
# built for another microarchitecture -- and a SIGILL there would be reported as
# a bad node, excluding a healthy one.
preflight_rc=0
bash .github/scripts/preflight.sh "$job_cluster" "$job_device" || preflight_rc=$?
if [ "$preflight_rc" -ne 0 ]; then
    exit "$preflight_rc"
fi

# --- Bench cluster flag ---
if [ "$job_cluster" = "phoenix" ]; then
    bench_cluster="phoenix-bench"
else
    bench_cluster="$job_cluster"
fi

# --- Frontier: keep Darshan out of the benchmark ---
# Since 2026-09-16 every Frontier bench case has taken ~3 min of solver time and
# ~17 min of wall time, and the job dies on the 2 h limit partway through the list.
# The gap is not in the solver: per-target exec times are identical to runs that
# passed before that date, and bench.py now prints how much of each case's wall
# time falls after the run printed its own End-time. That points at process exit,
# and Frontier preloads Darshan into every MPI job, which flushes its log in
# MPI_Finalize -- the same phase where its heatmap module asserted and killed a
# syscheck the same day. The benchmark has no use for I/O profiling, so turn it
# off here and let the timing say whether that was it.
if [ "$job_cluster" = "frontier" ] || [ "$job_cluster" = "frontier_amd" ]; then
    export DARSHAN_DISABLE=1
    echo "Darshan disabled for this benchmark (see bench.sh)."
fi

# --- Run benchmark ---
if [ "$job_device" = "gpu" ]; then
    ./mfc.sh bench --mem 4 -o "$job_slug.yaml" -- -c $bench_cluster $device_opts -n $n_ranks
else
    ./mfc.sh bench --mem 1 -o "$job_slug.yaml" -- -c $bench_cluster $device_opts -n $n_ranks
fi

# --- Phoenix cleanup (trap EXIT handles rm -rf "$currentdir") ---
if [ "$job_cluster" = "phoenix" ]; then
    sleep 10
    unset TMPDIR
fi
