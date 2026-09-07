#!/bin/bash
# Build and benchmark BOTH the master and PR trees inside ONE Phoenix SLURM job,
# on the same node and the same GPUs, so PR-vs-master is measured on identical
# hardware. Builds run first (compile only, no GPU contention); the two benchmark
# runs are then executed strictly back-to-back so they never share the GPUs.
#
# This replaces the old two-job model (one SLURM job per tree, both pinned to one
# partition via BENCH_GPU_PARTITION). That model needed two idle nodes in the
# SAME partition at once -- the main source of benchmark queue starvation. A
# single node needs only -G2, so the job submits to a partition LIST (see
# submit-slurm-job.sh) and backfills onto whichever partition frees first.
#
# Same-node is a strictly stronger fairness guarantee than the old "same GPU
# type": the two runs share the exact silicon, and by building first and
# benchmarking back-to-back they also run adjacently in time.
#
# The build/bench commands mirror .github/workflows/common/bench.sh. They are
# duplicated rather than folded in behind a phase flag because the pair
# orchestration (two trees, one node probe, shared GPUs) does not map cleanly
# onto that single-tree script.
#
# Runs inside the SLURM allocation, launched from the PR tree; master is a
# sibling at ../master. The module environment is already loaded by
# submit-slurm-job.sh (modules are cluster-level, so both trees build against
# them). Expects env: job_device, job_interface, job_cluster.

set -euo pipefail

pr_dir="$(pwd)"
master_dir="$(cd "${pr_dir}/../master" && pwd)"

# Output name the downstream comparison expects in each tree (see bench.yml
# "Generate & Post Comment": bench-<device>-<interface>.yaml).
bench_yaml="bench-${job_device}-${job_interface}.yaml"

# Cap parallel compile jobs (see common/bench.sh: GNR nodes have 192 cores but
# nproc is too aggressive for the build).
n_jobs=$(( $(nproc) > 64 ? 64 : $(nproc) ))

# $gpu_opts (e.g. "--gpu acc") from $job_device/$job_interface.
source "${pr_dir}/.github/scripts/gpu-opts.sh"

tmpbuild=/storage/project/r-sbryngelson3-0/sbryngelson3/mytmp_build

# Per-tree scratch dir, unique to this job so concurrent matrix jobs don't race.
tree_tmpdir() { echo "${tmpbuild}/run-$(basename "$1")-${SLURM_JOB_ID:-$$}"; }

build_tree() {                                  # <dir>
    local dir="$1"
    echo "===================================================="
    echo "BUILD: $dir"
    echo "===================================================="
    ( cd "$dir"
      export TMPDIR="$(tree_tmpdir "$dir")"; mkdir -p "$TMPDIR"
      # Always nuke stale builds: Phoenix compute nodes are heterogeneous, so a
      # binary left by another node risks a SIGILL microarchitecture mismatch.
      source .github/scripts/clean-build.sh; clean_build
      source .github/scripts/retry-build.sh
      retry_build ./mfc.sh build -j "$n_jobs" $gpu_opts )
}

bench_tree() {                                  # <dir>
    local dir="$1"
    echo "===================================================="
    echo "BENCH: $dir"
    echo "===================================================="
    ( cd "$dir"
      export TMPDIR="$(tree_tmpdir "$dir")"; mkdir -p "$TMPDIR"
      # $ngpus / $gpu_ids from the allocation (a -G2 job sees 2 GPUs: "0 1").
      source .github/scripts/detect-gpus.sh
      ./mfc.sh bench --mem 4 -o "$bench_yaml" \
          -- -c phoenix-bench $gpu_opts -g $gpu_ids -n "$ngpus" )
}

# --- Build both trees (GPUs idle; order/timing here does not affect fairness) ---
build_tree "$master_dir"
build_tree "$pr_dir"

# --- Probe the node once before spending the allocation on the benchmarks ---
# (exit 77 => submit-slurm-job.sh excludes this node and resubmits elsewhere.)
preflight_rc=0
bash "${pr_dir}/.github/scripts/preflight.sh" "$job_cluster" "$job_device" || preflight_rc=$?
if [ "$preflight_rc" -ne 0 ]; then
    exit "$preflight_rc"
fi

# --- Benchmarks: strictly back-to-back on the same GPUs (fair, no contention) ---
bench_tree "$master_dir"
bench_tree "$pr_dir"

# Best-effort scratch cleanup for this job's dirs.
sleep 10
rm -rf "${tmpbuild}/run-"*"-${SLURM_JOB_ID:-$$}" 2>/dev/null || true

echo "bench-pair complete: master and PR benchmarked on $(hostname -s)"
