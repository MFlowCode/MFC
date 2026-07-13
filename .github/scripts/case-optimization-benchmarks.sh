#!/bin/bash

# Single source of truth for the benchmark cases exercised by the
# case-optimization CI. The pre-build (compiles the case-optimized binaries)
# and the run (executes them) MUST iterate the same list in the same order:
# sharding partitions this list by index, so any drift silently pre-builds or
# runs the wrong cases. Sourced by prebuild-case-optimization.sh and
# run_case_optimization.sh.
benchmarks=(
    benchmarks/5eq_rk3_weno3_hllc/case.py
    benchmarks/viscous_weno5_sgb_acoustic/case.py
    benchmarks/hypo_hll/case.py
    benchmarks/ibm/case.py
    benchmarks/igr/case.py
)

# Parse an optional "$job_shard" of the form "i/N" (e.g. "2/3"). On success
# sets caseopt_shard_idx / caseopt_shard_count (both 1 when unset — a single
# shard covering every case). Aborts on a malformed value.
caseopt_parse_shard() {
    caseopt_shard_idx=1
    caseopt_shard_count=1
    [ -z "${job_shard:-}" ] && return 0
    caseopt_shard_idx="${job_shard%%/*}"
    caseopt_shard_count="${job_shard##*/}"
    case "$caseopt_shard_idx"   in ''|*[!0-9]*|0*) echo "ERROR: bad shard '$job_shard' (expected i/N)"; exit 1 ;; esac
    case "$caseopt_shard_count" in ''|*[!0-9]*|0*) echo "ERROR: bad shard '$job_shard' (expected i/N)"; exit 1 ;; esac
    if [ "$job_shard" != "$caseopt_shard_idx/$caseopt_shard_count" ] \
       || [ "$caseopt_shard_idx" -lt 1 ] || [ "$caseopt_shard_idx" -gt "$caseopt_shard_count" ]; then
        echo "ERROR: bad shard '$job_shard' (expected i/N with 1 <= i <= N)"; exit 1
    fi
}

# 0 (true) if the 1-based case index $1 belongs to this shard. Shard i owns
# every Nth case: indices i, i+N, i+2N, ...
caseopt_case_in_shard() {
    [ $((($1 - 1) % caseopt_shard_count)) -eq $((caseopt_shard_idx - 1)) ]
}
