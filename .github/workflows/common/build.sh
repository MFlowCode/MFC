#!/bin/bash
# Build-only script for all clusters.
# Runs inside a SLURM job via submit-slurm-job.sh.
# Builds MFC without running tests (--dry-run).
# Expects env vars: $job_device, $job_interface, $job_shard, $job_cluster

set -euo pipefail

source .github/scripts/gpu-opts.sh
build_opts="$gpu_opts"

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
# Phoenix: always start fresh to avoid SIGILL from stale binaries compiled
# on a different microarchitecture.
if [ "$job_cluster" = "phoenix" ]; then
    source .github/scripts/clean-build.sh
    clean_build
fi

source .github/scripts/retry-build.sh

# Phoenix: smoke-test the syscheck binary to catch architecture mismatches
# (SIGILL from binaries compiled on a different compute node).
#
# Launch under mpirun rather than directly. Phoenix's openmpi/4.1.5 predates the
# PMIx 4.2.6 that came with the Slurm 26.05.2 upgrade, so an MPI binary started
# bare inside an allocation misreads the PMIx environment as an srun launch and
# aborts in MPI_Init ("OPAL ERROR: Unreachable in file ext3x_client.c"). mpirun
# is unaffected, and it is how MFC launches every binary anyway. Output is left
# on stdout so a future failure is diagnosable from the CI log.
validate_cmd=""
if [ "$job_cluster" = "phoenix" ]; then
    validate_cmd='syscheck_bin=$(find build/install -name syscheck -type f 2>/dev/null | head -1); [ -z "$syscheck_bin" ] || mpirun -np 1 "$syscheck_bin"'
fi

# --- Variant selection ---
# The suite needs two simulation binaries: the base one and a chemistry one (the
# mechanism is compiled in). On amdflang each device link is ~1 hour, so building
# both serially exceeds the 2 h walltime. $job_variant lets the caller split them
# into concurrent jobs that write to disjoint staging directories:
#   base -> plain targets, what every non-chemistry test runs
#   chem -> the chemistry variant only
# Unset builds everything in this job (default for every other cluster).
case "${job_variant:-}" in
    base) build_cmd=(./mfc.sh build -j 8 $build_opts) ;;
    chem) build_cmd=(./mfc.sh test -v --dry-run -a -j 8 -o Chemistry $build_opts) ;;
    "")   build_cmd=(./mfc.sh test -v --dry-run -a -j 8 $build_opts) ;;
    *)    echo "ERROR: unknown job_variant '$job_variant'"; exit 1 ;;
esac

# Run a build step with its output teed, and classify a failure before giving up.
# Every ./mfc.sh call needs this, not just the solver build: the *first* one in
# the job is what bootstraps build/venv from PyPI, so it is the one that sees a
# package-index outage.
log_base="build-${job_slug:-${job_device}-${job_interface}}"

run_build_step() {
    local log="$1"
    shift
    set +e
    "$@" 2>&1 | tee "$log"
    local rc=${PIPESTATUS[0]}
    set -e
    if [ "$rc" -ne 0 ]; then
        exit "$rc"
    fi
}

# --- Probe this node before committing the solver build to it ---
# syscheck is a standalone target that links in 5-19 seconds, and it is already
# built second in the ordinary build order. Building it on its own first and
# running it here rejects an unusable node in about a minute, rather than after
# the ~40 minute solver build that used to precede the first GPU touch. In the
# Aug 2026 ECC failures that gap was a median of 38 minutes per job.
#
# Through retry_build so a transient blip still gets its nuke-and-retry, which a
# bare invocation here would have quietly dropped.
run_build_step "${log_base}-syscheck.log" retry_build ./mfc.sh build -t syscheck -j 8 $build_opts

preflight_rc=0
bash .github/scripts/preflight.sh "$job_cluster" "$job_device" || preflight_rc=$?
if [ "$preflight_rc" -ne 0 ]; then
    # 77 (bad node) travels back to submit-slurm-job.sh as this job's exit code,
    # which decides whether to requeue elsewhere.
    exit "$preflight_rc"
fi

# --- Solver build ---
# Output is teed so a failure can be classified afterwards. Some Frontier CCE
# and amdflang failures emit no compiler diagnostic at all, so whatever the
# build did print is the only evidence there is.
RETRY_VALIDATE_CMD="$validate_cmd" \
    run_build_step "${log_base}.log" retry_build "${build_cmd[@]}"
