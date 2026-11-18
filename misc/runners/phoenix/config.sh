#!/usr/bin/env bash
# Shared configuration for Phoenix GitHub Actions runner management.
#
# Sourced by all other scripts. Provides Phoenix constants and
# site-specific functions. Common functions live in ../common/runner-lib.sh.

# --- Phoenix constants ---
ORG="MFlowCode"
RUNNER_GROUP="phoenix"
RUNNER_LABEL="gt"
NODES=(login-phoenix-gnr-1 login-phoenix-gnr-2 login-phoenix-gnr-3)
CGROUP_LIMIT=4096  # per-user memory limit in MB on login nodes

SSH_OPTS="-o ConnectTimeout=5"

# Parent directories containing actions-runner-*/ installations on shared storage.
RUNNER_PARENT_DIRS=(
    /storage/project/r-sbryngelson3-0/sbryngelson3/mfc-runners-2
    /storage/cedar/cedar0/cedarp-sbryngelson3-0/mfc-runners-3
)

# Where create-runner installs new runners. Must be one of RUNNER_PARENT_DIRS,
# or find_runner_dirs() will not see them.
RUNNER_NEW_DIR=/storage/cedar/cedar0/cedarp-sbryngelson3-0/mfc-runners-3

source "$(dirname "${BASH_SOURCE[0]}")/../common/runner-lib.sh"

# --- Local filesystem ---

# No shared cache: each runner downloads its own tarball independently.
TARBALL_CACHE_DIR=""

# Return the directory where a named runner should be installed.
# Creates RUNNER_NEW_DIR/actions-runner-N, with N one past the highest number
# in use across all RUNNER_PARENT_DIRS so directory numbers stay unique.
# Args: $1 = runner name (unused; directory is numbered, not named), $2 = optional override dir
runner_install_dir() {
    local override="${2:-}"
    [ -n "$override" ] && echo "$override" && return
    local existing next_num
    existing=$(for parent in "${RUNNER_PARENT_DIRS[@]}"; do
                   ls -d "$parent"/actions-runner-* 2>/dev/null
               done | sed 's/.*actions-runner-//' | sort -n | tail -1)
    next_num=$(( ${existing:-0} + 1 ))
    echo "$RUNNER_NEW_DIR/actions-runner-$next_num"
}

# Find all runner directories on shared storage.
# Prints: one directory path per line.
find_runner_dirs() {
    for parent in "${RUNNER_PARENT_DIRS[@]}"; do
        for conf in "$parent"/actions-runner-*/.runner; do
            [ -f "$conf" ] && dirname "$conf"
        done
    done
}
