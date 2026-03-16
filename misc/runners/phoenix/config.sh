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
    /storage/scratch1/6/sbryngelson3/mfc-runners
    /storage/project/r-sbryngelson3-0/sbryngelson3/mfc-runners-2
)

source "$(dirname "${BASH_SOURCE[0]}")/../common/runner-lib.sh"

# --- Local filesystem ---

# No shared cache: each runner downloads its own tarball independently.
TARBALL_CACHE_DIR=""

# Return the directory where a named runner should be installed.
# Auto-increments the actions-runner-N suffix within RUNNER_PARENT_DIRS[0].
# Args: $1 = runner name (unused; directory is numbered, not named), $2 = optional override dir
runner_install_dir() {
    local override="${2:-}"
    [ -n "$override" ] && echo "$override" && return
    local parent="${RUNNER_PARENT_DIRS[0]}"
    local existing next_num
    existing=$(ls -d "$parent"/actions-runner-* 2>/dev/null | sed 's/.*actions-runner-//' | sort -n | tail -1)
    next_num=$(( ${existing:-0} + 1 ))
    echo "$parent/actions-runner-$next_num"
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
