#!/bin/bash
# Shared configuration for Phoenix GitHub Actions runner management.
#
# Sourced by all other scripts. Provides Phoenix constants, GitHub API
# helpers, and login-node process management functions.

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

# --- GitHub API ---

# List Phoenix runners from the GitHub API.
# Prints: id name status busy (one runner per line)
gh_list_runners() {
    gh api "orgs/$ORG/actions/runners" --paginate \
        --jq ".runners[]
              | select(.labels | map(.name) | index(\"$RUNNER_LABEL\"))
              | \"\(.id) \(.name) \(.status) \(.busy)\""
}

# --- Local filesystem ---

# Find all runner directories on shared storage.
# Prints: one directory path per line.
find_runner_dirs() {
    for parent in "${RUNNER_PARENT_DIRS[@]}"; do
        for conf in "$parent"/actions-runner-*/.runner; do
            [ -f "$conf" ] && dirname "$conf"
        done
    done
}

# --- Login-node process management (Phoenix-specific) ---

# Check if a runner process has /opt/slurm in PATH.
# Args: $1 = node, $2 = PID
has_slurm() {
    local count
    count=$(ssh $SSH_OPTS "$1" \
        "cat /proc/${2%% *}/environ 2>/dev/null | tr '\0' '\n' | grep -c /opt/slurm" \
        2>/dev/null || echo 0)
    [ "$count" -gt 0 ]
}
