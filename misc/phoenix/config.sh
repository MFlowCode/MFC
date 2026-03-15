#!/bin/bash
# Shared configuration for Phoenix GitHub Actions runner management.
#
# Sourced by all other scripts. Provides constants, GitHub API helpers,
# and login-node process management functions.

# --- Phoenix constants ---
ORG="MFlowCode"
RUNNER_GROUP="phoenix"
RUNNER_LABEL="gt"
NODES=(login-phoenix-gnr-1 login-phoenix-gnr-2 login-phoenix-gnr-3)
CGROUP_LIMIT=4096  # per-user memory limit in MB on login nodes

# Parent directories containing actions-runner-*/ installations on shared storage.
RUNNER_PARENT_DIRS=(
    /storage/scratch1/6/sbryngelson3/mfc-runners
    /storage/project/r-sbryngelson3-0/sbryngelson3/mfc-runners-2
)

# --- GitHub API ---

# List Phoenix runners from the GitHub API.
# Prints: id name status busy (one runner per line)
gh_list_runners() {
    gh api "orgs/$ORG/actions/runners" --paginate \
        --jq ".runners[]
              | select(.labels | map(.name) | index(\"$RUNNER_LABEL\"))
              | \"\(.id) \(.name) \(.status) \(.busy)\""
}

# Get a registration token for new runners.
gh_registration_token() {
    gh api "orgs/$ORG/actions/runners/registration-token" --jq .token
}

# Get the latest runner binary version.
gh_latest_runner_version() {
    gh api repos/actions/runner/releases/latest --jq '.tag_name | ltrimstr("v")'
}

# Remove a runner registration from GitHub.
# Args: $1 = runner ID (numeric, from API)
gh_remove_runner() {
    gh api "orgs/$ORG/actions/runners/$1" -X DELETE
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

# Get the GitHub runner name from a .runner config file.
# Args: $1 = runner directory
get_runner_name() {
    python3 -c "
import json
d = json.loads(open('$1/.runner').read().lstrip('\ufeff'))
print(d.get('agentName', ''))
" 2>/dev/null
}

# --- Login-node process management ---

# Find PIDs of a runner on a node by matching its CWD.
# (Runner.Listener's command line is just "Runner.Listener run" — no path.)
# Args: $1 = node, $2 = runner directory
# Prints: space-separated PIDs, or empty.
find_pids() {
    ssh -o ConnectTimeout=5 "$1" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            [ "$cwd" = "'"$2"'" ] && echo -n "$p "
        done
    ' 2>/dev/null || true
}

# Find which login node a runner is on.
# Args: $1 = runner directory
# Prints: node hostname, or "offline".
find_node() {
    for node in "${NODES[@]}"; do
        [ -n "$(find_pids "$node" "$1")" ] && echo "$node" && return
    done
    echo "offline"
}

# Start a runner on a node with login shell (for /opt/slurm PATH).
# Args: $1 = node, $2 = runner directory
# Returns: 0 if running after start, 1 otherwise.
start_runner() {
    local node="$1" dir="$2"
    ssh -o ConnectTimeout=5 "$node" \
        "setsid bash -lc 'cd $dir && nohup ./run.sh >> runner-nohup.log 2>&1 &'" \
        </dev/null 2>/dev/null &
    local ssh_pid=$!
    local i; for i in $(seq 1 10); do
        kill -0 $ssh_pid 2>/dev/null || break; sleep 1
    done
    kill $ssh_pid 2>/dev/null || true
    wait $ssh_pid 2>/dev/null || true
    sleep 3
    [ -n "$(find_pids "$node" "$dir")" ]
}

# Stop a runner on a node (SIGTERM then SIGKILL).
# Args: $1 = node, $2 = runner directory
stop_runner() {
    local node="$1" dir="$2" pids
    pids=$(find_pids "$node" "$dir")
    [ -z "$pids" ] && return 0
    for pid in $pids; do
        ssh -o ConnectTimeout=5 "$node" "kill $pid" 2>/dev/null || true
    done
    sleep 3
    pids=$(find_pids "$node" "$dir")
    for pid in $pids; do
        ssh -o ConnectTimeout=5 "$node" "kill -9 $pid" 2>/dev/null || true
    done
    sleep 1
}

# Check if a runner process has /opt/slurm in PATH.
# Args: $1 = node, $2 = PID
has_slurm() {
    local count
    count=$(ssh -o ConnectTimeout=5 "$1" \
        "cat /proc/${2%% *}/environ 2>/dev/null | tr '\0' '\n' | grep -c /opt/slurm" \
        2>/dev/null || echo 0)
    [ "$count" -gt 0 ]
}
