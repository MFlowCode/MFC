#!/bin/bash
# Shared configuration for Phoenix GitHub Actions runner management scripts.
#
# Sources this file to get:
#   NODES         — array of physical login node hostnames
#   CGROUP_LIMIT  — per-user memory limit in MB
#   discover_runners() — populates RUNNER_DIRS, RUNNER_NAMES, RUNNER_POOLS arrays
#   find_runner_node() — returns which node a runner is currently running on

# Physical login nodes (gnr-1 = login-1 = login-4, etc.)
NODES=(login-phoenix-gnr-1 login-phoenix-gnr-2 login-phoenix-gnr-3)

# Per-user cgroup memory limit on Phoenix login nodes (MB)
CGROUP_LIMIT=4096

# Parent directories that contain runner installations.
# Each may have actions-runner-N/ subdirectories with a .runner config file.
RUNNER_PARENT_DIRS=(
    /storage/scratch1/6/sbryngelson3/mfc-runners
    /storage/project/r-sbryngelson3-0/sbryngelson3/mfc-runners-2
    /storage/scratch1/6/sbryngelson3/cfdnn-runners
)

# Discover all registered runners on the shared filesystem.
# Populates parallel arrays: RUNNER_DIRS, RUNNER_NAMES, RUNNER_POOLS, RUNNER_ORGS
declare -a RUNNER_DIRS=()
declare -a RUNNER_NAMES=()
declare -a RUNNER_POOLS=()
declare -a RUNNER_ORGS=()

discover_runners() {
    RUNNER_DIRS=()
    RUNNER_NAMES=()
    RUNNER_POOLS=()
    RUNNER_ORGS=()

    for parent in "${RUNNER_PARENT_DIRS[@]}"; do
        for conf in "$parent"/actions-runner-*/.runner; do
            [ -f "$conf" ] || continue
            local dir=$(dirname "$conf")
            local info
            info=$(python3 -c "
import json, sys
d = json.loads(open('$conf').read().lstrip('\ufeff'))
print(d.get('agentName', ''))
print(d.get('poolName', ''))
print(d.get('gitHubUrl', ''))
" 2>/dev/null)
            local name pool org
            name=$(echo "$info" | sed -n '1p')
            pool=$(echo "$info" | sed -n '2p')
            org=$(echo "$info" | sed -n '3p')

            RUNNER_DIRS+=("$dir")
            RUNNER_NAMES+=("${name:-$(basename "$dir")}")
            RUNNER_POOLS+=("${pool:-unknown}")
            RUNNER_ORGS+=("${org:-unknown}")
        done
    done
}

# Find which physical node a runner is currently running on.
# Args: $1 = runner directory
# Prints: node hostname, or "offline" if not found on any node.
# Note: uses CWD matching (readlink /proc/PID/cwd) because the runner
# command line is just "Runner.Listener run" without the full path.
find_runner_node() {
    local dir="$1"
    for node in "${NODES[@]}"; do
        local found
        found=$(ssh -o ConnectTimeout=5 "$node" '
            for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
                cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
                [ "$cwd" = "'"$dir"'" ] && echo "$p" && break
            done
        ' 2>/dev/null || true)
        if [ -n "$found" ]; then
            echo "$node"
            return
        fi
    done
    echo "offline"
}

# Find PIDs of a runner on a node by matching CWD.
# Args: $1 = node, $2 = runner directory
# Prints: space-separated PIDs, or empty if not found.
_find_runner_pids() {
    local node="$1" dir="$2"
    ssh -o ConnectTimeout=5 "$node" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            cwd=$(readlink -f /proc/$p/cwd 2>/dev/null || true)
            [ "$cwd" = "'"$dir"'" ] && echo -n "$p "
        done
    ' 2>/dev/null || true
}

# Start a runner on a specific node with proper PATH and detachment.
# Args: $1 = node, $2 = runner directory
# Returns 0 on success, 1 on failure.
start_runner() {
    local node="$1" dir="$2"
    ssh -o ConnectTimeout=5 "$node" \
        "setsid bash -lc 'cd $dir && nohup ./run.sh >> runner-nohup.log 2>&1 &'" \
        </dev/null 2>/dev/null &
    local ssh_pid=$!
    # Wait up to 10s for SSH to exit
    local i
    for i in $(seq 1 10); do
        kill -0 $ssh_pid 2>/dev/null || break
        sleep 1
    done
    kill $ssh_pid 2>/dev/null || true
    wait $ssh_pid 2>/dev/null || true
    sleep 3
    # Verify
    local new_pid
    new_pid=$(_find_runner_pids "$node" "$dir")
    [ -n "$new_pid" ]
}

# Stop a runner on a specific node.
# Args: $1 = node, $2 = runner directory
stop_runner() {
    local node="$1" dir="$2"
    local pids
    pids=$(_find_runner_pids "$node" "$dir")
    if [ -z "$pids" ]; then
        return 0
    fi
    for pid in $pids; do
        ssh -o ConnectTimeout=5 "$node" "kill $pid" 2>/dev/null || true
    done
    sleep 3
    # Force kill survivors
    pids=$(_find_runner_pids "$node" "$dir")
    for pid in $pids; do
        ssh -o ConnectTimeout=5 "$node" "kill -9 $pid" 2>/dev/null || true
    done
    sleep 1
}

# Check if a runner has slurm in its PATH.
# Args: $1 = node, $2 = PID (whitespace trimmed)
# Returns 0 if slurm is in PATH, 1 otherwise.
check_slurm_path() {
    local node="$1" pid="${2// /}"
    local count
    count=$(ssh -o ConnectTimeout=5 "$node" \
        "cat /proc/$pid/environ 2>/dev/null | tr '\0' '\n' | grep -c /opt/slurm" \
        2>/dev/null || echo 0)
    [ "$count" -gt 0 ]
}
