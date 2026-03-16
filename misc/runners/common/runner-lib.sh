#!/usr/bin/env bash
# Shared GitHub Actions runner management library.
#
# Sourced by site-specific config.sh files (misc/frontier/config.sh,
# misc/phoenix/config.sh). Callers must define ORG, NODES, and SSH_OPTS
# before sourcing this file.

# Default: no cgroup memory limit displayed. Override in site config (e.g. CGROUP_LIMIT=4096).
CGROUP_LIMIT=${CGROUP_LIMIT:-0}

# --- GitHub API ---

# List runners from the GitHub API, filtered to this site's label.
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

# Get the GitHub runner name from a .runner config file.
# Uses sys.argv to avoid path injection into Python source code.
# Prints the agentName, or empty string if the file is missing or unparsable
# (with a warning to stderr).
# Args: $1 = runner directory
get_runner_name() {
    python3 -c "
import json, sys
d = json.loads(open(sys.argv[1] + '/.runner').read().lstrip('\ufeff'))
print(d.get('agentName', ''))
" "$1" 2>/dev/null \
        || echo "WARNING: could not read runner name from '$1/.runner'" >&2
}

# --- Login-node process management ---

# Find PIDs of a runner on a node by matching its executable path.
# Candidate PIDs are found by grepping ps for Runner.Listener; each
# candidate's /proc/$p/exe is then resolved and matched against
# $dir/bin/Runner.Listener to confirm identity. This makes the confirmation
# step independent of CWD or process arguments.
# Output is filtered to numeric lines only to strip SSH MOTD noise.
# Args: $1 = node, $2 = runner directory
# Prints: space-separated PIDs, or empty.
find_pids() {
    ssh $SSH_OPTS "$1" '
        for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
            exe=$(readlink -f /proc/$p/exe 2>/dev/null || true)
            [ "$exe" = "'"$2"'/bin/Runner.Listener" ] && echo "$p"
        done
    ' 2>/dev/null | grep -E '^[0-9]+$' | tr '\n' ' ' || true
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

# Check if a runner process has a slurm directory in its PATH.
# Works across sites regardless of the specific slurm installation path.
# Returns non-zero if slurm is absent OR if the SSH check itself fails
# (callers should treat non-zero as "could not confirm slurm present").
# Args: $1 = node, $2 = PID
has_slurm() {
    local node="$1" pid="$2"
    ssh $SSH_OPTS "$node" \
        "tr '\0' '\n' < /proc/$pid/environ 2>/dev/null | grep -q '^PATH=.*slurm'" \
        2>/dev/null
}

# Sweep all nodes in parallel, writing per-node result files to tmpdir.
# Each output line: RUNNER <node> <dir> <rss_mb> <slurm_ok>
#   dir      = runner directory derived from the Runner.Listener exe path
#   slurm_ok = "ok" if slurm appears in the process PATH, "MISSING" otherwise
# Warns to stderr for any node whose output file is empty (SSH likely failed).
# Caller must create tmpdir and parse the output files.
# Args: $1 = tmpdir
sweep_all_nodes() {
    local tmpdir="$1" node
    for node in "${NODES[@]}"; do
        ssh $SSH_OPTS "$node" '
            for p in $(ps aux | grep Runner.Listener | grep -v grep | awk "{print \$2}"); do
                exe=$(readlink -f /proc/$p/exe 2>/dev/null || true)
                [ -z "$exe" ] && continue
                dir=$(dirname "$(dirname "$exe")")
                rss=$(ps -p $p -o rss= 2>/dev/null | awk "{printf \"%.0f\", \$1/1024}" || echo 0)
                slurm=$(tr "\0" "\n" < /proc/$p/environ 2>/dev/null | grep -c "^PATH=.*slurm" || echo 0)
                [ "$slurm" -gt 0 ] && slurm_ok="ok" || slurm_ok="MISSING"
                echo "RUNNER '"$node"' $dir $rss $slurm_ok"
            done
        ' 2>/dev/null > "$tmpdir/$node.out" &
    done
    wait
    for node in "${NODES[@]}"; do
        if [ ! -s "$tmpdir/$node.out" ]; then
            echo "WARNING: no runner data from $node (SSH may have failed)" >&2
        fi
    done
}

# Start a runner on a node.
# Uses a login shell (bash -lc) so site PATH (e.g. SLURM) is available.
# SSH launch failures are logged to stderr but do not prevent the function
# from checking whether the process appeared (find_pids after 3s).
# Args: $1 = node, $2 = runner directory
# Returns: 0 if a Runner.Listener process is found after start, 1 otherwise.
start_runner() {
    local node="$1" dir="$2"
    timeout 15 ssh $SSH_OPTS "$node" \
        "cd $dir && setsid bash -lc 'nohup ./run.sh >> runner.log 2>&1 < /dev/null &'" \
        </dev/null 2>/dev/null \
        || echo "WARNING: SSH launch to $node failed; checking for process anyway..." >&2
    sleep 3
    [ -n "$(find_pids "$node" "$dir")" ]
}

# Stop a runner on a node (SIGTERM, 3s grace, then SIGKILL).
# After SIGKILL, waits 1s then re-checks with find_pids. If the process is
# still alive (e.g. SSH kill failed, wrong UID, kernel stuck), emits a warning
# to stderr and returns 1 so callers can react. Returns 0 if the runner is
# confirmed stopped or was never running.
# Args: $1 = node, $2 = runner directory
stop_runner() {
    local node="$1" dir="$2" pids
    pids=$(find_pids "$node" "$dir")
    [ -z "$pids" ] && return 0
    for pid in $pids; do
        ssh $SSH_OPTS "$node" "kill $pid" 2>/dev/null || true
    done
    sleep 3
    pids=$(find_pids "$node" "$dir")
    for pid in $pids; do
        ssh $SSH_OPTS "$node" "kill -9 $pid" 2>/dev/null || true
    done
    sleep 1
    pids=$(find_pids "$node" "$dir")
    if [ -n "$pids" ]; then
        echo "WARNING: process(es) $pids on $node survived SIGKILL; runner may still be running" >&2
        return 1
    fi
}
