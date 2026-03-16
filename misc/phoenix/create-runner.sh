#!/usr/bin/env bash
# Create and register a new GitHub Actions runner on Phoenix.
#
# Downloads the runner binary, registers with MFlowCode org, and starts
# on the specified login node. Uses config.sh for org/group/label defaults.
#
# Prerequisites: gh CLI with admin:org scope (gh auth refresh -s admin:org)
#
# Usage: bash create-runner.sh <runner-name> <node> [parent-dir]
#
# Examples:
#   bash create-runner.sh phoenix-11 login-phoenix-gnr-2
#   bash create-runner.sh phoenix-12 login-phoenix-gnr-3 /storage/project/.../mfc-runners-2

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <runner-name> <node> [parent-dir]"
    echo ""
    echo "  runner-name  Name for the runner (e.g. phoenix-11)"
    echo "  node         Login node (${NODES[*]})"
    echo "  parent-dir   Parent directory (default: ${RUNNER_PARENT_DIRS[0]})"
    exit 1
fi

runner_name="$1"
node="$2"
parent_dir="${3:-${RUNNER_PARENT_DIRS[0]}}"

# Determine next available runner directory number
existing=$(ls -d "$parent_dir"/actions-runner-* 2>/dev/null | sed 's/.*actions-runner-//' | sort -n | tail -1)
next_num=$(( ${existing:-0} + 1 ))
runner_dir="$parent_dir/actions-runner-$next_num"

echo "=== Creating Phoenix runner ==="
echo "  Name:      $runner_name"
echo "  Node:      $node"
echo "  Directory: $runner_dir"
echo "  Org:       $ORG"
echo "  Group:     $RUNNER_GROUP"
echo "  Label:     $RUNNER_LABEL"
echo ""

if [ -d "$runner_dir" ]; then
    echo "ERROR: Directory already exists: $runner_dir"
    exit 1
fi

# Registration token
echo "Getting registration token..."
token=$(gh_registration_token)
if [ -z "$token" ]; then
    echo "ERROR: Failed to get token. Run: gh auth refresh -h github.com -s admin:org"
    exit 1
fi

# Download runner
echo "Downloading latest runner binary..."
version=$(gh_latest_runner_version)
url="https://github.com/actions/runner/releases/download/v${version}/actions-runner-linux-x64-${version}.tar.gz"
echo "  Version: $version"

mkdir -p "$runner_dir"
cd "$runner_dir"
tmp="runner-download.tmp.$$"
curl -fsSL "$url" -o "$tmp"
tar xz < "$tmp"
rm -f "$tmp"
echo "  Extracted."

# Configure
echo "Configuring..."
./config.sh \
    --url "https://github.com/$ORG" \
    --token "$token" \
    --name "$runner_name" \
    --runnergroup "$RUNNER_GROUP" \
    --labels "$RUNNER_LABEL" \
    --work "_work" \
    --unattended \
    --replace
echo "  Configured."

# Start
echo "Starting on $node..."
if start_runner "$node" "$runner_dir"; then
    echo "$node" > "$runner_dir/runner.node"
    pids=$(find_pids "$node" "$runner_dir")
    pid=${pids%% *}
    if has_slurm "$node" "$pid"; then
        echo "  OK: PID $pid, slurm in PATH"
    else
        echo "  WARNING: PID $pid but slurm MISSING from PATH"
    fi
else
    echo "  ERROR: Failed to start."
    echo "  Try: ssh $node 'cd $runner_dir && setsid bash -lc \"nohup ./run.sh >> runner.log 2>&1 < /dev/null &\"'"
fi

echo ""
echo "Created $runner_name at $runner_dir"
