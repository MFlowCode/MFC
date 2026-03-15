#!/bin/bash
# Create and register a new GitHub Actions runner for Phoenix.
#
# Downloads the runner binary, configures it with a registration token,
# and starts it on the specified login node with proper PATH.
#
# Prerequisites:
#   - gh CLI authenticated with admin access to the target org
#   - The parent directory must exist and be on shared storage
#
# Usage: bash create-runner.sh <runner-name> <node> <parent-dir> [org] [runner-group]
#
# Examples:
#   bash create-runner.sh phoenix-11 login-phoenix-gnr-2 /storage/scratch1/6/sbryngelson3/mfc-runners
#   bash create-runner.sh phoenix-12 login-phoenix-gnr-3 /storage/project/.../mfc-runners-2 MFlowCode phoenix

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <runner-name> <node> <parent-dir> [org] [runner-group]"
    echo ""
    echo "  runner-name   Name for the runner (e.g. phoenix-11)"
    echo "  node          Login node to run on (e.g. login-phoenix-gnr-2)"
    echo "  parent-dir    Parent directory for the runner installation"
    echo "  org           GitHub org (default: MFlowCode)"
    echo "  runner-group  Runner group/pool (default: phoenix)"
    exit 1
fi

runner_name="$1"
node="$2"
parent_dir="$3"
org="${4:-MFlowCode}"
runner_group="${5:-phoenix}"

# Determine next available runner directory
existing=$(ls -d "$parent_dir"/actions-runner-* 2>/dev/null | sed 's/.*actions-runner-//' | sort -n | tail -1)
next_num=$(( ${existing:-0} + 1 ))
runner_dir="$parent_dir/actions-runner-$next_num"

echo "=== Creating runner ==="
echo "  Name:      $runner_name"
echo "  Node:      $node"
echo "  Directory: $runner_dir"
echo "  Org:       $org"
echo "  Group:     $runner_group"
echo ""

if [ -d "$runner_dir" ]; then
    echo "ERROR: Directory already exists: $runner_dir"
    exit 1
fi

# Get registration token
echo "Getting registration token from GitHub..."
token=$(gh api "orgs/$org/actions/runners/registration-token" --jq .token 2>/dev/null)
if [ -z "$token" ]; then
    echo "ERROR: Failed to get registration token. Check 'gh auth status' and org admin permissions."
    exit 1
fi
echo "  Token acquired."

# Get latest runner version
echo "Downloading runner..."
latest_version=$(gh api repos/actions/runner/releases/latest --jq .tag_name 2>/dev/null | sed 's/^v//')
if [ -z "$latest_version" ]; then
    echo "ERROR: Failed to determine latest runner version."
    exit 1
fi
runner_url="https://github.com/actions/runner/releases/download/v${latest_version}/actions-runner-linux-x64-${latest_version}.tar.gz"
echo "  Version: $latest_version"

mkdir -p "$runner_dir"
cd "$runner_dir"

curl -sL "$runner_url" | tar xz
echo "  Downloaded and extracted."

# Configure
echo "Configuring runner..."
./config.sh \
    --url "https://github.com/$org" \
    --token "$token" \
    --name "$runner_name" \
    --runnergroup "$runner_group" \
    --labels "gt" \
    --work "_work" \
    --unattended \
    --replace
echo "  Configured."

# Start on the target node
echo "Starting runner on $node..."
if start_runner "$node" "$runner_dir"; then
    pid=$(ssh -o ConnectTimeout=5 "$node" "pgrep -f 'Runner.Listener.*$runner_dir'" 2>/dev/null || true)
    if check_slurm_path "$node" "$pid"; then
        echo "  OK: Running as PID $pid on $node, slurm in PATH"
    else
        echo "  WARNING: Running as PID $pid but slurm NOT in PATH"
    fi
else
    echo "  ERROR: Failed to start. Try manually:"
    echo "  ssh $node 'cd $runner_dir && setsid bash -lc \"nohup ./run.sh >> runner-nohup.log 2>&1 &\"'"
fi

echo ""
echo "Runner $runner_name created at $runner_dir"
echo "Verify with: bash $SCRIPT_DIR/list-runners.sh"
