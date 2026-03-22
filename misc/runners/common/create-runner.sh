#!/usr/bin/env bash
# Create, register, and start a GitHub Actions runner.
#
# Sourced by misc/runners/runner.sh after config is loaded.
# Config must define runner_install_dir() and may set TARBALL_CACHE_DIR.
#
# runner_install_dir <name> [override-dir]
#   Returns the directory where the runner should be installed.
#   If override-dir is given it is used directly; otherwise the site
#   computes the path (e.g. SHARED_DIR/<name> on Frontier, or an
#   auto-numbered actions-runner-N/ directory on Phoenix).
#
# TARBALL_CACHE_DIR
#   If non-empty, the runner tarball is cached here and reused across
#   installs (useful on Frontier where shared Lustre is visible from all
#   login nodes). If empty or unset, a fresh download is made for each
#   runner and the temporary file is removed after extraction.
#
# Usage: runner.sh <site> create-runner <name> <node> [install-dir]
#   name         Runner name (e.g. frontier-23, phoenix-11)
#   node         Login node to start the runner on
#   install-dir  Optional: override the computed installation directory
set -euo pipefail

RUNNER_NAME="${1:?Usage: create-runner <name> <node> [install-dir]}"
TARGET_NODE="${2:?Usage: create-runner <name> <node> [install-dir]}"
INSTALL_DIR_OVERRIDE="${3:-}"

RUNNER_DIR=$(runner_install_dir "$RUNNER_NAME" "$INSTALL_DIR_OVERRIDE")
RUNNER_VERSION="${RUNNER_VERSION:-$(gh_latest_runner_version 2>/dev/null || echo "2.332.0")}"
TARBALL="actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
TARBALL_URL="https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${TARBALL}"

echo "=== Creating runner ==="
echo "  Name:      $RUNNER_NAME"
echo "  Node:      $TARGET_NODE"
echo "  Directory: $RUNNER_DIR"
echo "  Org:       $ORG"
echo "  Group:     $RUNNER_GROUP"
echo "  Label:     $RUNNER_LABEL"
echo "  Version:   $RUNNER_VERSION"
echo ""

if [ -d "$RUNNER_DIR" ]; then
    echo "ERROR: Directory already exists: $RUNNER_DIR" >&2
    exit 1
fi

# --- Download tarball ---
if [ -n "${TARBALL_CACHE_DIR:-}" ]; then
    if [ ! -f "$TARBALL_CACHE_DIR/$TARBALL" ]; then
        echo "==> Downloading runner v${RUNNER_VERSION} to cache..."
        tmp="$TARBALL_CACHE_DIR/$TARBALL.tmp.$$"
        curl -fsSL "$TARBALL_URL" -o "$tmp"
        mv "$tmp" "$TARBALL_CACHE_DIR/$TARBALL"
    fi
    tarball_path="$TARBALL_CACHE_DIR/$TARBALL"
else
    echo "==> Downloading runner v${RUNNER_VERSION}..."
    mkdir -p "$RUNNER_DIR"
    tarball_path="$RUNNER_DIR/runner-download.tmp.$$"
    curl -fsSL "$TARBALL_URL" -o "$tarball_path"
fi

# --- Extract ---
mkdir -p "$RUNNER_DIR"
echo "==> Extracting into $RUNNER_DIR..."
tar xzf "$tarball_path" -C "$RUNNER_DIR"
[ -z "${TARBALL_CACHE_DIR:-}" ] && rm -f "$tarball_path"

if [ ! -f "$RUNNER_DIR/run.sh" ]; then
    echo "ERROR: Extraction failed — run.sh not found in $RUNNER_DIR" >&2
    exit 1
fi

# --- Register ---
echo "==> Fetching registration token..."
token=$(gh_registration_token)
if [ -z "$token" ]; then
    echo "ERROR: Failed to get registration token." >&2
    echo "       Run: gh auth refresh -h github.com -s admin:org" >&2
    exit 1
fi

echo "==> Configuring runner..."
"$RUNNER_DIR/config.sh" \
    --url "https://github.com/$ORG" \
    --token "$token" \
    --name "$RUNNER_NAME" \
    --runnergroup "$RUNNER_GROUP" \
    --labels "$RUNNER_LABEL" \
    --work "_work" \
    --unattended \
    --replace
echo "==> Configured."

# --- Start ---
echo "==> Starting on $TARGET_NODE..."
if start_runner "$TARGET_NODE" "$RUNNER_DIR"; then
    echo "$TARGET_NODE" > "$RUNNER_DIR/runner.node"
    pids=$(find_pids "$TARGET_NODE" "$RUNNER_DIR")
    pid=${pids%% *}
    if has_slurm "$TARGET_NODE" "$pid"; then
        echo "==> OK: $RUNNER_NAME running on $TARGET_NODE (PID $pid, slurm in PATH)"
    else
        echo "==> WARNING: $RUNNER_NAME running on $TARGET_NODE (PID $pid) but slurm MISSING from PATH"
    fi
else
    echo "ERROR: $RUNNER_NAME did not start on $TARGET_NODE" >&2
    exit 1
fi

echo ""
echo "==> Log: $RUNNER_DIR/runner.log"
