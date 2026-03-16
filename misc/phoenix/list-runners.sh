#!/usr/bin/env bash
# List all Phoenix runners combining GitHub API status with live node process info.
# Thin wrapper — see misc/common/list-runners.sh for the implementation.
#
# Usage: bash list-runners.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/list-runners.sh"
