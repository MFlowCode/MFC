#!/usr/bin/env bash
# Restart a single Frontier runner on a given node.
# Thin wrapper — see misc/common/restart-runner.sh for the implementation.
#
# Usage: bash restart-runner.sh <node> <runner-dir>
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/restart-runner.sh"
