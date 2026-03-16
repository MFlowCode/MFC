#!/usr/bin/env bash
# Move a Phoenix runner to a different login node.
# Thin wrapper — see misc/common/move-runner.sh for the implementation.
#
# Usage: bash move-runner.sh <runner-name> <target-node>
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/../common/move-runner.sh"
