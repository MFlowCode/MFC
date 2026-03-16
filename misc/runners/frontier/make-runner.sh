#!/usr/bin/env bash
# Frontier convenience wrapper: derives runner name "frontier-<N>" from a number.
# For full name control, use: runner.sh frontier create-runner <name> <node>
#
# Usage: runner.sh frontier make-runner <num> [node]
#   num   Runner number (e.g. 23 creates "frontier-23")
#   node  Login node to start on (default: current host)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

NUM="${1:?Usage: $0 <runner-number> [node]}"
NODE="${2:-$(hostname -s)}"

source "$SCRIPT_DIR/../common/create-runner.sh" "frontier-$NUM" "$NODE"
