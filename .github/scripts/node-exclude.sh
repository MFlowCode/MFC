#!/bin/bash
# Bookkeeping for the sbatch --exclude list used when a node fails preflight.
#
# The in-allocation preflight prints "MFC_FAULT_NODE=<name>" into the job's
# output file when it finds the node unusable (dead GPU, SIGILL from a binary
# built for another microarchitecture). The submit wrapper reads that back, adds
# the node to --exclude and resubmits, so the retry lands elsewhere.
#
# Bad nodes are concentrated rather than scattered -- over 2026-08-18..31 one
# Phoenix node accounted for 25 of 29 ECC failures -- which is why excluding the
# offender is worth doing and why it was previously a hand-edited constant.
#
# Usage:
#   node-exclude.sh node-from <output-file>   print the faulted node, if any
#   node-exclude.sh merge <csv> <node>        print <csv> with <node> added once

set -uo pipefail

usage() {
    echo "Usage: $0 {node-from <output-file>|merge <csv> <node>}" >&2
}

case "${1:-}" in
    node-from)
        file="${2:-}"
        [ -f "$file" ] || exit 0
        # Last marker wins: one output path is reused across resubmits, so an
        # earlier attempt's marker can still be sitting above the current one.
        sed -n 's/.*MFC_FAULT_NODE=\([A-Za-z0-9._-]\{1,\}\).*/\1/p' "$file" | tail -n1
        ;;

    merge)
        csv="${2-}"
        node="${3-}"
        if [ -z "$node" ]; then
            printf '%s\n' "$csv"
            exit 0
        fi
        if [ -z "$csv" ]; then
            printf '%s\n' "$node"
            exit 0
        fi
        # Wrapping both sides in commas compares whole fields, so a shorter name
        # that happens to be a prefix of the new one is not mistaken for a match.
        case ",$csv," in
            *",$node,"*) printf '%s\n' "$csv" ;;
            *)           printf '%s,%s\n' "$csv" "$node" ;;
        esac
        ;;

    *)
        usage
        exit 2
        ;;
esac
