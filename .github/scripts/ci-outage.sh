#!/bin/bash
# Per-cluster circuit breaker for environment-wide outages.
#
# Some failures are not the node's fault and not the code's fault: pypi.org
# unreachable from a login node, a module tree mid-upgrade, a full project
# filesystem. Requeuing elsewhere cannot help, and every job that starts pays
# the same discovery cost -- on 2026-08-28, 17 Frontier jobs each spent ~33
# minutes learning that PyPI was down.
#
# The first job to notice records a marker on the shared filesystem (every
# self-hosted runner for a cluster shares $HOME); later jobs check it and exit
# immediately instead of submitting a SLURM job that is going to fail.
#
# The breaker is deliberately self-healing. A marker expires after
# MFC_CI_OUTAGE_TTL_SECONDS, and a marker that cannot be parsed is ignored, so
# neither a stale file nor a truncated write can wedge CI. Only jobs that
# actually observe the outage re-mark it, so once the outage clears the breaker
# closes on its own.
#
# Usage:
#   ci-outage.sh mark  <cluster> <reason>   record an outage
#   ci-outage.sh check <cluster>            exit 0 = clear, 1 = outage active
#   ci-outage.sh clear <cluster>            reset the breaker
#
# Env:
#   MFC_CI_STATE_DIR             where markers live (default ~/.mfc-ci-state)
#   MFC_CI_OUTAGE_TTL_SECONDS    marker lifetime in seconds (default 1200)

set -uo pipefail

STATE_DIR="${MFC_CI_STATE_DIR:-$HOME/.mfc-ci-state}"
TTL="${MFC_CI_OUTAGE_TTL_SECONDS:-1200}"

EXIT_CLEAR=0
EXIT_TRIPPED=1
EXIT_USAGE=2

usage() {
    echo "Usage: $0 {mark <cluster> <reason>|check <cluster>|clear <cluster>}" >&2
}

# Keep the marker name filesystem-safe regardless of what the caller passes.
marker_for() {
    local cluster
    cluster=$(printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_')
    printf '%s/outage-%s' "$STATE_DIR" "$cluster"
}

cmd="${1:-}"
cluster="${2:-}"

if [ -z "$cmd" ] || [ -z "$cluster" ]; then
    usage
    exit $EXIT_USAGE
fi

marker=$(marker_for "$cluster")

case "$cmd" in
    mark)
        reason="${3:-unspecified}"
        mkdir -p "$STATE_DIR" || exit $EXIT_USAGE
        # Write to a temporary file and rename so a concurrent `check` never
        # observes a half-written marker.
        tmp="${marker}.$$.tmp"
        {
            date +%s
            printf '%s\n' "$reason"
        } > "$tmp" && mv -f "$tmp" "$marker"
        echo "Recorded $cluster outage: $reason"
        echo "  marker: $marker (expires after ${TTL}s)"
        ;;

    check)
        [ -f "$marker" ] || exit $EXIT_CLEAR

        stamp=$(head -n1 "$marker" 2>/dev/null)
        reason=$(tail -n +2 "$marker" 2>/dev/null)

        # A marker we cannot parse is treated as absent: an unreadable breaker
        # must never be an un-clearable one.
        case "$stamp" in
            ''|*[!0-9]*)
                echo "Ignoring unparseable outage marker $marker"
                exit $EXIT_CLEAR
                ;;
        esac

        age=$(( $(date +%s) - stamp ))
        if [ "$age" -ge "$TTL" ] || [ "$age" -lt 0 ]; then
            exit $EXIT_CLEAR
        fi

        echo "::warning::Skipping: known $cluster outage recorded ${age}s ago: ${reason:-unspecified}"
        echo "Clear it early by deleting $marker"
        exit $EXIT_TRIPPED
        ;;

    clear)
        rm -f "$marker"
        echo "Cleared any $cluster outage marker ($marker)"
        ;;

    *)
        usage
        exit $EXIT_USAGE
        ;;
esac
