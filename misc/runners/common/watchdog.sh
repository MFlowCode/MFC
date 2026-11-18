#!/usr/bin/env bash
# Revive runners that GitHub reports offline, and do nothing else.
#
# Sourced by misc/runners/runner.sh after the site config is loaded.
#
# Unlike rebalance-runners, this never moves a healthy runner, so it is safe to
# run unattended on a timer. A login-node reboot or OLCF's firewall dropping the
# GitHub connection otherwise leaves runners dead until someone notices.
#
# For each runner GitHub reports offline (still offline after a grace period):
#   - listener found on a reachable node, idle, and up for a while
#       -> it lost its GitHub session; restart it in place
#   - listener found but busy, just started, or on more than one node
#       -> leave it (dedupe-runners handles duplicates)
#   - no listener on any reachable node
#       -> start it on the least-loaded reachable node, unless its runner.node
#          is a node that cannot be reached right now: it may still be running
#          there, and starting it elsewhere would create a duplicate listener
#
# Several watchdogs may run at once (one per login node, so a reboot cannot
# take them all out); a lock directory on shared storage serializes their
# passes. mkdir is used rather than flock because it is atomic on NFS, GPFS
# and Lustre alike, whereas flock across nodes depends on mount options.
#
# Usage: runner.sh <site> watchdog [interval-seconds]
#   no interval   one pass, then exit
#   interval      loop forever, one pass every interval seconds
#
# Environment:
#   APPLY=1                  act; otherwise print what would be done (dry run)
#   WATCHDOG_GRACE=60        seconds a runner must stay offline before acting
#   WATCHDOG_MIN_UPTIME=300  seconds a live-but-offline listener must have been
#                            up before it is restarted in place
#   WATCHDOG_LOCK=<path>     lock directory (default: next to the runner dirs)
#   WATCHDOG_LOCK_STALE=1800 seconds after which a held lock counts as stale
set -euo pipefail

WATCHDOG_INTERVAL="${1:-}"
WATCHDOG_GRACE="${WATCHDOG_GRACE:-60}"
WATCHDOG_MIN_UPTIME="${WATCHDOG_MIN_UPTIME:-300}"
WATCHDOG_LOCK_STALE="${WATCHDOG_LOCK_STALE:-1800}"

if [ -n "$WATCHDOG_INTERVAL" ] && ! [[ "$WATCHDOG_INTERVAL" =~ ^[0-9]+$ ]]; then
    echo "Usage: runner.sh <site> watchdog [interval-seconds]" >&2
    exit 1
fi

wd_log() { echo "[$(date '+%Y-%m-%d %H:%M:%S') $(hostname -s)] $*"; }

# Print "name status" for every runner of this site that GitHub knows about.
wd_github_status() {
    gh_list_runners | awk '{print $2, $3}'
}

# Print the runner dirs GitHub reports offline, as "dir name" lines.
wd_offline_runners() {
    local status_file="$1" dir name
    while IFS= read -r dir; do
        name=$(get_runner_name "$dir")
        [ -z "$name" ] && continue
        if grep -qxF "$name offline" "$status_file"; then
            echo "$dir $name"
        fi
    done < <(find_runner_dirs)
}

wd_lock_acquire() {
    local lock="$1" age
    if mkdir "$lock" 2>/dev/null; then
        echo "$(hostname -s) $$ $(date +%s)" > "$lock/owner"
        return 0
    fi
    age=$(( $(date +%s) - $(stat -c %Y "$lock" 2>/dev/null || date +%s) ))
    if [ "$age" -gt "$WATCHDOG_LOCK_STALE" ]; then
        wd_log "Removing stale lock ($age s old, owner: $(cat "$lock/owner" 2>/dev/null || echo unknown))"
        rm -rf "$lock"
        mkdir "$lock" 2>/dev/null || return 1
        echo "$(hostname -s) $$ $(date +%s)" > "$lock/owner"
        return 0
    fi
    return 1
}

wd_pass() {
    local tmpdir lock first_dir
    tmpdir=$(mktemp -d)

    if ! wd_github_status > "$tmpdir/status" 2>/dev/null || [ ! -s "$tmpdir/status" ]; then
        wd_log "Could not query GitHub runner status; skipping this pass"
        rm -rf "$tmpdir"
        return 0
    fi
    wd_offline_runners "$tmpdir/status" > "$tmpdir/offline"
    if [ ! -s "$tmpdir/offline" ]; then
        rm -rf "$tmpdir"
        return 0
    fi

    wd_log "Offline on GitHub: $(awk '{print $2}' "$tmpdir/offline" | tr '\n' ' ')- rechecking in ${WATCHDOG_GRACE}s"
    sleep "$WATCHDOG_GRACE"

    # `|| true`: head exits after one line, so under pipefail the SIGPIPE
    # find_runner_dirs takes would otherwise abort the pass.
    first_dir=$(find_runner_dirs | head -1 || true)
    lock="${WATCHDOG_LOCK:-$(dirname "$first_dir")/.runner-watchdog.lock}"
    if ! wd_lock_acquire "$lock"; then
        wd_log "Another watchdog holds $lock ($(cat "$lock/owner" 2>/dev/null || echo unknown)); skipping this pass"
        rm -rf "$tmpdir"
        return 0
    fi

    # Re-query under the lock: another watchdog, or a runner's own reconnect or
    # self-update, may have brought some back during the grace period.
    if ! wd_github_status > "$tmpdir/status" 2>/dev/null || [ ! -s "$tmpdir/status" ]; then
        wd_log "Could not re-query GitHub runner status; skipping this pass"
        rm -rf "$lock" "$tmpdir"
        return 0
    fi
    wd_offline_runners "$tmpdir/status" > "$tmpdir/offline"
    if [ ! -s "$tmpdir/offline" ]; then
        wd_log "All runners back online"
        rm -rf "$lock" "$tmpdir"
        return 0
    fi

    # Which nodes answer SSH right now?
    local node
    declare -A reachable=() load=()
    for node in "${NODES[@]}"; do
        ( ssh $SSH_OPTS "$node" true </dev/null >/dev/null 2>&1 && touch "$tmpdir/up.$node" ) &
    done
    wait
    for node in "${NODES[@]}"; do
        if [ -f "$tmpdir/up.$node" ]; then
            reachable[$node]=1
            load[$node]=0
        fi
    done
    if [ ${#reachable[@]} -eq 0 ]; then
        wd_log "No login node reachable; skipping this pass"
        rm -rf "$lock" "$tmpdir"
        return 0
    fi

    sweep_all_listeners "$tmpdir" > "$tmpdir/listeners"
    while read -r _l node _pid _dir _busy; do
        [ -n "${load[$node]+x}" ] && load[$node]=$(( ${load[$node]} + 1 ))
    done < "$tmpdir/listeners"

    local dir name hits count lnode lpid lbusy uptime recorded best
    while read -r dir name; do
        hits=$(awk -v d="$dir" '$4 == d' "$tmpdir/listeners")
        count=$(printf '%s' "$hits" | grep -c . || true)

        if [ "$count" -gt 1 ]; then
            wd_log "$name: listening on more than one node; leaving it for dedupe-runners"
            continue
        fi

        if [ "$count" -eq 1 ]; then
            read -r _l lnode lpid _d lbusy <<< "$hits"
            if [ "$lbusy" = "busy" ]; then
                wd_log "$name: offline on GitHub but running a job on $lnode; leaving it"
                continue
            fi
            uptime=$(ssh $SSH_OPTS "$lnode" "ps -o etimes= -p $lpid" 2>/dev/null | tr -d ' ' || true)
            if ! [[ "$uptime" =~ ^[0-9]+$ ]] || [ "$uptime" -lt "$WATCHDOG_MIN_UPTIME" ]; then
                wd_log "$name: listener on $lnode started ${uptime:-?}s ago; giving it time to connect"
                continue
            fi
            if [ "${APPLY:-0}" != "1" ]; then
                wd_log "$name: idle listener on $lnode lost its GitHub session; would restart it in place"
                continue
            fi
            wd_log "$name: idle listener on $lnode lost its GitHub session; restarting in place"
            if ! stop_runner "$lnode" "$dir"; then
                wd_log "  ERROR: could not stop $name on $lnode; leaving it"
                continue
            fi
            best="$lnode"
        else
            recorded=$(cat "$dir/runner.node" 2>/dev/null || true)
            if [ -n "$recorded" ] && [ -z "${reachable[$recorded]+x}" ]; then
                wd_log "$name: last ran on $recorded, which is unreachable; skipping (it may still be running there)"
                continue
            fi
            best=""
            for node in "${NODES[@]}"; do
                [ -n "${reachable[$node]+x}" ] || continue
                if [ -z "$best" ] || [ "${load[$node]}" -lt "${load[$best]}" ]; then
                    best="$node"
                fi
            done
            if [ "${APPLY:-0}" != "1" ]; then
                wd_log "$name: not running on any reachable node; would start it on $best"
                load[$best]=$(( ${load[$best]} + 1 ))
                continue
            fi
            wd_log "$name: not running on any reachable node; starting it on $best"
            load[$best]=$(( ${load[$best]} + 1 ))
        fi

        if start_runner "$best" "$dir"; then
            echo "$best" > "$dir/runner.node"
            lpid=$(find_pids "$best" "$dir" || true)
            lpid=${lpid%% *}
            if [ -n "$lpid" ] && has_slurm "$best" "$lpid"; then
                wd_log "  OK: $name started on $best (PID $lpid, slurm in PATH)"
            else
                wd_log "  WARNING: $name started on $best but slurm not confirmed in its PATH"
            fi
        else
            wd_log "  ERROR: $name did not start on $best"
        fi
    done < "$tmpdir/offline"

    rm -rf "$lock" "$tmpdir"
}

if [ -z "$WATCHDOG_INTERVAL" ]; then
    wd_pass
else
    wd_log "Watchdog started (every ${WATCHDOG_INTERVAL}s, APPLY=${APPLY:-0})"
    while :; do
        wd_pass || wd_log "Pass failed; will retry"
        sleep "$WATCHDOG_INTERVAL"
    done
fi
