#!/bin/bash
# Open an issue when a cluster's runners have stopped picking up jobs, and close it when they
# resume. Read from the jobs themselves, so it needs no token that can list runners: a cluster
# is stuck when a job has waited over two hours for one of its runners while none of its
# runners has a job in progress or started one in that time. Runners that are all busy are
# not stuck, and neither is a cluster nobody is waiting for.
# Env: GH_TOKEN, REPO.
set -euo pipefail

limit=$((2 * 3600)) now=$(date -u +%s)
jobs=$(for status in queued in_progress; do
    gh api "repos/$REPO/actions/runs?status=$status&per_page=100" --jq '.workflow_runs[].id'
done | sort -u | while read -r run; do
    gh api --paginate "repos/$REPO/actions/runs/$run/jobs?per_page=100" \
        --jq '.jobs[] | select(.labels | any(. == "frontier" or . == "gt")) | [(.labels | map(select(. == "frontier" or . == "gt"))[0]), .status, .created_at, (.started_at // "")] | @tsv' || true
done)

for pair in "frontier Frontier" "gt Phoenix"; do
    read -r label cluster <<<"$pair"
    title="No $cluster runner is picking up jobs"
    mine=$(awk -F'\t' -v l="$label" '$1 == l' <<<"$jobs")
    oldest=0 active=0
    while IFS=$'\t' read -r _ status created started; do
        [ -n "$status" ] || continue
        if [ "$status" = "queued" ]; then
            waited=$((now - $(date -u -d "$created" +%s)))
            if [ "$waited" -gt "$oldest" ]; then oldest=$waited; fi
        elif [ "$status" = "in_progress" ] || { [ -n "$started" ] && [ $((now - $(date -u -d "$started" +%s))) -lt "$limit" ]; }; then
            active=$((active + 1))
        fi
    done <<<"$mine"
    echo "$cluster: longest wait $((oldest / 60)) min, $active job(s) running or recently started"
    if [ "$oldest" -gt "$limit" ] && [ "$active" -eq 0 ]; then
        bash "$(dirname "$0")/ci-issue.sh" open "$title" "A job has waited $((oldest / 60)) minutes for a $cluster runner, and no $cluster runner has a job in progress or has started one in the last two hours. The runners are probably down; see \`misc/runners\` for how to restart them."
    elif [ "$active" -gt 0 ]; then
        bash "$(dirname "$0")/ci-issue.sh" close "$title" "A $cluster runner is running jobs again."
    fi
done
