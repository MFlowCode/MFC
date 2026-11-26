#!/bin/bash
# Keep one issue open per cluster while it fails in the scheduled master run, and close it on
# the next green one. A scheduled run that found master unchanged skips the cluster jobs and
# changes nothing here.
# A job that ends cancelled or timed out (preempted until the job timeout, or never
# scheduled) gave no answer, so it counts too, labelled "no result".
#
# Env: GH_TOKEN, REPO, RUN_ID, RUN_ATTEMPT, HEAD_SHA, CLUSTER (Phoenix or Frontier; matched
# against job names). DRY_RUN=1 prints instead of writing; RUN_EVENT (default schedule) is the
# kind of master run compared against.
set -euo pipefail

CLUSTER="${CLUSTER:-Phoenix}"
RUN_EVENT="${RUN_EVENT:-schedule}"
TITLE="$CLUSTER CI failing on master"
api() { gh api "$@"; }
body=$(mktemp)
trap 'rm -f "$body"' EXIT

# The cluster's job conclusions of one run attempt, one per line: "<conclusion>\t<name>\t<url>".
phoenix_jobs() {
    api --paginate "repos/$REPO/actions/runs/$1/attempts/$2/jobs?per_page=100" \
        --jq ".jobs[] | select(.name | test(\"$CLUSTER\")) | [(.conclusion // \"none\"), .name, .html_url] | @tsv"
}

jobs=$(phoenix_jobs "$RUN_ID" "$RUN_ATTEMPT")
failed=$(awk -F'\t' '$1 == "failure" || $1 == "cancelled" || $1 == "timed_out"' <<<"$jobs")
passed=$(awk -F'\t' '$1 == "success"' <<<"$jobs")

issue=$(gh issue list -R "$REPO" --state open --search "in:title \"$TITLE\"" --json number,title \
    --jq ".[] | select(.title == \"$TITLE\") | .number" | head -1)
run_url="https://github.com/$REPO/actions/runs/$RUN_ID"

write() {  # write <description> <gh args...>
    if [ "${DRY_RUN:-0}" = "1" ]; then echo "DRY_RUN: $1"; case "$1" in bisect*) ;; *) cat "$body" ;; esac; return; fi
    shift; "$@"
}

if [ -z "$failed" ]; then
    # Close only on a run where every job of this cluster passed.
    if [ -n "$issue" ] && [ -n "$passed" ] && ! awk -F'\t' '$1 != "success"' <<<"$jobs" | grep -q .; then
        printf '%s passed on %s (%s).\n' "$CLUSTER" "$HEAD_SHA" "$run_url" > "$body"
        write "close #$issue" gh issue close "$issue" -R "$REPO" --comment "$(cat "$body")"
    fi
    exit 0
fi

# Last scheduled master run whose jobs on this cluster all passed, to list the merges in between.
last_green=""
created=$(api "repos/$REPO/actions/runs/$RUN_ID" --jq .created_at)
while IFS=$'\t' read -r rid attempt sha; do
    j=$(phoenix_jobs "$rid" "$attempt")
    if [ -n "$j" ] && ! awk -F'\t' '$1 != "success"' <<<"$j" | grep -q .; then last_green=$sha; break; fi
done < <(api "repos/$REPO/actions/workflows/test.yml/runs?branch=master&event=$RUN_EVENT&status=completed&created=<$created&per_page=100" \
    --jq '.workflow_runs[] | [.id, .run_attempt, .head_sha] | @tsv')

{
    printf '%s did not pass on master at %s: %s\n\n' "$CLUSTER" "$HEAD_SHA" "$run_url"
    awk -F'\t' '{printf "- [%s](%s): %s\n", $2, $3, ($1 == "failure" ? "failed" : "no result (" $1 ")")}' <<<"$failed"
    if [ -n "$last_green" ]; then
        printf '\nMerges since the last green %s run (%s):\n\n' "$CLUSTER" "${last_green:0:9}"
        api "repos/$REPO/compare/$last_green...$HEAD_SHA" \
            --jq '.commits[] | "- " + .sha[0:9] + " " + (.commit.message | split("\n")[0])' > "$body.merges"
        tail -30 "$body.merges"
        # One scheduled run covers every merge since the previous one.
        if [ "$(wc -l < "$body.merges")" -gt 1 ]; then
            # shellcheck disable=SC2016 # backticks are Markdown
            printf '\nTo test one commit: `git push origin <sha>:refs/heads/bisect-<sha>` and `gh workflow run test.yml --ref bisect-<sha>`, then delete the branch.\n'
        fi
        rm -f "$body.merges"
    else
        printf '\nNo green %s run in the previous 100 scheduled master runs.\n' "$CLUSTER"
    fi
} > "$body"

if [ -n "$issue" ]; then
    write "comment on #$issue" gh issue comment "$issue" -R "$REPO" --body-file "$body"
else
    url=$(write "open issue" gh issue create -R "$REPO" --title "$TITLE" --label continuous-integration --body-file "$body")
    echo "$url"
    issue=${url##*/}
fi

# Bisect each test lane that has just gone red: it failed here and passed the last time a
# scheduled run ran it, so the merges between the two are the suspects, whatever else is already
# red. A lane that was already red is not bisected again, and a lane with no result (cancelled,
# timed out) has nothing to bisect.
prev_runs=$(api "repos/$REPO/actions/workflows/test.yml/runs?branch=master&event=$RUN_EVENT&status=completed&created=<$created&per_page=30" \
    --jq '.workflow_runs[] | [.id, .run_attempt, .head_sha] | @tsv')
declare -A prev_jobs
# The lane's test jobs with a verdict, in one run's jobs: a skipped job keeps its unexpanded
# matrix name and does not match, and a cancelled one says nothing.
lane_jobs() { awk -F'\t' -v c="$compiler" -v i="gpu-$interface" 'index($2, c) && index($2, i) && $2 !~ /^Case Opt/ && ($1 == "success" || $1 == "failure" || $1 == "timed_out")'; }
# curl, not `gh api`: gh 2.97 and later refuses to print a log, which always holds colour codes.
job_log() { curl -fsSL -H "Authorization: Bearer $GH_TOKEN" "https://api.github.com/repos/$REPO/actions/jobs/$1/logs"; }

# What the lane's failed jobs among <job rows> give a bisect to reproduce: sets tests (the ones
# cluster-lane.sh printed in the job log) and build_failed (the Build step failed).
reproducible() {
    local u
    tests="" build_failed=""
    while read -r u; do
        [ -n "$u" ] || continue
        tests="$tests $(job_log "${u##*/}" | sed -n 's/.*FAILED_TESTS: //p' | tr '\n' ' ' || true)"
        if [ -n "$(api "repos/$REPO/actions/jobs/${u##*/}" --jq '.steps[] | select(.conclusion == "failure" and .name == "Build") | .name')" ]; then build_failed=1; fi
    done < <(awk -F'\t' -v c="$compiler" -v i="gpu-$interface" '$1 == "failure" && index($2, c) && index($2, i) && $2 !~ /^Case Opt/ {print $3}' <<<"$1")
    # At most 15: from 20 tests on, the runner aborts a run in which 30% fail and names none of
    # them, and a bisect step that fails without a named test cannot judge its commit.
    tests=$(tr ' ' '\n' <<<"$tests" | { grep -E '^[0-9A-Fa-f]{8}$' || true; } | sort -u | head -15 | xargs)
}

earlier_attempt=""
if [ "$RUN_ATTEMPT" -gt 1 ]; then earlier_attempt=$(phoenix_jobs "$RUN_ID" "$((RUN_ATTEMPT - 1))"); fi
done_lanes=" "
while IFS=$'\t' read -r conclusion name job_url; do
    [ "$conclusion" = "failure" ] || continue
    case "$name" in
        "Case Opt"*) continue ;;
        *"Frontier (CCE)"*) lane=frontier     compiler="Frontier (CCE)" ;;
        *"Frontier (AMD)"*) lane=frontier_amd compiler="Frontier (AMD)" ;;
        *Phoenix*)          lane=phoenix      compiler="Phoenix" ;;
        *) continue ;;
    esac
    interface=$(sed -n 's/.*(gpu-\(acc\|omp\)[ )].*/\1/p' <<<"$name")
    [ -n "$interface" ] || continue
    case "$done_lanes" in *" $lane-$interface "*) continue ;; esac     # one bisect per lane, not per shard
    done_lanes="$done_lanes$lane-$interface "

    # A re-run in which the lane fails again is bisected once: not if the attempt before had a
    # failure a bisect could reproduce, which that attempt's alert then dispatched.
    reproducible "$earlier_attempt"
    if [ -n "$tests" ] || [ -n "$build_failed" ]; then continue; fi

    # The last earlier run that ran this lane: a run that found master unchanged skipped it.
    good=""
    while IFS=$'\t' read -r rid attempt sha; do
        [ -n "$rid" ] || continue
        [ -n "${prev_jobs[$rid]+x}" ] || prev_jobs[$rid]=$(phoenix_jobs "$rid" "$attempt")
        before=$(lane_jobs <<<"${prev_jobs[$rid]}")
        [ -n "$before" ] || continue
        if [ -z "$(awk -F'\t' '$1 != "success"' <<<"$before")" ]; then good=$sha; fi
        break
    done <<<"$prev_runs"
    if [ -z "$good" ] || [ "$good" = "$HEAD_SHA" ]; then continue; fi     # already red, or never seen green
    if [ "$(api "repos/$REPO/compare/$good...$HEAD_SHA" --jq '.commits | length')" -le 1 ]; then continue; fi   # one merge: named above

    reproducible "$failed"
    # Neither a failed test nor a failed build: nothing a bisect could reproduce.
    if [ -z "$tests" ] && [ -z "$build_failed" ]; then
        echo "$lane gpu-$interface failed without a failed test or a failed build; not bisecting."
        continue
    fi
    write "bisect $lane gpu-$interface on: ${tests:-the build}" gh workflow run bisect.yml -R "$REPO" --ref master \
        -f cluster="$lane" -f interface="$interface" -f good="$good" -f bad="$HEAD_SHA" -f tests="$tests" -f issue="$issue"
done <<<"$failed"
