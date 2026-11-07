#!/bin/bash
# Keep one issue open while Phoenix fails on master, and close it on the next green run.
# Phoenix runs only after merge (see the phoenix job in test.yml), so this is how a
# failure gets noticed. Preemption is resubmitted inside the job, so a job that ends in
# "failure" is a build or test failure.
#
# Env: GH_TOKEN, REPO, RUN_ID, RUN_ATTEMPT, HEAD_SHA. DRY_RUN=1 prints instead of writing.
set -euo pipefail

TITLE="Phoenix CI failing on master"
api() { gh api "$@"; }
body=$(mktemp)
trap 'rm -f "$body"' EXIT

# Phoenix job conclusions of one run attempt, one per line: "<conclusion>\t<name>\t<url>".
phoenix_jobs() {
    api --paginate "repos/$REPO/actions/runs/$1/attempts/$2/jobs?per_page=100" \
        --jq '.jobs[] | select(.name | test("Phoenix")) | [(.conclusion // "none"), .name, .html_url] | @tsv'
}

jobs=$(phoenix_jobs "$RUN_ID" "$RUN_ATTEMPT")
failed=$(awk -F'\t' '$1 == "failure"' <<<"$jobs")
passed=$(awk -F'\t' '$1 == "success"' <<<"$jobs")

issue=$(gh issue list -R "$REPO" --state open --search "in:title \"$TITLE\"" --json number,title \
    --jq ".[] | select(.title == \"$TITLE\") | .number" | head -1)
run_url="https://github.com/$REPO/actions/runs/$RUN_ID"

write() {  # write <description> <gh args...>
    if [ "${DRY_RUN:-0}" = "1" ]; then echo "DRY_RUN: $1"; cat "$body"; return; fi
    shift; "$@"
}

if [ -z "$failed" ]; then
    # Close only on a run where every Phoenix job passed; skipped or cancelled says nothing.
    if [ -n "$issue" ] && [ -n "$passed" ] && ! awk -F'\t' '$1 != "success"' <<<"$jobs" | grep -q .; then
        printf 'Phoenix passed on %s (%s).\n' "$HEAD_SHA" "$run_url" > "$body"
        write "close #$issue" gh issue close "$issue" -R "$REPO" --comment "$(cat "$body")"
    fi
    exit 0
fi

# Last master push whose Phoenix jobs all passed, to list the merges in between.
last_green=""
created=$(api "repos/$REPO/actions/runs/$RUN_ID" --jq .created_at)
while IFS=$'\t' read -r rid attempt sha; do
    j=$(phoenix_jobs "$rid" "$attempt")
    if [ -n "$j" ] && ! awk -F'\t' '$1 != "success"' <<<"$j" | grep -q .; then last_green=$sha; break; fi
done < <(api "repos/$REPO/actions/workflows/test.yml/runs?branch=master&event=push&status=completed&created=<$created&per_page=100" \
    --jq '.workflow_runs[] | [.id, .run_attempt, .head_sha] | @tsv')

{
    printf 'Phoenix failed on master at %s: %s\n\n' "$HEAD_SHA" "$run_url"
    awk -F'\t' '{printf "- [%s](%s)\n", $2, $3}' <<<"$failed"
    if [ -n "$last_green" ]; then
        printf '\nMerges since the last green Phoenix run (%s):\n\n' "${last_green:0:9}"
        api "repos/$REPO/compare/$last_green...$HEAD_SHA" \
            --jq '.commits[] | "- " + .sha[0:9] + " " + (.commit.message | split("\n")[0])' > "$body.merges"
        tail -30 "$body.merges"
        # Master runs coalesce, so one red run can cover several merges.
        if [ "$(wc -l < "$body.merges")" -gt 1 ]; then
            # shellcheck disable=SC2016 # backticks are Markdown
            printf '\nTo test one commit: `git push origin <sha>:refs/heads/bisect-<sha>` and `gh workflow run test.yml --ref bisect-<sha>`, then delete the branch.\n'
        fi
        rm -f "$body.merges"
    else
        printf '\nNo green Phoenix run in the previous 100 master pushes.\n'
    fi
} > "$body"

if [ -n "$issue" ]; then
    write "comment on #$issue" gh issue comment "$issue" -R "$REPO" --body-file "$body"
else
    write "open issue" gh issue create -R "$REPO" --title "$TITLE" --label continuous-integration --body-file "$body"
fi
