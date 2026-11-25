#!/bin/bash
# Keep one issue open while a named condition holds, and close it when it clears.
#   ci-issue.sh open   "<title>" "<what happened>"     once: an open issue is left as it is
#   ci-issue.sh report "<title>" "<what happened>"     every time: an open issue gets a comment
#   ci-issue.sh close  "<title>" "<what cleared it>"
# Env: GH_TOKEN, REPO.
set -euo pipefail

action=$1 title=$2 text=${3:-}
issue=$(gh issue list -R "$REPO" --state open --search "in:title \"$title\"" --json number,title \
    --jq "[.[] | select(.title == \"$title\")][0].number // \"\"")

case "$action" in
    open|report)
        if [ -n "$issue" ] && [ "$action" = report ]; then gh issue comment "$issue" -R "$REPO" --body "$text"; exit 0; fi
        if [ -n "$issue" ]; then echo "Issue #$issue is already open."; exit 0; fi
        gh issue create -R "$REPO" --title "$title" --label continuous-integration --body "$text"
        ;;
    close)
        if [ -z "$issue" ]; then exit 0; fi
        gh issue close "$issue" -R "$REPO" --comment "$text"
        ;;
esac
