#!/bin/bash
# Find the commit that broke one cluster lane on master: git bisect between the last green
# scheduled run and the red one, building each candidate and running only the failed tests.
# Dispatched by phoenix-master-alert.sh through bisect.yml.
#
# Env: CLUSTER (frontier|frontier_amd|phoenix), INTERFACE (acc|omp), GOOD, BAD, TESTS (UUIDs;
#      empty bisects the build), ISSUE (optional: issue to report on), REPO, GH_TOKEN.
#   cluster-bisect.sh        run the bisection
#   cluster-bisect.sh step   test the checked-out commit (called by git bisect run)
set -uo pipefail

if [ "${1:-}" = "step" ]; then
    lane=(gpu "$INTERFACE" "$CLUSTER")
    # 125 tells git bisect this commit cannot be judged.
    [ -f .github/scripts/cluster-lane.sh ] || exit 125
    rm -f ./*.out
    if [ "$CLUSTER" != "phoenix" ]; then
        bash ".github/workflows/$CLUSTER/build.sh" gpu "$INTERFACE" || exit 125
    fi
    rc=0
    bash .github/scripts/cluster-lane.sh build "${lane[@]}" || rc=$?
    # 78: the cluster is down, which says nothing about this commit.
    if [ "$rc" -eq 78 ]; then exit 125; fi
    if [ "$rc" -ne 0 ]; then
        if [ -n "$TESTS" ]; then exit 125; else exit 1; fi
    fi
    [ -n "$TESTS" ] || exit 0
    rc=0
    MFC_TEST_EXTRA_OPTS="--only $TESTS" bash .github/scripts/cluster-lane.sh test "${lane[@]}" || rc=$?
    if [ "$rc" -eq 78 ]; then exit 125; fi
    [ "$rc" -eq 0 ] || exit 1
    exit 0
fi

# The step script is run from a copy: bisect checks out other commits underneath it.
self=$(mktemp)
log=$(mktemp)
trap 'rm -f "$self" "$log"; git bisect reset >/dev/null 2>&1 || true' EXIT
cp "$0" "$self"

git bisect start --first-parent "$BAD" "$GOOD" || exit 1
git bisect run bash "$self" step | tee "$log"
culprit=$(sed -n 's/^\([0-9a-f]\{40\}\) is the first bad commit$/\1/p' "$log" | head -1)

{
    printf 'Bisected %s gpu-%s between %s (green) and %s (red)' "$CLUSTER" "$INTERFACE" "${GOOD:0:9}" "${BAD:0:9}"
    if [ -n "$TESTS" ]; then printf ', on tests %s' "$TESTS"; else printf ', on the build'; fi
    if [ -n "$culprit" ]; then
        printf '.\n\nFirst bad commit: %s %s\n' "$culprit" "$(git log -1 --format=%s "$culprit")"
    else
        printf '.\n\nNo single commit found: the failure did not reproduce, or a candidate could not be built.\n'
    fi
} > "$log.msg"
cat "$log.msg"
if [ -n "${ISSUE:-}" ]; then gh issue comment "$ISSUE" -R "$REPO" --body-file "$log.msg"; fi
rm -f "$log.msg"
[ -n "$culprit" ]
