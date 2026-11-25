#!/bin/bash
# Append lines to a file of MFlowCode/MFC-ci-data and push.
#
#   ci-data-append.sh <file in the data repository> <file of new lines> <commit message>
#
# DEPLOY_KEY is a deploy key with write access there. Other recorders push to the same
# repository: on a rejected push, start again from a fresh clone.
set -euo pipefail

target=$1; new=$2; message=$3
if [ -z "${DEPLOY_KEY:-}" ]; then echo "::error::The CI_DATA_DEPLOY_KEY secret is not set and $target cannot be recorded."; exit 1; fi

install -m 600 /dev/null ~/deploy_key && printf '%s\n' "$DEPLOY_KEY" > ~/deploy_key
export GIT_SSH_COMMAND="ssh -i ~/deploy_key -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new"
history=$(mktemp -d)
for attempt in 1 2 3 4 5; do
    rm -rf "$history"
    git clone --depth=1 git@github.com:MFlowCode/MFC-ci-data.git "$history"
    # Not the lines already there: a re-run records again what its first attempt recorded.
    touch "$history/$target"
    grep -vxFf "$history/$target" "$new" >> "$history/$target" || true
    git -C "$history" config user.name  'MFC Action'
    git -C "$history" config user.email '<>'
    git -C "$history" add "$target"
    git -C "$history" commit -m "$message" || { echo "Already recorded."; exit 0; }
    if git -C "$history" push; then exit 0; fi
    sleep $((attempt * 7))
done
echo "::error::Could not push $target after 5 attempts."; exit 1
