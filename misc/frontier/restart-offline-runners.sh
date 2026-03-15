#!/usr/bin/env bash
# Restart all offline frontier runners.
# Queries GitHub for offline runners, then SSHes to each runner's node
# and restarts run.sh in the background. All restarts happen in parallel.
set -euo pipefail

SSH_OPTS="-o StrictHostKeyChecking=no -o ConnectTimeout=30 -o ServerAliveInterval=10 -o ServerAliveCountMax=3"
SHARED_DIR="/lustre/orion/cfd154/proj-shared/runners"
ORG="MFlowCode"

echo "==> Checking for offline frontier runners..."
OFFLINE=$(gh api orgs/${ORG}/actions/runners \
    --jq '.runners[] | select(.name | startswith("frontier")) | select(.status == "offline") | .name')

if [ -z "${OFFLINE}" ]; then
    echo "==> All frontier runners are online. Nothing to do."
    exit 0
fi

echo "==> Offline runners: $(echo ${OFFLINE} | tr '\n' ' ')"

for RUNNER_NAME in ${OFFLINE}; do
    RUNNER_DIR="${SHARED_DIR}/${RUNNER_NAME}"
    NODE_FILE="${RUNNER_DIR}/runner.node"

    if [ ! -d "${RUNNER_DIR}" ]; then
        echo "WARN: No directory for ${RUNNER_NAME}, skipping."
        continue
    fi

    if [ ! -f "${NODE_FILE}" ]; then
        echo "WARN: No runner.node file for ${RUNNER_NAME}, skipping."
        continue
    fi

    TARGET_NODE=$(cat "${NODE_FILE}")
    CURRENT_NODE=$(hostname -s)

    echo "==> Restarting ${RUNNER_NAME} on ${TARGET_NODE}..."
    if [ "${TARGET_NODE}" = "${CURRENT_NODE}" ]; then
        nohup "${RUNNER_DIR}/run.sh" >> "${RUNNER_DIR}/runner.log" 2>&1 < /dev/null &
        echo $! > "${RUNNER_DIR}/runner.pid"
        echo "    PID: $(cat ${RUNNER_DIR}/runner.pid)"
    else
        (
            PID=$(ssh ${SSH_OPTS} "${TARGET_NODE}" \
                "nohup ${RUNNER_DIR}/run.sh >> ${RUNNER_DIR}/runner.log 2>&1 < /dev/null & echo \$!")
            echo "${PID}" > "${RUNNER_DIR}/runner.pid"
            echo "    ${RUNNER_NAME} PID: ${PID}"
        ) &
    fi
done

wait
echo "==> Done. Waiting 5s for runners to register..."
sleep 5
gh api orgs/${ORG}/actions/runners \
    --jq '.runners[] | select(.name | startswith("frontier")) | {name, status}'
