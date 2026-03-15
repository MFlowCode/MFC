# Frontier Runner Management Scripts

Scripts for managing GitHub Actions self-hosted runners on OLCF Frontier login
nodes. The runners submit SLURM jobs to Frontier compute nodes for MFC's CI/CD
pipeline using the `service` partition with the `develop` QOS.

## Background

Frontier has 11 login nodes (`login01`–`login11`). Runner binaries live on
shared Lustre storage (`/lustre/orion/cfd154/proj-shared/runners/`), so moving
a runner between nodes only requires stopping the process on one node and
starting it on another — no binary copying needed.

Each runner directory contains a `runner.node` file recording which login node
it runs on, and a `runner.pid` file with its process ID.

Runners occasionally die due to OLCF's firewall/proxy dropping long-lived TCP
connections to GitHub's broker. The `restart-offline-runners.sh` script handles
recovery. Login nodes vary in stability — if a runner keeps dying on a
particular node, move it to a quieter one (login01 tends to have low load).

## Quick Reference

```bash
# Deploy a new runner on a specific node
bash make-runner.sh 23 login01

# Deploy multiple runners across nodes (e.g. runners 23, 24, 25)
bash deploy-runners.sh 23 login01 login02 login03

# Restart all offline runners in place
bash restart-offline-runners.sh

# Stop and deregister a runner
bash stop-runner.sh frontier-12
```

## Scripts

| Script | Purpose |
|---|---|
| `make-runner.sh` | Download runner binary, register with GitHub via API, start on target node. Usage: `make-runner.sh <num> [login-node]` |
| `deploy-runners.sh` | Deploy multiple runners across login nodes in parallel. Usage: `deploy-runners.sh <start-num> <node1> [node2 ...]` |
| `restart-offline-runners.sh` | Query GitHub for offline frontier runners, SSH to their recorded node, and restart in parallel. Prints final online/offline status. |
| `stop-runner.sh` | Kill the runner process and deregister it from GitHub. Usage: `stop-runner.sh <runner-name>` |

## Troubleshooting

**Runner goes OFFLINE repeatedly on the same node** — That login node may have
process culling or high memory pressure. Move it:
```bash
# Kill on current node, start on a quieter one
ssh login10 "kill $(cat /lustre/orion/cfd154/proj-shared/runners/frontier-1/runner.pid)"
echo "login01" > /lustre/orion/cfd154/proj-shared/runners/frontier-1/runner.node
ssh login01 "nohup /lustre/orion/cfd154/proj-shared/runners/frontier-1/run.sh \
    >> /lustre/orion/cfd154/proj-shared/runners/frontier-1/runner.log 2>&1 < /dev/null &"
```

**Multiple runners OFFLINE at once** — Usually a transient OLCF network blip
to GitHub. Run `restart-offline-runners.sh` to recover all at once.

**runner.pid is wrong after a restart** — This was a known bug (now fixed in
`restart-offline-runners.sh`) where the PID was written before the SSH
command completed. The fix wraps the SSH call in a subshell so the PID is
captured synchronously.
