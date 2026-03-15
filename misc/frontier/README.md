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
it was last started on. This is used as a fallback hint when restarting offline
runners. The authoritative source of truth for whether a runner is running (and
on which node) is CWD-based process discovery — not any PID file.

`runner.node` is self-healing: `restart-offline-runners.sh` detects when a
runner is actually running on a different node than `runner.node` records (e.g.
after a manual restart) and corrects the file automatically.

Runners occasionally die due to OLCF's firewall/proxy dropping long-lived TCP
connections to GitHub's broker. The `restart-offline-runners.sh` script handles
recovery. Login nodes vary in stability — if a runner keeps dying on a
particular node, move it to a quieter one (login01 tends to have low load).

## Quick Reference

```bash
# List all runners with GitHub status, node, and memory usage
bash list-runners.sh

# Check runner health across all login nodes
bash check-runners.sh

# Deploy a new runner on a specific node
bash make-runner.sh 23 login01

# Deploy multiple runners across nodes (e.g. runners 23, 24, 25)
bash deploy-runners.sh 23 login01 login02 login03

# Restart all offline runners in place
bash restart-offline-runners.sh

# Move a runner to a different login node
bash move-runner.sh frontier-1 login01

# Stop and deregister a runner
bash stop-runner.sh frontier-12

# Rerun failed CI workflows (site-agnostic, also available at misc/common/)
bash ../common/rerun-failed.sh
APPLY=1 bash ../common/rerun-failed.sh
```

## Scripts

| Script | Purpose |
|---|---|
| `config.sh` | Shared configuration, constants, GitHub API helpers, and CWD-based process management functions. Sourced by all other scripts. |
| `check-runners.sh` | SSH to each login node, show Runner.Listener processes with name, status (idle/BUSY), and RSS memory. |
| `list-runners.sh` | List all runners with GitHub API status, actual node (from parallel SSH sweep), and RSS memory. Flags stale `runner.node` entries. |
| `make-runner.sh` | Download runner binary, register with GitHub via API, start on target node. Usage: `make-runner.sh <num> [login-node]` |
| `move-runner.sh` | Move a runner to a different login node: stop on current node, update `runner.node`, start on target. Usage: `move-runner.sh <runner-name> <target-node>` |
| `deploy-runners.sh` | Deploy multiple runners across login nodes in parallel. Usage: `deploy-runners.sh <start-num> <node1> [node2 ...]` |
| `restart-offline-runners.sh` | Query GitHub for offline frontier runners, locate via CWD-based discovery, stop stale processes, then restart in parallel. Self-heals stale `runner.node` files. Prints final status. |
| `stop-runner.sh` | Locate runner via CWD-based discovery, stop the process, and deregister from GitHub. Usage: `stop-runner.sh <runner-name>` |
| `../common/rerun-failed.sh` | Rerun failed GitHub Actions workflows on open PRs and master. No site-specific code. |

## Troubleshooting

**Runner goes OFFLINE repeatedly on the same node** — That login node may have
process culling or high memory pressure. Move it to a different node:
```bash
bash move-runner.sh frontier-1 login01
```

**Multiple runners OFFLINE at once** — Usually a transient OLCF network blip
to GitHub. Run `restart-offline-runners.sh` to recover all at once.

**Runner appears offline on GitHub but process is running** — GitHub status can
lag. `restart-offline-runners.sh` uses CWD-based process discovery first: if a
process is found running, it will stop it before restarting, preventing
duplicate runner processes.
