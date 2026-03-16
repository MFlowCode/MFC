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
on which node) is EXE-based process discovery via `/proc/$pid/exe` — not any
PID file.

`runner.node` is self-healing: `rebalance-runners.sh` calls `sync_runner_nodes`
at startup, which sweeps all nodes in parallel and corrects any stale
`runner.node` files automatically.

Runners occasionally die due to OLCF's firewall/proxy dropping long-lived TCP
connections to GitHub's broker. Run `rebalance-runners.sh` to restart and
redistribute them. Login nodes vary in stability — if a runner keeps dying on a
particular node, move it to a quieter one (login01 tends to have low load).

## Quick Reference

```bash
# List all runners with GitHub status, node, slurm, and memory usage
bash list-runners.sh

# Check runner health across all login nodes
bash check-runners.sh

# Rebalance runners across all 11 nodes (also restarts offline runners)
bash rebalance-runners.sh              # dry run
APPLY=1 bash rebalance-runners.sh      # execute
APPLY=1 FORCE=1 bash rebalance-runners.sh  # move busy runners too

# Restart all runners in place (e.g. after a node reboot)
APPLY=1 bash restart-all.sh

# Restart one specific runner
bash restart-runner.sh login01 /path/to/runner-dir

# Move a runner to a different login node
bash move-runner.sh frontier-1 login01

# Stop and deregister a runner
bash stop-runner.sh frontier-12

# Deploy a new runner on a specific node
bash make-runner.sh 23 login01

# Deploy multiple runners across nodes (e.g. runners 23, 24, 25)
bash deploy-runners.sh 23 login01 login02 login03

# Rerun failed CI workflows
bash ../common/rerun-failed.sh
APPLY=1 bash ../common/rerun-failed.sh
```

## Scripts

| Script | Purpose |
|---|---|
| `config.sh` | Shared configuration: Frontier constants, `find_runner_dirs()`, and `sync_runner_nodes()`. Sources `../common/runner-lib.sh` for shared functions. |
| `check-runners.sh` | SSH to each login node, show Runner.Listener processes with name, status (idle/BUSY), slurm PATH, and RSS memory. |
| `list-runners.sh` | List all runners with GitHub API status, actual node (parallel SSH sweep), slurm status, and RSS. Flags stale `runner.node`. |
| `rebalance-runners.sh` | Sync node locations, compute optimal distribution, move runners across all 11 nodes. Handles offline runners. Dry run by default. |
| `restart-all.sh` | Restart all runners in place. Skips busy unless `FORCE=1`. Dry run by default. |
| `restart-runner.sh` | Stop and restart one runner on a given node. Usage: `restart-runner.sh <node> <runner-dir>` |
| `move-runner.sh` | Move a runner to a different login node by name. Usage: `move-runner.sh <runner-name> <target-node>` |
| `stop-runner.sh` | Stop the runner process and deregister from GitHub. Usage: `stop-runner.sh <runner-name>` |
| `make-runner.sh` | Download runner binary, register with GitHub, start on target node. Usage: `make-runner.sh <num> [node]` |
| `deploy-runners.sh` | Deploy multiple runners across nodes in parallel. Usage: `deploy-runners.sh <start-num> <node1> [node2 ...]` |
| `../common/rerun-failed.sh` | Rerun failed GitHub Actions workflows on open PRs and master. |

## Troubleshooting

**Runner goes OFFLINE repeatedly on the same node** — That login node may have
process culling or high memory pressure. Move it to a different node:
```bash
bash move-runner.sh frontier-1 login01
```

**Multiple runners OFFLINE at once** — Usually a transient OLCF network blip
to GitHub. Run `rebalance-runners.sh` to recover and redistribute all at once.

**Runner appears offline on GitHub but process is running** — GitHub status can
lag. `rebalance-runners.sh` uses EXE-based process discovery first: if a
process is found running, it will stop it before restarting, preventing
duplicate runner processes.
