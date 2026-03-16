# Phoenix Runner Management Scripts

Scripts for managing GitHub Actions self-hosted runners on Georgia Tech Phoenix
login nodes. The runners submit SLURM jobs to Phoenix compute nodes for MFC's
CI/CD pipeline.

## Background

Phoenix has 3 physical login nodes (`login-phoenix-gnr-{1,2,3}`), each with a
**4 GB per-user cgroup memory limit**. Each runner uses ~60-130 MB, so
distributing them evenly (~3-4 per node) is important to avoid OOM kills.

Runner binaries live on shared storage, so moving a runner between nodes only
requires stopping the process on one node and starting it on another.

Runners must be started with a **login shell** (`bash -l`) so they inherit
`/opt/slurm/current/bin` in PATH (required for `sbatch`, `squeue`, `sacct`).

## Quick Reference

```bash
# Check health (quick, one SSH per node)
bash check-runners.sh

# Detailed table with GitHub API status
bash list-runners.sh

# Auto-rebalance across nodes (also restarts offline runners)
bash rebalance-runners.sh              # dry run
APPLY=1 bash rebalance-runners.sh      # execute

# Restart all runners in place (e.g. after a node reboot)
APPLY=1 bash restart-all.sh

# Restart one specific runner
bash restart-runner.sh login-phoenix-gnr-2 /path/to/actions-runner-3

# Move a runner to a different login node
bash move-runner.sh phoenix-3 login-phoenix-gnr-1

# Stop and deregister a runner
bash stop-runner.sh phoenix-3

# Create a new runner (needs gh CLI with admin:org scope)
bash create-runner.sh phoenix-11 login-phoenix-gnr-2

# Rerun failed CI on open PRs
bash rerun-failed.sh                   # dry run
APPLY=1 bash rerun-failed.sh           # execute
```

## Scripts

| Script | Purpose |
|---|---|
| `config.sh` | Shared config: Phoenix constants (`ORG`, `RUNNER_GROUP`, `RUNNER_LABEL`, `NODES`, `CGROUP_LIMIT`, `RUNNER_PARENT_DIRS`) and `find_runner_dirs()`. Sources `../common/runner-lib.sh` for shared functions. |
| `check-runners.sh` | Quick per-node health check. One SSH per node. Shows runner names, idle/BUSY status, slurm PATH, RSS, and total cgroup memory. |
| `list-runners.sh` | Table combining GitHub API status with live node info from a parallel SSH sweep. Shows slurm status and flags stale `runner.node` entries. |
| `rebalance-runners.sh` | Auto-compute optimal distribution and move runners. Prefers idle runners. Also places OFFLINE runners. Dry run by default. |
| `restart-all.sh` | Restart all runners in place. Skips BUSY runners unless `FORCE=1`. Dry run by default. |
| `restart-runner.sh` | Stop and restart one runner on a given node. Usage: `restart-runner.sh <node> <runner-dir>` |
| `move-runner.sh` | Move a runner to a different login node by name. Usage: `move-runner.sh <runner-name> <target-node>` |
| `stop-runner.sh` | Stop the runner process and deregister from GitHub. Usage: `stop-runner.sh <runner-name>` |
| `create-runner.sh` | Download runner binary, register with GitHub via API, start on target node. Usage: `create-runner.sh <name> <node> [parent-dir]` |
| `rerun-failed.sh` | Scan open non-draft PRs and master for failed workflows, rerun failed jobs only. |

## Safety

- **Dry run by default**: `rebalance-runners.sh`, `restart-all.sh`, and
  `rerun-failed.sh` show what they would do unless `APPLY=1` is set.
- **Busy runner protection**: Scripts skip BUSY runners unless `FORCE=1`.
- **Slurm PATH verification**: After starting, scripts verify `slurm` appears
  in the runner's PATH and warn if missing.

## Configuration

Edit `config.sh` to change:

- `ORG` / `RUNNER_GROUP` / `RUNNER_LABEL` — GitHub org and runner registration
- `NODES` — physical login node hostnames
- `CGROUP_LIMIT` — per-user memory limit in MB
- `RUNNER_PARENT_DIRS` — directories containing `actions-runner-*/` installations

## Troubleshooting

**"sbatch: command not found"** — Runner started without login shell.
Fix: `bash restart-runner.sh <node> <runner-dir>`

**OOM kills** — Too many runners on one node.
Fix: `bash check-runners.sh` then `APPLY=1 bash rebalance-runners.sh`

**Runner OFFLINE** — Process died or node rebooted.
Fix: `APPLY=1 bash rebalance-runners.sh` (auto-places on least-loaded node)

**All runners down** — Node maintenance.
Fix: `APPLY=1 bash restart-all.sh`
