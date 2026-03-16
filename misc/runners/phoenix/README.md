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

All commands are run via the dispatcher at `misc/runner.sh`:

```bash
bash misc/runner.sh phoenix <command> [args...]
```

## Quick Reference

```bash
R="bash misc/runner.sh phoenix"

# Check health (quick, one SSH per node)
$R check-runners

# Detailed table with GitHub API status
$R list-runners

# Auto-rebalance across nodes (also restarts offline runners)
$R rebalance-runners              # dry run
APPLY=1 $R rebalance-runners     # execute

# Restart all runners in place (e.g. after a node reboot)
APPLY=1 $R restart-all

# Restart one specific runner
$R restart-runner login-phoenix-gnr-2 /path/to/actions-runner-3

# Move a runner to a different login node
$R move-runner phoenix-3 login-phoenix-gnr-1

# Stop and deregister a runner
$R stop-runner phoenix-3

# Create a new runner (needs gh CLI with admin:org scope)
$R create-runner phoenix-11 login-phoenix-gnr-2

# Rerun failed CI on open PRs
$R rerun-failed                   # dry run
APPLY=1 $R rerun-failed           # execute
```

## Scripts

| Script | Purpose |
|---|---|
| `config.sh` | Shared config: Phoenix constants (`ORG`, `RUNNER_GROUP`, `RUNNER_LABEL`, `NODES`, `CGROUP_LIMIT`, `RUNNER_PARENT_DIRS`) and `find_runner_dirs()`. Sources `../common/runner-lib.sh` for shared functions. |
| `create-runner.sh` | Download runner binary, register with GitHub via API, start on target node. Usage: `create-runner <name> <node> [parent-dir]` |
| `../common/` | All other commands (`check-runners`, `list-runners`, `rebalance-runners`, etc.) live here and are dispatched via `misc/runner.sh`. |

## Safety

- **Dry run by default**: `rebalance-runners`, `restart-all`, and
  `rerun-failed` show what they would do unless `APPLY=1` is set.
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
Fix: `bash misc/runner.sh phoenix restart-runner <node> <runner-dir>`

**OOM kills** — Too many runners on one node.
Fix: `bash misc/runner.sh phoenix check-runners` then `APPLY=1 bash misc/runner.sh phoenix rebalance-runners`

**Runner OFFLINE** — Process died or node rebooted.
Fix: `APPLY=1 bash misc/runner.sh phoenix rebalance-runners` (auto-places on least-loaded node)

**All runners down** — Node maintenance.
Fix: `APPLY=1 bash misc/runner.sh phoenix restart-all`
