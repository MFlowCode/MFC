# Phoenix Runner Management Scripts

Scripts for managing GitHub Actions self-hosted runners on Georgia Tech Phoenix
login nodes. The runners submit SLURM jobs to Phoenix compute nodes for MFC's
CI/CD pipeline.

## Background

Phoenix has 3 physical login nodes (`login-phoenix-gnr-{1,2,3}`), each with a
**4 GB per-user cgroup memory limit**. Each runner process uses ~60–130 MB, so
distributing runners evenly (~3–4 per node) is important to avoid OOM kills.

Runner binaries live on shared storage (`/storage/scratch1/` or
`/storage/project/`), so moving a runner between nodes only requires stopping
the process on one node and starting it on another — no files are copied.

Runners must be started with a **login shell** (`bash -l`) so they inherit the
full PATH, which includes `/opt/slurm/current/bin` (required for `sbatch`,
`squeue`, `sacct`). Starting without a login shell causes "sbatch: command not
found" errors.

## Scripts

All scripts live in this directory and share configuration via `config.sh`.

### Quick Reference

```bash
# Check health of all runners
bash check-runners.sh

# Detailed table of all registered runners
bash list-runners.sh

# Restart a single runner
bash restart-runner.sh login-phoenix-gnr-2 /path/to/actions-runner-3

# Restart all runners in place (e.g. after a node reboot)
APPLY=1 bash restart-all.sh

# Auto-rebalance runners across nodes
bash rebalance-runners.sh              # dry run
APPLY=1 bash rebalance-runners.sh      # execute

# Create a new runner
bash create-runner.sh phoenix-11 login-phoenix-gnr-2 /path/to/parent-dir

# Rerun failed CI workflows on open PRs
bash rerun-failed.sh                   # dry run
APPLY=1 bash rerun-failed.sh           # execute
```

### Script Details

| Script | Purpose |
|---|---|
| `config.sh` | Shared configuration sourced by all other scripts. Defines login node list, cgroup limit, runner parent directories, and helper functions (`discover_runners`, `find_runner_node`, `start_runner`, `stop_runner`, `check_slurm_path`). |
| `check-runners.sh` | Quick per-node health check. Shows each runner's name, idle/BUSY status, whether slurm is in PATH, RSS memory, and total user memory per node. Lightweight (one SSH call per node). |
| `list-runners.sh` | Detailed table of all registered runners discovered from `.runner` config files on shared storage. Shows GitHub runner name, which node it's on (or OFFLINE), status, slurm PATH, RSS, pool, and directory. Slower than `check-runners.sh` (one SSH per runner). |
| `restart-runner.sh` | Stop and restart a single runner on a given node. Uses `setsid bash -lc` for proper PATH and `</dev/null` for SSH detachment. Verifies slurm is in PATH after start. |
| `restart-all.sh` | Restart all runners on their current nodes. Skips BUSY runners unless `FORCE=1`. Useful after node reboots or to pick up environment changes. |
| `rebalance-runners.sh` | Auto-computes the optimal runner distribution (equal per node), identifies which runners to move, and executes the moves. Prefers moving idle runners over busy ones. Also places OFFLINE runners on underloaded nodes. |
| `create-runner.sh` | Downloads the latest GitHub Actions runner binary, registers it with the GitHub org using `gh` CLI, and starts it on the target node. Requires `gh` CLI with org admin access. |
| `rerun-failed.sh` | Scans open non-draft PRs and the master branch for failed GitHub Actions workflows and reruns the failed jobs (not the entire workflow). Requires `gh` CLI. |

### Safety Features

- **Dry run by default**: `rebalance-runners.sh`, `restart-all.sh`, and
  `rerun-failed.sh` only show what they would do unless `APPLY=1` is set.
- **Busy runner protection**: Scripts that stop runners skip BUSY runners
  (those with an active `Runner.Worker` process) unless `FORCE=1` is set.
- **Slurm PATH verification**: After starting a runner, scripts verify that
  `/opt/slurm` is in the process's PATH and warn if it's missing.

## Configuration

Edit `config.sh` to update:

- **`NODES`** — list of physical login node hostnames
- **`CGROUP_LIMIT`** — per-user memory limit in MB (currently 4096)
- **`RUNNER_PARENT_DIRS`** — directories containing `actions-runner-*/` installations

## Troubleshooting

**Runner shows "sbatch: command not found"**
The runner was started without a login shell. Restart it:
```bash
bash restart-runner.sh <node> <runner-dir>
```

**OOM kills during CI jobs**
Check memory distribution with `bash check-runners.sh`. If one node is
overloaded, rebalance: `APPLY=1 bash rebalance-runners.sh`.

**Runner shows OFFLINE in `list-runners.sh`**
The runner process is not running on any node. Either the node rebooted or the
process was killed. Restart it on an appropriate node:
```bash
bash restart-runner.sh login-phoenix-gnr-2 /path/to/runner-dir
```
Or let `rebalance-runners.sh` place it automatically:
```bash
APPLY=1 bash rebalance-runners.sh
```

**All runners down after node maintenance**
```bash
APPLY=1 bash restart-all.sh
```
