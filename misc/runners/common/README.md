# Common Runner Management Scripts

Site-agnostic scripts shared between the Frontier and Phoenix runner setups.
All shared logic lives here; site directories contain only site-specific files
(`config.sh` and scripts unique to that cluster).

Scripts are invoked via the dispatcher at `misc/runner.sh`:
```bash
bash misc/runner.sh <site> <command> [args...]
```

## Scripts

| Script | Purpose |
|---|---|
| `runner-lib.sh` | Shared library: GitHub API helpers, EXE-based process discovery, parallel node sweep, start/stop primitives. Sourced by site `config.sh` files. |
| `check-runners.sh` | Per-node health check: Runner.Listener processes with name, idle/BUSY, slurm PATH, RSS. Optional cgroup memory footer. |
| `list-runners.sh` | Full table: GitHub API status × parallel node sweep. Shows slurm status, flags stale `runner.node`. |
| `rebalance-runners.sh` | Compute optimal distribution and move runners across nodes. Handles offline runners. Writes `runner.node`. Dry run by default. |
| `restart-runner.sh` | Stop and restart one runner on a given node. Verifies slurm in PATH. Writes `runner.node`. |
| `restart-all.sh` | Restart all runners in place. Skips busy unless `FORCE=1`. Dry run by default. |
| `move-runner.sh` | Move a runner to a different login node by name. Stops on current node, starts on target. Writes `runner.node`. |
| `stop-runner.sh` | Stop a runner process and remove its GitHub registration. |
| `rerun-failed.sh` | Rerun failed GitHub Actions workflows on open non-draft PRs and master. Dry run by default. |
