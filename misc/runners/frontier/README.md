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

`runner.node` is self-healing: `rebalance-runners` calls `sync_runner_nodes`
at startup, which sweeps all nodes in parallel and corrects any stale
`runner.node` files automatically.

Runners occasionally die due to OLCF's firewall/proxy dropping long-lived TCP
connections to GitHub's broker. Run `rebalance-runners` to restart and
redistribute them. Login nodes vary in stability — if a runner keeps dying on a
particular node, move it to a quieter one (login01 tends to have low load).

All commands are run via the dispatcher at `misc/runners/runner.sh`:

```bash
bash misc/runners/runner.sh frontier <command> [args...]
```

## Quick Reference

```bash
R="bash misc/runners/runner.sh frontier"

# List all runners with GitHub status, node, slurm, and memory usage
$R list-runners

# Check runner health across all login nodes
$R check-runners

# Rebalance runners across all 11 nodes (also restarts offline runners)
$R rebalance-runners              # dry run
APPLY=1 $R rebalance-runners     # execute
APPLY=1 FORCE=1 $R rebalance-runners  # move busy runners too

# Restart all runners in place (e.g. after a node reboot)
APPLY=1 $R restart-all

# Restart one specific runner
$R restart-runner login01 /path/to/runner-dir

# Move a runner to a different login node
$R move-runner frontier-1 login01

# Stop and deregister a runner
$R stop-runner frontier-12

# Deploy a new runner on a specific node
$R make-runner 23 login01

# Deploy multiple runners across nodes (e.g. runners 23, 24, 25)
$R deploy-runners 23 login01 login02 login03

# Rerun failed CI workflows
$R rerun-failed
APPLY=1 $R rerun-failed
```

## Scripts

| Script | Purpose |
|---|---|
| `config.sh` | Shared configuration: Frontier constants, `find_runner_dirs()`, and `sync_runner_nodes()`. Sources `../common/runner-lib.sh` for shared functions. |
| `make-runner.sh` | Download runner binary, register with GitHub, start on target node. Usage: `make-runner <num> [node]` |
| `deploy-runners.sh` | Deploy multiple runners across nodes in parallel. Usage: `deploy-runners <start-num> <node1> [node2 ...]` |
| `../common/` | All other commands (`check-runners`, `list-runners`, `rebalance-runners`, etc.) live here and are dispatched via `misc/runners/runner.sh`. |

## Troubleshooting

**Runner goes OFFLINE repeatedly on the same node** — That login node may have
process culling or high memory pressure. Move it to a different node:
```bash
bash misc/runners/runner.sh frontier move-runner frontier-1 login01
```

**Multiple runners OFFLINE at once** — Usually a transient OLCF network blip
to GitHub. Run `rebalance-runners` to recover and redistribute all at once.

**Runner appears offline on GitHub but process is running** — GitHub status can
lag. `rebalance-runners` uses EXE-based process discovery first: if a
process is found running, it will stop it before restarting, preventing
duplicate runner processes.
