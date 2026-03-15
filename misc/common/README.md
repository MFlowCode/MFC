# Common Runner Management Scripts

This directory contains site-agnostic scripts shared between the Phoenix and
Frontier runner management setups. Scripts here have no site-specific logic and
can be invoked directly or via thin site-specific wrappers.

## Scripts

| Script | Purpose |
|---|---|
| `rerun-failed.sh` | Rerun failed GitHub Actions workflows on open non-draft MFC PRs and master. Dry-run by default; set `APPLY=1` to actually trigger reruns. |

## Usage

```bash
# Dry run — show which failed workflows would be rerun
bash misc/common/rerun-failed.sh

# Actually rerun failed workflows
APPLY=1 bash misc/common/rerun-failed.sh
```

## Site wrappers

`misc/phoenix/rerun-failed.sh` is a thin wrapper that delegates to this
script, so both `bash misc/phoenix/rerun-failed.sh` and
`bash misc/common/rerun-failed.sh` invoke the same logic.
