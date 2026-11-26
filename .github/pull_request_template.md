## Contribution Policy

All contributions are expected to demonstrate:
- A clear understanding of the codebase
- Alignment with product direction
- Thoughtful reasoning behind changes
- Evidence of real-world usage or hands-on experience with the problem

If these expectations are not met, we would prefer to implement the changes ourselves rather than spend time reviewing low-effort submissions.

## What to fill in

| Your change touches | 1. Problem | 2. Tests to run | 3. Speed |
|---|---|---|---|
| Docs only | N/A | N/A | N/A |
| Anything else (toolchain, CI, tests, examples) | required | N/A | N/A |
| `src/` | required | `./mfc.sh test --select-enforce` | required |
| Build flags or compiler versions (`CMakeLists.txt`, `cmake/`, `toolchain/cmake/`, `toolchain/modules`) | required | `./mfc.sh test` | required |

If more than one row applies, do everything each of them requires. Write N/A with a reason; don't leave blanks.

## 1. Does it solve the problem?

Command that shows it:

| Quantity | master | this PR | expected |
|---|---|---|---|
|  |  |  |  |

## 2. Does it create another one?

| Check | Result |
|---|---|
| Tests: command, number passed, number failed |  |
| Golden files changed: count, largest difference, why |  |
| Not tested |  |

## 3. Is it slower?

Required for `src/` and build changes: benchmark your base commit and your PR on the same GPU.
Any supported GPU will do.

GPU and compiler:

Base commit:

PR commit:

<!-- Paste build/speed.md here -->

Explain any grind speedup below 0.95 that shows up twice.

<details>
<summary>Commands</summary>

```bash
# Run inside a GPU job with your compiler loaded (supported clusters: source ./mfc.sh load -c <slug> -m g).
(
set -e    # stop at the first failure
rm -f build/speed.md
git diff --quiet HEAD || { echo "Commit your changes first."; exit 1; }
git fetch https://github.com/MFlowCode/MFC.git master
base=$(git merge-base FETCH_HEAD HEAD) || { echo "No commit shared with master. Shallow clone? Run: git fetch --unshallow"; exit 1; }
base=$(git rev-parse --short $base)
[ $base != $(git rev-parse --short HEAD) ] || { echo "Nothing to compare: your commit is already in master."; exit 1; }
[ -d ../mfc-$base ] || git worktree add --detach ../mfc-$base $base
source .github/scripts/detect-gpus.sh

# NVIDIA or Cray: --gpu acc. AMD flang: --gpu mp. On a cluster, add -c <computer> after "--".
for tree in ../mfc-$base .; do
    (cd $tree &&
        ./mfc.sh build --gpu acc -j 8 &&
        ./mfc.sh bench --mem 4 -o build/bench.yaml -- --gpu acc -g $gpu_ids -n $ngpus)
done

GITHUB_STEP_SUMMARY=build/speed.md ./mfc.sh bench_diff ../mfc-$base/build/bench.yaml build/bench.yaml
echo "Base commit: $base"; echo "PR commit: $(git rev-parse --short HEAD)"
)
```

</details>
