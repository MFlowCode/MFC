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

When required, run the benchmarks on a GPU yourself. A pull request without these results is incomplete.

GPU and compiler:

Base commit:

PR commit:

<!-- Replace this line with the contents of build/speed.md -->

If a case is below 0.95, run the commands again from the `for` loop, paste both tables, and explain it only if it stays below. If no benchmark case uses your change, add master and PR grind times for a case that does.

<details>
<summary>Commands (run from your PR checkout, inside a GPU allocation)</summary>

```bash
# Commit your changes first, so the PR commit id matches what you measure.
git fetch https://github.com/MFlowCode/MFC.git master
base=$(git merge-base FETCH_HEAD HEAD)
# Creates ../mfc-base, or moves one left from an earlier run to this base.
if [ -d ../mfc-base ]; then
    git -C ../mfc-base checkout --detach $base
else
    git worktree add --detach ../mfc-base $base
fi

# <slug> is from toolchain/modules; <computer> names a file in toolchain/templates. On Frontier they are f and frontier.
# Elsewhere, load your own modules and drop both -c options.
source ./mfc.sh load -c <slug> -m g

# One rank on every GPU in the allocation, as CI does: sets $gpu_ids and $ngpus.
source .github/scripts/detect-gpus.sh

# Use --gpu mp for OpenMP offload.
for tree in ../mfc-base .; do
    (cd $tree &&
        rm -f build/bench.yaml &&
        ./mfc.sh build --gpu acc -j 8 &&
        ./mfc.sh bench --mem 4 -o build/bench.yaml -- -c <computer> --gpu acc -g $gpu_ids -n $ngpus)
done

rm -f build/speed.md
GITHUB_STEP_SUMMARY=build/speed.md ./mfc.sh bench_diff ../mfc-base/build/bench.yaml build/bench.yaml

echo "Base commit: $(git -C ../mfc-base rev-parse --short HEAD)"
echo "PR commit: $(git rev-parse --short HEAD)"
```

</details>
