## Description

Summarize your changes and the motivation behind them.

Fixes #(issue)

### Type of change

- [ ] Bug fix
- [ ] New feature
- [ ] Refactor
- [ ] Documentation
- [ ] Other: _describe_

## Testing

Describe how you tested your changes. List the compilers and platforms you used.

## Checklist

- [ ] New and existing tests pass locally (`./mfc.sh test`)
- [ ] I added or updated tests for new behavior
- [ ] I updated documentation if user-facing behavior changed

See the [developer guide](https://mflowcode.github.io/documentation/md_contributing.html) for full coding standards.

<details>
<summary><strong>GPU changes</strong> (expand if you modified code in <code>src/</code>)</summary>

- [ ] Code compiles with NVHPC (nvfortran)
- [ ] Code compiles with Cray (ftn)
- [ ] GPU results match CPU results
- [ ] Tested on NVIDIA GPU (V100/A100/H100) or AMD GPU (MI200+)
- [ ] Consider attaching an Nsight Systems or rocprof-systems profile if performance-sensitive
- [ ] Consider testing with multiple GPUs (e.g. 1, 2, 8) to verify scaling

</details>
