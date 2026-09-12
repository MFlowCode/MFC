@page testing Testing

## Testing

To run MFC's test suite, run
```shell
./mfc.sh test -j <thread count>
```

It will generate and run test cases, comparing their output to previous runs from versions of MFC considered accurate.
*golden files*, stored in the `tests/` directory contain this data, aggregating `.dat` files generated when running MFC.
A test is considered passing when our error tolerances are met in order to maintain a high level of stability and accuracy.
`./mfc.sh test` has the following unique options:
- `-l` outputs the full list of tests
- `--from` (`-f)` and `--to` (`t`) restrict testing to a range of contiguous slugs
- `--only` (`-o`) restricts testing to a non-contiguous range of tests whose trace contains a given whole trace element (see [Selection and Execution Pitfalls](#selection-and-execution-pitfalls) for the exact matching rules)
- `--test-all` (`a`) test post process and ensure the Silo database files are correct
- `--percent` (`%`) to specify a percentage of the test suite to select at random and test
- `--max-attempts` (`-m`) the maximum number of attempts to make on a test before considering it failed
- `--no-examples` skips the testing of cases in the examples folder
- `--rdma-mpi` runs additional tests where RDMA MPI is enabled.

To specify a computer, pass the `-c` flag to `./mfc.sh run` like so:
```shell
./mfc.sh test -j <thread count> -- -c <computer name>
```
where `<computer name>` could be `phoenix` or any of the others in the [templates](https://github.com/MFlowCode/MFC/tree/master/toolchain/templates)).
You can create new templates with the appropriate run commands or omit this option.
The use of `--` in the above command passes options to the `./mfc.sh run` command underlying the `./mfc.sh test`.

### Creating Tests

Creating and updating test cases can be done with the following command line arguments:
- `--generate` to generate golden files for a new test case
- `--add-new-variables` to similar to `--generate`, but rather than generating a golden file from scratch, it generates a gold file with new variables for an updated test without changing the original golden file values.
- `--remove-old-tests` to remove the directories of tests that no longer exist

It is recommended that a range be specified when generating golden files for new test cases, as described in the previous section, in an effort not to regenerate the golden files of existing test cases.

Adding a new test case can be done by modifying [cases.py](https://github.com/MFlowCode/MFC/tree/master/toolchain/mfc/test/cases.py).
The function `list_cases` is responsible for generating the list of test cases.
Loops and conditionals are used to vary parameters, whose defaults can be found in the `BASE_CFG` case object within [case.py](https://github.com/MFlowCode/MFC/tree/master/toolchain/mfc/test/case.py).
The function operates on two variables:

- `stack`: A stack that holds the variations to the default case parameters.
By pushing and popping the stack inside loops and conditionals, it is easier to nest test case descriptions, as it holds the variations that are common to all future test cases within the same indentation level (in most scenarios).

- `cases`: A list that holds fully-formed `Case` objects, that will be returned at the end of the function.

Internally a test case is described as:
```python
@dataclasses.dataclass(init=False)
class Case:
    trace:  str
    params: dict
    ppn:    int
```

where:
- The `trace` is a string that contains a human-readable description of what parameters were varied, or more generally what the case is meant to test.
**Each `trace` must be distinct.**
- `params` is the fully resolved case dictionary, as would appear in a Python case input file.
- `ppn` is the number of processes per node to use when running the case.

To illustrate, consider the following excerpt from `list_cases`:

```python
for weno_order in [3, 5]:
  stack.push(f"weno_order={weno_order}", {'weno_order': weno_order})

  for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
      stack.push(f"(mapped_weno={mapped_weno},mp_weno={mp_weno})", {
          'mapped_weno': mapped_weno,
          'mp_weno':     mp_weno
      })

      if not (mp_weno == 'T' and weno_order != 5):
          cases.append(define_case_d(stack, '', {}))

      stack.pop()

  stack.pop()
```

When pushing to the stack or creating a new case with the `define_case_d` function, you must specify:
- `stack`: The current stack.
- `trace`: A human-readable string describing what you are currently varying.
- `variations`: A Python dictionary with case parameter variations.
- (Optional) `ppn`: The number of processes per node to use (default is 1).

If a trace is empty (that is, the empty string `""`), it will not appear in the final trace, but any case parameter variations associated with it will still be applied.

Finally, the case is appended to the `cases` list, which will be returned by the `list_cases` function.

### Selection and Execution Pitfalls

Each of these fails quietly rather than loudly.

- **`--only` matches whole trace elements, not substrings**, and it ANDs labels while ORing
  UUIDs (`_filter_only` in `toolchain/mfc/test/test.py`). `--only bubbles` matches nothing,
  because the trace element is `Bubbles`; `--only low_Mach=1 low_Mach=2` asks for cases
  carrying both labels at once and also matches nothing. An empty selection then exits
  **143**, which reads like an external kill rather than an empty filter. Pass UUIDs when
  you want the union of several groups.
- **Sibling `define_case_d` calls at the same stack level are never combined.** Two switches
  that only matter together therefore get no effective coverage unless one is pushed onto
  the stack and the other defined beneath it — `avg_state=1`, for instance, is only read
  when `wave_speeds=2`. Check reachability before trusting that a flag is tested.
- **`--no-build` silently runs whatever binary is already on disk**, including one built for
  a different configuration. Chemistry has its own configuration that a plain `./mfc.sh
  build` never produces, so a `--no-build` run can report failures from stale binaries and
  hide real compile breaks. Run chemistry-touching sets without it.
- **Identify the newest binary by the binary's own mtime**, not by its install directory's:
  a stale configuration's directory can be newer than a fresh build's.
- **The pre-commit hook lives in the main repository's `.git/hooks/`**, and git exports
  `GIT_DIR` there during a commit, so from a worktree the toolchain lint enumerates the
  other checkout and fails. Run `./mfc.sh precheck` by hand and commit with `--no-verify`.
- **`/tmp` is node-local.** Scratch does not survive a compute-node change, and its absence
  is silence rather than an error. Keep patches and resource baselines on a shared
  filesystem.
- **An unexplained golden-file difference is a bug report, not noise to be regenerated
  away.** Regenerate only the affected tests.

Tests are generated programmatically in `toolchain/mfc/test/cases.py`; a test's UUID is the
CRC32 of its trace string, and `./mfc.sh test -l` lists every one.

### Testing Post Process

To test the post-processing code, append the `-a` or `--test-all` option:
```shell
./mfc.sh test -a -j 8
```

This argument will re-run the test stack with ``parallel_io='T'``, which generates silo_hdf5 files.
It will also turn most write parameters (`*_wrt`) on.
Then, it searches through the silo files using `h5dump` to ensure that there are no `NaN`s or `Infinity`s.
Although adding this option does not guarantee that accurate `.silo` files are generated, it does ensure that the post-process code does not fail or produce malformed data.



<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-02-15</div>
