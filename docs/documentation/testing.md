## Testing
 
To run MFC's test suite, run
```shell
./mfc.sh test -j <thread count>
```

It will generate and run test cases, comparing their output to that of previous runs from versions of MFC considered to be accurate.
*golden files*, stored in the `tests/` directory contain this data, by aggregating `.dat` files generated when running MFC.
A test is considered passing when our error tolerances are met, in order to maintain a high level of stability and accuracy.
Run `./mfc.sh test -h` for a full list of accepted arguments.

Most notably, you can consult the full list of tests by running
```
./mfc.sh test -l
```

To restrict to a given range, use the `--from` (`-f`) and `--to` (`-t`) options.
To run a (non-contiguous) subset of tests, use the `--only` (`-o`) option instead.

### Creating Tests

To (re)generate *golden files*, append the `--generate` option:
```shell
./mfc.sh test --generate -j 8
```

It is recommended that a range be specified when generating golden files for new test cases, as described in the previous section, in an effort not to regenerate the golden files of existing test cases.

**Note:** If you output new variables and want to update the golden files to include these without modifying the original data, use the `--add-new-variables` option instead.

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

When pushing to the stack, or creating a new case with the `define_case_d` function, you must specify:
- `stack`: The current stack.
- `trace`: A human-readable string describing what you are currently varying.
- `variations`: A Python dictionary with case parameter variations.
- (Optional) `ppn`: The number of processes per node to use (default is 1).

If a trace is empty (that is, the empty string `""`), it will not appear in the final trace, but any case parameter variations associated with it will still be applied.

Finally, the case is appended to the `cases` list, which will be returned by the `list_cases` function.

### Testing Post Process

To test updated post process code, append the `-a` or `--test-all` option: 
```shell
./mfc.sh test -a -j 8
```

This argument will re-run the test stack with `parallel_io=True`, which generates silo_hdf5 files.
It will also turn most write parameters (`*_wrt`) on.
Then, it searches through the silo files using `h5dump` to ensure that there are no NaNs or Infinitys.
Although adding this option does not guarantee that accurate silo files are generated, it does ensure that post process does not fail or produce malformed data. 

