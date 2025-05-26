## Description

Please include a summary of the changes and the related issue(s) if they exist.
Please also include relevant motivation and context.

Fixes #(issue) [optional]

### Type of change

Please delete options that are not relevant.

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Something else

### Scope

- [ ] This PR comprises a set of related changes with a common goal

If you cannot check the above box, please split your PR into multiple PRs that each have a common goal.

## How Has This Been Tested?

Please describe the tests that you ran to verify your changes. 
Provide instructions so we can reproduce. 
Please also list any relevant details for your test configuration

- [ ] Test A
- [ ] Test B

**Test Configuration**:

* What computers and compilers did you use to test this:

## Checklist

- [ ] I have added comments for the new code
- [ ] I added Doxygen docstrings to the new code
- [ ] I have made corresponding changes to the documentation (`docs/`)
- [ ] I have added regression tests to the test suite so that people can verify in the future that the feature is behaving as expected
- [ ] I have added example cases in `examples/` that demonstrate my new feature performing as expected.
They run to completion and demonstrate "interesting physics"
- [ ] I ran `./mfc.sh format` before committing my code
- [ ] New and existing tests pass locally with my changes, including with GPU capability enabled (both NVIDIA hardware with NVHPC compilers and AMD hardware with CRAY compilers) and disabled
- [ ] This PR does not introduce any repeated code (it follows the [DRY](https://en.wikipedia.org/wiki/Don%27t_repeat_yourself) principle)
- [ ] I cannot think of a way to condense this code and reduce any introduced additional line count

### If your code changes any code source files (anything in `src/simulation`)

To make sure the code is performing as expected on GPU devices, I have:
- [ ] Checked that the code compiles using NVHPC compilers
- [ ] Checked that the code compiles using CRAY compilers
- [ ] Ran the code on either V100, A100, or H100 GPUs and ensured the new feature performed as expected (the GPU results match the CPU results)
- [ ] Ran the code on MI200+ GPUs and ensure the new features performed as expected (the GPU results match the CPU results)
- [ ] Enclosed the new feature via `nvtx` ranges so that they can be identified in profiles
- [ ] Ran a Nsight Systems profile using `./mfc.sh run XXXX --gpu -t simulation --nsys`, and have attached the output file (`.nsys-rep`) and plain text results to this PR
- [ ] Ran a Rocprof Systems profile using `./mfc.sh run XXXX --gpu -t simulation --rsys --hip-trace`, and have attached the output file and plain text results to this PR.
- [ ] Ran my code using various numbers of different GPUs (1, 2, and 8, for example) in parallel and made sure that the results scale similarly to what happens if you run without the new code/feature
