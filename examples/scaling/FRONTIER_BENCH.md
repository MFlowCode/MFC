# Description

The scripts and case file in this directory are set up to benchmarking strong
and weak scaling performance as well as single device absolute performance on
OLCF Frontier. The case file is for a three dimensional, two fluid liquid--gas
problem without viscosity or surface tension. The scripts contained here have
been tested for the default node counts and problem sizes in the scripts. The
reference data in `reference.dat` also makes use of the default node counts and
problem sizes and will need to be regenerated if either changes. The benchmarks
can be run with the following steps:

## Getting the code

The code is hosted on GitHub and can be cloned with the following command:

```bash
git clone git@github.com:MFlowCode/MFC.git; cd MFC; chmod u+x examples/scaling/*.sh;
```

The above command clones the repository, changes directories in the repository
root, and makes the benchmark scripts executable.

## Running the benchmarks

### Step 1: Building

The code for the benchmarks is built with the following command
```
./examples/scaling/build.sh
```

### Step 2: Running

The benchmarks can be run in their default configuration with the following
```
./examples/scaling/submit_all.sh --account <account_name>
```
By default this will submit the following jobs for benchmarking

| Job                | Nodes | Description                                                         |
| ------------------ | ----- | ------------------------------------------------------------------- |
| `MFC-W-16-64`      | 16    | Weak scaling calculation with a ~64GB problem per GCD on 16 nodes   |
| `MFC-W-128-64`     | 128   | Weak scaling calculation with a ~64GB problem per GCD on 128 nodes  |
| `MFC-W-1024-64`    | 1024  | Weak scaling calculation with a ~64GB problem per GCD on 1024 nodes |
| `MFC-W-8192-64`    | 8192  | Weak scaling calculation with a ~64GB problem per GCD on 8192 nodes |
| `MFC-S-8-4096`     | 8     | Strong scaling calculation with a ~4096GB problem on 8 nodes        |
| `MFC-S-64-4096`    | 64    | Strong scaling calculation with a ~4096GB problem on 64 nodes       |
| `MFC-S-512-4096`   | 512   | Strong scaling calculation with a ~4096GB problem on 512 nodes      |
| `MFC-S-4096-4096`  | 4096  | Strong scaling calculation with a ~4096GB problem on 4096 nodes     |
| `MFC-G-8`          | 1     | Single device grind time calculation with ~8GB per GCD              |
| `MFC-G-16`         | 1     | Single device grind time calculation with ~16GB per GCD             |
| `MFC-G-32`         | 1     | Single device grind time calculation with ~32GB per GCD             |
| `MFC-G-64`         | 1     | Single device grind time calculation with ~64GB per GCD             |
Strong and weak scaling cases run `pre_process` once and then run `simulation`
with and without GPU-aware MPI in a single job. Individual benchmarks can be run
by calling the `submit_[strong,weak,grind].sh` scripts directly, or modifying
the `submit_all.sh` script to fit your needs.

#### Modifying the benchmarks
The submitted jobs can be modified by appending options to the `submit_all.sh`
script. For examples, appending
```
--nodes "1,2,4,8"
```
to the `submit_strong.sh` and `submit_weak.sh` scripts will run the strong and
weak scaling benchmarks on 1, 2, 4, and 8 nodes. Appending
```
--mem "x,y"
```
will modify the approximate problem size in terms of GB of memory
(see the `submit_[strong,weak,grind].sh` for details on what this number refers
to for the different types of tests).

### Step 3: Post processing

The log files can be post processed into a more human readable format with
```
python3 examples/scaling/analyze.py
```
This Python script generates a table of results in the command line with
comparison to the reference data in `reference.dat`. The `rel_perf` column
compares the raw run times of the current results to the reference data.
Relative performance numbers small than 1.0 indicate a speedup and numbers larger
than one indicate a slowdown relative to the reference data. The selected problem
sizes are intended to be comparable to the tiny, small, medium, and large labels
used by the SpecHPC benchmark.

## Common errors

The only common failure point identified during testing were "text file busy"
errors causing job failures. These errors are intermittent and are usually
resolved by resubmitting the test.
