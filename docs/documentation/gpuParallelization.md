# GPU Parallelization

MFC compiles GPU code via OpenACC and in the future OpenMP as well.

In order to swap between OpenACC and OpenMP, custom GPU macros are used that translate to equivalent OpenACC and OpenMP directives.
FYPP is used to process the GPU macros.

[OpenACC Quick start Guide](https://openacc-best-practices-guide.readthedocs.io/en/latest/01-Introduction.html)

[OpenACC API Documentation](https://www.openacc.org/sites/default/files/inline-files/API%20Guide%202.7.pdf)

------------------------------------------------------------------------------------------

## Macro API Documentation

Note: Ordering is not guaranteed or stable, so use key-value pairing when using macros

### Data Type Meanings

- Integer is a number

- Boolean is a pythonic boolean - Valid options: `True` or `False`

- String List is given as a comma separated list surrounding by brackets and inside quotations
  - Ex: ``'[hello, world, Fortran]'``

- 2-level string list is given as a comma separated list of string lists surrounding by brackets and inside quotations
  - Ex: ``'[[hello, world], [Fortran, MFC]]'`` or ``'[[hello]]'``

### Data Flow

- Data on the GPU has a reference counter
- When data is referred to being allocated, it means that GPU memory is allocated if it is not already present in GPU memory. If a variable is already present, the reference counter is just incremented.
- When data is referred to being dellocated, it means that the reference counter is decremented. If the reference counter is zero, then the data is actually deallocated from GPU memory
- When data is referred to being attached, it means that the device pointer attaches to target if it not already attached. If pointer is already attached, then the attachment counter is just incremented
- When data is referred to being detached, it means that the attachment counter is decremented. If attachment counter is zero, then actually detached

------------------------------------------------------------------------------------------

### Computation Macros

<details>
  <summary><code>GPU_PARALLEL_LOOP</code> -- <code>(Execute the following loop on the GPU in parallel)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_PARALLEL_LOOP(...)`

**Parameters**

| name             | data type           | Default Value     | description                                                                               |
|------------------|---------------------|-------------------|-------------------------------------------------------------------------------------------|
| `collapse`       | integer             | None              | Number of loops to combine into 1 loop                                                    |
| `parallelism`    | string list         | '\[gang,vector\]' | Parallelism granularity to use for this loop                                              |
| `default`        | string              | 'present'         | Implicit assumptions compiler should make                                                 |
| `private`        | string list         | None              | Variables that are private to each iteration/thread                                       |
| `firstprivate`   | string list         | None              | Initialized variables that are private to each iteration/thread                           |
| `reduction`      | 2-level string list | None              | Variables unique to each iteration and reduced at the end                                 |
| `reductionOp`    | string list         | None              | Operator that each list of reduction will reduce with                                     |
| `copy`           | string list         | None              | Allocates and copies data to GPU on entrance, then deallocated and copies to CPU on exit  |
| `copyin`         | string list         | None              | Allocates and copies data to GPU on entrance and then deallocated on exit                 |
| `copyinReadOnly` | string list         | None              | Allocates and copies readonly data to GPU and then deallocated on exit                    |
| `copyout`        | string list         | None              | Allocates data on GPU on entrance and then deallocates and copies to CPU on exit          |
| `create`         | string list         | None              | Allocates data on GPU on entrance and then deallocates on exit                            |
| `no_create`      | string list         | None              | Use data in CPU memory unless data is already in GPU memory                               |
| `present`        | string list         | None              | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit |
| `deviceptr`      | string list         | None              | Pointer variables that are already allocated on GPU memory                                |
| `attach`         | string list         | None              | Attaches device pointer to device targets on entrance, then detach on exit                |
| `extraAccArgs`   | string              | None              | String of any extra arguments added to the OpenACC directive                              |

**Parameter Restrictions**

| name          | Restricted range                                  |
|---------------|---------------------------------------------------|
| `collapse`    | Must be greater than 1                            |
| `parallelism` | Valid elements: 'gang', 'worker', 'vector', 'seq' |
| `default`     | 'present' or 'none'                               |

**Additional information**

- default present means that the any non-scalar data in assumed to be present on the GPU
- default none means that the compiler should not implicitly determine the data attributes for any variable
- reduction and reductionOp must match in length
- With ``reduction='[[sum1, sum2], [largest]]'`` and ``reductionOp='[+, max]'``, `sum1` and `sum2` will be the sum of sum1/sum2 in each loop iteration, and `largest` will the maximum value of `largest` all the loop iterations
- A reduction implies a copy, so it does not need to be added for both

**Example**

```python
 $:GPU_PARALLEL_LOOP(collapse=3, private='[tmp, r]', reduction='[[vol, avg], [max_val]]', reductionOp='[+, MAX]')
 $:GPU_PARALLEL_LOOP(collapse=2, private='[sum_holder]', copyin='[starting_sum]', copyout='[eigenval,C]')
```

</details>

<details>
  <summary><code>GPU_LOOP</code> -- <code>(Execute loop on GPU)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_LOOP(...)`

**Parameters**

| name              | data type           | Default Value | description                                                                                      |
|-------------------|---------------------|---------------|--------------------------------------------------------------------------------------------------|
| `collapse`        | integer             | None          | Number of loops to combine into 1 loop                                                           |
| `parallelism`     | string list         | None          | Parallelism granularity to use for this loop                                                     |
| `data_dependency` | string              | None          | 'independent'-> assert loop iterations are independent, 'auto->let compiler analyze dependencies |
| `private`         | string list         | None          | Variables that are private to each iteration/thread                                              |
| `reduction`       | 2-level string list | None          | Variables unique to each iteration and reduced at the end                                        |
| `reductionOp`     | string list         | None          | Operator that each list of reduction will reduce with                                            |
| `extraAccArgs`    | string              | None          | String of any extra arguments added to the OpenACC directive                                     |

**Parameter Restrictions**

| name              | Restricted range                                  |
|-------------------|---------------------------------------------------|
| `collapse`        | Must be greater than 1                            |
| `parallelism`     | Valid elements: 'gang', 'worker', 'vector', 'seq' |
| `data_dependency` | 'auto' or 'independent'                           |

**Additional information**

- Loop parallelism is most commonly ``'[seq]'``
- reduction and reductionOp must match in length
- With ``reduction='[[sum1, sum2], [largest]]'`` and ``reductionOp='[+, max]'``, `sum1` and `sum2` will be the sum of sum1/sum2 in each loop iteration, and `largest` will the maximum value of `largest` all the loop iterations

**Example**

```python
 $:GPU_LOOP(parallelism='[seq]')
 $:GPU_LOOP(collapse=3, parallelism='[seq]',private='[tmp, r]')
```

</details>

<details>
  <summary><code>GPU_PARALLEL</code> -- <code>(Execute the following on the GPU in parallel)</code></summary>

**Macro Invocation**

Uses FYPP call directive using `#:call`

```C
#:call GPU_PARALLEL(...)
   {code}
#:endcall GPU_PARALLEL 
```

**Parameters**

| name             | data type           | Default Value     | description                                                                               |
|------------------|---------------------|-------------------|-------------------------------------------------------------------------------------------|
| `default`        | string              | 'present'         | Implicit assumptions compiler should make                                                 |
| `private`        | string list         | None              | Variables that are private to each iteration/thread                                       |
| `firstprivate`   | string list         | None              | Initialized variables that are private to each iteration/thread                           |
| `reduction`      | 2-level string list | None              | Variables unique to each iteration and reduced at the end                                 |
| `reductionOp`    | string list         | None              | Operator that each list of reduction will reduce with                                     |
| `copy`           | string list         | None              | Allocates and copies data to GPU on entrance, then deallocated and copies to CPU on exit  |
| `copyin`         | string list         | None              | Allocates and copies data to GPU on entrance and then deallocated on exit                 |
| `copyinReadOnly` | string list         | None              | Allocates and copies readonly data to GPU and then deallocated on exit                    |
| `copyout`        | string list         | None              | Allocates data on GPU on entrance and then deallocates and copies to CPU on exit          |
| `create`         | string list         | None              | Allocates data on GPU on entrance and then deallocates on exit                            |
| `no_create`      | string list         | None              | Use data in CPU memory unless data is already in GPU memory                               |
| `present`        | string list         | None              | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit |
| `deviceptr`      | string list         | None              | Pointer variables that are already allocated on GPU memory                                |
| `attach`         | string list         | None              | Attaches device pointer to device targets on entrance, then detach on exit                |
| `extraAccArgs`   | string              | None              | String of any extra arguments added to the OpenACC directive                              |

**Parameter Restrictions**

| name          | Restricted range                                  |
|---------------|---------------------------------------------------|
| `default`     | 'present' or 'none'                               |

**Additional information**

- default present means that the any non-scalar data in assumed to be present on the GPU
- default none means that the compiler should not implicitly determine the data attributes for any variable
- reduction and reductionOp must match in length
- With ``reduction='[[sum1, sum2], [largest]]'`` and ``reductionOp='[+, max]'``, `sum1` and `sum2` will be the sum of sum1/sum2 in each loop iteration, and `largest` will the maximum value of `largest` all the loop iterations
- A reduction implies a copy, so it does not need to be added for both

**Example**

```C
 #:call GPU_PARALLEL()
      {code}
      ...
 #:endcall GPU_PARALLEL
 #:call GPU_PARALLEL(create='[pixel_arr]', copyin='[initial_index]')
      {code}
      ...
 #:endcall
```

</details>

------------------------------------------------------------------------------------------

### Data Control Macros

<details>
 <summary><code>GPU_DATA</code> -- <code>(Make data accessible on GPU in specified region)</code></summary>

**Macro Invocation**

Uses FYPP call directive using `#:call`

```C
#:call GPU_DATA(...)
   {code}
#:endcall GPU_DATA 
```

**Parameters**

| name             | data type   | Default Value | description                                                                                  |
|------------------|-------------|---------------|----------------------------------------------------------------------------------------------|
| `code`           | code        | Required      | Region of code where defined data is accessible                                              |
| `copy`           | string list | None          | Allocates and copies variable to GPU on entrance, then deallocated and copies to CPU on exit |
| `copyin`         | string list | None          | Allocates and copies data to GPU on entrance and then deallocated on exit                    |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU and then deallocated on exit                 |
| `copyout`        | string list | None          | Allocates data on GPU on entrance and then deallocates and copies to CPU on exit             |
| `create`         | string list | None          | Allocates data on GPU on entrance and then deallocates on exit                               |
| `no_create`      | string list | None          | Use data in CPU memory unless data is already in GPU memory                                  |
| `present`        | string list | None          | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit    |
| `deviceptr`      | string list | None          | Pointer variables that are already allocated on GPU memory                                   |
| `attach`         | string list | None          | Attaches device pointer to device targets on entrance, then detach on exit                   |
| `default`        | string      | None          | Implicit assumptions compiler should make                                                    |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive                                 |

**Parameter Restrictions**

| name   | Restricted range                                 |
|--------|--------------------------------------------------|
| `code` | Do not assign it manually with key-value pairing |

**Example**

```C
 #:call GPU_DATA(copy='[pixel_arr]', copyin='[starting_pixels, initial_index]',attach='[p_real, p_cmplx, p_fltr_cmplx]')
      {code}
      ...
 #:endcall GPU_DATA
 #:call GPU_DATA(create='[pixel_arr]', copyin='[initial_index]')
      {code}
      ...
 #:endcall
```

</details>

<details>
 <summary><code>GPU_ENTER_DATA</code> -- <code>(Allocate/move data to GPU until matching GPU_EXIT_DATA or program termination)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_ENTER_DATA(...)`

**Parameter**

| name             | data type   | Default Value | description                                                  |
|------------------|-------------|---------------|--------------------------------------------------------------|
| `copyin`         | string list | None          | Allocates and copies data to GPU on entrance                 |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU on entrance  |
| `create`         | string list | None          | Allocates data on GPU on entrance                            |
| `attach`         | string list | None          | Attaches device pointer to device targets on entrance        |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive |

**Example**

```python
 $:GPU_ENTER_DATA(copyin='[pixels_arr]', copyinReadOnly='[starting_pixels, initial_index]')
 $:GPU_ENTER_DATA(create='[bc_buffers(1:num_dims, -1:1)]', copyin='[initial_index]')
```

</details>

<details>
 <summary><code>GPU_EXIT_DATA</code> -- <code>(Deallocate/move data from GPU created by GPU_ENTER_DATA)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_EXIT_DATA(...)`

**Parameters**

| name           | data type   | Default Value | description                                                  |
|----------------|-------------|---------------|--------------------------------------------------------------|
| `copyout`      | string list | None          | Deallocates and copies data from GPU to CPU on exit          |
| `delete`       | string list | None          | Deallocates data on GPU on exit                              |
| `detach`       | string list | None          | Detach device pointer from device targets on exit            |
| `extraAccArgs` | string      | None          | String of any extra arguments added to the OpenACC directive |

**Example**

```python
 $:GPU_EXIT_DATA(copyout='[pixels_arr]', delete='[starting_pixels, initial_index]')
 $:GPU_EXIT_DATA(delete='[bc_buffers(1:num_dims, -1:1)]', copyout='[initial_index]')
```

</details>

<details>
 <summary><code>GPU_DECLARE</code> -- <code>(Allocate module variables on GPU or for implicit data region )</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_DECLARE(...)`

**Parameters**

| name             | data type   | Default Value | description                                                                               |
|------------------|-------------|---------------|-------------------------------------------------------------------------------------------|
| `copy`           | string list | None          | Allocates and copies data to GPU on entrance, then deallocated and copies to CPU on exit  |
| `copyin`         | string list | None          | Allocates and copies data to GPU on entrance and then deallocated on exit                 |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU and then deallocated on exit              |
| `copyout`        | string list | None          | Allocates data on GPU on entrance and then deallocates and copies to CPU on exit          |
| `create`         | string list | None          | Allocates data on GPU on entrance and then deallocates on exit                            |
| `present`        | string list | None          | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit |
| `deviceptr`      | string list | None          | Pointer variables that are already allocated on GPU memory                                |
| `link`           | string list | None          | Declare global link, and only allocate when variable used in data clause.                 |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive                              |

**Additional information**

- An implicit data region is created at the start of each procedure and ends after the last executable statement in that procedure.
- Use only create, copyin, device_resident or link clauses for module variables
- GPU_DECLARE exit is the end of the implicit data region
- Link is useful for large global static data objects

**Example**

```python
 $:GPU_DECLARE(create='[x_cb,y_cb,z_cb,x_cc,y_cc,z_cc,dx,dy,dz,dt,m,n,p]')
 $:GPU_DECLARE(create='[x_cb,y_cb,z_cb]', copyin='[x_cc,y_cc,z_cc]', link='[dx,dy,dz,dt,m,n,p]')
```

</details>

<details>
 <summary><code>GPU_UPDATE</code> -- <code>(Updates data from CPU to GPU or GPU to CPU)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_UPDATE(...)`

**Parameters**

| name           | data type   | Default Value | description                                                  |
|----------------|-------------|---------------|--------------------------------------------------------------|
| `host`         | string list | None          | Updates data from GPU to CPU                                 |
| `device`       | string list | None          | Updates data from CPU to GPU                                 |
| `extraAccArgs` | string      | None          | String of any extra arguments added to the OpenACC directive |

**Example**

```python
 $:GPU_UPDATE(host='[arr1, arr2]')
 $:GPU_UPDATE(host='[updated_gpu_val]', device='[updated_cpu_val]')
```

</details>

<details>
 <summary><code>GPU_HOST_DATA</code> -- <code>(Make GPU memory address available on CPU)</code></summary>

**Macro Invocation**

Uses FYPP call directive using `#:call`

```C
 #:call GPU_HOST_DATA(...)
    {code}
 #:endcall GPU_HOST_DATA 
```

**Parameters**

| name           | data type   | Default Value | description                                                      |
|----------------|-------------|---------------|------------------------------------------------------------------|
| `code`         | code        | Required      | Region of code where GPU memory addresses is accessible          |
| `use_device`   | string list | None          | Use GPU memory address of variable instead of CPU memory address |
| `extraAccArgs` | string      | None          | String of any extra arguments added to the OpenACC directive     |

**Parameter Restrictions**

| name   | Restricted range                                 |
|--------|--------------------------------------------------|
| `code` | Do not assign it manually with key-value pairing |

**Example**

```C
 #:call GPU_HOST_DATA(use_device='[addr1, addr2]')
      {code}
      ...
 #:endcall GPU_HOST_DATA
 #:call GPU_HOST_DATA(use_device='[display_arr]')
      {code}
      ...
  #:endcall
```

</details>

------------------------------------------------------------------------------------------

### Synchronization Macros

<details>
 <summary><code>GPU_WAIT</code> -- <code>(Makes CPU wait for async GPU activities)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_WAIT(...)`

**Parameters**

| name           | data type | Default Value | description                                                  |
|----------------|-----------|---------------|--------------------------------------------------------------|
| `extraAccArgs` | string    | None          | String of any extra arguments added to the OpenACC directive |

**Example**

```python
 $:GPU_WAIT()
```

</details>

<details>
 <summary><code>GPU_ATOMIC</code> -- <code>(Do an atomic operation on the GPU)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_ATOMIC(...)`

**Parameters**

| name           | data type | Default Value | description                                                  |
|----------------|-----------|---------------|--------------------------------------------------------------|
| `atomic`       | string    | Required      | Which atomic operation is performed                          |
| `extraAccArgs` | string    | None          | String of any extra arguments added to the OpenACC directive |

**Parameter Restrictions**

| name     | Restricted range                        |
|----------|-----------------------------------------|
| `atomic` | 'read', 'write', 'update', or 'capture' |

**Additional information**

- read atomic is reading in a value
  - Ex: `v=x`
- write atomic is writing a value to a variable
  - Ex:`x=square(tmp)`
- update atomic is updating a variable in-place
  - Ex:`x= x .and. 1`
- Capture is a pair of read/write/update operations with one dependent on the other
  - Ex:

  ```Fortran
      x=x .and. 1
      v=x
  ```

**Example**

```python
 $:GPU_ATOMIC(atomic='update')
 x = square(x)
 $:GPU_ATOMIC(atomic='capture')
 x = square(x)
 v = x
```

</details>

------------------------------------------------------------------------------------------

### Miscellaneous Macros

<details>
 <summary><code>GPU_ROUTINE</code> -- <code>(Compile a procedure for the GPU)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_ROUTINE(...)`

**Parameters**

| name            | data type   | Default Value | description                                                  |
|-----------------|-------------|---------------|--------------------------------------------------------------|
| `function_name` | string      | None          | Name of subroutine/function                                  |
| `parallelism`   | string list | None          | Parallelism granularity to use for this routine              |
| `nohost`        | boolean     | False         | Do not compile procedure code for CPU                        |
| `cray_inline`   | boolean     | False         | Inline procedure on cray compiler                            |
| `extraAccArgs`  | string      | None          | String of any extra arguments added to the OpenACC directive |

**Parameter Restrictions**

| name          | Restricted range                                  |
|---------------|---------------------------------------------------|
| `parallelism` | Valid elements: 'gang', 'worker', 'vector', 'seq' |

**Additional information**

- Function name only needs to be given when cray_inline is True
- Future capability is to parse function header for function name
- Routine parallelism is most commonly ``'[seq]'``

**Example**

```python
 $:GPU_ROUTINE(parallelism='[seq]')
 $:GPU_ROUTINE(function_name='s_matmult', parallelism='[seq]', cray_inline=True)
```

</details>

<details>
 <summary><code>GPU_CACHE</code> -- <code>(Data to be cache in software-managed cache)</code></summary>

**Macro Invocation**

Uses FYPP eval directive using `$:`

`$:GPU_CACHE(...)`

**Parameters**

| name             | data type   | Default Value | description                                                  |
|------------------|-------------|---------------|--------------------------------------------------------------|
| `cache`          | string list | Required      | Data that should to stored in cache                          |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive |

**Example**

```python
 $:GPU_CACHE(cache='[pixels_arr]')
```

</details>

------------------------------------------------------------------------------------------

# Debugging Tools and Tips for GPUs

## Compiler agnostic tools

## OpenMP tools
```bash
OMP_DISPLAY_ENV=true | false | verbose
```
- Prints out the internal control values and environment variables at the beginning of the program if `true` or `verbose`
- `verbose` will also print out vendor-specific internal control values and environment variables

```bash
OMP_TARGET_OFFLOAD = MANDATORY | DISABLED | DEFAULT
```
- Quick way to turn off off-load (`DISABLED`) or make it abort if a GPU isn't found (`MANDATORY`)
- Great first test: does the problem disappear when you drop back to the CPU? 

```bash
OMP_THREAD_LIMIT=<positive_integer>
```
- Sets the maximum number of OpenMP threads to use in a contention group
- Might be useful in checking for issues with contention or race conditions

```bash
OMP_DISPLAY_AFFINITY=TRUE
```
- Will display affinity bindings for each OpenMP thread, containing hostname, process identifier, OS thread identifier, OpenMP thread identifier, and affinity binding.

## Cray Compiler Tools

### Cray General Options

```bash
CRAY_ACC_DEBUG: 0 (off), 1, 2, 3 (very noisy)
```

- Dumps a time-stamped log line (`ACC: ...`) for every allocation, data transfer, kernel launch, wait, etc. Great first stop when "nothing seems to run on the GPU".
- Outputs on STDERR by default. Can be changed by setting `CRAY_ACC_DEBUG_FILE`.
  - Recognizes `stderr`, `stdout`, and `process`.
  - `process` automatically generates a new file based on `pid` (each MPI process will have a different file)
- While this environment variable specifies ACC, it can be used for both OpenACC and OpenMP

```bash
CRAY_ACC_FORCE_EARLY_INIT=1
```

- Force full GPU initialization at program start so you can see start-up hangs immediately
- Default behavior without an environment variable is to defer initialization on first use
- Device initialization includes initializing the GPU vendor’s low-level device runtime library (e.g., libcuda for NVIDIA GPUs) and establishing all necessary software contexts for interacting with the device

### Cray OpenACC Options

```bash
CRAY_ACC_PRESENT_DUMP_SAVE_NAMES=1
```
- Will cause `acc_present_dump()` to output variable names and file locations in addition to variable mappings
- Add `acc_present_dump()` around hotspots to help find problems with data movements
  - Helps more if adding `CRAY_ACC_DEBUG` environment variable

## NVHPC Compiler Options

### NVHPC General Options

```bash
STATIC_RANDOM_SEED=1
```
- Forces the seed returned by `RANDOM_SEED` to be constant, so it generates the same sequence of random numbers 
- Useful for testing issues with randomized data

```bash
NVCOMPILER_TERM=option[,option]
```
- `[no]debug`: Enables/disables just-in-time debugging (debugging invoked on error)
- `[no]trace`: Enables/disables stack traceback on error

### NVHPC OpenACC Options

```bash
NVCOMPILER_ACC_NOTIFY= <bitmask>
```
- Assign the environment variable to a bitmask to print out information to stderr for the following
  - kernel launches: 1
  - data transfers: 2
  - region entry/exit: 4
  - wait operation of synchronizations with the device: 8
  - device memory allocations and deallocations: 16
- 1 (kernels only) is the usual first step.3 (kernels + copies) is great for "why is it so slow?" 

```bash
NVCOMPILER_ACC_TIME=1
```
- Lightweight profiler
- prints a tidy end-of-run table with per-region and per-kernel times and bytes moved
- Do not use with CUDA profiler at the same time

```bash
NVCOMPILER_ACC_DEBUG=1
```
- Spews everything the runtime sees: host/device addresses, mapping events, present-table look-ups, etc.
- Great for "partially present" or "pointer went missing" errors. 
- [Doc for NVCOMPILER_ACC_DEBUG](https://docs.nvidia.com/hpc-sdk/archive/20.9/pdf/hpc209openacc_gs.pdf)
  - Ctrl+F for `NVCOMPILER_ACC_DEBUG`

### NVHPC OpenMP Options

```bash
LIBOMPTARGET_PROFILE=run.json
```
- Emits a Chrome-trace (JSON) timeline you can open in chrome://tracing or Speedscope
- Great lightweight profiler when Nsight is overkill.
- Granularity in µs via `LIBOMPTARGET_PROFILE_GRANULARITY` (default 500).

```bash
LIBOMPTARGET_INFO=<bitmask>
```
- Prints out different types of runtime information
- Human-readable log of data-mapping inserts/updates, kernel launches, copies, waits.
- Perfect first stop for "why is nothing copied?" 
- Flags
  - Print all data arguments upon entering an OpenMP device kernel: 0x01
  - Indicate when a mapped address already exists in the device mapping table: 0x02
  - Dump the contents of the device pointer map at kernel exit: 0x04
  - Indicate when an entry is changed in the device mapping table: 0x08
  - Print OpenMP kernel information from device plugins: 0x10
  - Indicate when data is copied to and from the device: 0x20

```bash
LIBOMPTARGET_DEBUG=1
```
- Developer-level trace (host-side)
- Much noisier than `INFO`
- Only works if the runtime was built with `-DOMPTARGET_DEBUG`.

```bash
LIBOMPTARGET_JIT_OPT_LEVEL=-O{0,1,2,3}
```
- This environment variable can be used to change the optimization pipeline used to optimize the embedded device code as part of the device JIT. 
- The value corresponds to the `-O{0,1,2,3}` command line argument passed to clang.

```bash
LIBOMPTARGET_JIT_SKIP_OPT=1
```
- This environment variable can be used to skip the optimization pipeline during JIT compilation.
- If set, the image will only be passed through the backend.
- The backend is invoked with the `LIBOMPTARGET_JIT_OPT_LEVEL` flag.

## Compiler Documentation

- [Cray & OpenMP Docs](https://cpe.ext.hpe.com/docs/24.11/cce/man7/intro_openmp.7.html#environment-variables)
- [Cray & OpenACC Docs](https://cpe.ext.hpe.com/docs/24.11/cce/man7/intro_openacc.7.html#environment-variables)
- [NVHPC & OpenACC Docs](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html?highlight=NVCOMPILER_#environment-variables)
- [NVHPC & OpenMP Docs](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html?highlight=NVCOMPILER_#id2)
- [LLVM & OpenMP Docs](https://openmp.llvm.org/design/Runtimes.html)
    - NVHPC is built on top of LLVM
- [OpenMP Docs](https://www.openmp.org/spec-html/5.1/openmp.html)
- [OpenACC Docs](https://www.openacc.org/sites/default/files/inline-files/OpenACC.2.7.pdf)
