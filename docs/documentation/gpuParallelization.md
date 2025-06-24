# Parallelization via GPUs

MFC compiles GPU code via OpenACC and in the future OpenMP as well.

In order to swap between OpenACC and OpenMP, ACC and MP directives are implemented via FYPP macros.

[OpenACC Quick start Guide](https://openacc-best-practices-guide.readthedocs.io/en/latest/01-Introduction.html)
[OpenACC API Documentation](https://www.openacc.org/sites/default/files/inline-files/API%20Guide%202.7.pdf)

------------------------------------------------------------------------------------------

## Macro API Documentation

Note: Ordering is not guarennteed or stable, so use key-value pairing when using macros

### Data Type Meanings

- Integer is a number

- Boolean is a pythonic boolean - Valid options: `True` or `False`

- String List is given as a comma separated list surrounding by brackets and inside quotations
  - Ex: `'[hello, world, Fortran]'`

- 2-level string list is given as a comma separated list of string lists surrounding by brackets and inside quotations
  - Ex: `'[[hello, world], [Fortran, MFC]]'` or `[[hello]]`

### Data Flow

- Data on the GPU has a reference counter
- When a variable is referred to being allocated, it means that GPU memory is allocated if it is not already present in GPU memory. If a variable is already present, the reference counter is just incremented.
- When a variable is referred to being dellocated, it means that the reference counter is decremented. If the reference counter is zero, then the data is actually deallocated from GPU memory
- When a variable is referred to being attached, it means that the device pointer attaches to target if it not already attached. If pointer is already attached, then the attachment counter is just incremented
- When a variable is referred to being detached, it means that the attachment counter is decremented. If attachment counter is zero, then actually detached

### Creating new/overwriting existing stubs & proxy configs

<details>
 <summary><code>GPU_PARALLEL_LOOP</code> <code>(Execute the following loop on the GPU in parallel)</code></summary>

#### Macro Invocation

Uses FYPP eval directive using `$:`
> `$:GPU_PARALLEL_LOOP(...)`

#### Parameters

> | name             | data type           | Default Value     | description                                                                                  |
> |------------------|---------------------|-------------------|----------------------------------------------------------------------------------------------|
> | `collapse`       | integer             | None              | Number of loops to combine into 1 loop                                                       |
> | `parallelism`    | string list         | '\[gang,vector\]' | Parallelism granularity to use for this loop                                                 |
> | `default`        | string              | 'present'         | Implicit assumptions compiler should make                                                    |
> | `private`        | string list         | None              | Variables that are private to each iteration/thread                                          |
> | `firstprivate`   | string list         | None              | Initialized variables that are private to each iteration/thread                              |
> | `reduction`      | 2-level string list | None              | Variables unique to each iteration and reduced at the end                                    |
> | `reductionOp`    | string list         | None              | Operator that each list of reduction will reduce with                                        |
> | `copy`           | string list         | None              | Allocates and copies variable to GPU on entrance, then deallocated and copies to CPU on exit |
> | `copyin`         | string list         | None              | Allocates and copies a variable to GPU on entrance and then deallocated on exit              |
> | `copyinReadOnly` | string list         | None              | Allocates and copies a readonly variable to GPU and then deallocated on exit                 |
> | `copyout`        | string list         | None              | Allocates a variable on GPU on entrance and then deallocates and copies to CPU on exit       |
> | `create`         | string list         | None              | Allocates a variable on GPU on entrance and then deallocates on exit                         |
> | `no_create`      | string list         | None              | Use data in CPU memory unless data is already in GPU memory                                  |
> | `present`        | string list         | None              | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit    |
> | `deviceptr`      | string list         | None              | Pointer variables that are already allocated on GPU memory                                   |
> | `attach`         | string list         | None              | Attaches device pointer to device targets on entrance, then detach on exit                   |
> | `extraAccArgs`   | string              | None              | String of any extra arguments added to the OpenACC directive                                 |

#### Parameter Restrictions

| name        | Restricted range                                  |
|-------------|---------------------------------------------------|
| collapse    | Must be greater than 1                            |
| parallelism | Valid elements: 'gang', 'worker', 'vector', 'seq' |
| default     | 'present' or 'none'                               |

#### Additional information

- default present means that the any non-scalar data in assumed to be present on the GPU
- default none means that the compiler should not implicitly determine the data attributes for any variable
- reduction and reductionOp must match in length
- With `reduction='[[sum1, sum2], [largest]]'` and `reductionOp='[+, max]'`, `sum1` and `sum2` will be the sum of sum1/sum2 in each loop iteration, and `largest` will the maximum value of `largest` all the loop iterations
- A reduction implies a copy, so it does not need to be added for both

#### Example

> ```python
>  $:GPU_PARALLEL_LOOP(collapse=3, private='[tmp, r]', reduction='[[vol, avg], [max_val]]', reductionOp='[+, MAX]')
>  $:GPU_PARALLEL_LOOP(collapse=2, private='[sum_holder]', copyin='[starting_sum]', copyout='[eigenval,C]')
> ```

</details>

------------------------------------------------------------------------------------------

### Listing existing stubs & proxy configs as YAML string

<details>
 <summary><code>GPU_ROUTINE</code> <code>(Compile a procedure for the GPU)</code></summary>

#### Macro Invocation

Uses FYPP eval directive using `$:`
> `$:GPU_ROUTINE(...)`


#### Parameters

| name            | data type   | Default Value | description                                                  |
|-----------------|-------------|---------------|--------------------------------------------------------------|
| `function_name` | string      | None          | Name of subroutine/function                                  |
| `parallelism`   | string list | None          | Parallelism granularity to use for this routine              |
| `nohost`        | boolean     | False         | Do not compile procedure code for CPU                        |
| `cray_inline`   | boolean     | False         | Inline procedure on cray compiler                            |
| `extraAccArgs`  | string      | None          | String of any extra arguments added to the OpenACC directive |

#### Parameter Restrictions

| name        | Restricted range                                  |
|-------------|---------------------------------------------------|
| parallelism | Valid elements: 'gang', 'worker', 'vector', 'seq' |

#### Additional information

- Function name only needs to be given when cray_inline is True
- Future capability is to parse function header for function name
- Routine parallelism is most commonly `'[seq]'`

#### Example

> ```python
>  $:GPU_ROUTINE(parallelism='[seq]')
>  $:GPU_ROUTINE(function_name='s_matmult', parallelism='[seq]', cray_inline=True)
> ```

</details>

<details>
 <summary><code>GPU_DECLARE</code> <code>(Allocate module variables on GPU or for implicit data region )</code></summary>

#### Macro Invocation

Uses FYPP eval directive using `$:`
> `$:GPU_DECLARE(...)`

#### Parameters

| name             | data type   | Default Value | description                                                                                  |
|------------------|-------------|---------------|----------------------------------------------------------------------------------------------|
| `copy`           | string list | None          | Allocates and copies variable to GPU on entrance, then deallocated and copies to CPU on exit |
| `copyin`         | string list | None          | Allocates and copies a variable to GPU on entrance and then deallocated on exit              |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU and then deallocated on exit                 |
| `copyout`        | string list | None          | Allocates a variable on GPU on entrance and then deallocates and copies to CPU on exit       |
| `create`         | string list | None          | Allocates a variable on GPU on entrance and then deallocates on exit                         |
| `present`        | string list | None          | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit    |
| `deviceptr`      | string list | None          | Pointer variables that are already allocated on GPU memory                                   |
| `link`           | string list | None          | Declare global link, and only allocate when variable used in data clause.                    |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive                                 |

#### Additional information

- An implicit data region is created at the start of each procedure
and ends after the last executable statement in that procedure.
- Use only create, copyin, device_resident or link clauses for module variables
- GPU_DECLARE exit is the end of the implicit data region
- Link is useful for large global static data objects

#### Example

> ```python
>  $:GPU_DECLARE(create='[x_cb,y_cb,z_cb,x_cc,y_cc,z_cc,dx,dy,dz,dt,m,n,p]')
>  $:GPU_DECLARE(create='[x_cb,y_cb,z_cb]', copyin='[x_cc,y_cc,z_cc]', link='[dx,dy,dz,dt,m,n,p]')
> ```

</details>

<details>
  <summary><code>GPU_LOOP</code> <code>(Execute loop on GPU)</code></summary>

#### Macro Invocation

Uses FYPP eval directive using `$:`
> `$:GPU_LOOP(...)`

#### Parameters

| name              | data type           | Default Value | description                                                                                      |
|-------------------|---------------------|---------------|--------------------------------------------------------------------------------------------------|
| `collapse`        | integer             | None          | Number of loops to combine into 1 loop                                                           |
| `parallelism`     | string list         | None          | Parallelism granularity to use for this loop                                                     |
| `data_dependency` | string              | None          | 'independent'-> assert loop iterations are independent, 'auto->let compiler analyze dependencies |
| `private`         | string list         | None          | Variables that are private to each iteration/thread                                              |
| `reduction`       | 2-level string list | None          | Variables unique to each iteration and reduced at the end                                        |
| `reductionOp`     | string list         | None          | Operator that each list of reduction will reduce with                                            |
| `extraAccArgs`    | string              | None          | String of any extra arguments added to the OpenACC directive                                     |

#### Parameter Restrictions

| name            | Restricted range                                  |
|-----------------|---------------------------------------------------|
| collapse        | Must be greater than 1                            |
| parallelism     | Valid elements: 'gang', 'worker', 'vector', 'seq' |
| data_dependency | 'auto' or 'independent'                           |

#### Additional information

- Loop parallelism is most commonly `'[seq]'`
- reduction and reductionOp must match in length
- With `reduction='[[sum1, sum2], [largest]]'` and `reductionOp='[+, max]'`, `sum1` and `sum2` will be the sum of sum1/sum2 in each loop iteration, and `largest` will the maximum value of `largest` all the loop iterations

#### Example

> ```python
>  $:GPU_PARALLEL_LOOP(parallelism='[seq]')
>  $:GPU_PARALLEL_LOOP(collapse=3, parallelism='[seq]',private='[tmp, r]')
> ```

</details>

<details>
 <summary><code>GPU_DATA</code> <code>(Allocate module variables on GPU or for implicit data region )</code></summary>

#### Macro Invocation

Uses FYPP call directive using `#:call`
> ```C
>
> #:call GPU_DATA(...)
>    {code}
> #:endcall GPU_DATA 
>```
> 

#### Parameters

| name             | data type   | Default Value | description                                                                                  |
|------------------|-------------|---------------|----------------------------------------------------------------------------------------------|
| `code`           | code        | Required      | Region of code where defined data is accessible                                              |
| `copy`           | string list | None          | Allocates and copies variable to GPU on entrance, then deallocated and copies to CPU on exit |
| `copyin`         | string list | None          | Allocates and copies a variable to GPU on entrance and then deallocated on exit              |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU and then deallocated on exit                 |
| `copyout`        | string list | None          | Allocates a variable on GPU on entrance and then deallocates and copies to CPU on exit       |
| `create`         | string list | None          | Allocates a variable on GPU on entrance and then deallocates on exit                         |
| `no_create`      | string list | None          | Use data in CPU memory unless data is already in GPU memory                                  |
| `present`        | string list | None          | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit    |
| `deviceptr`      | string list | None          | Pointer variables that are already allocated on GPU memory                                   |
| `attach`         | string list | None          | Attaches device pointer to device targets on entrance, then detach on exit                   |
| `default`        | string      | None          | Implicit assumptions compiler should make                                                    |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive                                 |

#### Parameter Restrictions

| name            | Restricted range                                  |
|-----------------|---------------------------------------------------|
| code        | Do not assign it manually with key-value pairing                            |

#### Additional information

#### Example

> ```C
>  #:call GPU_DATA(copy='[pixel_arr]', copyin='[starting_pixels, inital_index]',attach='[p_real, p_cmplx, p_fltr_cmplx]')
>       {code}
>       ...
>  #:endcall
>  #:call GPU_DATA(create='[pixel_arr]', copyin='[inital_index]')
>       {code}
>       ...
>   #:endcall
> ```

</details>

<details>
 <summary><code>GPU_DATA</code> <code>(Allocate module variables on GPU or for implicit data region )</code></summary>

#### Macro Invocation

Uses FYPP call directive using `#:call`
> ```C
>
> #:call GPU_DATA(...)
>    {code}
> #:endcall GPU_DATA 
>```
> 

#### Parameters

| name             | data type   | Default Value | description                                                                                  |
|------------------|-------------|---------------|----------------------------------------------------------------------------------------------|
| `code`           | code        | Required      | Region of code where defined data is accessible                                              |
| `copy`           | string list | None          | Allocates and copies variable to GPU on entrance, then deallocated and copies to CPU on exit |
| `copyin`         | string list | None          | Allocates and copies a variable to GPU on entrance and then deallocated on exit              |
| `copyinReadOnly` | string list | None          | Allocates and copies a readonly variable to GPU and then deallocated on exit                 |
| `copyout`        | string list | None          | Allocates a variable on GPU on entrance and then deallocates and copies to CPU on exit       |
| `create`         | string list | None          | Allocates a variable on GPU on entrance and then deallocates on exit                         |
| `no_create`      | string list | None          | Use data in CPU memory unless data is already in GPU memory                                  |
| `present`        | string list | None          | Data that must be present in GPU memory. Increment counter on entrance, decrement on exit    |
| `deviceptr`      | string list | None          | Pointer variables that are already allocated on GPU memory                                   |
| `attach`         | string list | None          | Attaches device pointer to device targets on entrance, then detach on exit                   |
| `default`        | string      | None          | Implicit assumptions compiler should make                                                    |
| `extraAccArgs`   | string      | None          | String of any extra arguments added to the OpenACC directive                                 |

#### Parameter Restrictions

| name            | Restricted range                                  |
|-----------------|---------------------------------------------------|
| code        | Do not assign it manually with key-value pairing                            |

#### Additional information

#### Example

> ```C
>  #:call GPU_DATA(copy='[pixel_arr]', copyin='[starting_pixels, inital_index]',attach='[p_real, p_cmplx, p_fltr_cmplx]')
>       {code}
>       ...
>  #:endcall
>  #:call GPU_DATA(create='[pixel_arr]', copyin='[inital_index]')
>       {code}
>       ...
>   #:endcall
> ```

</details>


------------------------------------------------------------------------------------------
