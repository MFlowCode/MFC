# Case Files

Example Python case files, also referred to as *input files*, can be found in the [examples/](https://github.com/MFlowCode/MFC/tree/master/examples) directory. They print a Python dictionary containing input parameters for MFC. Their contents, and a guide to filling them out, are documented in the user manual. A commented, tutorial script
can also be found in [examples/3d_sphbubcollapse/](https://github.com/MFlowCode/MFC/blob/master/examples/3D_sphbubcollapse/case.py).

## Basic Skeleton

The skeleton for an input file may look like the following:

```python
#!/usr/bin/env python3

import json

# Configuring case dictionary
print(json.dumps({
  # Insert case parameters here
  ...
}))
```

Thus, you can run your case file with Python to view the computed case dictionary that will be processed by MFC when you run:

```shell
python3 my_case_file.py
```

This is particularly useful when computations are done in Python to generate the case.

## (Optional) Accepting command line arguments

Input files can accept command line arguments, forwarded by `mfc.sh run`.
Consider this example from the `scaling` case:

```python
import json, argparse

parser = argparse.ArgumentParser(
    prog="scaling",
    description="Weak- and strong-scaling benchmark case.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("-s", "--scaling", type=str, metavar="SCALING", choices=["weak", "strong"],
                    help="Whether weak- or strong-scaling is being exercised.")

# Your parsed arguments are here
args = parser.parse_args()
```

The `--mfc` argument is a JSON string representing `mfc.sh run`'s internal
state, passed in when MFC runs your input file.
It contains all the runtime information you might want from the build/run system.
You can add as many additional arguments and options as you may need.

To run such a case, use the following format:

```shell
./mfc.sh run <path/to/case.py> <mfc.sh run arguments> -- <case arguments>
```

For example, to run the `scaling` case in "weak-scaling" mode:

```shell
./mfc.sh run examples/scaling/case.py -t pre_process -j 8 -- --scaling weak
```

## Parameters

There are multiple sets of parameters that must be specified in the python input file:
1. [Runtime Parameters](#1-runtime)
2. [Computational Domain Parameters](#2-computational-domain)
3. [Patch Parameters](#3-patches)
4. [Immersed Boundary Patches](#4-immersed-boundary-patches)
5. [Fluid Material's Parameters](#5-fluid-materials)
6. [Simulation Algorithm Parameters](#6-simulation-algorithm)
7. [Formatted Database and Structure Parameters](#7-formatted-output)
8. [(Optional) Acoustic Source Parameters](#8-acoustic-source)
9. [(Optional) Ensemble-Averaged Bubble Model Parameters](#9-ensemble-averaged-bubble-model)
10. [(Optional) Velocity Field Setup Parameters](#10-velocity-field-setup)
11. [(Optional) Phase Change Parameters](#11-Phase-Change-Model)
12. [(Optional) Artificial Mach Number Parameters](#12-artificial-Mach-number)

Items 8, 9, 10, 11 and 12 are optional sets of parameters that activate the acoustic source model, ensemble-averaged bubble model, initial velocity field setup, phase change, artificial Mach number respectively.
Definition of the parameters is described in the following subsections.

### 1. Runtime

| Parameter        | Type           | Description                               |
| ---:             |    :----:      |          :---                             |
| `run_time_info`  | Logical        | Output run-time information               |
| `rdma_mpi`       | Logical        | (GPUs) Enable RDMA for MPI communication. |

- `run_time_info` generates a text file that includes run-time information including the CFL number(s) at each time-step.
- `rdma_mpi` optimizes data transfers between GPUs using Remote Direct Memory Access (RDMA).
The underlying MPI implementation and communication infrastructure must support this
feature, detecting GPU pointers and performing RDMA accordingly.

### 2. Computational Domain

| Parameter                | Type    | Description                      |
| ---:                     | :----:  |          :---                    |
| `x[y,z]_domain%%beg[end]`| Real    | Beginning [ending] of the $x$[y,z]-direction domain    |
| `stretch_x[y,z]`         | Logical | Stretching of the mesh in the $x$[y,z]-direction |
| `a_x[y,z]`               | Real    | Rate at which the grid is stretched in the $x$[y,z]-direction |
| `x[y,z]_a`               | Real    | Beginning of the stretching in the negative $x$[y,z]-direction |
| `x[y,z]_b`               | Real    | Beginning of the stretching in the positive $x$[y,z]-direction |
| `loops_x[y,z]`           | Integer | Number of times to recursively apply grid stretching |
| `cyl_coord`              | Logical | Cylindrical coordinates (2D: Axisymmetric, 3D: Cylindrical) |
| `m`                      | Integer | Number of grid cells in the $x$-coordinate direction |
| `n`                      | Integer | Number of grid cells in the $y$-coordinate direction |
| `p`                      | Integer | Number of grid cells in the $z$-coordinate direction |

The parameters define the boundaries of the spatial and temporal domains, and their discretization that are used in simulation.

- `[x,y,z]_domain%[beg,end]` define the spatial domain in $x$, $y$, and $z$ Cartesian coordinates:

$$ x \in \left[ x \\_ domain \\% beg, x \\_ domain \\% end \right], y \in \left[ y \\_ domain \\% beg, y \\_ domain \\% end \right], z \in \left[ z \\_ domain \\% beg, z \\_ domain \\% end \right] $$

- $m$, $n$, and $p$ define the number of finite volume cells that uniformly discretize the domain along the $x$, $y$, and $z$ axes, respectively.
Note that the actual number of cells in each coordinate axis is given as $[m,n,p]+1$.
For example, $m=n=p=499$ discretizes the domain into $500^3$ cells.
When the simulation is 2D/axi-symmetric or 1D, it requires that $p=0$ or $p=n=0$, respectively.

- `stretch_[x,y,z]` activates grid stretching in the $[x,y,z]$ directions.
The grid is gradually stretched such that the domain boundaries are pushed away from the origin along a specified axis.

- `a_[x,y,z]`, `[x,y,z]_a`, and `[x,y,z]_b` are parameters that define the grid stretching function. When grid stretching along the $x$ axis is considered, the stretching function is given as:

$$ x_{cb,stretch} = x_{cb} + \frac{x_{cb}}{a_x} \Bigg[ \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_a)}{L} \right) \right] + \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_b)}{L} \right) \right] -2 \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_b-x_a)}{2L} \right) \right]  \Bigg] $$

where `x_cb` and `x_[cb,stretch]` are the coordinates of a cell boundary at the original and stretched domains, respectively.
`L` is the domain length along the `x` axis: `L`=`x_domain%%end`-`x_domain%%beg`.
Crudely speaking, `x_a` and `x_b` define the coordinates at which the grid begins to get stretched in the negative and positive directions along the $x$ axis, respectively.
$a_x$ defines the smoothness of the stretching.
Stretching along the $y$ and $z$ axes follows the same logistics.
Optimal choice of the parameters for grid stretching is case-dependent and left to the user.
`loops_x[y,z]` defines the number of times the grid stretching function is applied and has a default value of one.

- `cyl_coord` activates cylindrical coordinates.
The domain is defined in $x$-$y$-$z$ cylindrical coordinates, instead of Cartesian coordinates.
Domain discretization is accordingly conducted along the axes of cylindrical coordinates.
When $p=0$, the domain is defined on $x$-$y$ axi-symmetric coordinates.
In both Coordinates, mesh stretching can be defined along the $x$- and $y$-axes.
MPI topology is automatically optimized to maximize the parallel efficiency for given choice of coordinate systems.

### 3. Patches

| Parameter            | Type    | Analytical Definition | Description                                                  |
| ---:                 | :----:  | :----:                | :---                                                         |
| `num_patches`        | Integer | Not Supported         | Number of initial condition geometric patches.               |
| `num_fluids`	       | Integer | Not Supported         | Number of fluids/components present in the flow.             |
| `geometry` *         | Integer | Not Supported         | Geometry configuration of the patch.                         |
| `alter_patch(i)` *   | Logical | Not Supported         | Alter the $i$-th patch.                                      |
| `x[y,z]_centroid` *  | Real    | Not Supported         | Centroid of the applied geometry in the $[x,y,z]$-direction. |
| `length_x[y,z]` *    | Real    | Not Supported         | Length, if applicable, in the $[x,y,z]$-direction.           |
| `radius` *           | Real    | Not Supported         | Radius, if applicable, of the applied geometry.              |
| `smoothen` *         | Logical | Not Supported         | Smoothen the applied patch.                                  |
| `smooth_patch_id` *  | Integer | Not Supported         | A patch with which the applied patch is smoothened.          |
| `smooth_coeff` *     | Real    | Not Supported         | Smoothen coefficient.                                        |
| `alpha(i)` *         | Real    | Supported             | Volume fraction of fluid $i$.                                |
| `alpha_rho(i)` *     | Real    | Supported             | Partial density of fluid $i$.                                |
| `pres` *             | Real    | Supported             | Pressure.                                                    |
| `vel(i)` *           | Real    | Supported             | Velocity in direction $i$.                                   |
| `tau_e(i)` *         | Real    | Supported             | Elastic stresses.                                            |
| `hcid` *             | Integer | N/A                   | Hard coded patch id                                          |
| `cf_val` *           | Real    | Supported             | Surface tension color function value                         |
| `model_filepath`     | String  | Not Supported         | Path to an STL or OBJ file (not all OBJs are supported).     |
| `model_scale(i)`     | Real    | Not Supported         | Model's (applied) scaling factor for component $i$.          |
| `model_rotate(i)`    | Real    | Not Supported         | Model's (applied) angle of rotation about axis $i$.          |
| `model_translate(i)` | Real    | Not Supported         | Model's $i$-th component of (applied) translation.           |
| `model_spc`          | Integer | Not Supported         | Number of samples per cell when discretizing the model into the grid. |
| `model_threshold`    | Real    | Not Supported         | Ray fraction inside the model patch above which the fraction is set to one.|

*: These parameters should be prepended with `patch_icpp(j)%` where $j$ is the patch index.

The Table lists the patch parameters.
The parameters define the geometries and physical parameters of fluid components (patch) in the domain at initial condition.
Note that the domain must be fully filled with patche(s).
The code outputs error messages when an empty region is left in the domain.

- `tau_e(i)` is the `i`-th component of the elastic stress tensor, ordered as `tau_xx`, `tau_xy`, `tau_yy`, `tau_xz`, `tau_yz`, and `tau_zz`. 1D simulation requires `tau(1)`, 2D `tau(1:3)`, and 3D `tau(1:6)`.

#### Analytical Definition of Primitive Variables

Some parameters, as described above, can be defined by analytical functions in the input file. For example, one can define the following patch:

```shell
'patch_icpp(2)%geometry'    : 1,
'patch_icpp(2)%x_centroid'  : 0.25,
'patch_icpp(2)%length_x'    : 9.5,
'patch_icpp(2)%vel(1)'      : 0.,
'patch_icpp(2)%pres'        : 1.,
'patch_icpp(2)%alpha_rho(1)': '1 + 0.1*sin(20*x*pi)',
'patch_icpp(2)%alpha(1)'    : 1.,
```

where `alpha_rho` is defined with the `1 + 0.1*sin(20*x*pi)` function.

The expressions are interpreted as Fortran code, in the execution context and scope of the function that defines the patch.
Additionally, the following variables are made available as shorthand:

| Shorthand | Expands To               | Shorthand | Expands To                | Shorthand | Expands To               |
| --------- | ------------------------ | --------- | ------------------------- | --------- | ------------------------ |
| `x`       | `x_cc(i)`                | `lx`      | The patch's `length_x`    | `xc`      | The patch's `x_centroid` |
| `y`       | `y_cc(j)`                | `ly`      | The patch's `length_y`    | `yc`      | The patch's `y_centroid` |
| `z`       | `z_cc(k)`                | `lz`      | The patch's `length_z`    | `zc`      | The patch's `z_centroid` |
| `eps`     | The patch's `epsilon`    | `beta`    | The patch's `beta`        | `radii`   | The patch's `radii`      |
| `tau_e`   | The patch's `tau_e`      | `r`       | The patch's `radius`      | `pi` and `e` | $\pi$ and $e$         |

where $(i,j,k)$ are the grid-indices of the current cell in each coordinate direction.

In the example above, the following code is generated:

```f90
if (patch_id == 2) then
    q_prim_vf(contxb)%sf(i, 0, 0) = 1 + 0.1*sin(20*x_cc(i)*3.141592653589793)
end if
```

#### Hard Coded Patches

Some patch configurations are not adequately handled with the above analytic variable definitions.
In this case, a hard coded patch can be used.
Hard coded patches can be added by adding additional hard coded patch identifiers to `src/pre_process/include/1[2,3]dHardcodedIC.fpp`.
When using a hard coded patch, the `patch_icpp(patch_id)%%hcid` must be set to the hard-coded patch id.
For example, to add a 2D Hardcoded patch with an id of 200, one would add the following to `src/pre_process/include/2dHardcodedIC.fpp`

```f90
    case(200)
        ! Primitive variables assignment
```

and use `patch_icpp(i)%%hcid = 200` in the input file.
Additional variables can be declared in `Hardcoded1[2,3]DVariables` and used in `hardcoded1[2,3]D`.
As a convention, any hard coded patches that are part of the MFC master branch should be identified as 1[2,3]xx where the first digit indicates the number of dimensions.

The code provides three pre-built patches for dimensional extrusion of initial conditions:

- `case(170)`: Load 1D profile from data files
- `case(270)`: Extrude 1D data to 2D domain
- `case(370)`: Extrude 2D data to 3D domain

Setup: Only requires specifying `init_dir` and filename pattern via `zeros_default`. Grid dimensions are automatically detected from the data files.
Implementation: All variables and file handling are managed in `src/pre_process/include/ExtrusionHardcodedIC.fpp` with no manual grid configuration needed.
Usage: Ideal for initializing simulations from lower-dimensional solutions, enabling users to add perturbations or modifications to the base extruded fields for flow instability studies.

#### Parameter Descriptions

- `num_patches` defines the total number of patches defined in the domain.
The number has to be a positive integer.

- `num_fluids` defines the total number of fluids defined in each of the patches.
The number has to be a positive integer.

- `patch_icpp(j)%%geometry` defines the type of geometry of $j$-th patch by using an integer from 1 to 13.
Definition of the patch type for each integer is listed in table [Patch Types](#patch-types).

- `[x,y,z]_centroid`, `length_[x,y,z]`, and/or `radius` are used to uniquely define the geometry of the patch with given type.
Requisite combinations of the parameters for each type can be found in is listed in table [Patch types](#patch-types).

- `patch_icpp(j)%%alter_patch(i)` activates alternation of `patch(i)` with `patch(j)`.
For instance, in a 2D simulation, when a cylindrical `patch(2)` is immersed in a rectangular `patch(1)`:
  - `patch_icpp(1)%%geometry = 3`
  - `patch_icpp(2)%%geometry = 2`
  - ``patch_icpp(2)%%alter_patch(1) = 'T'``

- `smoothen` activates smoothening of the boundary of the patch that alters the existing patch.
When smoothening occurs, fluids of the two patches are mixed in the region of the boundary.
For instance, in the aforementioned case of the cylindrical patch immersed in the rectangular patch, smoothening occurs when ``patch_icpp(2)smoothen = 'T'``.
`smooth_coeff` controls the thickness of the region of smoothening (sharpness of the mixture region).
The default value of `smooth_coeff` is unity. The region of smoothening is thickened with decreasing the value.
Optimal choice of the value of `smooth_coeff` is case-dependent and left to the user.

- `patch_icpp(j)alpha(i)`, `patch_icpp(j)alpha_rho(i)`, `patch_icpp(j)pres`, and `patch_icpp(j)vel(i)` define for $j$-th patch the void fraction of `fluid(i)`, partial density of `fluid(i)`, the pressure, and the velocity in the $i$-th coordinate direction.
These physical parameters must be consistent with fluid material's parameters defined in the next subsection.

- `model_filepath` defines the root directory of the STL or OBJ model file.

- `model_scale`, `model_rotate` and `model_translate` define how the model should be transformed to domain-space by first scaling by `model_scale`, then rotating about the Z, X, and Y axes (using `model_rotate`), and finally translating by `model_translate`.

- `model_spc` and `model_threshold` are ray-tracing parameters. `model_spc` defines the number of rays per cell to render the model. `model_threshold` defines the ray-tracing threshold at which the cell is marked as the model.

#### Elliptic Smoothing

Initial conditions in which not all patches support the `patch_icpp(j)%smoothen` parameter can still be smoothed by applying iterations of the heat equation to the initial condition.
This is enabled by adding `'elliptic_smoothing': "T",` and `'elliptic_smoothing_iters': N,` to the case dictionary, where `N` is the number of smoothing iterations to apply.

### 4. Immersed Boundary Patches

| Parameter            | Type    | Description |
| ---:                 | :----:  | :---                |
| `geometry`             | Integer | Geometry configuration of the patch.|
| `x[y,z]_centroid`      | Real    | Centroid of the applied geometry in the [x,y,z]-direction. |
| `length_x[y,z]`        | Real    | Length, if applicable, in the [x,y,z]-direction. |
| `radius`               | Real    | Radius, if applicable, of the applied geometry. |
| `theta`                | Real    | Angle of attach applied to airfoil IB patches |
| `c`                    | Real    | NACA airfoil parameters (see below) |
| `t`                    | Real    | NACA airfoil parameters (see below) |
| `m`                    | Real    | NACA airfoil parameters (see below) |
| `p`                    | Real    | NACA airfoil parameters (see below) |
| `slip`                 | Logical | Apply a slip boundary |
| `model_filepath`      | String  | Path to an STL or OBJ file (not all OBJs are supported).     |
| `model_scale(i)`      | Real    | Model's (applied) scaling factor for component $i$.          |
| `model_rotate(i)`     | Real    | Model's (applied) angle of rotation about axis $i$.          |
| `model_translate(i)`  | Real    | Model's $i$-th component of (applied) translation.           |
| `model_spc`           | Integer | Number of samples per cell when discretizing the model into the grid. |
| `model_threshold`     | Real    | Ray fraction inside the model patch above which the fraction is set to one.|

These parameters should be prepended with `patch_ib(j)%` where $j$ is the patch index.

#### Parameter Descriptions

- `geometry` defines the type of geometry of an immersed boundary patch with an integer number.
Definitions for currently implemented immersed boundary patch types are listed in table [Immersed Boundary Patch Type](#immersed-boundary-patch-types).

- `x[y,z]_centroid` is the centroid location of the patch in the x[y,z]-direction

- `length_x[y,z]` is the length of the patch in the x[y,z]-direction.

- `radius` is the radius to be used for circular patches.

- `theta` allows for the angle of attach of airfoil patches to be changed.

- `c`, `t`, `p`, and `m` specify the parameters for a NACA airfoil.
`m` is the maximum camber, `p` is the location of maximum camber, `c` is the coord length, and `t` is the thickness.
Additional details on this specification can be found in [The Naca Airfoil Series](https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf)

- `slip` applies a slip boundary to the surface of the patch if true and a no-slip boundary condition to the surface if false.

- Please see [Patch Parameters](#3-patches) for the descriptions of `model_filepath`, `model_scale`, `model_rotate`, `model_translate`, `model_spc`, and `model_threshold`.

### 5. Fluid Material’s

| Parameter | Type   | Description                                    |
| ---:      | :----: |          :---                                  |
| `gamma`   | Real   | Stiffened-gas parameter $\Gamma$ of fluid.     |
| `pi_inf`  | Real   | Stiffened-gas parameter $\Pi_\infty$ of fluid. |
| `Re(1)` * | Real   | Shear viscosity of fluid.                      |
| `Re(2)` * | Real   | Volume viscosity of fluid.                     |
| `cv`   ** | Real   | Sffened-gas parameter $c_v$ of fluid.          |
| `qv`   ** | Real   | Stiffened-gas parameter $q$ of fluid.          |
| `qvp`  ** | Real   | Stiffened-gas parameter $q'$ of fluid.         |
| `sigma`   | Real   | Surface tension coefficient                    |
| `G`       | Real   | Shear modulus of solid.                        |

Fluid material's parameters. All parameters except for sigma should be prepended with `fluid_pp(i)` where $i$ is the fluid index.

*: Parameters that work only with `model_eqns = 2`.

**: Parameters that work only with `model_eqns = 3`.

The table lists the fluid material's parameters.
The parameters define material's property of compressible fluids that are used in simulation.

- `fluid_pp(i)%%gamma` and `fluid_pp(i)%%pi_inf` define $\Gamma$ and $\Pi$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

- `fluid_pp(i)%%Re(1)` and `fluid_pp(i)%%Re(2)` define the shear and volume viscosities of $i$-th fluid, respectively.

When these parameters are undefined, fluids are treated as inviscid.
Details of implementation of viscosity in MFC can be found in [Coralic (2015)](references.md).

- `fluid_pp(i)%%cv`, `fluid_pp(i)%%qv`, and `fluid_pp(i)%%qvp` define $c_v$, $q$, and $q'$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

- `fluid_pp(i)%%G` is required for `hypoelasticity`.

### 6. Simulation Algorithm

| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `bc_[x,y,z]%%beg[end]` | Integer | Beginning [ending] boundary condition in the $[x,y,z]$-direction (negative integer, see table [Boundary Conditions](#boundary-conditions)) |
| `bc_[x,y,z]%%vb[1,2,3]`‡| Real   | Velocity in the (x,1), (y, 2), (z,3) direction applied to `bc_[x,y,z]%%beg` |
| `bc_[x,y,z]%%ve[1,2,3]`‡| Real   | Velocity in the (x,1), (y, 2), (z,3) direction applied to `bc_[x,y,z]%%end` |
| `model_eqns`           | Integer | Multicomponent model: [1] $\Gamma/\Pi_\infty$; [2] 5-equation; [3] 6-equation; [4] 4-equation |
| `alt_soundspeed` *     | Logical | Alternate sound speed and $K \nabla \cdot u$ for 5-equation model |
| `adv_n`   	         | Logical | Solving directly for the number density (in the method of classes) and compute void fraction from the number density |
| `mpp_lim`	             | Logical | Mixture physical parameters limits |
| `mixture_err`          | Logical | Mixture properties correction |
| `time_stepper`         | Integer | Runge--Kutta order [1-3] |
| `adap_dt`              | Logical | Strang splitting scheme with adaptive time stepping |
| `recon_type`           | Integer | Reconstruction Type: [1] WENO; [2] MUSCL |
| `adap_dt_tol`          | Real    | Tolerance for adaptive time stepping in Strang splitting scheme|
| `adap_dt_max_iters`    | Integer | Max iteration for adaptive time stepping in Strang splitting scheme |
| `weno_order`	         | Integer | WENO order [1,3,5] |
| `weno_eps`	         | Real    | WENO perturbation (avoid division by zero) |
| `mapped_weno`	         | Logical | WENO-M (WENO with mapping of nonlinear weights) |
| `wenoz`	             | Logical | WENO-Z |
| `wenoz_q`              | Real    | WENO-Z power parameter q (only for WENO7) |
| `teno`                 | Logical | TENO (Targeted ENO) |
| `teno_CT`              | Real    | TENO threshold for smoothness detection |
| `null_weights`         | Logical | Null WENO weights at boundaries |
| `mp_weno`              | Logical | Monotonicity preserving WENO |
| `muscl_order`          | Integer | MUSCL order [1,2] |
| `muscl_lim`            | Integer | MUSCL Slope Limiter: [1] minmod; [2] monotonized central; [3] Van Albada; [4] Van Leer; [5] SUPERBEE |
| `int_comp`             | Logical | THINC Interface Compression |
| `ic_eps`               | Real | Interface compression threshold (default: 1e-4) |
| `ic_beta`              | Real | Interface compression sharpness parameter (default: 1.6) |
| `riemann_solver`       | Integer | Riemann solver algorithm: [1] HLL*; [2] HLLC; [3] Exact*; [4] HLLD	(only for MHD) |
| `low_Mach`             | Integer | Low Mach number correction for HLLC Riemann solver: [0] None; [1] Pressure (Chen et al. 2022); [2] Velocity (Thornber et al. 2008)	 |
| `avg_state`	         | Integer | Averaged state evaluation method: [1] Roe average*; [2] Arithmetic mean  |
| `wave_speeds`          | Integer | Wave-speed estimation: [1] Direct (Batten et al. 1997); [2] Pressure-velocity* (Toro 1999)	 |
| `weno_Re_flux`         | Logical | Compute velocity gradient using scalar divergence theorem	 |
| `weno_avg`          	 | Logical | Arithmetic mean of left and right, WENO-reconstructed, cell-boundary values |
| `dt`                   | Real    | Time step size |
| `t_step_start`         | Integer | Simulation starting time step |
| `t_step_stop`          | Integer | Simulation stopping time step |
| `t_step_save`          | Integer | Frequency to output data |
| `t_step_print`         | Integer | Frequency to print the current step number to standard output (default 1) |
| `cfl_adap_dt`          | Logical | CFL based adaptive time-stepping |
| `cfl_const_dt`         | Logical | CFL based non-adaptive time-stepping |
| `cfl_target`           | Real    | Specified CFL value |
| `n_start`              | Integer | Save file from which to start simulation |
| `t_save`               | Real    | Time duration between data output |
| `t_stop`               | Real    | Simulation stop time |
| `surface_tension`      | Logical | Activate surface tension |
| `viscous`              | Logical | Activate viscosity |
| `hypoelasticity`       | Logical | Activate hypoelasticity* |
| `igr`                  | Logical | Enable solution via information geometric regularization (IGR) [Cao (2024)](references.md) |
| `igr_order`            | Integer | Order of reconstruction for IGR [3,5] |
| `alf_factor`           | Real    | Alpha factor for IGR entropic pressure (default 10) |
| `igr_pres_lim`         | Logical | Limit IGR pressure to avoid negative values (default F) |
| `igr_iter_solver`      | Integer | Solution method for IGR elliptic solve [1] Jacobi [2] Gauss-Seidel |
| `num_igr_iters`        | Integer | Number of iterations for for the IGR elliptic solve (default 2) |
| `num_igr_warm_start_iters` | Integer | Number of iterations for the IGR elliptic solve at the first time step (default 50) |

- \* Options that work only with `model_eqns = 2`.
- † Options that work only with ``cyl_coord = 'F'``.
- ‡ Options that work only with `bc_[x,y,z]%[beg,end] = -15` and/or `bc_[x,y,z]%[beg,end] = -16`.

The table lists simulation algorithm parameters.
The parameters are used to specify options in algorithms that are used to integrate the governing equations of the multi-component flow based on the initial condition.
Models and assumptions that are used to formulate and discritize the governing equations are described in [Bryngelson et al. (2019)](references.md).
Details of the simulation algorithms and implementation of the WENO scheme can be found in [Coralic (2015)](references.md).

- `bc_[x,y,z]%[beg,end]` specifies the boundary conditions at the beginning and the end of domain boundaries in each coordinate direction by a negative integer from -1 through -16.
See table [Boundary Conditions](#boundary-conditions) for details.
Boundary condition patches can be used with non-characteristic boundary conditions.
Their use is detailed in [Boundary Condition Patches](#boundary-condition-patches).

- `bc_[x,y,z]%%vb[1,2,3]` specifies the velocity in the (x,1), (y,2), (z,3) direction applied to `bc_[x,y,z]%%beg` when using `bc_[x,y,z]%%beg = -16`.
Tangential velocities require viscosity, `weno_avg = T`, and `bc_[x,y,z]%%beg = -16` to work properly. Normal velocities require `bc_[x,y,z]%%end = -15` or `\bc_[x,y,z]%%end = -16` to work properly.

- `bc_[x,y,z]%%ve[1,2,3]` specifies the velocity in the (x,1), (y,2), (z,3) direction applied to `bc_[x,y,z]%%beg` when using `bc_[x,y,z]%%end = -16`.
Tangential velocities require viscosity, `weno_avg = T`, and `bc_[x,y,z]%%end = 16` to work properly. Normal velocities require `bc_[x,y,z]%%end = -15` or `\bc_[x,y,z]%%end = -16` to work properly.

- `model_eqns` specifies the choice of the multi-component model that is used to formulate the dynamics of the flow using integers from 1 through 3.
`model_eqns = 1`, `2`, and `3` correspond to $\Gamma$-$\Pi_\infty$ model ([Johnsen, 2008](references.md)), 5-equation model ([Allaire et al., 2002](references.md)), and 6-equation model ([Saurel et al., 2009](references.md)), respectively.
The difference of the two models is assessed by ([Schmidmayer et al., 2019](references.md)).
Note that some code parameters are only compatible with 5-equation model.

- `alt_soundspeed` activates the source term in the advection equations for the volume fractions, $K\nabla\cdot \underline{u}$, that regularizes the speed of sound in the mixture region when the 5-equation model is used.
The effect and use of the source term are assessed by [Schmidmayer et al., 2019](references.md).

- `adv_n` activates the direct computation of number density by the Riemann solver instead of computing number density from the void fraction in the method of classes.

- `mpp_lim` activates correction of solutions to avoid a negative void fraction of each component in each grid cell, such that $\alpha_i>\varepsilon$ is satisfied at each time step.

- `mixture_err` activates correction of solutions to avoid imaginary speed of sound at each grid cell.

- `time_stepper` specifies the order of the Runge-Kutta (RK) time integration scheme that is used for temporal integration in simulation, from the 1st to 5th order by corresponding integer.
Note that `time_stepper = 3` specifies the total variation diminishing (TVD), third order RK scheme ([Gottlieb and Shu, 1998](references.md)).

- `adap_dt` activates the Strang operator splitting scheme which splits flux and source terms in time marching, and an adaptive time stepping strategy is implemented for the source term. It requires ``bubbles_euler = 'T'``, ``polytropic = 'T'``, ``adv_n = 'T'`` and `time_stepper = 3`. Additionally, it can be used with ``bubbles_lagrange = 'T'`` and `time_stepper = 3`. `adap_dt_tol` and `adap_dt_max_iters` are 1e-4 and 100, respectively, by default.

- `weno_order` specifies the order of WENO scheme that is used for spatial reconstruction of variables by an integer of 1, 3, 5, and 7, that correspond to the 1st, 3rd, 5th, and 7th order, respectively.

- `weno_eps` specifies the lower bound of the WENO nonlinear weights.
It is recommended to set `weno_eps` to $10^{-6}$ for WENO-JS, and to $10^{-40}$ for other WENO variants.

- `mapped_weno` activates the WENO-M scheme in place of the default WENO-JS scheme ([Henrick et al., 2005](references.md)). WENO-M a variant of the WENO scheme that remaps the nonlinear WENO-JS weights by assigning larger weights to non-smooth stencils, reducing dissipation compared to the default WENO-JS scheme, at the expense of higher computational cost. Only one of `mapped_weno`, `wenoz`, and `teno` can be activated.

- `wenoz` activates the WENO-Z scheme in place of the default WENO-JS scheme ([Borges et al., 2008](references.md)). WENO-Z is a variant of the WENO scheme that further reduces the dissipation compared to the WENO-M scheme. It has similar computational cost to the WENO-JS scheme.

- `wenoz_q` specifies the power parameter `q` used in the WENO-Z scheme. It controls how aggressively the smoothness coefficients scale the weights. A higher value of `wenoz_q` increases the sensitivity to smoothness, improving stability but worsening numerical dissipation. For WENO3 and WENO5, `q=1` is fixed, so `wenoz_q` must not be set. For WENO7, `wenoz_q` can be set to 2, 3, or 4.

- `teno` activates the TENO scheme in place of the default WENO-JS scheme ([Fu et al., 2016](references.md)). TENO is a variant of the ENO scheme that is the least dissipative, but could be less robust for extreme cases. It uses a threshold to identify smooth and non-smooth stencils, and applies optimal weights to the smooth stencils. Only available for `weno_order = 5` and `7`. Requires `teno_CT` to be set. Does not support grid stretching.

- `teno_CT` specifies the threshold for the TENO scheme. This dimensionless constant, also known as $C_T$, sets a threshold to identify smooth and non-smooth stencils. Larger values make the scheme more robust but also more dissipative. A recommended value for teno_CT is `1e-6`. When adjusting this parameter, it is recommended to try values like `1e-5` or `1e-7` for TENO5. A smaller value can be used for TENO7.

- `null_weights` activates nullification of the nonlinear WENO weights at the buffer regions outside the domain boundaries when the Riemann extrapolation boundary condition is specified (`bc_[x,y,z]%%beg[end]} = -4`).

- `mp_weno` activates monotonicity preservation in the WENO reconstruction (MPWENO) such that the values of reconstructed variables do not reside outside the range spanned by WENO stencil ([Balsara and Shu, 2000](references.md); [Suresh and Huynh, 1997](references.md)).

- `muscl_order` specifies the order of the MUSCL scheme that is used for spatial reconstruction of variables by an integer of 1, or 2, that corresponds to the 1st, and 2nd order respectively. When using `muscl_order = 2`, `muscl_lim` must be defined. 

- `muscl_lim` specifies the slope limiter that is used in 2nd order MUSCL Reconstruction by an integer from 1 through 5. 
`muscl_lim = 1`, `2`, `3`, `4`, and `5` correspond to minmod, monotonized central, Van Albada, Van Leer, and SUPERBEE, respectively.

- `int_comp` activates interface compression using THINC used in MUSCL Reconstruction, with control parameters (`ic_eps`, and `ic_beta`).

- `riemann_solver` specifies the choice of the Riemann solver that is used in simulation by an integer from 1 through 4.
`riemann_solver = 1`, `2`, and `3` correspond to HLL, HLLC, and Exact Riemann solver, respectively ([Toro, 2013](references.md)).
`riemann_solver = 4` is only for MHD simulations. It resolves 5 of the full seven-wave structure of the MHD equations ([Miyoshi and Kusano, 2005](references.md)).

- `low_Mach` specifies the choice of the low Mach number correction scheme for the HLLC Riemann solver. `low_Mach = 0` is default value and does not apply any correction scheme. `low_Mach = 1` and `2` apply the anti-dissipation pressure correction method ([Chen et al., 2022](references.md)) and the improved velocity reconstruction method ([Thornber et al., 2008](references.md)). This feature requires `model_eqns = 2` or `3`. `low_Mach = 1` works for `riemann_solver = 1` and `2`, but `low_Mach = 2` only works for `riemann_solver = 2`.

- `avg_state` specifies the choice of the method to compute averaged variables at the cell-boundaries from the left and the right states in the Riemann solver by an integer of 1 or 2.
`avg_state = 1` and `2` correspond to Roe- and arithmetic averages, respectively.

- `wave_speeds` specifies the choice of the method to compute the left, right, and middle wave speeds in the Riemann solver by an integer of 1 and 2.
`wave_speeds = 1` and `2` correspond to the direct method ([Batten et al., 1997](references.md)), and indirect method that approximates the pressures and velocity ([Toro, 2013](references.md)), respectively.

- `weno_Re_flux` activates the scalar divergence theorem in computing the velocity gradients using WENO-reconstructed cell boundary values.
If this option is false, velocity gradient is computed using finite difference scheme of order 2 which is independent of the WENO order.

- `weno_avg` it activates the arithmetic average of the left and right, WENO-reconstructed, cell-boundary values.
This option requires `weno_Re_flux` to be true because cell boundary values are only utilized when employing the scalar divergence method in the computation of velocity gradients.

- `surface_tension` activates surface tension when set to ``'T'``. Requires `sigma` to be set and `num_fluids`. The color function in each patch should be assigned such that `patch_icpp(i)%cf_val = 1` in patches where `patch_icpp(i)%alpha = 1 - eps` and `patch_icpp(i)%cf_val = 0` in patches where `patch_icpp(i)%alpha = eps`.

- `viscous` activates viscosity when set to ``'T'``. Requires `Re(1)` and `Re(2)` to be set.

- `hypoelasticity` activates elastic stress calculations for fluid-solid interactions. Requires `G` to be set in the fluid material's parameters.

#### Boundary Condition Patches

| Parameter              | Type    | Description                                                     |
| ---:                   | :----:  |          :---                                                   |
| `num_bc_patches`       | Integer | Number of boundary condition patches                            |
| `dir`*                 | Integer | Direction of the boundary patch. [1]: x; [2]: y; [3]: z         |
| `loc`*                 | Integer | Location of the patch in the domain                             |
| `type`*                | Integer | The geometry of the patch. [1]: Line [2]: Circle [3]: Rectangle |
| `x[y,z]_centroid`*     | Real    | Centroid of the boundary patch in the x[y,z]-direction          |
| `length_x[y,z]`*       | Real    | Length of the boundary patch in the x[y,z]-direction            |
| `radiue`*              | Real    | Radius of the boundary patch                                    |
*: These parameters should be prepended with `patch_bc(j)%` where $j$ is the patch index.

Boundary condition patches can be used with the following boundary condition types:
- `-2` reflective
- `-3` ghost cell extrapolation
- `-15` slip wall
- `-16` no-slip wall
- `-17` Dirichlet

Line segments along each domain edge are supported for 2D simulations.
Squares and circles on each face are supported for 3D simulations.

#### Constant Time-Stepping

- `dt` specifies the constant time step size used in the simulation.
The value of `dt` needs to be sufficiently small to satisfy the Courant-Friedrichs-Lewy (CFL) condition.

- `t_step_start` and `t_step_end` define the time steps at which the simulation starts and ends.

`t_step_save` is the time step interval for data output during simulation.
To newly start the simulation, set `t_step_start = 0`.
To restart the simulation from $k$-th time step, set `t_step_start = k`; see [Restarting Cases](running.md).

##### Adaptive Time-Stepping

- `cfl_adap_dt` enables adaptive time stepping with a constant CFL when true

- `cfl_const_dt` enables constant `dt` time-stepping where `dt` results in a specified CFL for the initial condition

- `cfl_target` specifies the target CFL value

- `n_start` specifies the save file to start at

- `t_save` specifies the time interval between data output during the simulation

- `t_stop` specifies at what time the simulation should stop

To newly start the simulation, set `n_start = 0`.
To restart the simulation from $k$-th time step, see [Restarting Cases](running.md).

### 7. Formatted Output

| Parameter            | Type    | Description                                    |
| ---:                 | :----:  |          :---                                  |
| `format`             | Integer | Output format. [1]: Silo-HDF5; [2] Binary	|
| `precision`          | Integer | [1] Single; [2] Double	 |
| `parallel_io`        | Logical | Parallel I/O	|
| `file_per_process`   | Logical | Whether or not to write one IO file per process |
| `cons_vars_wrt`      | Logical | Write conservative variables |
| `prim_vars_wrt`      | Logical | Write primitive variables	|
| `alpha_rho_wrt(i)`   | Logical | Add the partial density of the fluid $i$ to the database \|
| `rho_wrt`            | Logical | Add the mixture density to the database	 |
| `mom_wrt(i)`         | Logical | Add the $i$-direction momentum to the database	 |
| `vel_wrt(i)`         | Logical | Add the $i$-direction velocity to the database	  |
| `E_wrt`              | Logical | Add the total energy to the database	 |
| `pres_wrt`           | Logical | Add the pressure to the database	|
| `alpha_wrt(i)`       | Logical | Add the volume fraction of fluid $i$ to the database	|
| `gamma_wrt`          | Logical | Add the specific heat ratio function to the database	|
| `heat_ratio_wrt`     | Logical | Add the specific heat ratio to the database	|
| `pi_inf_wrt`         | Logical | Add the liquid stiffness function to the database |
| `pres_inf_wrt`       | Logical | Add the liquid stiffness to the formatted database	 |
| `c_wrt`              | Logical | Add the sound speed to the database	 |
| `omega_wrt(i)`       | Logical | Add the $i$-direction vorticity to the database	 |
| `schlieren_wrt`      | Logical | Add the numerical schlieren to the database|
| `qm_wrt`             | Logical | Add the Q-criterion to the database|
| `tau_wrt`            | Logical | Add the elastic stress components to the database|
| `fd_order`           | Integer | Order of finite differences for computing the vorticity and the numerical Schlieren function [1,2,4] |
| `schlieren_alpha(i)` | Real    | Intensity of the numerical Schlieren computed via `alpha(i)` |
| `probe_wrt`          | Logical | Write the flow chosen probes data files for each time step	|
| `num_probes`         | Integer | Number of probes	|
| `probe(i)%[x,y,z]`   | Real	 | Coordinates of probe $i$	|
| `output_partial_domain` | Logical | Output part of the domain |
| `[x,y,z]_output%beg` | Real    | Beginning of the output domain in the [x,y,z]-direction |
| `[x,y,z]_output%end` | Real    | End of the output domain in the [x,y,z]-direction |

The table lists formatted database output parameters. The parameters define variables that are outputted from simulation and file types and formats of data as well as options for post-processing.

- `format` specifies the choice of the file format of data file outputted by MFC by an integer of 1 and 2. `format = 1` and `2` correspond to Silo-HDF5 format and binary format, respectively.

- `precision` specifies the choice of the floating-point format of the data file outputted by MFC by an integer of 1 and 2. `precision = 1` and `2` correspond to single-precision and double-precision formats, respectively.

- `parallel_io` activates parallel input/output (I/O) of data files. It is highly recommended to activate this option in a parallel environment.
With parallel I/O, MFC inputs and outputs a single file throughout pre-process, simulation, and post-process, regardless of the number of processors used.
Parallel I/O enables the use of different numbers of processors in each process (e.g., simulation data generated using 1000 processors can be post-processed using a single processor).

- `file_per_process` deactivates shared file MPI-IO and activates file per process MPI-IO.
The default behavior is to use a shared file.
File per process is useful when running on >10K ranks.
If `file_per_process` is true, then pre_process, simulation, and post_process must be run with the same number of ranks.

- `cons_vars_wrt` and `prim_vars_wrt` activate the output of conservative and primitive state variables into the database.

- `[variable's name]_wrt` activates the output of each specified variable into the database.

- `schlieren_alpha(i)` specifies the intensity of the numerical Schlieren of $i$-th component.

- `fd_order` specifies the order of the finite difference scheme used to compute the vorticity from the velocity field and the numerical schlieren from the density field using an integer of 1, 2, and 4.
`fd_order = 1`, `2`, and `4` correspond to the first, second, and fourth-order finite difference schemes.

- `probe_wrt` activates the output of state variables at coordinates specified by `probe(i)%[x;y,z]`.

- `output_partial_domain` activates the output of part of the domain specified by `[x,y,z]_output%beg` and `[x,y,z]_output%end`.
This is useful for large domains where only a portion of the domain is of interest.
It is not supported when `precision = 1` and `format = 1`. 
It also cannot be enabled with `flux_wrt`, `heat_ratio_wrt`, `pres_inf_wrt`, `c_wrt`, `omega_wrt`, `ib`, `schlieren_wrt`, or `qm_wrt`.

### 8. Acoustic Source {#acoustic-source}

| Parameter                             | Type    | Description |
| ---:                                  | :----:  | :--- |
| `acoustic_source`                     | Logical | Acoustic source module activation |
| `num_source`                          | Integer | Number of acoustic sources |
| `acoustic(i)%%support`                | Integer | Geometry of spatial support for the acoustic source |
| `acoustic(i)%%dipole`                 | Logical | Dipole source activation (optional; default = false -> monopole) |
| `acoustic(i)%%loc(j)`                 | Real    | $j$-th coordinate of the point that defines the acoustic source location |
| `acoustic(i)%%pulse`                  | Integer | Acoustic wave form: [1] Sine [2] Gaussian [3] Square [4] Broadband  |
| `acoustic(i)%%npulse`                 | Real    | Number of pulse cycles |
| `acoustic(i)%%mag`                    | Real    | Pulse magnitude	|
| `acoustic(i)%%frequency`              | Real    | Sine/Square - Frequency of the acoustic wave  (exclusive) |
| `acoustic(i)%%wavelength`             | Real    | Sine/Square - Wavelength of the acoustic wave (exclusive) |
| `acoustic(i)%%gauss_sigma_time`       | Real    | Gaussian - Gaussian pulse time width in terms of sigma  (exclusive) |
| `acoustic(i)%%gauss_sigma_dist`       | Real    | Gaussian - Gaussian pulse spatial width in terms of sigma (exclusive) |
| `acoustic(i)%%delay`                  | Real    | Time delay of the acoustic wave (optional for `%%pulse = 1` or `3`; default = 0) |
| `acoustic(i)%%dir`                    | Real    | Planer - Direction of acoustic propagation |
| `acoustic(i)%%length`                 | Real    | 2D/3D Planer - Spatial pulse length |
| `acoustic(i)%%height`                 | Real    | 3D Planer - Spatial pulse height |
| `acoustic(i)%%foc_length`             | Real    | Transducer - Focal length of the transducer |
| `acoustic(i)%%aperture`               | Real    | Transducer - Aperture of the transducer |
| `acoustic(i)%%num_elements`           | Integer | Transducer array - Number of transducer elements in a transducer array |
| `acoustic(i)%%element_on`             | Integer | Transducer array - Element number that is on (optional; default = 0 -> all elements) |
| `acoustic(i)%%element_spacing_angle`  | Real    | 2D Transducer array - Spacing angle (in rad) between adjacent transducer elements |
| `acoustic(i)%%element_polygon_ratio`  | Real    | 3D Transducer array - Ratio of polygon side length to transducer element radius |
| `acoustic(i)%%rotate_angle`           | Real    | 3D Transducer array - Rotation angle of the transducer array (optional; default = 0) |
| `acoustic(i)%%bb_num_freq`            | integer | Number of frequencies in broadband wave |
| `acoustic(i)%%bb_bandwidth`           | Real    | The bandwidth of each frequency in the broadband wave |
| `acoustic(i)%%bb_lowest_freq`         | Real    | The lower frequency bound of the broadband wave |

Details of the transducer acoustic source model can be found in [Maeda and Colonius (2017)](references.md).

- `acoustic_source` activates the acoustic source module.

- `num_source` defines the total number of source planes by an integer.

- `%%support` specifies the choice of the geometry of acoustic source distribution. See table [Acoustic Supports](#acoustic-supports) for details.

- `%%dipole` changes the default monopole (one-sided) source to a dipole source. It is only available for planar waves.

- `%%loc(j)` specifies the location of the acoustic source in the $j$-th coordinate direction. For planer support, the location defines midpoint of the source plane. For transducer arrays, the location defines the center of the transducer or transducer array (not the focal point; for 3D it's the tip of the spherical cap, for 2D it's the tip of the arc).

- `%%pulse` specifies the acoustic wave form. `%%pulse = 1`, `2`, `3` and `4` correspond to sinusoidal wave, Gaussian wave, square wave and broadband wave, respectively. The implementation of the broadband wave is based on [Tam (2005)](references.md)

- `%%npulse` specifies the number of cycles of the acoustic wave generated. Only applies to `%%pulse = 1 and 3` (sine and square waves), and must be an integer for non-planar waves.

- `%%mag` specifies the peak amplitude of the acoustic wave.

- `%%frequency` and `%%wavelength` specify the frequency and wavelength of the acoustic wave, respectively. These parameters are exclusive and exactly one of them must be specified for `%%pulse = 1` or `3` (sine or square waves). They are related by the speed of sound in the medium: `frequency = speed_of_sound / wavelength`.

- `%%gauss_sigma_time` and `%%gauss_sigma_dist` specify the time and spatial widths of the Gaussian pulse in terms of sigma, respectively. In particular, `%%gauss_sigma_time` is the standard deviation in the Gaussian equation. These parameters are exclusive and exactly one of them must be specified for `%%pulse = 2` (Gaussian wave). They are related by the speed of sound in the medium: `gauss_sigma_dist = speed_of_sound * gauss_sigma_time`.

- `%%delay` specifies the time delay of the acoustic wave. This parameter is optional for `%%pulse = 1` or `3` (sine or square waves) and defaults to 0. It must be specified for `%%pulse = 2` (Gaussian wave). It is important to note that setting the delay to 0 for a Gaussian pulse results in a half-Gaussian pulse, and delays that are too small may result in the pulse being cut off at the start of the simulation. `4*gauss_sigma_time` is a typical value for the delay of a Gaussian pulse.

- `%%dir` specifies the direction of acoustic wave propagation for planar waves. The direction is defined by the angle in degrees from the x-axis in the x-y plane. It applies to both 2D and 3D simulation of planar waves (support is infinite in z-direction for 3D).

- `%%length` specifies the spatial length of the 2D or 3D planar wave. It is the length of the source plane perpendicular to the direction of wave propagation.

- `%%height` specifies the spatial height of the planar wave. Since `%%dir` is in the x-y plane, the height is perpendicular to the direction of wave propagation.

- `%%foc_length` specifies the focal length of the transducer for transducer waves. It is the distance from the transducer to the focal point.

- `%%aperture` specifies the aperture of the transducer. It is the diameter of the projection of the transducer arc onto the y-axis (2D) or spherical cap onto the y-z plane (3D). Set the aperture to double the focal length to simulate a transducer enclosing half of the circle/sphere. For the transducer array, it is the total aperture of the array.

- `%%num_elements` specifies the number of transducer elements in a transducer array.

- `%%element_on` specifies the element number of the transducer array that is on. The element number starts from 1, if all elements are on, set `%%element_on` to 0.

- `%%element_spacing_angle` specifies the spacing angle between adjacent transducers in radians. The total aperture (`%%aperture`) is set, so each transducer element is smaller if `%%element_spacing_angle` is larger.

- `%%element_polygon_ratio` specifies the ratio of the polygon side length to the aperture diameter of each transducer element in a circular 3D transducer array. The polygon side length is calculated by using the total aperture (`%%aperture`) as the circumcircle diameter and `%%num_elements` as the number of sides of the polygon. The ratio is used to specify the aperture size of each transducer element in the array as a ratio of the total aperture.

- `%%rotate_angle` specifies the rotation angle of the 3D circular transducer array along the x-axis (principal axis). It is optional and defaults to 0.

- `%%bb_num_freq` specifies the number discretized frequencies in the broadband acoustic wave. If `%%bb_num_freq` is 1, the acoustic wave will be a discrete tone (i.e. single frequency sine wave).

- `%%bb_bandwidth` specifies the bandwidth of the discretized frequencies.

- `%%bb_lowest_freq` specifies the lower frequency bound of the broadband acoustic wave. The upper frequency bound will be calculated as `%%bb_lowest_freq + %%bb_num_freq * %%bb_bandwidth`. The wave is no longer broadband below the lower bound and above the upper bound.

### 9. Sub-grid Bubble Models

| Parameter         | Type    | Description                                    |
| ---:              | :----:  |          :---                                  |
| `bubbles_euler`   | Logical	| Ensemble-averaged bubble modeling	|
| `bubbles_lagrange`| Logical	| Volume-averaged bubble modeling	|
| `bubble_model`    | Integer	| [1] Gilmore; [2] Keller--Miksis |
| `mu_l0` *	        | Real 		|	Liquid viscosity (only specify in liquid phase)  |
| `ss` *		        | Real 		|	Surface tension (only specify in liquid phase) |
| `pv` *		        | Real 		|	Vapor pressure (only specify in liquid phase) |
| `gamma_v` † 	    | Real 	  |	Specific heat ratio |
| `M_v` †     	    | Real 		| Molecular weight |
| `mu_v` †	        | Real 		|	Viscosity |
| `k_v` †	          | Real 		|	Thermal conductivity |
| `cp_v` †	        | Real 		|	Specific heat capacity |

These options work only for gas-liquid two-component flows.
Component indexes are required to be 1 for liquid and 2 for gas.

- \* These parameters should be prepended with patch index $1$ that is filled with liquid: `fluid_pp(1)%`.
- †  These parameters should be prepended with patch indexes filled with liquid and gas: `fluid_pp(1)%` and `fluid_pp(2)%`.

This table lists the sub-grid bubble model parameters, which can be utilized in both the ensemble-averaged and volume-averaged bubble models.

- `bubbles_euler` activates the ensemble-averaged bubble model.

- `bubbles_lagrange` activates the volume-averaged bubble model.

- `bubble_model` specified a model for spherical bubble dynamics by an integer of 1 and 2.
`bubble_model = 1`, `2`, and `3` correspond to the Gilmore, Keller-Miksis, and Rayleigh-Plesset models.

- `mu_l0`, `ss`, and `pv`, `gamma_v`, `M_v`, `mu_v`, `k_v`, and `cp_v` specify simulation parameters for the non-polytropic gas compression model.
`mu_l0`, `ss`, and `pv` correspond to the liquid viscosity, surface tension, and vapor pressure, respectively.
`gamma_v`, `M_v`, `mu_v`, `k_v`, and `cp_v` specify the specific heat ratio, molecular weight, viscosity, thermal conductivity and specific heat capacity of a chosen component (`cp_v` only for ``bubbles_lagrange = 'T'``).
Implementation of the parameters into the model follow [Ando (2010)](references.md).

#### 9.1 Ensemble-Averaged Bubble Model

| Parameter      | Type    | Description                                    |
| ---:           | :----:  |          :---                                  |
| `bubbles_euler`| Logical	| Ensemble-averaged bubble modeling	|
| `polytropic`      | Logical	| Polytropic gas compression |
| `thermal` 		    | Integer	| Thermal model: [1] Adiabatic; [2] Isothermal; [3] Transfer |
| `R0ref` 			 | Real		  | Reference bubble radius |
| `polydisperse`   | Logical	| Polydispersity in equilibrium bubble radius R0 |
| `nb` 			     | Integer	| Number of bins: [1] Monodisperse; [$>1$] Polydisperse |
| `poly_sigma` 	       | Real 		|	Standard deviation for probability density function of polydisperse bubble populations |
| `Ca` 			     | Real		  | Cavitation number |
| `Web` 			   | Real		  | Weber number |
| `Re_inv` 		   | Real		  | Inverse Reynolds number |
| `qbmm` 	       | Logical 		|	Quadrature by  method of moments |
| `dist_type` 	       | Integer 		|	Joint probability density function for bubble radius and velocity (only for ``qbmm = 'T'``) |
| `sigR` 	       | Real 		|	Standard deviation for the probability density function of bubble radius (only for ``qbmm = 'T'``) (EE)|
| `sigV` 	       | Real 		|	Standard deviation for the probability density function of bubble velocity (only for ``qbmm = 'T'``) (EE)|
| `rhoRV`	       | Real 		|	Correlation coefficient for the joint probability density function of bubble radius and velocity (only for ``qbmm = 'T'``) |

This table lists the ensemble-averaged bubble model parameters.

- `polytropic` activates polytropic gas compression in the bubble.
When ``polytropic = 'F'``, the gas compression is modeled as non-polytropic due to heat and mass transfer across the bubble wall with constant heat and mass transfer coefficients based on ([Preston et al., 2007](references.md)).

- `thermal` specifies a model for heat transfer across the bubble interface by an integer from 1 through 3.
`thermal = 1`, `2`, and `3` correspond to no heat transfer (adiabatic gas compression), isothermal heat transfer, and heat transfer with a constant heat transfer coefficient based on [Preston et al., 2007](references.md), respectively.

- `polydisperse` activates polydispersity in the bubble model through a probability density function (PDF) of the equilibrium bubble radius. Simpson's rule is used for integrating the log-normal PDF of equilibrium bubble radius for polydisperse populations.

- `R0ref` specifies the reference bubble radius.

- `nb` specifies the number of discrete bins that define the probability density function (PDF) of the equilibrium bubble radius.

- `poly_sigma` specifies the standard deviation of the log-normal PDF of equilibrium bubble radius for polydisperse populations.

- `Ca`, `Web`, and `Re_inv` respectively specify the Cavitation number, Weber number, and the inverse Reynolds number that characterize the offset of the gas pressure from the vapor pressure, surface tension, and liquid viscosity when the polytropic gas compression model is used.

- `qbmm` activates quadrature by method of moments, which assumes a PDF for bubble radius and velocity.

- `dist_type` specifies the initial joint PDF of initial bubble radius and bubble velocity required in qbmm. `dist_type = 1`  and `2` correspond to binormal and lognormal-normal distributions respectively.

- `sigR` specifies the standard deviation of the PDF of bubble radius required in the QBMM feature.

- `sigV` specifies the standard deviation of the PDF of bubble velocity required in the QBMM feature.

- `rhoRV` specifies the correlation coefficient of the joint PDF of bubble radius and bubble velocity required in the QBMM feature.

#### 9.2 Volume-Averaged Bubble Model

| Parameter             | Type    | Description                                               |
| ---:                  | :---:   | :---                                                      |
| `bubbles_lagrange`    | Logical | Lagrangian subgrid bubble model switch                    |
| `nBubs_glb`           | Integer | Global number of bubbles                                  |
| `solver_approach`     | Integer | 1: One-way coupling, 2: two-way coupling                  |
| `cluster_type`        | Integer | Method to find p_inf                                      |
| `pressure_corrector`  | Logical | Cell pressure correction term                             |
| `smooth_type`         | Integer | Smoothing function. 1: Gaussian, 2:Delta 3x3              |
| `heatTransfer_model`  | Logical | Activates the interface heat transfer model               |
| `massTransfer_model`  | Logical | Activates the interface mass transfer model               |
| `write_bubbles`       | Logical | Write files to track the bubble evolution each time step  |
| `write_bubbles_stats` | Logical | Write the maximum and minimum radius of each bubble       |
| `epsilonb`            | Real    | Standard deviation scaling for the gaussian function      |
| `charwidth`           | Real    | Domain virtual depth (z direction, for 2D simulations)    |
| `valmaxvoid`          | Real    | Maximum void fraction permitted                           |
| `c0`                  | Real    | Reference speed                                           |
| `rho0`                | Real    | Reference density                                         |
| `T0`                  | Real    | Reference temperature                                     |
| `x0`                  | Real    | Reference length                                          |
| `Thost`               | Real    | Temperature of the surrounding liquid (host)              |
| `diffcoefvap`         | Real    | Vapor diffusivity in the gas                              |

- `nBubs_glb` Total number of bubbles. Their initial conditions need to be specified in the ./input/lag_bubbles.dat file. See the example cases for additional information.

- `solver_approach` Specifies the Euler-Lagrange coupling method: [1] enables a one-way coupling approach, where the bubbles do not influence the Eulerian field. [2] activates the two-way coupling approach based on [Maeda and Colonius (2018)](references.md), where the effect of the bubbles is added in the Eulerian field as source terms.

- `cluster_type` Specifies method to find p_inf (pressure that drives the bubble dynamics): [1] activates the bilinear interpolation of the pressure field, while [2] enables the bubble dynamic closure based on [Maeda and Colonius (2018)](references.md), the full model is obtained when `pressure_corrector` is true.

- `smooth_type` Specifies the smoothening method of projecting the lagrangian bubbles in the Eulerian field: [1] activates the gaussian kernel function described in  [Maeda and Colonius (2018)](references.md), while [2] activates the delta kernel function where the effect of the bubble is only seen in the specific bubble location cell.

- `heatTransfer_model` Activates the heat transfer model at the bubble's interface based on ([Preston et al., 2007](references.md)).

- `massTransfer_model` Activates the mass transfer model at the bubble's interface based on ([Preston et al., 2007](references.md)).

### 10. Velocity Field Setup

| Parameter              | Type    | Description |
| ---:                   | :----:  | :--- |
| `perturb_flow`         | Logical | Perturb the initial velocity field by random noise |
| `perturb_flow_fluid`   | Integer | Fluid density whose flow is to be perturbed |
| `perturb_flow_mag`     | Real    | Set the magnitude of flow perturbations |
| `perturb_sph`          | Logical | Perturb the initial partial density by random noise |
| `perturb_sph_fluid`    | Integer | Fluid component whose partial density is to be perturbed |
| `mixlayer_vel_profile` | Logical | Set the mean streamwise velocity to hyperbolic tangent profile |
| `mixlayer_vel_coef`    | Real    | Coefficient for the hyperbolic tangent profile of a mixing layer |
| `mixlayer_perturb`     | Logical | Perturb the initial velocity field using a spectrum-based synthetic turbulence generation method |

The table lists velocity field parameters.
The parameters are optionally used to define initial velocity profiles and perturbations.

- `perturb_flow` activates the perturbation of initial velocity by random noise.

- `perturb_flow_fluid` specifies the fluid component whose flow will be perturbed.

- `perturb_flow` activates the perturbation of initial velocity by random noise.

- `perturb_sph` activates the perturbation of initial partial density by random noise.

- `perturb_sph_fluid` specifies the fluid component whose partial density is to be perturbed.

- `mixlayer_vel_profile` activates setting the mean streamwise velocity to a hyperbolic tangent profile. This option works only for `n > 0`.

- `mixlayer_vel_coef` is a parameter for the hyperbolic tangent profile of a mixing layer when `mixlayer_vel_profile = 'T'`. The mean streamwise velocity profile is given as:

$$ u = \mbox{patch\_icpp(1)\%vel(1)} * \tanh( y_{cc} * \mbox{mixlayer\_vel\_coef}) $$

- `mixlayer_perturb` activates the velocity perturbation for a temporal mixing layer with hyperbolic tangent mean streamwise velocity profile, using an inverter version of the spectrum-based synthetic turbulence generation method proposed by Guo et al. (2023, JFM). This option only works for `p > 0` and `mixlayer_vel_profile = 'T'`.

### 11. Phase Change Model
| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `relax`                | Logical | Activates Phase Change model |
| `relax_model`          | Integer | Phase change model: [5] pT-equilibrium; [6] pTg-equilibrium |
| `palpha_eps`           | Real    | tolerance of the Newton Solver to activate pT-equilibrium  |
| `ptgalpha_eps`	     | Real    | tolerance of the Newton Solver to activate pTg-equilibrium |

- `relax` Activates the Phase Change model.

- `relax_model` Specifies the phase change model to be used: [5] enables pT-equilibrium, while [6] activates pTg-equilibrium (if criteria are met).

- `palpha_eps` Specifies the tolerance for the Newton Solvers used in the pT-equilibrium model.

- `ptgalpha_eps` Specifies the tolerance used for the Newton Solvers used in the pTg-equilibrium model.

### 12. Artificial Mach Number
| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `pi_fac`               | Real    | Ratio of artificial and true `pi_\infty` values|

- `pi_fac` specifies the ratio of artificial and true `pi_\infty` values (`=` artificial `pi_\infty` / true `pi_\infty`).
This parameter enables the use of true `pi_\infty` in bubble dynamics models when the `pi_\infty` given in the `case.py` file is an artificial value.

### 13. Body Forces

| Parameter         | Type    | Description                                |
| ---:              | :---:   | :---                                       |
| `bf_x[y,z]`       | Logical | Enable body forces in the x[y,z] direction |
| `k_x[y,y]`        | Real    | Magnitude of oscillating acceleration      |
| `w_x[y,z]`        | Real    | Frequency of oscillating acceleration      |
| `p_x[y,z]`        | Real    | Phase shift of oscillating acceleration    |
| `g_x[y,z]`        | Real    | Magnitude of background acceleration       |

`k_x[y,z]`, `w_x[y,z]`, `p_x[y,z]`, and `g_x[y,z]` define an oscillating acceleration in the `x[y,z]` direction with the form

$$ a_{x[y,z]} = g_{x[y,z]} + k_{x[y,z]}\sin\left(w_{x[y,z]}t + p_{x[y,z]}\right). $$

By convention, positive accelerations in the `x[y,z]` direction are in the positive `x[y,z]` direction.

### 14. Magnetohydrodynamics (MHD)

| Parameter         | Type    | Description                                         |
| ---:              | :---:   | :---                                                |
| `mhd`             | Logical | Enable ideal MHD simulation                         |
| `relativity`      | Logical | Enable relativistic MHD simulation                  |
| `powell`          | Logical | Enable Powell's method for solenoidal constraint    |
| `fd_order`        | Integer | Finite difference order for Powell's method         |
| `Bx[y,z]`         | Real    | Initial magnetic field in the x[y,z] direction      |
| `Bx0`             | Real    | Constant magnetic field in the x direction (1D only)|

- `mhd` is currently only available for single-component flows and 5-equation model. Its compatibility with most other features is work in progress.

- `relativity` only works for `mhd` enabled and activates relativistic MHD (RMHD) simulation.

- `powell` activates Powell's eight-wave method to impose the solenoidal constraint in the MHD simulation [Powell (1994)](references.md). It should not be used in conjunction with HLLD (`riemann_solver = 4`).

- `fd_order` specifies the finite difference order for computing the RHS of the Powell's method. `fd_order = 1`, `2`, and `4` are allowed.

- `Bx0` is only used in 1D simulations to specify the constant magnetic field in the x direction. It must be specified in 1D simulations. `Bx` must not be used in 1D simulations.

- `Bx`, `By`, and `Bz` are used to specify the initial magnetic field in the x, y, and z directions, respectively. They must be specified in all 1D/2D/3D MHD simulations, with the exception of `Bx` in 1D simulations.

Note: In 1D/2D/3D simulations, all three velocity components are treated as state variables and must be specified in the case file.

Note: For relativistic flow, the conservative and primitive densities are different. `rho_wrt` outputs the primitive (rest mass) density.

### 15. Elasticity

| Parameter         | Type    | Description                                         |
| ---:              | :---:   | :---                                                |
| `hypoelasticity`  | Logical | Enable hypoelasticity simulation                    |
| `hyperelasticity` | Logical | Enable hyperelasticity simulation                   |
| `cont_damage`     | Logical | Enable continuum damage model                       |
| `tau_star`        | Real    | Threshold stress for continuum damage model         |
| `cont_damage_s`   | Real    | Power `s` for continuum damage model                |
| `alpha_bar`       | Real    | Damage factor (rate) for continuum damage model     |

- `cont_damage` activates continuum damage model for solid materials. Requires `tau_star`, `cont_damage_s`, and `alpha_bar` to be set (empirically determined) ([Cao et al., 2019](references.md)).

### 16. Cylindrical Coordinates

When ``cyl_coord = 'T'`` is set in 3D the following constraints must be met:

- `bc_y%beg = -14`  enables the axis boundary condition

- `bc_z%beg = bc_z%end = -1`  enables periodic boundary conditions in the azimuthal direction

- `z_domain%beg = 0`  sets the azimuthal starting point to 0

- `z_comain%end = 2*math.pi` to set the azimuthal ending point to $2\pi$ (note, requires `import math` in the case file)

When ``cyl_coord = 'T'`` is set in 2D the following constraints must be met:

- `bc_y%beg = -2` to enable reflective boundary conditions

## Enumerations

### Boundary conditions

| #    | Type           | Description |
| ---: | :----:         | :---        |
|  -1  | Normal         | Periodic |
|  -2  | Normal         | Reflective |
|  -3  | Normal         | Ghost cell extrapolation |
|  -4  | Normal         | Riemann extrapolation |
|  -5  | Characteristic | Slip wall |
|  -6  | Characteristic | Non-reflecting subsonic buffer |
|  -7  | Characteristic | Non-reflecting subsonic inflow |
|  -8  | Characteristic | Non-reflecting subsonic outflow |
|  -9  | Characteristic | Force-free subsonic outflow |
|  -10 | Characteristic | Constant pressure subsonic outflow |
|  -11 | Characteristic | Supersonic inflow |
|  -12 | Characteristic | Supersonic outflow |
|  -14 | Normal         | Axis * |
|  -15 | Normal         | Slip wall |
|  -16 | Normal         | No-slip wall |

*: This boundary condition is only used for `bc_y%%beg` when using cylindrical coordinates (``cyl_coord = 'T'`` and 3D). For axisymmetric problems, use `bc_y%%beg = -2` with ``cyl_coord = 'T'`` in 2D.

The boundary condition supported by the MFC are listed in table [Boundary Conditions](#boundary-conditions).
Their number (`#`) corresponds to the input value in `input.py` labeled `bc_[x,y,z]%[beg,end]` (see table [Simulation Algorithm Parameters](#5-simulation-algorithm)).
The entries labeled "Characteristic." are characteristic boundary conditions based on [Thompson (1987)](references.md) and [Thompson (1990)](references.md).

### Generalized Characteristic Boundary conditions

| Parameter                     | Type    | Description |
| ---:                          | :----:  | :--- |
| `bc_[x,y,z]%grcbc_in`         | Logical | Enable grcbc for subsonic inflow |
| `bc_[x,y,z]%grcbc_out`        | Logical | Enable grcbc for subsonic outflow (pressure)|
| `bc_[x,y,z]%grcbc_vel_out`    | Logical | Enable grcbc for subsonic outflow (pressure + normal velocity) |
| `bc_[x,y,z]%vel_in`           | Real Array | Inflow velocities in x, y and z directions |
| `bc_[x,y,z]%vel_out`          | Real Array | Outflow velocities in x, y and z directions |
| `bc_[x,y,z]%pres_in`          | Real    | Inflow pressure |
| `bc_[x,y,z]%pres_out`         | Real    | Outflow pressure |
| `bc_[x,y,z]%alpha_rho_in`     | Real Array | Inflow density |
| `bc_[x,y,z]%alpha_in`         | Real Array | Inflow void fraction |

This boundary condition can be used for subsonic inflow (`bc_[x,y,z]%[beg,end]` = -7) and subsonic outflow (`bc_[x,y,z]%[beg,end]` = -8) characteristic boundary conditions. These are based on [Pirozzoli (2013)](references.md). This enables to provide inflow and outflow conditions outside the computational domain.

### Patch types

| #    | Name               | Dim.  | Smooth | Description |
| ---: | :----:             | :---: |  :---: | :--- |
| 1    | Line segment 	    | 1     | N      | Requires `x_centroid` and `length_x`. |
| 2    | Circle 		    | 2     | Y      | Requires `[x,y]_centroid` and `radius`. |
| 3    | Rectangle 	        | 2     | N      | Coordinate-aligned. Requires `[x,y]_centroid` and `length_[x,y]`. |
| 4    | Sweep line 		| 2     | Y      | Not coordinate aligned. Requires `[x,y]_centroid` and `normal(i)`. |
| 5    | Ellipse 		    | 2     | Y      | Requires `[x,y]_centroid` and `radii(i)`. |
| 6    | N/A 		        | N/A   | N/A    | No longer exists. Empty. |
| 7    | N/A        	    | N/A   | N/A    | No longer exists. Empty. |
| 8    | Sphere 		    | 3     | Y      | Requires `[x,y,z]_centroid` and `radius` |
| 9    | Cuboid 		    | 3     | N      | Coordinate-aligned. Requires `[x,y,z]_centroid` and `length_[x,y,z]`. |
| 10   | Cylinder 		    | 3     | Y      | Requires `[x,y,z]_centroid`, `radius`, and `length_[x,y,z]`. |
| 11   | Sweep plane 	    | 3     | Y      | Not coordinate-aligned. Requires `x[y,z]_centroid` and `normal(i)`. |
| 12   | Ellipsoid 		    | 3     | Y      | Requires `[x,y,z]_centroid` and `radii(i)`. |
| 13   | N/A        	    | N/A   | N/A    | No longer exists. Empty. |
| 14   | Spherical Harmonic | 3     | N      | Requires `[x,y,z]_centroid`, `radius`, `epsilon`, `beta` |
| 15   | N/A                | N/A   | N/A    | No longer exists. Empty.  |
| 16   | 1D bubble pulse    | 1     | N      | Requires `x_centroid`, `length_x` |
| 17   | Spiral             | 2     | N      | Requires `[x,y]_centroid` |
| 18   | 2D Varcircle       | 2     | Y      | Requires `[x,y]_centroid`, `radius`, and `thickness` |
| 19   | 3D Varcircle       | 3     | Y      | Requires `[x,y,z]_centroid`, `length_z`, `radius`, and `thickness` |
| 20   | 2D Taylor-Green Vortex  | 2  | N     | Requires `[x,y]_centroid`, `length_x`, `length_y`, `vel(1)`, and `vel(2)` |
| 21   | Model              | 2 & 3 | Y      | Imports a Model (STL/OBJ). Requires `model%%filepath`. |

The patch types supported by the MFC are listed in table [Patch Types](#patch-types).
This includes types exclusive to one-, two-, and three-dimensional problems.
The patch type number (`#`) corresponds to the input value in `input.py` labeled  `patch_icpp(j)%%geometry` where $j$ is the patch index.
Each patch requires a different set of parameters, which are also listed in this table.

### Immersed Boundary Patch Types

| #    | Name               | Dim.   |
| ---: | :----:             | :---   |
| 2    | 2D Circle          | 2      |
| 3    | 2D Rectangle       | 2      |
| 4    | 2D Airfoil         | 2      |
| 8    | 3D Sphere          | 3      |
| 10   | 3D Cylinder        | 3      |
| 11   | 3D Airfoil         | 3      |

### Acoustic Supports {#acoustic-supports}

| #    | Name                         | Dim.      | Requirements                                                                            |
| ---: | :----:                       | :---:     | :---                                                                                    |
|  1   | Planar source                | 1D        | `%%loc(1)`, `%%pulse`, `%%npulse`, `%%mag`, and `%%dir`                                 |
|  2   | Planar source                | 2D        | #1 requirements, `%%loc(2)` and `%%length`                                              |
|  3   | Planar source                | 3D        | #2 requirements and `%%height`                                                          |
|  5   | Cylindrical Transducer       | 2D        | `%%loc(1)`, `%%loc(2)`, `%%pulse`, `%%npulse`, `%%mag`, `%%foc_length`, `%%aperture`    |
|  6   | Spherical Transducer         | 2D-Axisym | #5 requirements                                                                         |
|  7   | Spherical Transducer         | 3D        | #5 requirements and `%%loc(3)`                                                          |
|  9   | Arcuate Transducer Array     | 2D        | #5 requirements, `%%num_elements`, `%%element_on`, `%%element_spacing_angle`            |
| 10   | Annular Transducer Array     | 2D-Axisym | #9 requirements                                                                         |
| 11   | Circular Transducer Array    | 3D        | #7 requirements, `%%element_polygon_ratio`, and `%%rotate_angle`                        |

The required parameters for each acoustic support type are listed in [Acoustic Source](#acoustic-source).
The acoustic support number (`#`) corresponds to the acoustic support type `Acoustic(i)%%support`, where $i$ is the acoustic source index.
For each `%%parameter`, prepend the parameter with `acoustic(i)%`.

Additional requirements for all acoustic support types:
- ``acoustic_source = 'T'`` must be used to activate the acoustic source module.

- `num_source` must be set to the total number of acoustic sources.

- `%%support` must be set to the acoustic support number in the table.

- `%%dipole` is only supported for planar sources.

- `%%npulse = 1 or 3` requires exactly one of `%%frequency` or `%%wavelength` to be set. It accepts `%%delay` as an optional parameter (default = 0).

- `%%npulse = 2` requires exactly one of `%%gauss_sigma_time` or `%%gauss_sigma_space` to be set. It requires `%%delay` to be set.

Description of the acoustic support types:
- `%%support = 1` specifies an infinite source plane that is normal to the $x$-axis and intersects with the axis at $x=$ `%%loc(1)` in 1D simulation. `%%dir > 0` specifies a rightward propagating wave, and `%%dir < 0` specifies a leftward propagating wave. `%%dir = 0` is not allowed.

- `%%support = 2` specifies a semi-infinite source plane in 2D simulation.
The midplane location is [`%%loc(1)`, `%%loc(2)`] and the normal vector is [$\mathrm{cos}$(`%%dir`), $\mathrm{sin}$(`%%dir`)]. The length of the source plane is `%%length`, and the plane is perpendicular to the direction of wave propagation (defined by `%%dir`).

- `%%support = 3` specifies a semi-infinite source plane in 3D simulation.
The midplane location is [`%%loc(1)`, `%%loc(2)`] and the normal vector is [$\mathrm{cos}$(`%%dir`), $\mathrm{sin}$(`%%dir`)]. The length of the source plane is `%%length`, and the plane is perpendicular to the direction of wave propagation (defined by `%%dir`). Note that the plane is infinite in the $z$-direction, so `%%loc(3)` is not required.

- `%%support = 5` specifies a circular transducer in 2D simulation. The transducer is centered at [`%%loc(1)`, `%%loc(2)`] with a focal length of `%%foc_length` and an aperture of `%%aperture`. The center location is not the focal point; it is the tip of the circular arc (intersection of the arc and the x-axis). The aperture is the length of the projection of the circular arc onto the y-axis. If a semi-circle is desired, set the aperture to double the focal length. Note that this is physically a cylindrical transducer, and due to the complexity of Green's function for 2D wave, no closed-form solution is available for the 2D circular transducer, and an approximate is used (see [Maeda and Colonius (2017)](references.md) for details). For the mass source term correction factor, the theoretical approximation factor of -0.5 in ($r_{foc}^{-0.5}$) is replaced by an empirically determined factor of -0.85.

- `%%support = 6` specifies a spherical transducer in 2D axisymmetric simulation. It is identical to `%%support = 5` in terms of simulation parameters. Note that this is physically a spherical 3D transducer, so the equation is exact.

- `%%support = 7` specifies a spherical transducer in 3D simulation. The transducer is centered at [`%%loc(1)`, `%%loc(2)`, `%%loc(3)`] with a focal length of `%%foc_length` and an aperture of `%%aperture`. The center location is not the focal point; it is the tip of the spherical cap (intersection of the cap and the x-axis). The aperture is the diameter of the projection of the spherical cap onto the y-z plane. If a semi-sphere is desired, set the aperture to double the focal length. Again, the equation is exact.

- `%%support = 9` specifies an arcuate transducer array in 2D simulation. The total aperture of the array is `%%aperture`, which is similar to `%%support = 5`. The parameters `%%num_elements` and `%%element_spacing_angle` specify the number of transducer elements and the spacing angle. The spacing angle is the angle of the gap between adjacent transducer elements in the array. Because the total aperture is set, each transducer element is smaller if the spacing angle is larger. Physically it represents curved panels. Note that similar to `%%support = 5`, the mass source term correction factor is empirically determined to be -0.85.

- `%%support = 10` specifies an annular transducer array in 2D axisymmetric simulation. It is identical to `%%support = 9` in terms of simulation parameters. It physically represents the a annulus obtained by revolving the arc in `%%support = 9` around the x-axis.

- `%%support = 11` specifies a circular transducer array in 3D simulation. The total aperture of the array is `%%aperture`, which is similar to `%%support = 7`. The parameters `%%num_elements`, `%%element_polygon_ratio`, and `%%rotate_angle` specify the number of transducer elements, the ratio of the polygon side length to the transducer element radius, and the rotation angle of the array. The polygon side length is calculated by using the total aperture as the circumcicle diameter, and the number of sides of the polygon as `%%num_elements`. The ratio is used specify the aperture size of each transducer element in the array, as a ratio of the total aperture. The rotation angle is optional and defaults to 0. Physically it represents a circular ring of transducer elements.

### Conservative Variables Ordering

| 5-eqn                           | 6-eqn |
| ----                            |  ---- |
| num_fluids continuity variables | num_fluids continuity variables        |
| num_dims momentum variables     | num_dims momentum variables          |
| 1 energy variable               | 1 energy variable                            |
| num_fluids advection variables  | num_fluids advection variables              |
| N/A                             | num_fluids internal energy variables |

The above variables are used for all simulations.

| 5-eqn                     | 6-eqn |
| ----                      |  ---- |
| sub-grid bubble variables | N/A   |
| hypoelastic variables     | N/A   |

The above variables correspond to optional physics.

### Primitive Variables Ordering

| 5-eqn                         | 6-eqn |
| -----                         | ----  |
| num_fluids densities          | num_fluids densities          |
| num_dims velocities           | num_dims velocities           |
| 1 pressure                    | 1 pressure                    |
| num_fluids volume fractions   | num_fluids volume fractions   |
| N/A                           | num_fluids partial pressures  |

The above variables are used for all simulations.

| 5-eqn | 6-eqn |
| ----  |  ---- |
| sub-grid bubble variables | N/A |
| hypoelastic variables     | N/A |

The above variables correspond to optional physics.
