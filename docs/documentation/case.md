# Case Files

Example Python case files, also referred to as *input files*, can be found in the [examples/](https://github.com/MFlowCode/MFC/tree/master/examples) directory. They print a Python dictionary containing input parameters for MFC. Their contents, and a guide to filling them out, are documented
in the user manual. A commented, tutorial script
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
import argparse

parser = argparse.ArgumentParser(
    prog="scaling",
    description="Weak- and strong-scaling benchmark case.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("dict", type=str, metavar="DICT")
parser.add_argument("-s", "--scaling",  type=str, metavar="SCALING",  choices=["weak", "strong"], help="Whether weak- or strong-scaling is being exercised.")

# Your parsed arguments are here
args = parser.parse_args()
```

The first argument is always a JSON string representing `mfc.sh run`'s internal
state.
It contains all the runtime information you might want from the build/run system.
You can add as many additional arguments as you may need.

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
| `dt`                     | Real    | Time step size |
| `t_step_start`           | Integer | Simulation starting time step |
| `t_step_stop`            | Integer | Simulation stopping time step |
| `t_step_save`            | Integer | Frequency to output data |
| `t_step_print`           | Integer | Frequency to print the current step number to standard output (default 1) |

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

- `dt` specifies the constant time step size that is used in simulation.
The value of `dt` needs to be sufficiently small such that the Courant-Friedrichs-Lewy (CFL) condition is satisfied.

- `t_step_start` and `t_step_end` define the time steps at which simulation starts and ends, respectively.
`t_step_save` is the time step interval for data output during simulation.
To newly start the simulation, set `t_step_start = 0`.
To restart simulation from $k$-th time step, set `t_step_start = k`, do not run `pre_process`, and run `simulation` directly (`./mfc.sh run [...] -t simulation`).
Ensure the data for the $k$-th time step is stored in the `restart_data/` directory within the case repository.

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
| `hcid` *             | Integer | N/A                   | Hard coded patch id                                          |
| `model%%filepath`     | String  | Not Supported         | Path to an STL or OBJ file (not all OBJs are supported).     |
| `model%%scale(i)`     | Real    | Not Supported         | Model's (applied) scaling factor for component $i$.          |
| `model%%rotate(i)`    | Real    | Not Supported         | Model's (applied) angle of rotation about axis $i$.          |
| `model%%translate(i)` | Real    | Not Supported         | Model's $i$-th component of (applied) translation.           |
| `model%%spc`          | Integer | Not Supported         | Number of samples per cell when discretizing the model into the grid. |
| `model%%threshold`    | Real    | Not Supported         | Ray fraction inside the model patch above which the fraction is set to one.|

*: These parameters should be prepended with `patch_icpp(j)%` where $j$ is the patch index. 

The Table lists the patch parameters.
The parameters define the geometries and physical parameters of fluid components (patch) in the domain at initial condition.
Note that the domain must be fully filled with patche(s).
The code outputs error messages when an empty region is left in the domain.

#### Analytical Definition of Primitive Variables

Some parameters, as described above, can be defined by analytical functions in the input file. For example, one can define the following patch:

```shell
'patch_icpp(2)%geometry'    : 15,
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
For example, to add a 2D Hardcoded patch with an id of 200, one would add the following to `src/pre_process/include/2dHardcodedIC.fpp`

```f90
    case(200)
        ! Primitive variables assignment
```

and use `patch_icpp(i)%%geometry = 7` and `patch_icpp(i)%%hcid = 200` in the input file.
Additional variables can be declared in `Hardcoded1[2,3]DVariables` and used in `hardcoded1[2,3]D`.
As a convention, any hard coded patches that are part of the MFC master branch should be identified as 1[2,3]xx where the first digit indites the number of dimensions.

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

- `model%%scale`, `model%%rotate` and `model%%translate` define how the model should be transformed to domain-space by first scaling by `model%%scale`, then rotating about the Z, X, and Y axes (using `model%%rotate`), and finally translating by `model%%translate`.

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

These parameters should be prepended with `patch_ib(j)%` where $j$ is the patch index. 

#### Parameter Descriptions

- `geometry` defines the type of geometry of a patch with an integer number.
Definitions for currently implemented patch types are list in table [Immersed Boundary Patch Type](#immersed-boundary-patch-types)

- `x[y,z]_centroid` is the centroid location of the patch in the x[y,z]-direction

- `length_x[y,z]` is the length of the patch in the x[y,z]-direction.

- `radius` is the radius to be used for circular patches.

- `theta` allows for the angle of attach of airfoil patches to be changed.

- `c`, `t`, `p`, and `m` specify the parameters for a NACA airfoil.
`m` is the maximum camber, `p` is the location of maximum camber, `c` is the coord length, and `t` is the thickness.
Additional details on this specification can be found in [The Naca Airfoil Series](https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf)

- `slip` applies a slip boundary to the surface of the patch if true and a no-slip boundary condition to the surface if false.

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

Fluid material's parameters. All parameters except for sigma should be prepended with `fluid_pp(i)` where $i$ is the fluid index.

*: Parameters that work only with `model_eqns = 2`.

**: Parameters that work only with `model_eqns = 3`.

The table lists the fluid material's parameters.
The parameters define material's property of compressible fluids that are used in simulation.

- `fluid_pp(i)%%gamma` and `fluid_pp(i)%%pi_inf` define $\Gamma$ and $\Pi$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

- `fluid_pp(i)%%Re(1)` and `fluid_pp(i)%%Re(2)` define the shear and volume viscosities of $i$-th fluid, respectively.

When these parameters are undefined, fluids are treated as inviscid.
Details of implementation of viscosity in MFC can be found in [Coralic (2015)](references.md#Coralic15).

- `fluid_pp(i)%%cv`, `fluid_pp(i)%%qv`, and `fluid_pp(i)%%qvp` define $c_v$, $q$, and $q'$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

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
| `weno_order`	         | Integer | WENO order [1,3,5] |
| `weno_eps`	         | Real    | WENO perturbation (avoid division by zero) |
| `mapped_weno`	         | Logical | WENO-M (WENO with mapping of nonlinear weights) |
| `wenoz`	             | Logical | WENO-Z |
| `teno`                 | Logical | TENO (Targeted ENO) |
| `teno_CT`              | Real    | TENO threshold for smoothness detection | 
| `null_weights`         | Logical | Null WENO weights at boundaries |
| `mp_weno`              | Logical | Monotonicity preserving WENO |
| `riemann_solver`       | Integer | Riemann solver algorithm: [1] HLL*; [2] HLLC; [3] Exact*	 |
| `low_Mach`             | Integer | Low Mach number correction for HLLC Riemann solver: [0] None; [1] Pressure (Chen et al. 2022); [2] Velocity (Thornber et al. 2008)	 |
| `avg_state`	         | Integer | Averaged state evaluation method: [1] Roe averagen*; [2] Arithmetic mean  |
| `wave_speeds`          | Integer | Wave-speed estimation: [1] Direct (Batten et al. 1997); [2] Pressure-velocity* (Toro 1999)	 |
| `weno_Re_flux`         | Logical | Compute velocity gradient using scaler divergence theorem	 |
| `weno_avg`          	 | Logical | Arithmetic mean of left and right, WENO-reconstructed, cell-boundary values |

- \* Options that work only with `model_eqns = 2`.
- † Options that work only with ``cyl_coord = 'F'``.
- ‡ Options that work only with `bc_[x,y,z]%[beg,end] = -15` and/or `bc_[x,y,z]%[beg,end] = -16`.

The table lists simulation algorithm parameters.
The parameters are used to specify options in algorithms that are used to integrate the governing equations of the multi-component flow based on the initial condition.
Models and assumptions that are used to formulate and discritize the governing equations are described in [Bryngelson et al. (2019)](references.md#Bryngelson19).
Details of the simulation algorithms and implementation of the WENO scheme can be found in [Coralic (2015)](references.md#Coralic15).

- `bc_[x,y,z]%[beg,end]` specifies the boundary conditions at the beginning and the end of domain boundaries in each coordinate direction by a negative integer from -1 through -12.
See table [Boundary Conditions](#boundary-conditions) for details.

- `bc_[x,y,z]%%vb[1,2,3]` specifies the velocity in the (x,1), (y,2), (z,3) direction applied to `bc_[x,y,z]%%beg` when using `bc_[x,y,z]%%beg = -16`.
Tangential velocities require viscosity, `weno_avg = T`, and `bc_[x,y,z]%%beg = -16` to work properly. Normal velocities require `bc_[x,y,z]%%end = -15` or `\bc_[x,y,z]%%end = -16` to work properly.

- `bc_[x,y,z]%%ve[1,2,3]` specifies the velocity in the (x,1), (y,2), (z,3) direction applied to `bc_[x,y,z]%%beg` when using `bc_[x,y,z]%%end = -16`.
Tangential velocities require viscosity, `weno_avg = T`, and `bc_[x,y,z]%%end = 16` to work properly. Normal velocities require `bc_[x,y,z]%%end = -15` or `\bc_[x,y,z]%%end = -16` to work properly.

- `model_eqns` specifies the choice of the multi-component model that is used to formulate the dynamics of the flow using integers from 1 through 3. 
`model_eqns = 1`, `2`, and `3` correspond to $\Gamma$-$\Pi_\infty$ model ([Johnsen, 2008](references.md#Johnsen08)), 5-equation model ([Allaire et al., 2002](references.md#Allaire02)), and 6-equation model ([Saurel et al., 2009](references.md#Saurel09)), respectively.
The difference of the two models is assessed by ([Schmidmayer et al., 2019](references.md#Schmidmayer19)).
Note that some code parameters are only compatible with 5-equation model.

- `alt_soundspeed` activates the source term in the advection equations for the volume fractions, $K\nabla\cdot \underline{u}$, that regularizes the speed of sound in the mixture region when the 5-equation model is used.
The effect and use of the source term are assessed by [Schmidmayer et al., 2019](references.md#Schmidmayer19).

- `adv_n` activates the direct computation of number density by the Riemann solver instead of computing number density from the void fraction in the method of classes.

- `mpp_lim` activates correction of solutions to avoid a negative void fraction of each component in each grid cell, such that $\alpha_i>\varepsilon$ is satisfied at each time step.

- `mixture_err` activates correction of solutions to avoid imaginary speed of sound at each grid cell.

- `time_stepper` specifies the order of the Runge-Kutta (RK) time integration scheme that is used for temporal integration in simulation, from the 1st to 5th order by corresponding integer. 
Note that `time_stepper = 3` specifies the total variation diminishing (TVD), third order RK scheme ([Gottlieb and Shu, 1998](references.md#Gottlieb98)).

- `adap_dt` activates the Strang operator splitting scheme which splits flux and source terms in time marching, and an adaptive time stepping strategy is implemented for the source term. It requires ``bubbles = 'T'``, ``polytropic = 'T'``, ``adv_n = 'T'`` and `time_stepper = 3`.

- `weno_order` specifies the order of WENO scheme that is used for spatial reconstruction of variables by an integer of 1, 3, and 5, that correspond to the 1st, 3rd, and 5th order, respectively.

- `weno_eps` specifies the lower bound of the WENO nonlinear weights.
Practically, `weno_eps` $<10^{-6}$ is used.

- `mapped_weno` activates the WENO-M scheme in place of the default WENO-JS scheme ([Henrick et al., 2005](references.md#Henrick05)). WENO-M a variant of the WENO scheme that remaps the nonlinear WENO-JS weights by assigning larger weights to non-smooth stencils, reducing dissipation compared to the default WENO-JS scheme, at the expense of higher computational cost. Only one of `mapped_weno`, `wenoz`, and `teno` can be activated.

- `wenoz` activates the WENO-Z scheme in place of the default WENO-JS scheme ([Borges et al., 2008](references.md#Borges08)). WENO-Z is a variant of the WENO scheme that further reduces the dissipation compared to the WENO-M scheme. It has similar computational cost to the WENO-JS scheme.

- `teno` activates the TENO scheme in place of the default WENO-JS scheme ([Fu et al., 2016](references.md#Fu16)). TENO is a variant of the ENO scheme that is the least dissipative, but could be less robust for extreme cases. It uses a threshold to identify smooth and non-smooth stencils, and applies optimal weights to the smooth stencils. Only available for `weno_order = 5`. Requires `teno_CT` to be set.

- `teno_CT` specifies the threshold for the TENO scheme. This dimensionless constant, also known as $C_T$, sets a threshold to identify smooth and non-smooth stencils. Larger values make the scheme more robust but also more dissipative. A recommended value for teno_CT is `1e-6`. When adjusting this parameter, it is recommended to try values like `1e-5` or `1e-7`. 

- `null_weights` activates nullification of the nonlinear WENO weights at the buffer regions outside the domain boundaries when the Riemann extrapolation boundary condition is specified (`bc_[x,y,z]%%beg[end]} = -4`).

- `mp_weno` activates monotonicity preservation in the WENO reconstruction (MPWENO) such that the values of reconstructed variables do not reside outside the range spanned by WENO stencil ([Balsara and Shu, 2000](references.md#Balsara00); [Suresh and Huynh, 1997](references.md#Suresh97)).

- `riemann_solver` specifies the choice of the Riemann solver that is used in simulation by an integer from 1 through 3.
`riemann_solver = 1`, `2`, and `3` correspond to HLL, HLLC, and Exact Riemann solver, respectively ([Toro, 2013](references.md#Toro13)).

- `low_Mach` specifies the choice of the low Mach number correction scheme for the HLLC Riemann solver. `low_Mach = 0` is default value and does not apply any correction scheme. `low_Mach = 1` and `2` apply the anti-dissipation pressure correction method ([Chen et al., 2022](references.md#Chen22)) and the improved velocity reconstruction method ([Thornber et al., 2008](references.md#Thornber08)). This feature requires `riemann_solver = 2` and `model_eqns = 2`.

- `avg_state` specifies the choice of the method to compute averaged variables at the cell-boundaries from the left and the right states in the Riemann solver by an integer of 1 or 2.
`avg_state = 1` and `2` correspond to Roe- and arithmetic averages, respectively.

- `wave_speeds` specifies the choice of the method to compute the left, right, and middle wave speeds in the Riemann solver by an integer of 1 and 2.
`wave_speeds = 1` and `2` correspond to the direct method ([Batten et al., 1997](references.md#Batten97)), and indirect method that approximates the pressures and velocity ([Toro, 2013](references.md#Toro13)), respectively.

- `weno_Re_flux` activates the scaler divergence theorem in computing the velocity gradients using WENO-reconstructed cell boundary values.
If this option is false, velocity gradient is computed using finite difference scheme of order 2 which is independent of the WENO order.

- `weno_avg` it activates the arithmetic average of the left and right, WENO-reconstructed, cell-boundary values.
This option requires `weno_Re_flux` to be true because cell boundary values are only utilized when employing the scalar divergence method in the computation of velocity gradients.


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
| `fd_order`           | Integer | Order of finite differences for computing the vorticity and the numerical Schlieren function [1,2,4] |
| `schlieren_alpha(i)` | Real    | Intensity of the numerical Schlieren computed via `alpha(i)` |
| `probe_wrt`          | Logical | Write the flow chosen probes data files for each time step	|
| `num_probes`         | Integer | Number of probes	|
| `probe(i)%[x,y,z]`   | Real	 | Coordinates of probe $i$	|

The table lists formatted database output parameters. The parameters define variables that are outputted from simulation and file types and formats of data as well as options for post-processing.

- `format` specifies the choice of the file format of data file outputted by MFC by an integer of 1 and 2. `format = 1` and `2` correspond to Silo-HDF5 format and binary format, respectively.

- `precision` specifies the choice of the floating-point format of the data file outputted by MFC by an integer of 1 and 2. `precision = 1` and `2` correspond to single-precision and double-precision formats, respectively.

- `parallel_io` activates parallel input/output (I/O) of data files. It is highly recommended to activate this option in a parallel environment.
With parallel I/O, MFC inputs and outputs a single file throughout pre-process, simulation, and post-process, regardless of the number of processors used.
Parallel I/O enables the use of different number of processors in each of the processes (i.e., simulation data generated using 1000 processors can be post-processed using a single processor).

- `file_per_process` deactivates shared file MPI-IO and activates file per process MPI-IO.
The default behavior is to use a shared file.
File per process is useful when running on 10's of thousands of ranks.
If `file_per_process` is true, then pre_process, simulation, and post_process must be run with the same number of ranks.

- `cons_vars_wrt` and `prim_vars_wrt` activate output of conservative and primitive state variables into the database, respectively.

- `[variable's name]_wrt` activates output of the each specified variable into the database.

- `schlieren_alpha(i)` specifies the intensity of the numerical Schlieren of $i$-th component.

- `fd_order` specifies the order of the finite difference scheme that is used to compute the vorticity from the velocity field and the numerical schlieren from the density field by an integer of 1, 2, and 4.
`fd_order = 1`, `2`, and `4` correspond to the first, second, and fourth-order finite difference schemes, respectively.

- `probe_wrt` activates output of state variables at coordinates specified by `probe(i)%[x;y,z]`.


### 8. Acoustic Source {#acoustic-source}

| Parameter                             | Type    | Description |
| ---:                                  | :----:  | :--- |
| `acoustic_source`                     | Logical | Acoustic source module activation |
| `num_source`                          | Integer | Number of acoustic sources |
| `acoustic(i)%%support`                | Integer | Geometry of spatial support for the acoustic source |
| `acoustic(i)%%dipole`                 | Logical | Dipole source activation (optional; default = false for monopole) |
| `acoustic(i)%%loc(j)`                 | Real    | $j$-th coordinate of the point that defines the acoustic source location |
| `acoustic(i)%%pulse`                  | Integer | Acoustic wave form: [1] Sine [2] Gaussian [3] Square |
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
| `acoustic(i)%%element_on`             | Integer | Transducer array - Element number that is on (optional; default = 0 for all elements) |
| `acoustic(i)%%element_spacing_angle`  | Real    | 2D Transducer array - Spacing angle (in rad) between adjacent transducer elements |
| `acoustic(i)%%element_polygon_ratio`  | Real    | 3D Transducer array - Ratio of polygon side length to transducer element radius |
| `acoustic(i)%%rotate_angle`           | Real    | 3D Transducer array - Rotation angle of the transducer array (optional; default = 0) |

Details of the transducer acoustic source model can be found in [Maeda and Colonius (2017)](references.md#Maeda17).

- `acoustic_source` activates the acoustic source module.

- `num_source` defines the total number of source planes by an integer.

- `%%support` specifies the choice of the geometry of acoustic source distribution. See table [Acoustic Supports](#acoustic-supports) for details.

- `%%dipole` changes the default monopole (one-sided) source to a dipole source. It is only available for planar waves.

- `%%loc(j)` specifies the location of the acoustic source in the $j$-th coordinate direction. For planer support, the location defines midpoint of the source plane. For transducer arrays, the location defines the center of the transducer or transducer array (not the focal point; for 3D it's the tip of the spherical cap, for 2D it's the tip of the arc). 

- `%%pulse` specifies the acoustic wave form. `%%pulse = 1`, `2`, and `3` correspond to sinusoidal wave, Gaussian wave, and square wave, respectively.

- `%%npulse` specifies the number of cycles of the acoustic wave generated. Only applies to `%%pulse = 1 and 3` (sine and square waves), and must be an integer for non-planar waves.

- `%%mag` specifies the peak amplitude of the acoustic wave.

- `%%frequency` and `%%wavelength` specify the frequency and wavelength of the acoustic wave, respectively. These parameters are exclusive and exactly one of them must be specified for `%%pulse = 1` or `3` (sine or square waves). They are related by the speed of sound in the medium: `frequency = speed_of_sound / wavelength`.

- `%%gauss_sigma_time` and `%%gauss_sigma_dist` specify the time and spatial widths of the Gaussian pulse in terms of sigma, respectively. In particular, `%%gauss_sigma_time` is the standard deviation in the Gaussian equation. These parameters are exclusive and exactly one of them must be specified for `%%pulse = 2` (Gaussian wave). They are related by the speed of sound in the medium: `gauss_sigma_dist = speed_of_sound * gauss_sigma_time`.

- `%%delay` specifies the time delay of the acoustic wave. This parameter is optional for `%%pulse = 1` or `3` (sine or square waves) and defaults to 0. It must be specified for `%%pulse = 2` (Gaussian wave). It is important to note that setting the delay to 0 for a Gaussian pulse results in a half-Gaussian pulse, and delays that are too small may result in the pulse being cut off at the start of the simulation. `4*gauss_sigma_time` is a typical value for the delay of a Gaussian pulse.

- `%%dir` specifies the direction of acoustic wave propagation for planar waves. The direction is defined by the angle in degrees from the x-axis in the x-y plane. It applies to both 2D and 3D simulation of planar waves (support is infinite in z-direction for 3D).

- `%%length` specifies the spatial length of the 2D or 3D planar wave. It is the length of the source plane perpendicular to the direction of wave propagation.

- `%%height` specifies the spatial height of the planar wave. Since `%%dir` is in the x-y plane, the height is perpendicular to the direction of wave propagation.

- `%%foc_length` specifies the focal length of the transducer for transducer waves. It is the distance from the transducer to the focal point.

- `%%aperture` specifies the aperture of the transducer. It is the diameter of the projection of the transducer arc onto the y-axis (2D) or spherical cap onto the y-z plane (3D). To simulate a transducer enclosing half of the circle/sphere, set the aperture to double the focal length. For transducer array, it is the total aperture of the array.

- `%%num_elements` specifies the number of transducer elements in a transducer array. 

- `%%element_on` specifies the element number of the transducer array that is on. The element number starts from 1. If all elements are on, set `%%element_on` to 0.

- `%%element_spacing_angle` specifies the spacing angle between adjacent transducer in radian. The total aperture (`%%aperture`) is set, so each transducer element is smaller if `%%element_spacing_angle` is larger.

- `%%element_polygon_ratio` specifies the ratio of the polygon side length to the aperture diameter of each transducer element in a circular 3D transducer array. The polygon side length is calculated by using the total aperture (`%%aperture`) as the circumcicle diameter, and `%%num_elements` as the number of sides of the polygon. The ratio is used specify the aperture size of each transducer element in the array, as a ratio of the total aperture. 

- `%%rotate_angle` specifies the rotation angle of the 3D circular transducer array along the x-axis (principal axis). It is optional and defaults to 0.

### 9. Ensemble-Averaged Bubble Model

| Parameter      | Type    | Description                                    |
| ---:           | :----:  |          :---                                  |
| `bubbles` 		 | Logical	| Ensemble-averaged bubble modeling	 |
| `bubble_model` | Integer	| [1] Gilmore; [2] Keller--Miksis |
| `polytropic`   | Logical	| Polytropic gas compression |
| `thermal` 		 | Integer	| Thermal model: [1] Adiabatic; [2] Isothermal; [3] Transfer |
| `R0ref` 			 | Real		  | Reference bubble radius |
| `polydisperse`   | Logical	| Polydispersity in equilibrium bubble radius R0|
| `nb` 			     | Integer	| Number of bins: [1] Monodisperse; [$>1$] Polydisperse |
| `poly_sigma` 	       | Real 		|	Standard deviation for probability density function of polydisperse bubble populations |
| `R0_type` 	       | Integer 		|	Quadrature rule for probability density function of polydisperse bubble populations |
| `Ca` 			     | Real		  | Cavitation number |
| `Web` 			   | Real		  | Weber number |
| `Re_inv` 		   | Real		  | Inverse Reynolds number |
| `mu_l0` *	     | Real 		|	Liquid viscosity (only specify in liquid phase)  |
| `ss` *		     | Real 		|	Surface tension (only specify in liquid phase) |
| `pv` *		     | Real 		|	Vapor pressure (only specify in liquid phase) | 
| `gamma_v` † 	 | Real 	  |	Specific heat ratio |
| `M_v` †     	 | Real 		| Molecular weight |
| `mu_v` †	     | Real 		|	Viscosity |
| `k_v` †	       | Real 		|	Thermal conductivity |
| `qbmm` 	       | Logical 		|	Quadrature by  method of moments|
| `dist_type` 	       | Integer 		|	Joint probability density function for bubble radius and velocity (only when qbmm is true)|
| `sigR` 	       | Real 		|	Standard deviation for probability density function of bubble radius (only when qbmm is true) |
| `sigV` 	       | Real 		|	Standard deviation for probability density function of bubble velocity (only when qbmm is true) |
| `rhoRV`	       | Real 		|	Correlation coefficient for joint probability density function of bubble radius and velocity (only when qbmm is true) |

These options work only for gas-liquid two component flows.
Component indexes are required to be 1 for liquid and 2 for gas.

- \* These parameters should be pretended with patch index $1$ that is filled with liquid: `fluid_pp(1)%`.
- †  These parameters should be pretended with patch indexes that are respectively filled with liquid and gas: `fluid_pp(1)%` and `fluid_pp(2)%`.

This table lists the ensemble-averaged bubble model parameters.

- `bubbles` activates the ensemble-averaged bubble model.

- `bubble_model` specified a model for spherical bubble dynamics by an integer of 1 and 2.
`bubble_model = 1`, `2`, and `3` correspond to the Gilmore, Keller-Miksis, and Rayleigh-Plesset models.

- `polytropic` activates polytropic gas compression in the bubble.
When `polytropic` is set `False`, the gas compression is modeled as non-polytropic due to heat and mass transfer across the bubble wall with constant heat and mass transfer coefficients based on ([Preston et al., 2007](references.md#Preston07)).

- `polydisperse` activates polydispersity in the bubble model by means of a probability density function (PDF) of the equiliibrium bubble radius. 

- `thermal` specifies a model for heat transfer across the bubble interface by an integer from 1 through 3.
`thermal = 1`, `2`, and `3` correspond to no heat transfer (adiabatic gas compression), isothermal heat transfer, and heat transfer with a constant heat transfer coefficient based on [Preston et al., 2007](references.md#Preston07), respectively.

- `R0ref` specifies the reference bubble radius.

- `nb` specifies the number of discrete bins that define the probability density function (PDF) of the equilibrium bubble radius.

- `R0_type` specifies the quadrature rule for integrating the log-normal PDF of equilibrium bubble radius for polydisperse populations.
`R0_type = 1` corresponds to Simpson's rule. 

- `poly_sigma` specifies the standard deviation of the log-normal PDF of equilibrium bubble radius for polydisperse populations. 

- `Ca`, `Web`, and `Re_inv` respectively specify the Cavitation number, Weber number, and the inverse Reynolds number that characterize the offset of the gas pressure from the vapor pressure, surface tension, and liquid viscosity when the polytropic gas compression model is used.

- `mu_l0`, `ss`, and `pv`, `gamma_v`, `M_v`, `mu_v`, and `k_v` specify simulation parameters for the non-polytropic gas compression model.
`mu_l0`, `ss`, and `pv` correspond to the liquid viscosity, surface tension, and vapor pressure, respectively. 
`gamma_v`, `M_v`, `mu_v`, and `k_v` specify the specific heat ratio, molecular weight, viscosity, and thermal conductivity of a chosen component.
Implementation of the parameters into the model follow [Ando (2010)](references.md#Ando10).

- `qbmm` activates quadrature by method of moments, which assumes a PDF for bubble radius and velocity. 

- `dist_type` specifies the initial joint PDF of initial bubble radius and bubble velocity required in qbmm. `dist_type = 1`  and `2` correspond to binormal and lognormal-normal distributions respectively. 

- `sigR` specifies the standard deviation of the PDF of bubble radius required in qbmm.  

- `sigV` specifies the standard deviation of the PDF of bubble velocity required in qbmm.  

- `rhoRV` specifies the correlation coefficient of the joint PDF of bubble radius and bubble velocity required in qbmm.  

### 10. Velocity Field Setup

| Parameter              | Type    | Description |
| ---:                   | :----:  | :--- |
| `perturb_flow`         | Logical | Perturb the initlal velocity field by random noise |
| `perturb_flow_fluid`   | Integer | Fluid density whose flow is to be perturbed |
| `perturb_flow_mag`     | Real    | Set the magnitude of flow perturbations |
| `perturb_sph`          | Logical | Perturb the initial partial density by random noise |
| `perturb_sph_fluid`    | Integer | Fluid component whose partial density is to be perturbed |
| `mixlayer_vel_profile` | Logical | Set the mean streamwise velocity to hyperbolic tangent profile |
| `mixlayer_vel_coef`    | Real    | Coefficient for the hyperbolic tangent profile of a mixing layer |
| `mixlayer_perturb`     | Logical | Perturb the initial velocity field by instability waves |
| `mixlayer_domain`      | Real    | Domain size of a mixing layer for the linear stability analysis |

The table lists velocity field parameters.
The parameters are optionally used to define initial velocity profiles and perturbations.

- `perturb_flow` activates the perturbation of initial velocity by random noise.

- `perturb_flow_fluid` specifies the fluid component whose flow is to be perturbed.

- `perturb_flow` activates the perturbation of initial velocity by random noise.

- `perturb_sph` activates the perturbation of initial partial density by random noise. 

- `perturb_sph_fluid` specifies the fluid component whose partial density is to be perturbed.

- `mixlayer_vel_profile` activates setting of the mean streamwise velocity to hyperbolic tangent profile. This option works only for 2D and 3D cases.

- `mixlayer_vel_coef` is a parameter for the hyperbolic tangent profile of a mixing layer when `mixlayer_vel_profile = 'T'`. The mean streamwise velocity profile is given as:

$$ u = patch\_icpp(1)\%vel(1) * tanh(y\_cc * mixlayer\_vel\_profile) $$

- `mixlayer_perturb` activates the perturbation of initial velocity by instability waves obtained from linear stability analysis for a mixing layer with hyperbolic tangent mean streamwise velocity profile. This option only works for `n > 0`, `bc_y%[beg,end] = -6`, `num_fluids = 1`, `model_eqns = 2` and `mixlayer_vel_profile = 'T'`.

- `mixlayer_domain` defines the domain size to compute spatial eigenvalues of the linear instability analysis when `mixlayer_perturb = 'T'`. For example, the spatial eigenvalue in `x` direction in 2D problem will be $2 \pi \alpha / (mixlayer\_domain*patch\_icpp(1)\%length\_y)$ for $\alpha = 1$, $2$ and $4$.

### 11. Phase Change Model
| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `relax`                | Logical | Activates Phase Change model |
| `relax_model`          | Integer | Phase change model: [5] pT-equilibrium; [6] pTg-equilibrium |
| `palpha_eps`           | Real    | tolerance of the Newton Solver to activate pT-equilibrium  |
| `ptgalpha_eps`	     | Real    | tolerance of the Newton Solver to activate pTg-equilibrium |

- `relax` Activates the Phase Change model.

- `relax_model` Specifies the phase change model to be used: [5] enables pT-equilibrium, while [6] activates pTg-equilibrium (if criteria are met).

- `palpha_eps` Specifies the tolerance used for the Newton Solvers used in the pT-equilibrium model. 

- `ptgalpha_eps` Specifies the tolerance used for the Newton Solvers used in the pTg-equilibrium model. 

### 12. Artificial Mach Number
| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `pi_fac`               | Real    | Ratio of artificial and true `pi_\infty` values|

- `pi_fac` specifies the ratio of artificial and true `pi_\infty` values (`=` artificial `pi_\infty` / true `pi_\infty`).
This parameter enables the use of true `pi_\infty` in bubble dynamics models, when the `pi_\infty` given in the `case.py` file is an artificial value.

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

Positive accelerations are in the `x[y,z]` direction are in the positive `x[y,z]` direction by convention.

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
The entries labeled "Characteristic." are characteristic boundary conditions based on [Thompson (1987)](references.md#Thompson87) and [Thompson (1990)](references.md#Thompson90).

### Patch types

| #    | Name               | Dim.  | Smooth | Description |
| ---: | :----:             | :---: |  :---: | :--- |
| 1    | Line segment 	    | 1     | N      | Requires `x_centroid` and `length_x`. |
| 2    | Circle 		    | 2     | Y      | Requires `[x,y]_centroid` and `radius`. |
| 3    | Rectangle 	        | 2     | N      | Coordinate-aligned. Requires `[x,y]_centroid` and `length_[x,y]`. |
| 4    | Sweep line 		| 2     | Y      | Not coordinate aligned. Requires `[x,y]_centroid` and `normal(i)`. |
| 5    | Ellipse 		    | 2     | Y      | Requires `[x,y]_centroid` and `radii(i)`. |
| 6    | N/A 		        | 2     | N      | No longer exists. Empty. |
| 7    | 2D analytical 	    | 2     | N      | Assigns the primitive variables as analytical functions. |
| 8    | Sphere 		    | 3     | Y      | Requires `[x,y,z]_centroid` and `radius` |
| 9    | Cuboid 		    | 3     | N      | Coordinate-aligned. Requires `[x,y,z]_centroid` and `length_[x,y,z]`. |
| 10   | Cylinder 		    | 3     | Y      | Requires `[x,y,z]_centroid`, `radius`, and `length_[x,y,z]`. |
| 11   | Sweep plane 	    | 3     | Y      | Not coordinate-aligned. Requires `x[y,z]_centroid` and `normal(i)`. |
| 12   | Ellipsoid 		    | 3     | Y      | Requires `[x,y,z]_centroid` and `radii(i)`. |
| 13   | 3D analytical 	    | 3     | N      | Assigns the primitive variables as analytical functions |
| 14   | Spherical Harmonic | 3     | N      | Requires `[x,y,z]_centroid`, `radius`, `epsilon`, `beta` |   
| 15   | 1D analytical      | 1     | N      | Assigns the primitive variables as analytical functions  |
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
| 11   | Circular Transducer Array    | 3D        | #7 requirements, `%%element_polygon_ratio`, and `%%rotate_angle`(optional; default = 0) |

Details of the required parameters for each acoustic support type are listed in [Acoustic Source](#acoustic-source).
The acoustic support number (`#`) corresponds to the acoustic support type `Acoustic(i)%%support`, where $i$ is the acoustic source index.
For each `%%parameter`, prepend the parameter with `acoustic(i)%`.

Additional requirements for all acoustic support types:
- ``acoustic_source = 'T'`` must be used to activate the acoustic source module.

- `num_source` must be set to the total number of acoustic sources.

- `%%support` must be set to the acoustic support number listed in the table.

- `%%dipole` is only supported for planar sources.

- `%%npulse = 1 or 3` requires exactly one of `%%frequency` or `%%wavelength` to be set. It accepts `%%delay` as an optional parameter (default = 0).

- `%%npulse = 2` requires exactly one of `%%gauss_sigma_time` or `%%gauss_sigma_space` to be set. It requires `%%delay` to be set.

Description of the acoustic support types:
- `%%support = 1` specifies an infinite source plane that is normal to the $x$-axis and intersects with the axis at $x=$ `%%loc(1)` in 1D simulation. `%%dir > 0` specifies a rightward propagating wave, and `%%dir < 0` specifies a leftward propagating wave. `%%dir = 0` is not allowed.

- `%%support = 2` specifies a semi-infinite source plane in 2D simulation.
The midplane location is [`%%loc(1)`, `%%loc(2)`] and the normal vector is [$\mathrm{cos}$(`%%dir`), $\mathrm{sin}$(`%%dir`)]. The length of the source plane is `%%length`, and the plane is perpendicular to the direction of wave propagation (defined by `%%dir`).

- `%%support = 3` specifies a semi-infinite source plane in 3D simulation.
The midplane location is [`%%loc(1)`, `%%loc(2)`] and the normal vector is [$\mathrm{cos}$(`%%dir`), $\mathrm{sin}$(`%%dir`)]. The length of the source plane is `%%length`, and the plane is perpendicular to the direction of wave propagation (defined by `%%dir`). Note that the plane is infinite in the $z$-direction, so `%%loc(3)` is not required.

- `%%support = 5` specifies a circular transducer in 2D simulation. The transducer is centered at [`%%loc(1)`, `%%loc(2)`] with a focal length of `%%foc_length` and an aperture of `%%aperture`. The center location is not the focal point; it is the tip of the circular arc (intersection of the arc and the x-axis). The aperture is the length of the projection of the circular arc onto the y-axis. If a semi-circle is desired, set the aperture to double the focal length. Note that this is physically a cylindrical transducer, and due to the complexity of Green's function for 2D wave, no closed-form solution is available for the 2D circular transducer, and an approximate is used (see [Maeda and Colonius (2017)](references.md#Maeda17) for details). For the mass source term correction factor, the theoretical approximation factor of -0.5 in ($r_{foc}^{-0.5}$) is replaced by an empirically determined factor of -0.85.

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
