# Case Files

Example Python case files, also referred to as *input files*, can be found in the [examples/](examples/) directory. They print a Python dictionary containing input parameters for MFC. Their contents, and a guide to filling them out, are documented
in the user manual. A commented, tutorial script
can also be found in [examples/3d_sphbubcollapse/](examples/3D_sphbubcollapse/case.py).

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

Thus, you can run your case file with Python to view the computed case dictionary
that will be processed by MFC when you run:

```console
$ python3 my_case_file.py
```

This is particularly useful when computations are done in Python to generate the case.

## (Optional) Accepting command line arguments

Input files can accept **positional** command line arguments, forwarded by `mfc.sh run`.
Consider this example from the 3D_weak_scaling case:

```python
import argparse

parser = argparse.ArgumentParser(
    prog="3D_weak_scaling",
    description="This MFC case was created for the purposes of weak scaling.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("dict", type=str, metavar="DICT", help=argparse.SUPPRESS)
parser.add_argument("gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory")

# Your parsed arguments are here
ARGS = vars(parser.parse_args())
```

The first argument is always a JSON string representing `mfc.sh run`'s internal
state. It contains all the runtime information you might want from the build/run system.
We hide it from the help menu with `help=argparse.SUPPRESS` since it is not meant
to be passed in by users. You can add as many additional positional arguments as
you may need.

To run such a case, use the following format:

```console
./mfc.sh run <path/to/case.py> <positional arguments> <regular mfc.sh run arguments>
```

For example, to run the 3D_weak_scaling case with `gbpp=2`:

```console
./mfc.sh run examples/3D_weak_scaling/case.py 2 -t pre_process -j 8
```

## Parameters

There are multiple sets of parameters that must be specified in the python input file:
1. [Runtime Parameters](#1-runtime)
2. [Computational Domain Parameters](#2-computational-domain)
3. [Patch Parameters](#3-patches)
4. [Fluid Material's Parameters](#4-fluid-materials)
5. [Simulation Algorithm Parameters](#5-simulation-algorithm)
6. [Formatted Database and Structure Parameters](#6-formatted-output)
7. [(Optional) Acoustic Source Parameters](#7-acoustic-source)
8. [(Optional) Ensemble-Averaged Bubble Model Parameters](#8-ensemble-averaged-bubble-model)
9. [(Optional) Velocity Field Setup Parameters](#9-velocity-field-setup)

Items 7, 8, and 9 are optional sets of parameters that activate the acoustic source model, ensemble-averaged bubble model, and initial velocity field setup, respectively.
Definition of the parameters is described in the following subsections.

### 1. Runtime

| Parameter        | Type           | Description                      |
| ---:             |    :----:      |          :---                    |
| `run_time_info`  | Logical        | Output run-time information      |

- `run_time_info` generates a text file that includes run-time information including the CFL number(s) at each time-step.

### 2. Computational Domain

| Parameter                | Type    | Description                      |
| ---:                     | :----:  |          :---                    |
| `x[y,z]_domain%beg[end]` | Real    | Beginning [ending] of the $x$[y,z]-direction domain    |
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

The parameters define the boundaries of the spatial and temporal domains, and their discritization that are used in simulation.

- `[x,y,z]_domain%[beg,end]` define the spatial domain in $x$, $y$, and $z$ Cartesian coordinates:

$$ x \in \left[ x \\_ domain \\% beg, x \\_ domain \\% end \right], y \in \left[ y \\_ domain \\% beg, y \\_ domain \\% end \right], z \in \left[ z \\_ domain \\% beg, z \\_ domain \\% end \right] $$

- $m$, $n$, and $p$ define the number of finite volume cells that uniformly discritize the domain along the $x$, $y$, and $z$ axes, respectively.
Note that the actual number of cells in each coordinate axis is given as $[m,n,p]+1$.
For example, $m=n=p=499$ discretizes the domain into $500^3$ cells. 
When the simulation is 2D/axi-symmetric or 1D, it requires that $p=0$ or $p=n=0$, respectively.

- `stretch_[x,y,z]` activates grid stretching in the $[x,y,z]$ directions.
The grid is gradually stretched such that the domain boundaries are pushed away from the origin along a specified axis.

- `a_[x,y,z]`, `[x,y,z]_a`, and `[x,y,z]_b` are parameters that define the grid stretching function. When grid stretching along the $x$ axis is considered, the stretching function is given as:

$$ x_{cb,stretch} = x_{cb} + \frac{x_{cb}}{a_x} \Bigg[ \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_a)}{L} \right) \right] + \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_{cb}-x_b)}{L} \right) \right] -2 \mathrm{log}\left[\mathrm{cosh} \left( \frac{a_x(x_b-x_a)}{2L} \right) \right]  \Bigg] $$

where `x_cb` and `x_[cb,stretch]` are the coordinates of a cell boundary at the original and stretched domains, respectively. `L` is the domain length along the `x` axis: `L`=`x_domain%end`-`x_domain%beg`. Crudely speaking, `x_a` and `x_b` define the coordinates at which the grid begins to get stretched in the negative and positive directions along the $x$ axis, respectively. $a_x$ defines the smoothness of the stretching. Stretching along the $y$ and $z$ axes follows the same logistics. Optimal choice of the parameters for grid stretching is case-dependent and left to the user. `loops_x[y,z]` defines the number of times
the grid stretching funciton is applied and has a default value of one.

- `cyl_coord` activates cylindrical coordinates. The domain is defined in $x$-$y$-$z$ cylindrical coordinates, instead of Cartesian coordinates. Domain discritization is accordingly conducted along the axes of cylindrical coordinates. Wnen $p=0$, the domain is defined on $x$-$y$ axi-symmetric coordinates. In both Coordinates, mesh stretching can be defined along the $x$- and $y$-axes. MPI topology is automatically optimized to maximize the parallel efficiency for given choice of coordinate systems.

- `dt` specifies the constant time step size that is used in simulation. The value of `dt` needs to be sufficiently small such that the Courant-Friedrichs-Lewy (CFL) condition is satisfied.

- `t_step_start` and `t_step_end` define the time steps at which simulation starts and ends, respectively. `t_step_save` is the time step interval for data output during simulation. To newly start simulation, set `t_step_start`=0. To restart simulation from $k$-th time step, set `t_step_start`=k.

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
| `model%filepath`     | String  | Not Supported         | Path to an STL or OBJ file (not all OBJs are supported).     |
| `model%scale(i)`     | Real    | Not Supported         | Model's (applied) scaling factor for component $i$.          |
| `model%rotate(i)`    | Real    | Not Supported         | Model's (applied) angle of rotation about axis $i$.          |
| `model%translate(i)` | Real    | Not Supported         | Model's $i$-th component of (applied) translation.           |
| `model%spc`          | Integer | Not Supported         | Number of samples per cell when discretizing the model into the grid. |

*: These parameters should be prepended with `patch_icpp(j)%` where $j$ is the patch index. 

The Table lists the patch parameters. The parameters define the geometries and physical parameters of fluid components (patch) in the domain at initial condition. Note that the domain must be fully filled with patche(s). The code outputs error messages when an empty region is left in the domain.

#### Analytical Definition of Primitive Variables

Some parameters, as described above, can be defined by analytical functions in the input file. For example, one can define the following patch:

```json
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

Some patch configurations are not adequatley handeled with the above analytic variable definitions. In this case, a hard coded patch can be used. Hard coded patches can be added by adding additional hard coded patch identifers to `src/pre_process/include/1[2,3]dHardcodedIC.fpp`. For example, to add a 2D Hardcoded patch with an id of 200, one would add the following to `src/pre_process/include/2dHardcodedIC.fpp`

```f90
    case(200)
        ! Primitive variables assignment
```

and use `patch_icpp(i)%geometry = 7` and `patch_icpp(i)%hcid = 200` in the input file. Additional variables can be declared in `Hardcoded1[2,3]DVariables` and used in `hardcoded1[2,3]D`. As a convention, any hard coded patches that are part of the MFC master branch should be identied as 1[2,3]xx where the first digit indiates the number of dimensions.

#### Parameter Descriptions

- `num_patches` defines the total number of patches defined in the domain. The number has to be a positive integer.

- `num_fluids` defines the total number of fluids defined in each of the patches. The number has to be a positive integer.

- `patch_icpp(j)%geometry` defines the type of geometry of $j$-th patch by using an integer from 1 to 13. Definition of the patch type for each integer is listed in table [Patch Types](#patch-types).

- `[x,y,z]_centroid`, `length_[x,y,z]`, and/or `radius` are used to uniquely define the geometry of the patch with given type. Requisite combinations of the parameters for each type can be found in is listed in table [Patch types](#patch-types).

- `patch_icpp(j)%alter_patch(i)` activates alternation of `patch(i)` with `patch(j)`.
For instance, in a 2D simulation, when a cylindrical `patch(2)` is immersed in a rectangular `patch(1)`:
  - `patch_icpp(1)%geometry`=3
  - `patch_icpp(2)%geometry`=2
  - `patch_icpp(2)%alter_patch(1)`=TRUE

- `smoothen` activates smoothening of the boundary of the patch that alters the existing patch.
When smoothening occurs, fluids of the two patches are mixed in the region of the boundary.
For instance, in the aforementioned case of the cylindrical patch immersed in the rectangular patch, smoothening occurs when `patch_icpp(2)smoothen`=TRUE. `smooth_coeff` controls the thickness of the region of smoothening (sharpness of the mixture region). The default value of `smooth_coeff` is unity. The region of smoothening is thickened with decreasing the value.
Optimal choice of the value of `smooth_coeff` is case-dependent and left to the user.

- `patch_icpp(j)alpha(i)`, `patch_icpp(j)alpha_rho(i)`, `patch_icpp(j)pres`, and `texttt{patch_icpp(j)vel(i)` define for $j$-th patch the void fraction of `fluid(i)`, partial density of `fluid(i)`, the pressure, and the velocity in the $i$-th coordinate direction. These physical parameters must be consistent with fluid material's parameters defined in the next subsection.
See also `adv_alphan` in table [Simulation Algorithm Parameters](#5-simulation-algorithm).

- 'model%scale', 'model%rotate` and `model%translate` define how the model should be transformed to domain-space by first
scaling by `model%scale`, then rotating about the Z, X, and Y axes (using `model%rotate`), and finally translating by `model%translate`.

### 4. Fluid Material’s

| Parameter | Type   | Description                                    |
| ---:      | :----: |          :---                                  |
| `gamma`   | Real   | Stiffened-gas parameter $\Gamma$ of fluid.     |
| `pi_inf`  | Real   | Stiffened-gas parameter $\Pi_\infty$ of fluid. |
| `Re(1)` * | Real   | Shear viscosity of fluid.                      |
| `Re(2)` * | Real   | Volume viscosity of fluid.                     |

Fluid material's parameters. All parameters should be prepended with `fluid_pp(i)` where $i$ is the fluid index.

*: Parameters that work only with `model_eqns`=2.

The table lists the fluid material's parameters.
The parameters define material's property of compressible fluids that are used in simulation.

- `fluid_pp(i)%gamma` and `fluid_pp(i)%pi_inf` define $\Gamma$ and $\Pi$ as parameters of $i$-th fluid that are used in stiffened gas equation of state.

- `fluid_pp(i)%Re(1)` and `fluid_pp(i)%Re(2)` define the shear and volume viscosities of $i$-th fluid, respectively.
When these parameters are undefined, fluids are treated as inviscid.
Details of implementation of viscosity in MFC can be found in [Coralic (2015)](references.md#Coralic15).

### 5. Simulation Algorithm

| Parameter              | Type    | Description                                    |
| ---:                   | :----:  |          :---                                  |
| `bc_[x,y,z]\%beg[end]` | Integer | Beginning [ending] boundary condition in the $[x,y,z]$-direction (negative integer, see table [Boundary Conditions](#boundary-conditions)) |
| `model_eqns`           | Integer | Multicomponent model: [1] $\Gamma/\Pi_\infty$; [2] 5-equation; [3] 6-equation\\%;%[4] 4-equation |
| `alt_soundspeed` *     | Logical | Alternate sound speed and $K \nabla \cdot u$ for 5-equation model |
| `adv_alphan`	         | Logical | Equations for all $N$ volume fractions (instead of $N-1$) |
| `mpp_lim`	             | Logical | Mixture physical parameters limits |
| `mixture_err`          | Logical | Mixture properties correction |
| `time_stepper`         | Integer | Runge--Kutta order [1--3] |
| `weno_order`	         | Integer | WENO order [1,3,5] |
| `weno_eps`	           | Real    | WENO perturbation (avoid division by zero) |
| `mapped_weno`	         | Logical | WENO with mapping of nonlinear weights |
| `null_weights`         | Logical | Null WENO weights at boundaries |
| `mp_weno`              | Logical | Monotonicity preserving WENO |
| `riemann_solver`       | Integer | Riemann solver algorithm: [1] HLL*; [2] HLLC; [3] Exact*	 |
| `avg_state`	           | Integer | Averaged state evaluation method: [1] Roe averagen*; [2] Arithmetic mean  |
| `wave_speeds`          | Integer | Wave-speed estimation: [1] Direct (Batten et al. 1997); [2] Pressure-velocity* (Toro 1999)	 |
| `weno_Re_flux`          | Logical | Compute velocity gradient using scaler divergence theorem	 |
| `weno_avg`          	 | Logical | Arithmetic mean of left and right, WENO-reconstructed, cell-boundary values |

- \* Options that work only with `model_eqns` $=2$.
- † Options that work only with `cyl_coord` $=$ `False`.

The table lists simulation algorithm parameters.
The parameters are used to specify options in algorithms that are used to integrate the governing equations of the multi-component flow based on the initial condition.
Models and assumptions that are used to formulate and discritize the governing equations are described in [Bryngelson et al. (2019)](references.md#Bryngelson19).
Details of the simulation algorithms and implementation of the WENO scheme can be found in [Coralic (2015)](references.md#Coralic15).

- `bc_[x,y,z]%[beg,end]` specifies the boundary conditions at the beginning and the end of domain boundaries in each coordinate direction by a negative integer from -1 through -12. See table [Boundary Conditions](#boundary-conditions) for details.

- `model_eqns` specifies the choice of the multi-component model that is used to formulate the dynamics of the flow using integers from 1 through 3. 
`model_eqns` $=$ 1, 2, and 3 correspond to $\Gamma$-$\Pi_\infty$ model ([Johnsen, 2008](references.md#Johnsen08)), 5-equation model ([Allaire et al., 2002](references.md#Allaire02)), and 6-equation model ([Saurel et al., 2009](references.md#Saurel09)), respectively.
The difference of the two models is assessed by ([Schmidmayer et al., 2019](references.md#Schmidmayer19)).
Note that some code parameters are only compatible with 5-equation model.

- `alt_soundspeed` activates the source term in the advection equations for the volume fractions, $K\nabla\cdot \underline{u}$, that regularizes the speed of sound in the mixture region when the 5-equation model is used. The effect and use of the source term are assessed by [Schmidmayer et al., 2019](references.md#Schmidmayer19).

- `adv_alphan` activates the advection equations of all the components of fluid. If this parameter is set false, the void fraction of $N$-th component is computed as the residual of the void fraction of the other components at each cell:

$$ \alpha_N=1-\sum^{N-1}_{i=1} \alpha_i $$

where $\alpha_i$ is the void fraction of $i$-th component. When a single-component flow is simulated, it requires that `adv_alphan` $=$ `True`.

- `mpp_lim` activates correction of solutions to avoid a negative void fraction of each component in each grid cell, such that $\alpha_i>\varepsilon$ is satisfied at each time step.

- `mixture_err` activates correction of solutions to avoid imaginary speed of sound at each grid cell.

- `time_stepper` specifies the order of the Runge-Kutta (RK) time integration scheme that is used for temporal integration in simulation, from the 1st to 5th order by corresponding integer. 
Note that `time_stepper` $=$ 3 specifies the total variation diminishing (TVD), third order RK scheme ([Gottlieb and Shu, 1998](references.md#Gottlieb98)).

- `weno_order` specifies the order of WENO scheme that is used for spatial reconstruction of variables by an integer of 1, 3, and 5, that correspond to the 1st, 3rd, and 5th order, respectively.

- `weno_eps` specifies the lower bound of the WENO nonlinear weights. Practically, `weno_eps` $<10^{-6}$ is used.

- `mapped_weno` activates mapping of the nonlinear WENO weights to the more accurate nonlinear weights in order to reinstate the optimal order of accuracy of the reconstruction in the proximity of critical points ([Henrick et al., 2005](references.md#Henrick05)).

- `null_weights` activates nullification of the nonlinear WENO weights at the buffer regions outside the domain boundaries when the Riemann extrapolation boundary condition is specified (`bc_[x,y,z]\%beg[end]}` $=-4$).

- `mp_weno` activates monotonicity preservation in the WENO reconstruction (MPWENO) such that the values of reconstructed variables do not reside outside the range spanned by WENO stencil ([Balsara and Shu, 2000](references.md#Balsara00); [Suresh and Huynh, 1997](references.md#Suresh97)).

- `riemann_solver` specifies the choice of the Riemann solver that is used in simulation by an integer from 1 through 3. `riemann_solver` $=$ 1,2, and 3 correspond to HLL, HLLC, and Exact Riemann solver, respectively ([Toro, 2013](references.md#Toro13)).

- `avg_state` specifies the choice of the method to compute averaged variables at the cell-boundaries from the left and the right states in the Riemann solver by an integer of 1 or 2. `avg_state` $=$ 1 and 2 correspond to Roe- and arithmetic averages, respectively.

- `wave_speeds` specifies the choice of the method to compute the left, right, and middle wave speeds in the Riemann solver by an integer of 1 and 2.
`wave_speeds` $=$ 1 and 2 correspond to the direct method ([Batten et al., 1997](references.md#Batten97)), and indirect method that approximates the pressures and velocity ([Toro, 2013](references.md#Toro13)), respectively.

- `weno_Re_flux` activates the scaler divergence theorem in computing the velocity gradients using WENO-reconstructed cell boundary values. If this option is false, velocity gradient is computed using finite difference scheme of order 2 which is independent of the WENO order.

- `weno_avg` it activates the arithmetic average of the left and right, WENO-reconstructed, cell-boundary values. This option requires `weno_Re_flux` to be true because cell boundary values are only utilized when employing the scalar divergence method in the computation of velocity gradients.


### 6. Formatted Output

| Parameter            | Type    | Description                                    |
| ---:                 | :----:  |          :---                                  |
| `format`             | Integer | Output format. [1]: Silo-HDF5; [2] Binary	|
| `precision`          | Integer | [1] Single; [2] Double	 |
| `parallel_io`        | Logical | Parallel I/O	|
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

- `format` specifies the choice of the file format of data file outputted by MFC by an integer of 1 and 2. `format` $=$ 1 and 2 correspond to Silo-HDF5 format and binary format, respectively.

- `precision` specifies the choice of the floating-point format of the data file outputted by MFC by an integer of 1 and 2. `precision` $=$ 1 and 2 correspond to single-precision and double-precision formats, respectively.

- `parallel_io` activates parallel input/output (I/O) of data files. It is highly recommended to activate this option in a parallel environment.
With parallel I/O, MFC inputs and outputs a single file throughout pre-process, simulation, and post-process, regardless of the number of processors used.
Parallel I/O enables the use of different number of processors in each of the processes (i.e. simulation data generated using 1000 processors can be post-processed using a single processor).

- `cons_vars_wrt` and `prim_vars_wrt} activate output of conservative and primitive state variables into the database, respectively.

- `[variable's name]_wrt` activates output of the each specified variable into the database.

- `schlieren_alpha(i)` specifies the intensity of the numerical Schlieren of $i$-th component.

- `fd_order` specifies the order of finite difference scheme that is used to compute the vorticity from the velocity field and the numerical schlieren from the density field by an integer of 1, 2, and 4. `fd_order` $=$ 1, 2, and 4 correspond to the first, second, and fourth order finite difference schemes, respectively.

- `probe_wrt` activates output of state variables at coordinates specified by `probe(i)%[x;y,z]`.


### 7. Acoustic Source

| Parameter         | Type    | Description |
| ---:              | :----:  | :--- |
| `Monopole`        | Logical | Acoustic source |
| `num_mono`        | Integer | Number of acoustic sources |
| `Mono(i)%pulse`   | Integer | Acoustic wave form: [1] Sine [2] Gaussian [3] Square |
| `Mono(i)%npulse`  | Integer | Number of pulse cycles |
| `Mono(i)%support` | Integer | Type of the spatial support of the acoustic source : [1] 1D [2] Finite width (2D) [3] Support for finite line/patch [4] General support for 3D simulation in cartesian systems [5] Support along monopole acoustic transducer [6] Support for cylindrical coordinate system along axial-dir |
| `Mono(i)%loc(j)`  | Real    | $j$-th coordinate of the point that consists of $i$-th source plane |
| `Mono(i)%dir`     | Real    | Direction of acoustic propagation	|
| `Mono(i)%mag`     | Real    | Pulse magnitude	|
| `Mono(i)%length`  | Real    | Spatial pulse length |

The table lists acoustic source parameters. The parameters are optionally used to define a source plane in the domain that generates an acoustic wave that propagates in a specified direction normal to the source plane (one-way acoustic source). Details of the acoustic source model can be found in [Maeda and Colonius (2017)](references.md#Maeda17).

- `Monopole` activates the acoustic source.

- `num_mono` defines the total number of source planes by an integer.

- `Mono(i)%pulse` specifies the choice of the acoustic wave form generated from $i$-th source plane by an integer.
`Mono(i)%pulse` $=$ 1, 2, and 3 correspond to sinusoidal wave, Gaussian wave, and square wave, respectively.

- `Mono(i)%npulse` defines the number of cycles of the acoustic wave generated from $i$-th source plane by an integer.

- `Mono(i)%mag` defines the peak amplitude of the acoustic wave generated from $i$-th source plane with a given wave form.

- `Mono(i)%length` defines the characteristic wavelength of the acoustic wave generated from $i$-th source plane.

- `Mono(i)%support` specifies the choice of the geometry of acoustic source distribution of $i$-th source plane by an integer from 1 through 3:\\
`Mono(i)%support` $=1$ specifies an infinite source plane that is normal to the $x$-axis and intersects with the axis at $x=$ `Mono(i)%loc(1)` in 1-D simulation.\\
`Mono(i)%support` $=2$ specifies a semi-infinite source plane in 2-D simulation.
The $i$-th source plane is determined by the point at [`Mono(i)%loc(1)`, `Mono(i)%loc(2)`] and the normal vector [$\mathrm{cos}$(`Mono(i)%dir`), $\mathrm{sin}$(`Mono(i)%dir`)] that consists of this point. 
The source plane is defined in the finite region of the domain: $x\in[-\infty,\infty]$ and $y\in$[-`mymono_length`/2, `mymono_length`/2].\\
`Mono(i)%support` $=3$ specifies a semi-infinite source plane in 3-D simulation.
The $i$-th source plane is determined by the point at [`Mono(i)%loc(1)`, `Mono(i)%loc(2)`, `Mono(i)%loc(3)`] and the normal vector [$\mathrm{cos}$(`Mono(i)%dir`), $\mathrm{sin}$(`Mono(i)%dir`), 1] that consists of this point.
The source plane is defined in the finite region of the domain: $x\in[-\infty,\infty]$ and $y,z\in$[-`mymono_length`/2, `mymono_length`/2]. There are a few additional spatial support types available for special source types and coordinate systems tabulated in [Monopole supports](#monopole-supports).

### 8. Ensemble-Averaged Bubble Model

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

These options work only for gas-liquid two component flows. Component indexes are required to be 1 for liquid and 2 for gas.

- \* These parameters should be pretended with patch index $1$ that is filled with liquid: `fluid_pp(1)%`.
- †  These parameters should be pretended with patch indexes that are respectively filled with liquid and gas: `fluid_pp(1)%` and `fluid_pp(2)%`.

This table lists the ensemble-averaged bubble model parameters.

- `bubbles` activates the ensemble-averaged bubble model.

- `bubble_model` specified a model for spherical bubble dynamics by an integer of 1 and 2.
`bubble_model` $=$ 1 and 2 correspond to the Gilmore and the Keller-Miksis equations, respectively.

- `polytropic` activates polytropic gas compression in the bubble.
When `polytropic` is set `False`, the gas compression is modeled as non-polytropic due to heat and mass transfer across the bubble wall with constant heat and mass transfer coefficients based on ([Preston et al., 2007](references.md#Preston07)).

- `polydisperse` activates polydispersity in the bubble model by means of a probability density function (PDF) of the equiliibrium bubble radius. 

- `thermal` specifies a model for heat transfer across the bubble interface by an integer from 1 through 3.
`thermal` $=$ 1, 2, and 3 correspond to no heat transfer (adiabatic gas compression), isothermal heat transfer, and heat transfer with a constant heat transfer coefficient based on [Preston et al., 2007](references.md#Preston07), respectively.

- `R0ref` specifies the reference bubble radius.

- `nb` specifies the number of discrete bins that define the probability density function (PDF) of the equilibirum bubble radius.

- `R0_type` specifies the quadrature rule for integrating the log-normal PDF of equilibrium bubble radius for polydisperse populations. `R0_type` $=$ 1 corresponds to simpson's rule. 

- `poly_sigma` specifies the standard deviation of the log-normal PDF of equilibirium bubble radius for polydisperse populations. 

- `Ca`, `Web`, and `Re_inv` respectively specify the Cavitation number, Weber number, and the inverse Reynolds number that characterize the offset of the gas pressure from the vapor pressure, surface tension, and liquid viscosity when the polytropic gas compression model is used.

- `mu_l0`, `ss`, and `pv`, `gamma_v`, `M_v`, `mu_v`, and `k_v` specify simulation parameters for the non-polytropic gas compression model.
`mu_l0`, `ss`, and `pv` correspond to the liquid viscosity, surface tension, and vapor pressure, respectively. 
`gamma_v`, `M_v`, `mu_v`, and `k_v` specify the specific heat ratio, molecular weight, viscosity, and thermal conductivity of a chosen component.
Implementation of the parameterse into the model follow [Ando (2010)](references.md#Ando10).

- `qbmm` activates quadrature by method of moments, which assumes a PDF for bubble radius and velocity. 

- `dist_type` specifies the initial joint PDF of initial bubble radius and bubble velocity required in qbmm. `dist_type` $=$ 1  and 2 correspond to binormal and lognormal-normal distributions respectively. 

- `sigR` specifies the standard deviation of the PDF of bubble radius required in qbmm.  

- `sigV` specifies the standard deviation of the PDF of bubble velocity required in qbmm.  

- `rhoRV` specifies the correlation coefficient of the joint PDF of bubble radius and bubble velocity required in qbmm.  

### 9. Velocity Field Setup

| Parameter           | Type    | Description |
| ---:                | :----:  | :--- |
| `perturb_flow`      | Logical | Perturb the initlal velocity field by random noise |
| `perturb_sph`       | Logical | Perturb the initial partial density by random noise |
| `perturb_sph_fluid` | Integer | Fluid component whose partial density to be perturbed |
| `vel_profile`       | Logical | Set the mean streamwise velocity to hyperbolic tangent profile |
| `instability_wave`  | Logical | Perturb the initial velocity field by instability waves |

The table lists velocity field parameters. The parameters are optionally used to define initial velocity profiles and perturbations.

- `perturb_flow` activates the perturbation of initial velocity by random noise.

- `perturb_sph` activates the perturbation of intial partial density by random noise. 

- `perturb_sph_fluid` specifies the fluid component whose the partial density to be perturbed.

- `vel_profile` activates setting the mean streamwise velocity to hyperbolic tangent profile. This option works only for 2D and 3D cases.

- `instability_wave` activates the perturbation of initial velocity by instability waves obtained from linear stability analysis for a mixing layer with hyperbolic tangent mean streamwise velocity profile. This option only works for 2D and 3D cases, together with `vel_profile`=TRUE.


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

*: This boundary condition is only used for `bc_y%beg` when using cylindrical coordinates (`cyl_coord = 'T'` and 3d). For axisymmetric problems, use `bc_y%beg = -2` with `cyl_coord = 'T'` in 2D.

The boundary condition supported by the MFC are listed in table [Boundary Conditions](#boundary-conditions). Their number (`#`)
corresponds to the input value in `input.py` labeled `bc_[x,y,z]%[beg,end]` (see table [Simulation Algorithm Parameters](#5-simulation-algorithm)).
The entries labeled "Characteristic." are characteristic boundary conditions based on [Thompson (1987)](references.md#Thompson87) and [Thompson (1990)](references.md#Thompson90).

### Patch types

| #    | Name               | Dim.  | Smooth | Description |
| ---: | :----:             | :---: |  :---: | :--- |
| 1    | Line segment 	    | 1     | N      | Requires `x_centroid` and `x_length`. |
| 2    | Circle 		    | 2     | Y      | Requires `[x,y]_centroid` and `radius`. |
| 3    | Rectangle 	        | 2     | N      | Coordinate-aligned. Requires `[x,y]_centroid` and `[x,y]_length`. |
| 4    | Sweep line 		| 2     | Y      | Not coordinate aligned. Requires `[x,y]_centroid` and `normal(i)`. |
| 5    | Ellipse 		    | 2     | Y      | Requires `[x,y]_centroid` and `radii(i)`. |
| 6    | N/A 		    | 2     | N      | No longer exists. Empty. |
| 7    | 2D analytical 	    | 2     | N      | Assigns the primitive variables as analytical functions. |
| 8    | Sphere 		    | 3     | Y      | Requires `[x,y,z]_centroid` and `radius` |
| 9    | Cuboid 		    | 3     | N      | Coordinate-aligned. Requires `[x,y,z]_centroid` and `[x,y,z]_length`. |
| 10   | Cylinder 		    | 3     | Y      | Requires `[x,y,z]_centroid`, `radius`, and `[x,y,z]_length`. |
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
| 21   | Model              | 2 & 3 | Y      | Imports a Model (STL/OBJ). Requires `model%filepath`. |

The patch types supported by the MFC are listed in table [Patch Types](#patch-types). This includes
types exclusive to one-, two-, and three-dimensional problems. The patch type number (`#`)
corresponds to the input value in `input.py` labeled  `patch_icpp(j)%geometry` where
$j$ is the patch index. Each patch requires a different set of parameters, which are 
also listed in this table.

### Monopole supports

| #    | Description |
| ---- | ----        |
|    1 | 1D normal to x-axis |
|    2 | 2D semi-infinite source plane |
|    3 | 3D semi-infinite source plane along some lines |
|    4 | 3D semi-infinite source plane |
|    5 | Transducer |
|    6 | Cyl_coord along axial-dir |

The monopole support types available in MFC are listed in table [Monopole supports](#monopole-supports). This includes
types exclusive to one-, two-, and three-dimensional problems with special souce geometry like transducers as well as coordinate systems such as cylindrical coordinates. The monopole support number (`#`) corresponds to the input value in `input.py` labeled  `Mono(i)%support` where
$i$ is the monopole source index.

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

The above variables are used for all simualtions.

| 5-eqn | 6-eqn |
| ----  |  ---- |
| sub-grid bubble variables | N/A |
| hypoelastic variables     | N/A |

The above variables correspond to optional physics.
