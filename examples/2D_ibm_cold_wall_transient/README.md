# 2D Cold Cylinder in Hot Gas: Heat Absorbed with the Wall Conduction Column

A cylinder of radius 50 um at 300 K (an immersed boundary with `thermal_bc = 1`) meets quiescent air at 1500 K and 1 atm, as a cold
particle meets hot post-shock gas. The gas starts in a near-step conduction layer, so for the first fraction of the heat-diffusion time
the wall layer is much thinner than a cell and the ghost-cell wall cannot carry its flux. `patch_ib(1)%wall_column = T` gives each
first-layer ghost a one-dimensional conduction column along the wall normal that does. `D/ib_heat.dat` (`ib_heat_wrt`) records the
heat the body takes every step.

`check_heat.py` integrates it and compares the heat absorbed, Q(t) per unit depth, with a body-fitted 1D reference (`ref1d.py`): conduction
only, at constant pressure, in the mass coordinate, with the Cantera conductivity and heat capacity of the same frozen air. The results,
the model limit and the references are on the @ref validation_wall_thermal page.

## Run

    ./mfc.sh run examples/2D_ibm_cold_wall_transient/case.py -t pre_process simulation -n 1 -- --cells_per_R 40 --half 4 --tend 0.1 --wall_column 1
    python3 examples/2D_ibm_cold_wall_transient/check_heat.py examples/2D_ibm_cold_wall_transient

`--wall_column 0` is the default ghost-cell wall; `--tend` is in units of R^2/alpha of the hot gas (7.4 us). Each run writes into the
example directory, so check it before the next one. Without arguments the case is the CI-sized run (12 cells per R, a box of
half width 3 R, to 0.03 R^2/alpha, about 100 steps); the suite runs it as `2D -> Chemistry -> IBM -> Cold Wall Transient -> Wall Column`
(the golden holds the conservative fields).

## Results

GNU 13, one rank, double precision, `master` at 0cdfdf02a. Error = Q_MFC / Q_ref - 1 at the time in R^2/alpha:

| cells per R | wall | 0.003 | 0.01 | 0.03 | 0.1 |
|---:|:---|---:|---:|---:|---:|
| 10 | ghost mirror | -51.9 % | -33.7 % | -19.7 % | -9.6 % |
| 20 | ghost mirror | -37.5 % | -23.1 % | -12.9 % | -6.1 % |
| 40 | ghost mirror | -20.8 % | -10.3 % | -4.2 % | -1.4 % |
| 10 | wall column | +12.7 % | +6.4 % | +2.5 % | +0.7 % |
| 20 | wall column | +7.1 % | +0.8 % | -1.0 % | -0.4 % |
| 40 | wall column | +6.0 % | -0.3 % | -1.3 % | +0.0 % |

The ghost-cell wall loses most of the heat in the first part of the transient and converges slowly; the column is within 1.3 % of the
reference from 0.01 R^2/alpha on at 20 and 40 cells per R. Its excess at the earliest time does not fall with refinement in this range
(the column is a low-Mach model; see the page for the likely cause).
