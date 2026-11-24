# 2D Inviscid Flow Around a Fixed Slip Cylinder: Mass Budget and Exact Vortex

An immersed-boundary cylinder with a slip wall (`patch_ib(1)%slip = T`) in inviscid gas. Two checks of the ghost-cell wall, selected
by `--mode`:

- `leak` (the default): uniform flow at Mach 0.5 in a fully periodic box. No gas can enter or leave, so the total gas mass must stay
  constant. A slip wall does not do that by itself: the slip ghost copies the tangential velocity, which crosses the oblique faces of
  the staircase, and the Riemann flux there exchanges mass with ghost cells that are rewritten every stage. `ib_mass_relax = 1` returns
  each step's intake to the fluid within two cells of the body.
- `vortex`: the steady isentropic vortex u_theta = M c0 R / r (stagnation state 1e5 Pa and 1.2 kg/m^3 at infinity) is an exact Euler
  solution tangent to the wall. The pressure and velocity errors in the cells next to the wall measure the ghost state: the curvature
  term of the ghost pressure, the image-point interpolation (`ib_ip_interp`) and the isentropic ghost density with a total-enthalpy
  tangential speed (`ib_slip_ccst`). The initial condition is an analytic patch, so this mode compiles its own `pre_process` and
  `simulation` and is not in the test suite.

`check_slip.py` reads the run and prints the numbers of the tables. The results, limits and the slip-wall cases that fail are on the
@ref validation_slip_walls page.

## Run

    ./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 -- --steps 0 --cells_per_R 16 --half 8 --mass_relax 1
    python3 examples/2D_ibm_slip_cylinder/check_slip.py examples/2D_ibm_slip_cylinder --case-args '--steps 0 --cells_per_R 16 --half 8 --mass_relax 1'

    ./mfc.sh run examples/2D_ibm_slip_cylinder/case.py -t pre_process simulation post_process -n 1 -- --mode vortex --mach 0.3 --cells_per_R 20 --ip_interp 2 --ccst
    python3 examples/2D_ibm_slip_cylinder/check_slip.py examples/2D_ibm_slip_cylinder --case-args '--mode vortex --mach 0.3 --cells_per_R 20 --ip_interp 2 --ccst'

Each run writes into the example directory, so check it before the next one. `--steps 0` runs one flow-through of the box in `leak`
mode; the default 100 fixed steps on 8 cells per radius is the CI-sized case. The suite runs it as
`2D -> IBM -> Slip Cylinder -> Mass Budget` (the golden holds the conservative fields; skipped in the single-precision lane).
On an NVIDIA A100 `--mass_relax 1` had no effect in this case (see the page); the CPU result is the one in the tables.

## Results

Gas mass of the fluid cells, change relative to the first save, one flow-through of a box of side 16 R (GNU 13, one rank, `master` at 0cdfdf02a):

| cells per R | ib_mass_relax | half a flow-through | one flow-through |
| ---: | :--- | ---: | ---: |
| 16 | 0 (what -1 gives for a fixed body) | +1.69e-3 | +3.61e-3 |
| 16 | 30 | +4.3e-5 | +5.8e-5 |
| 16 | 1 | -5.2e-14 | -1.2e-13 |
| 32 | 0 | +1.00e-3 | +2.08e-3 |
| 32 | 1 | -1.0e-13 | -2.3e-13 |

See the page for the exact-vortex tables (M = 0.3 and 0.6, 10 to 40 cells per R, four option sets). The first-band pressure error of the M = 0.3 vortex at 10 / 20 / 40 cells per R is 3.28e-2 / 1.58e-2 / 6.38e-3 with the defaults and 4.51e-3 / 1.22e-3 / 2.40e-4 with `ib_ip_interp = 2` and `ib_slip_ccst = T`.
