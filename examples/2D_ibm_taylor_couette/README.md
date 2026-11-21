# 2D Taylor-Couette Flow Between Two Immersed Walls

A circle of radius R = 1 turns at a prescribed rate (`moving_ibm = 1`) inside a fixed no-slip ring (an STL immersed boundary
with inner radius 2 R). The gas (rho = 1, sound speed 1, constant viscosity) starts from rest and is run to steady state at tip
Mach number 0.05 and Re = Omega R^2 / nu = 10. The steady incompressible solution is u_theta = A r + B / r, so the velocity
profile, the torque on each wall and the wall speed are known exactly. The example measures what `ib_wall_order = 1` (mirror
ghost velocity, default) and `ib_wall_order = 2` (two-probe quadratic) do at a curved, rotating wall.

`check_couette.py` reads the run directory (serial text output) and prints the errors:

    python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette [--tol X]

| cells per R | wall order | L2 / Omega R | azimuthal-mean part | inner torque error | outer torque error | wall speed error (fit) |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1 | 1.66e-2 | 5.3e-3 | -3.0e-3 | -1.5e-2 | -5.4e-2 |
| 10 | 2 | 1.23e-2 | 5.2e-3 | -1.7e-3 | -9.7e-3 | -3.9e-2 |
| 20 | 1 | 1.09e-2 | 8.3e-3 | -1.2e-2 | -1.9e-2 | -6.3e-2 |
| 20 | 2 | 6.5e-3 | 4.2e-3 | -3.7e-3 | -5.8e-3 | -3.8e-2 |
| 40 | 1 | 8.1e-3 | 8.5e-3 | -1.6e-2 | -1.9e-2 | -7.2e-2 |
| 40 | 2 | 4.7e-3 | 4.9e-3 | -6.6e-3 | -8.8e-3 | -3.4e-2 |

Order 2 lowers the velocity error by 26 to 42 % at every resolution and about halves the torque error from 20 cells per radius.
Neither order converges cleanly on current master, which evaluates a rotating body's ghost-cell wall velocity at the ghost point, not at the wall point (a fix is in preparation); the azimuthal-mean error stays near 4e-3 to 9e-3 and the fitted inner wall speed stays 3 to 7 %
low. The measurement, its limits and the rotating-wall variants are on the page @ref validation_wall_order_couette.

## Run

    ./mfc.sh run examples/2D_ibm_taylor_couette/case.py -t pre_process simulation -n 1 -- --cells_per_R 20 --ntau 12 --wall_order 2
    python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette

`case.py` defaults to a CI-sized run (12 cells per radius, 64 x 64 cells, 300 fixed steps, far from steady state). Options:
`--cells_per_R`, `--wall_order`, `--ip_interp` (1 inverse distance, 2 TSC), `--visc_compact` (compact viscous gradients next to a
body), `--mach`, `--re`, `--sign`, `--shift` (grid shift in cells), `--steps`, `--ntau` (run this many slowest viscous times
instead of `--steps`) and `--stl`. Each run writes into the example directory; remove `D/`, `restart_data/` and the `.inp` files
between runs. The output is serial text, so use one rank. `ring_stl.py` writes the ring STL (`ring.stl` is
`ring_stl.py ring.stl 2 5 512`). The suite registers the CI-sized run as `2D -> IBM -> Taylor-Couette -> ib_wall_order=2`
(a regression golden; the accuracy is what the checker measures).

## Limits

2D only, one Mach number and Reynolds number, one machine and compiler for the table (NVIDIA A100, NVHPC 25.11, double
precision). The wall shear read from the near-wall cells does not converge for either order. Mean density rises by up to 2 %
over the run because the immersed boundary does not conserve mass; this is part of the plateau and is not quantified here.
