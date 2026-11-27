# 2D Plane Poiseuille Flow Between Two Flat Immersed Walls

Two static rectangle immersed boundaries at y = 0 and y = H = 1 bound a channel periodic in x and driven by a body
acceleration. The gas (rho = 1, sound speed 1, constant viscosity) starts from rest and is run to steady state at peak Mach
number 0.05 and Re = u_max H / nu = 1. The exact solution is u = g y (H - y) / (2 nu), with flow rate g H^3 / (12 nu), so the
profile and the flow rate are known. The example measures `ib_wall_order = 1` (wall value) and `ib_wall_order = 2` at a flat wall,
on the cell face and 0.37 cell inside it.

`check_poiseuille.py` reads the run directory (serial text output):

    python3 examples/2D_ibm_plane_poiseuille/check_poiseuille.py examples/2D_ibm_plane_poiseuille [--tol X]

Wall on a cell face, 12 slowest viscous times (error of the rms velocity over u_max; flow-rate error):

| cells across gap | L2, order 1 | L2, order 2 | flow rate, order 1 | flow rate, order 2 |
|---:|---:|---:|---:|---:|
| 8 | 0.455 | 5.4e-2 | +0.682 | +8.7e-2 |
| 16 | 0.216 | 1.8e-2 | +0.322 | +2.8e-2 |
| 32 | 0.113 | 6.7e-3 | +0.159 | +9.7e-3 |

With the wall 0.37 cell inside a face: 8 cells 0.474 and 7.2e-2, 16 cells 0.223 and 2.0e-2 (orders 1 and 2). Order 1 acts
as if the no-slip sat about 0.75 cell outside the wall; order 2 removes most of it. The wall force from the face fluxes is
a body-force balance and agrees within 1.2e-3 on a cell face for both orders; with the wall inside a cell the order 2 force splits
unequally between the two walls (+-7 % at 8 cells, +-4 % at 16). The page @ref validation_wall_order_couette has the full tables.

## Run

    ./mfc.sh run examples/2D_ibm_plane_poiseuille/case.py -t pre_process simulation -n 1 -- --cells 16 --ntau 12 --wall_order 2
    python3 examples/2D_ibm_plane_poiseuille/check_poiseuille.py examples/2D_ibm_plane_poiseuille

`case.py` defaults to a CI-sized run (8 cells across the gap, 32 x 28 cells, 300 fixed steps, far from steady state). Options:
`--cells`, `--wall_order`, `--mach`, `--re`, `--offset` (wall offset from a cell face, in cells), `--nx`, `--steps`, `--ntau`. Each
run writes into the example directory; remove `D/`, `restart_data/` and the `.inp` files between runs. Use one rank (serial text
output). The suite registers the CI-sized run as `2D -> IBM -> Plane Poiseuille -> ib_wall_order=2` (a regression golden).

## Limits

Static walls, one Mach number and Reynolds number, 2D, one machine and compiler (NVIDIA A100, NVHPC 25.11, double precision). The
shear read from the three cells nearest a wall does not converge for either order. No high-Reynolds-number boundary layer and no
moving flat wall was measured.
