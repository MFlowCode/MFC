@page validation_wall_order_couette Wall-Normal Ghost Velocity at Viscous Walls

# Wall-Normal Ghost Velocity at Viscous Walls

**Status:** partially validated (the gain of `ib_wall_order = 2` is measured on a flat and a curved wall; two limits are not explained).
**Code:** `master` at dd9778b26 (contains #2176 `ib_wall_order`, #2169 centripetal ghost pressure, #2165 TSC image points, #2212 compact viscous gradients), unmodified. NVIDIA A100, NVHPC 25.11, OpenACC, double precision, one rank.
**Last run:** 2026-11-21. **Examples:** `examples/2D_ibm_plane_poiseuille`, `examples/2D_ibm_taylor_couette`.

#2176 added `ib_wall_order = 2`, a two-probe quadratic wall-normal ghost velocity at no-slip walls, and reported no accuracy gain against an exact solution. This page measures it against two exact viscous solutions, for `ib_wall_order` 1 (the default mirror) and 2, and records what happens at a rotating curved wall.

## Setup

| Item | Channel (flat wall) | Annulus (curved wall) |
| :--- | :--- | :--- |
| Geometry | two static rectangle walls at y = 0 and y = H = 1, periodic in x | inner circle R = 1 (rotating, `moving_ibm = 1`), outer wall an STL ring at r = 2R (fixed); both no slip |
| Gas | rho = 1, sound speed 1, constant mu, 1 fluid | same |
| Drive, regime | body acceleration, peak Mach 0.05, Re = u_max H / nu = 1 | tip Mach 0.05, Re = Omega R^2 / nu = 10 |
| Run | from rest, 12 slowest viscous times, fixed dt at acoustic CFL 0.3 | same |
| Grids | 8, 16, 32 cells across the gap; wall on a cell face or 0.37 cell inside it (8, 16) | 10, 20, 40 cells per R, grid centred on the body |
| Numerics | WENO5, HLLC, SSP-RK3, `fd_order = 4`, `weno_Re_flux = T`, `ib_wall_order` 1 or 2 | same |

## Reference data

1. Channel: \f$u = g y (H-y)/(2\nu)\f$, \f$u_\text{max} = g H^2/(8\nu)\f$, flow rate \f$g H^3/(12\nu)\f$, wall shear \f$\rho g H/2\f$.
2. Annulus (incompressible): \f$u_\theta = A r + B/r\f$ with \f$A = -\Omega R^2/(r_2^2-R^2)\f$, \f$B = \Omega R^2 r_2^2/(r_2^2-R^2)\f$; torque on the inner body per unit depth \f$-4\pi\mu B\f$ (the outer wall receives \f$+4\pi\mu B\f$).
3. Compressibility (O(M^2)), viscous heating and the density drift of the immersed boundary are not in the references. The mean density rises over the run (about 2 % at 10 and 0.4 % at 40 cells per R in an earlier measurement, not repeated here).

## Quantities

- **Velocity error (L2):** rms over the fluid cells (centre inside the gas and more than half a cell from a wall in the annulus), divided by u_max (channel) or Omega R (annulus). In the annulus it is split into the azimuthal-mean part (error averaged in bins of one cell width in r) and the rest.
- **Flow rate (channel):** sum of u dy over the fluid cells against the exact value.
- **Torque (annulus):** from the face fluxes between each body and the fluid (`ib_state_wrt`, `ib_force_wrt`), mean over the last 10 % of the run, relative to the exact value.
- **Wall speed:** cubic fit (annulus) or quadratic fit through the three nearest cells (channel) extrapolated to the wall; exact value Omega R (annulus) or 0 (channel).
- **Wall force (channel):** from the face fluxes against g H Lx / 2 per wall. In steady state this is the body-force balance, so it cannot separate the two orders.

The scripts are `check_poiseuille.py` and `check_couette.py`; the text output labels each cell with its right face, so cell centres are the labels minus dx/2.

## Results

### Flat channel, Mach 0.05

Wall on a cell face (offset 0) and 0.37 cell inside it. Errors are signed (MFC minus exact over exact) for the flow rate.

| Cells across | Offset | L2, order 1 | L2, order 2 | Flow rate, order 1 | Flow rate, order 2 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 0 | 0.455 | 0.0545 | +0.682 | +0.0866 |
| 8 | 0.37 | 0.474 | 0.0724 | +0.679 | +0.0872 |
| 16 | 0 | 0.216 | 0.0185 | +0.322 | +0.0279 |
| 16 | 0.37 | 0.223 | 0.0204 | +0.320 | +0.0250 |
| 32 | 0 | 0.113 | 0.00668 | +0.159 | +0.00971 |

- Order 1 behaves as if the no-slip sat about 0.75 cell outside the wall (flow rate 1.7 times the exact one at 8 cells across). Order 2 reduces the velocity error by 8 to 17 times from 8 to 32 cells and converges at about 1.5 orders (L2 0.0545, 0.0185, 0.00668).
- Wall velocity from the quadratic fit, in units of u_max (exact 0): order 1 0.77, 0.42, 0.46 at 8, 16, 32 cells; order 2 0.115, 0.051, 0.020.
- The wall force agrees with g H / 2 within 1.2e-3 for both orders when the wall is on a cell face. With the wall inside a cell the order 2 force splits unequally between the two walls (+-7 % at 8 cells, +-4 % at 16), because the fluxes are taken at the staircase faces, not at the wall.

### Taylor-Couette annulus, Mach 0.05

| Cells per R | Order | L2 / Omega R | Azimuthal mean | Rest | Inner torque | Outer torque | Wall speed (fit) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 1 | 1.66e-2 | 5.3e-3 | 1.60e-2 | -3.0e-3 | -1.48e-2 | -5.4e-2 |
| 10 | 2 | 1.23e-2 | 5.2e-3 | 1.14e-2 | -1.7e-3 | -9.7e-3 | -3.9e-2 |
| 20 | 1 | 1.09e-2 | 8.3e-3 | 8.5e-3 | -1.19e-2 | -1.87e-2 | -6.3e-2 |
| 20 | 2 | 6.5e-3 | 4.2e-3 | 5.6e-3 | -3.7e-3 | -5.8e-3 | -3.8e-2 |
| 40 | 1 | 8.1e-3 | 8.5e-3 | 3.8e-3 | -1.63e-2 | -1.94e-2 | -7.2e-2 |
| 40 | 2 | 4.7e-3 | 4.9e-3 | 2.0e-3 | -6.6e-3 | -8.8e-3 | -3.4e-2 |

- The velocity error is 26 to 42 % lower with order 2 at every resolution. Both converge at only about 0.4 to 0.5 order between 20 and 40 cells per R (order 1 1.09e-2 to 8.1e-3, order 2 6.5e-3 to 4.7e-3).
- The azimuthal-mean part does not converge for either order (order 1: 5.3e-3, 8.3e-3, 8.5e-3; order 2: 5.2e-3, 4.2e-3, 4.9e-3). It is a deficit of angular velocity near the inner wall: the fitted wall speed is 5 to 7 % low for order 1 and 3 to 4 % low for order 2, with no trend with resolution.
- The torque from the face fluxes is within 2 % of the exact value for both orders, with order 2 about 2 to 3 times closer at 20 and 40 cells per R. It does not converge monotonically, and the inner and outer torque differ by up to 1.2e-2, so the run is not exactly steady or the angular momentum is not exactly conserved.
- The change of the azimuthal velocity over the last save interval is below 1.3e-5 of Omega R in these runs.

### Rotating inner wall with the optional wall treatments, Mach 0.01

Same annulus at tip Mach 0.01 and Re = 11, run for 4 slowest viscous times, `ib_wall_order = 1` unless noted. The centripetal ghost pressure of #2169 is always on for a spinning body; `ib_ip_interp = 2` (TSC) and `ib_visc_compact = T` (compact viscous gradients next to a body) are opt-in. Shear is the inner-wall shear \f$S = r\,d(u_\theta/r)/dr\f$ from local cubic fits at 96 wall points (3.5 to 5 cells), rms relative error; the fit applied to the exact profile gives 3.2e-3 (20 cells per R) and 9.2e-4 (30).

| Cells per R | Options | L2 / Omega R | Azimuthal mean | Wall speed (fit) | Inner torque | Shear rms |
| ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 20 | default | 1.34e-2 | 8.5e-3 | -3.5e-2 | +3.8e-2 | 0.47 |
| 20 | TSC | 1.34e-2 | 8.4e-3 | -3.4e-2 | +3.7e-2 | 0.47 |
| 20 | TSC, compact | 1.35e-2 | 9.9e-3 | -2.0e-2 | +3.6e-2 | 0.47 |
| 20 | TSC, compact, order 2 | 9.6e-3 | 7.8e-3 | -6.6e-3 | +2.9e-2 | 0.38 |
| 30 | default | 7.7e-3 | 4.3e-3 | -4.4e-2 | +1.5e-2 | 0.47 |
| 30 | TSC | 7.8e-3 | 4.3e-3 | -4.4e-2 | +1.4e-2 | 0.47 |
| 30 | TSC, compact | 6.5e-3 | 3.3e-3 | -1.8e-2 | +1.5e-2 | 0.35 |
| 30 | TSC, compact, order 2 | 4.1e-3 | 2.5e-3 | -6.4e-3 | +1.4e-2 | 0.28 |

- TSC alone changes nothing here. Compact viscous gradients cut the wall-speed deficit by about half, and adding order 2 brings it to 0.6 %; at 30 cells per R the velocity error is then 1.9 times smaller than the default (4.1e-3 against 7.7e-3).
- The inner-wall shear from local fits stays at 28 to 47 % for every option and does not converge with resolution. This differs from an earlier measurement of the same quantity (0.2 to 2 % with these options, post-processed from restart files) that this run does not reproduce. Starting from the exact field instead of from rest (20 cells per R, default and the last option) gives the same shear (0.47 and 0.38), so the start-up is not the cause. The cause of the difference is not identified (not verified); treat the shear numbers as unresolved and the wall speed, velocity error and torque as the measured results.
- The torque is within 4 % of the exact value for every option, with no trend between options.

## Assessment

- **Flat wall:** order 2 is a large gain in the velocity profile and the flow rate (8 to 17 times at 8 to 32 cells across). The default order leaves a large first-order error because the no-slip acts about 0.75 cell outside the wall.
- **Curved rotating wall:** the gain is modest (26 to 42 % in the velocity, 2 to 3 times in the torque). Neither order converges cleanly, because of a wall-speed deficit at the rotating inner wall that order 2 halves but does not remove.
- Order 2 is not uniformly better: it splits the wall force unevenly between two walls when a wall lies inside a cell.

## Known limitations and open issues

- The default order acts about 0.75 cell outside a flat wall; the same code path acts less far out at the curved fixed wall. The cause is not identified (not verified).
- The rotating inner wall keeps a 3 to 7 % wall-speed deficit and a non-converging mean error for both orders. Whether the density drift of the immersed boundary causes it is not verified.
- The wall shear read from near-wall cells (a three-cell quadratic fit in the channel, a cubic fit in the annulus) does not converge for either order. Use the flow rate (channel) and the torque (annulus) as the integral measures.
- 2D, one Mach number and Reynolds number per problem, one machine and compiler. A third grid beyond 40 cells per R and 3D are not measured, so the status stays "partially validated". High-Reynolds-number boundary layers, moving flat walls and slip walls are not measured.
- The rigid-rotation check (an inviscid circle in gas already in rigid rotation) was not reproduced for this page: it needs an initial condition that no case-file patch provides.

## Reproduce

```
./mfc.sh run examples/2D_ibm_plane_poiseuille/case.py -t pre_process simulation -n 1 -- --cells 16 --ntau 12 --wall_order 2
python3 examples/2D_ibm_plane_poiseuille/check_poiseuille.py examples/2D_ibm_plane_poiseuille
./mfc.sh run examples/2D_ibm_taylor_couette/case.py -t pre_process simulation -n 1 -- --cells_per_R 20 --ntau 12 --wall_order 2
python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette
```

Each run writes into the example directory (clear `D/`, `restart_data/` and the `.inp` files between runs). Channel: `--offset 0.37`, `--wall_order 1`. Annulus: `--wall_order 1`, `--ip_interp 2`, `--visc_compact 1`, `--mach 0.01 --re 11 --ntau 4`.
The CI-sized runs are the tests `2D -> IBM -> Taylor-Couette -> ib_wall_order=2` and `2D -> IBM -> Plane Poiseuille -> ib_wall_order=2` (regression goldens at 300 steps, GNU 13 double; the tolerance is @@TOL@@).
