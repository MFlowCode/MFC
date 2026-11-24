@page validation_wall_order_couette Wall-Normal Ghost Velocity at Viscous Walls

# Wall-Normal Ghost Velocity at Viscous Walls

**Status:** partially validated (`ib_wall_order = 2` is measured on a flat wall and on a rotating curved wall; 2D only, two grids per regime beyond the first, and the Mach 0.05 wall shear is open).
**Code:** channel: `master` at dd9778b26; annulus: `master` at 13b4d8b (contains #2176 `ib_wall_order`, #2169 centripetal ghost pressure, #2165 TSC image points, #2212 compact viscous gradients and #2383 wall-point ghost velocity), unmodified. NVIDIA A100, NVHPC 25.11, OpenACC, double precision, one rank.
**Last run:** 2026-11-22. **Examples:** `examples/2D_ibm_plane_poiseuille`, `examples/2D_ibm_taylor_couette`.

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
| 10 | 1 | 1.48e-2 | 1.31e-2 | 9.3e-3 | -1.9e-2 | -3.2e-2 | -5.7e-2 |
| 10 | 2 | 1.36e-3 | 1.23e-3 | 7.7e-4 | -2.7e-4 | -2.3e-3 | -7.4e-3 |
| 20 | 1 | 8.0e-3 | 7.7e-3 | 4.7e-3 | -1.28e-2 | -1.90e-2 | -4.4e-2 |
| 20 | 2 | 3.3e-4 | 3.2e-4 | 1.9e-4 | +7.6e-4 | -3.8e-4 | -2.2e-3 |
| 40 | 1 | 5.8e-3 | 6.3e-3 | 2.4e-3 | -1.29e-2 | -1.56e-2 | -4.7e-2 |
| 40 | 2 | 9.6e-5 | 1.0e-4 | 4.1e-5 | +1.4e-4 | -4.4e-5 | -7.2e-4 |

- With order 2 the velocity error is 11, 24 and 60 times lower than with order 1 at 10, 20 and 40 cells per R, and converges at about second order (1.36e-3, 3.3e-4, 9.6e-5: observed orders 2.0 and 1.8). Order 1 converges at 0.9 and then 0.5 order (1.48e-2, 8.0e-3, 5.8e-3).
- The azimuthal-mean part dominates the error for both orders. With order 1 the fitted wall speed stays 4 to 6 % low and the torque 1.3 to 3 % off at every resolution; with order 2 the wall speed deficit falls from 0.7 % to 0.07 % and the torque error from 0.2 % to 0.01 % (inner) at 10 to 40 cells per R.
- The inner-wall shear from local cubic fits at 96 wall points, net of the fit's own error: order 2 gives 1.9e-2, 2.0e-2 and 1.9e-2 at 10, 20 and 40 cells per R (the fit's own error is 1.9e-2, 3.2e-3 and 4.5e-4 on the exact profile); order 1 gives 0.12, 0.20 and 0.98. The order 2 shear does not converge at this Mach and Reynolds number (about 2 %); the cause is not identified (not verified). At Mach 0.01 it is 0.4 to 0.6 % (next section).
- The change of the azimuthal velocity over the last save interval is below 1.1e-5 of Omega R in these runs.

### Rotating inner wall with the optional wall treatments, Mach 0.01

Same annulus at tip Mach 0.01 and Re = 11, run for 4 slowest viscous times from rest, `ib_wall_order = 1` unless noted. The centripetal ghost pressure of #2169 is always on for a spinning body; `ib_ip_interp = 2` (TSC) and `ib_visc_compact = T` (compact viscous gradients next to a body) are opt-in. Shear is the inner-wall shear \f$S = r\,d(u_\theta/r)/dr\f$ from local cubic fits at 96 wall points (3.5 to 5 cells), rms relative error; the same fit applied to the exact profile gives 3.2e-3 (20 cells per R) and 9.2e-4 (30). The runs are not fully steady (the velocity changes by 1.7e-3 to 1.9e-3 of Omega R over the last save interval).

| Cells per R | Options | L2 / Omega R | Azimuthal mean | Wall speed (fit) | Inner torque | Shear rms |
| ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 20 | default | 9.6e-3 | 7.7e-3 | -2.9e-2 | +1.3e-2 | 0.21 |
| 20 | TSC | 9.7e-3 | 7.8e-3 | -2.9e-2 | +1.2e-2 | 0.21 |
| 20 | TSC, compact | 8.5e-3 | 6.8e-3 | -1.3e-2 | +1.1e-2 | 0.18 |
| 20 | TSC, compact, order 2 | 4.0e-3 | 4.0e-3 | -5.7e-4 | +8.4e-3 | 3.8e-3 |
| 30 | default | 8.2e-3 | 7.8e-3 | -3.2e-2 | +5.0e-3 | 0.29 |
| 30 | TSC | 8.3e-3 | 7.8e-3 | -3.2e-2 | +4.6e-3 | 0.29 |
| 30 | TSC, compact | 7.1e-3 | 6.7e-3 | -1.2e-2 | +4.2e-3 | 0.13 |
| 30 | TSC, compact, order 2 | 3.9e-3 | 3.9e-3 | -2.9e-4 | +8.9e-3 | 6.4e-3 |

- TSC alone changes nothing here. Compact viscous gradients halve the wall-speed deficit and the shear error, and adding order 2 removes the deficit (0.06 % and 0.03 %) and takes the shear to 0.4 % (20 cells per R) and 0.6 % (30), which agrees with the 0.2 to 2 % range measured earlier on a build that already evaluated the wall velocity at the wall point. The shear does not fall from 20 to 30 cells per R here (3.8e-3 to 6.4e-3); the runs are short and not fully steady, so this is not a convergence measurement.
- The inner torque stays within 1.3 % with order 1 and 0.9 % with order 2; the outer torque is 1.8 to 3.3 % low for every option.
- Before #2383, `master` evaluated a rotating body's wall velocity at the ghost point instead of the wall point (in `s_compute_ghost_point_velocity`), so the imposed wall speed was \f$\Omega(R-|\phi|)\f$, not \f$\Omega R\f$. The same runs on that build gave an inner shear of 28 to 47 % for every option, which is why #2383 was made. The numbers above are with #2383.

## Assessment

- **Flat wall:** order 2 is a large gain in the velocity profile and the flow rate (8 to 17 times at 8 to 32 cells across). The default order leaves a large first-order error because the no-slip acts about 0.75 cell outside the wall.
- **Curved rotating wall (with #2383):** order 2 gives second-order convergence of the velocity (1.4e-3, 3.3e-4, 9.6e-5), 11 to 60 times below order 1, and removes the wall-speed deficit. Order 1 converges at 0.5 to 0.9 order and keeps a 4 to 6 % wall-speed deficit.
- **Wall shear:** from near-wall cells it does not converge for order 1, and for order 2 it stays near 2 % at Mach 0.05 and Re 10. At Mach 0.01 with order 2, TSC and compact gradients it is 0.4 to 0.6 %. The flow rate (channel) and the torque (annulus) are the integral measures that converge.
- Order 2 is not uniformly better: it splits the wall force unevenly between two walls when a wall lies inside a cell.

## Known limitations and open issues

- The default order acts about 0.75 cell outside a flat wall; the same code path acts less far out at the curved fixed wall. The cause is not identified (not verified).
- The order 2 inner-wall shear at Mach 0.05 and Re 10 stays near 2 % from 10 to 40 cells per R. The cause is not identified (not verified).
- The wall shear read from near-wall cells (a three-cell quadratic fit in the channel, a cubic fit in the annulus) does not converge for the channel with either order.
- 2D, two Mach numbers and Reynolds numbers, one machine and compiler. A third grid beyond 40 cells per R and 3D are not measured, so the status stays "partially validated". High-Reynolds-number boundary layers, moving flat walls and slip walls are not measured. The Mach 0.01 runs are not fully steady.
- The rigid-rotation check (an inviscid circle in gas already in rigid rotation) was not reproduced for this page: it needs an initial condition that no case-file patch provides.
- The registered tests agree across the CPU lanes to 1e-13 and between an NVIDIA A100 and the CPU golden to 1e-14, so their tolerance is 1e-10. An A100 difference of 5e-6 reported earlier came from an annulus golden generated in single precision, not from a CPU/GPU gap.

## Reproduce

```
./mfc.sh run examples/2D_ibm_plane_poiseuille/case.py -t pre_process simulation -n 1 -- --cells 16 --ntau 12 --wall_order 2
python3 examples/2D_ibm_plane_poiseuille/check_poiseuille.py examples/2D_ibm_plane_poiseuille
./mfc.sh run examples/2D_ibm_taylor_couette/case.py -t pre_process simulation -n 1 -- --cells_per_R 20 --ntau 12 --wall_order 2
python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette
```

Each run writes into the example directory (clear `D/`, `restart_data/` and the `.inp` files between runs). Channel: `--offset 0.37`, `--wall_order 1`. Annulus: `--wall_order 1`, `--ip_interp 2`, `--visc_compact 1`, `--mach 0.01 --re 11 --ntau 4`.
The CI-sized runs are the tests `2D -> IBM -> Taylor-Couette -> ib_wall_order=2` and `2D -> IBM -> Plane Poiseuille -> ib_wall_order=2` (regression goldens, GNU CPU, double: 100 steps for the annulus, 300 for the channel; the tolerance is 1e-10: the GNU, Intel, NVHPC and macOS CI lanes agree to 1e-13 and an A100 to 1e-14).
