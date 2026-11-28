@page ib_accuracy_limits Accuracy Limits of the Ghost-Cell Immersed Boundary

# Accuracy limits of the ghost-cell immersed boundary

**Status:** open discrepancy (every item below is a limitation by construction; none is a bug report). **Code:** `master` at 0cdfdf02a.
**Last update:** 2026-11-27. The wall-temperature rows reproduce with `examples/2D_ibm_accuracy_limits` (`order_study.py`, re-measured on `master` after #2486, same orders). Its examples are `examples/2D_ibm_conduction_annulus` (wall temperature, image-point weights) and `examples/2D_ibm_taylor_couette` (wall velocity, torque and wall shear, `ib_wall_order`); both are registered tests, and their READMEs hold the `master` tables.

## Purpose

A sharp ghost-cell immersed boundary (IB) enforces a wall value (velocity or temperature) through the state of its ghost cells.
Field-wide errors and global integrals (drag, total wall heat) converge faster than the two quantities the wall actually sets: the wall value
and the wall-normal gradient (wall shear, heat flux). This page states what is known about their order, which options help, and what does not.

## Measured orders

Orders are observed orders between successive grids (cells per radius, R). "Pre-merge" means measured earlier on a development build
with the options named, not repeated on `master`; the options exist on `master` under the same names but their merged form may differ in detail.

| Quantity | Setup | Order (finest pairs) | Provenance |
| :--- | :--- | :--- | :--- |
| Wall temperature, inverse-distance image point (`ib_ip_interp = 1`, the default before TSC) | annulus, exact conduction solution | 0.6 then 1.2 (idw), error 5.0e-3 of \f$T_2-T_1\f$ at 40 cells/R | `master`, @ref validation_tsc_annulus |
| Wall temperature, TSC (`ib_ip_interp = 2`) | same | 2.3 then 2.0 (inner wall), 1.1 then 1.0 (outer wall) | `master`, @ref validation_tsc_annulus |
| Wall velocity, `ib_wall_order = 1` | Taylor-Couette, 10-30 cells/R | about 1 to 2, not monotonic at the outer wall | pre-merge |
| Wall shear, static curved wall, `ib_wall_order = 1` | same | 0.7 to 1 or worse, 1.4 % error at 30 cells/R | pre-merge |
| Wall shear, static curved wall, `ib_wall_order = 2` | same | 2.8 then 2.4, 0.065 % at 30 cells/R | pre-merge |
| Wall shear, rotating wall, any option | same | below 0.5, plateau near 2 % | pre-merge |
| Thermal wall gradient (fit through the first cells), inverse distance | annulus, 10-40 cells/R | 0.4 to 1.4 | pre-merge |
| Thermal wall gradient, TSC | same | 1.4 to 1.6, error about 0.1 % at 40 cells/R | pre-merge |
| Curved isothermal wall, temperature \f$L_\infty\f$, \f$L_1\f$, total heat, inverse-distance weights | annulus, 600/300 K, 10-80 cells/R | 0.85, 1.1-1.4, 0.9-1.3 between 40 and 80 cells/R | pre-merge, before #2248 |
| Flat isothermal wall at 0 degrees to the grid | slab, 10-80 cells per gap | 2.0 in all norms | pre-merge, before #2248 |
| Flat isothermal wall at 26.6 and 45 degrees | same | 1.6 to 1.75 on average, scattered from grid to grid | pre-merge, before #2248 |

Reading the table:

- **The wall value is about second order** once the grid is fine enough with TSC or `ib_wall_order = 2`. With inverse-distance weights (`ib_ip_interp = 1`) it is first order on the annulus.
- **The wall-normal gradient is first order or worse** with every option tried. TSC lowers its error about fourfold in the thermal annulus
  (order 1.4 to 1.6); none of the options reaches second order at 40 cells/R.
- **A rotating wall keeps a shear error of about 2 % for every option** (rms over angle, 10-30 cells/R, orders below 0.5). The wall value there is right to 0.1 %, and the
  azimuthal-mean velocity error is 0.05 % of the wall speed, so the cause is not an interpolation effect. Not explained.
- **Curved isothermal walls**: the max-norm temperature error was first order and the total wall heat about order 1.3 (0.06 % of \f$T_1-T_2\f$ and 0.01 % of the exact heat at 80 cells/R).
  These were measured before the unbounded mirror temperature of [#2248](https://github.com/MFlowCode/MFC/pull/2248) and with inverse-distance weights; they were not repeated on `master`.
- Wall values that look converged do not imply a converged wall flux: the wall heat of a closed conduction problem is fixed by conservation and agrees between the two bodies to five digits,
  so it is a global integral, not a measure of the local gradient.

## What to set

| Goal | Setting | Evidence |
| :--- | :--- | :--- |
| Isothermal wall temperature of a conduction or heated wall | `ib_ip_interp = 2` | wall error 5 to 130 times lower, whole-gap error 6 to 14 times lower: @ref validation_tsc_annulus |
| No-slip velocity at a flat or curved static wall | `ib_wall_order = 2` | measured gain and where it fails (4 to 6 cells across a gap): the `ib_wall_order` entry of @ref case "the case documentation" |
| Force and heat of a cylinder or sphere | face-flux force `ib_force_wrt`, heat `ib_heat_wrt` | agree with references and with a momentum balance: @ref validation_cyl_drag_nu |
| Gas mass of a closed box with a body | `ib_flux_wall = T` (reacting or inert fixed no-slip body); `ib_mass_relax` for slip or moving bodies | `examples/2D_ibm_wall_mass_budget` |
| A body in contact with a no-slip or moving wall differing between CPU and GPU builds (about 1e-5) | no setting closes it: the central tangential velocity of the HLLC viscous work (always on) removed one cause, but the measured A100-against-CPU difference stays ([#2419](https://github.com/MFlowCode/MFC/issues/2419)); budget the tolerance | the `riemann_solver` entry of the case documentation |
| Wall shear or local heat flux | read an integrated force or total wall heat instead | a local fit through the near-wall cells does not converge for any option above |

The options are described in @ref ibm "the immersed-boundary overview". Resolution is the remaining lever: the gradient errors above fall by about a factor of 2 per grid doubling.

## What does not help

- **A larger minimum image-point distance** (`ib_ip_min_dist = 1.5`) is no better than TSC for the wall value (errors 2 to 3 times those of TSC in the annulus) and 5 to 7 % worse for the rotating-wall shear.
- **The quadratic wall fill for a rotating wall**: the shear plateau near 2 % is unchanged (pre-merge).
- **The mirror-weight cap, as a cure.** The pre-merge ghost-temperature rule capped the mirror weight at 1. A steady discrete reference solver of the same equations (a prototype
  in Python, not shipped) that removed the cap converged at second order on the curved wall, which would explain the first-order curved-wall error. This was *not verified in MFC*:
  the ghost rule changed with [#2248](https://github.com/MFlowCode/MFC/pull/2248), and no `master` run isolates the cap. Treat the cause of the curved-wall first-order error as open.
- **A conservation-only volume-fraction approach for moving bodies** (cut-cell volume fractions that book the displaced gas once, after Seo and Mittal \cite SeoMittal2011) was investigated and is not part of MFC.
  The variants tried did not keep a uniform free stream exact and the claimed heat gain on a moving twin
  (a cylinder translating through still gas against the same cylinder fixed in the same relative flow) was not reproduced. Exact free stream together with local conservation with level-set
  fractions would need face apertures, i.e. a full cut-cell method, which this IB does not have. No code ships.

## Known limitations of the measurements

- The Taylor-Couette runs stop at 30 cells/R: at 40 cells/R every scheme, order 1 with inverse distance included, aborted with an image-point bounds error (cause not investigated). TSC has no 30 cells/R point for the same reason.
- Orders come from three grids; the measurement floor is not negligible on the coarsest (inner-wall gradient floor 1.4 % at 10 cells/R). Coarse-pair orders (10 to 20 cells/R) were flattered by being pre-asymptotic and fell on refinement, so only the finest pairs are quoted.
- The wall gradient here is a cubic or quadratic fit through the fluid cells near the wall, not the ghost-cell gradient (ghost cells are never used); it measures what the solution holds at the wall.
- Slab orders at 26.6 and 45 degrees scatter because the wall-to-cell-centre offsets change with each grid; no single order is claimed.
- Thermal runs used quiescent moist air with mixture-averaged transport and an exact Kirchhoff solution; velocity and temperature were never coupled in these measurements.

## Reproduce

The two wall-temperature rows, with the page check, in one command (about 9 minutes on one GPU):

```
python3 examples/2D_ibm_accuracy_limits/order_study.py -- -n 1 --gpu mp
```

No script is shipped for the pre-merge tables: they came from analysis tools tied to run directories of a development build. The `master` rows reproduce with

```
./mfc.sh run examples/2D_ibm_conduction_annulus/case.py -n 1 -- --cells_per_R 20 --interp 2
python3 examples/2D_ibm_conduction_annulus/check_annulus.py examples/2D_ibm_conduction_annulus
```

The Taylor-Couette rows reproduce with

```
./mfc.sh run examples/2D_ibm_taylor_couette/case.py -t pre_process simulation -n 1 -- --cells_per_R 20 --ntau 12 --wall_order 2
python3 examples/2D_ibm_taylor_couette/check_couette.py examples/2D_ibm_taylor_couette
```

Its README measures on `master` (with #2383) the wall velocity of `ib_wall_order = 1` and `2` at 10, 20 and 40 cells/R: order 2 converges at
about second order and order 1 at 0.5 to 0.9; TSC image points with `ib_visc_compact` bring the near-wall shear to 0.1 % at 40 cells/R, below
the 2 % plateau of the pre-merge rows above.

Open items to repeat on `master`: the annulus \f$L_\infty\f$, \f$L_1\f$ and total heat orders with `ib_ip_interp` 1 and 2 at 10 to 80 cells/R; the Taylor-Couette shear with `ib_wall_order` 1 and 2
(and the 40 cells/R abort).
