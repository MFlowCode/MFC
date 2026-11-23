@page validation_standoff Validation: Bow-Shock Stand-Off of a Cylinder and a Sphere

# Supersonic bluff bodies: bow-shock stand-off at Mach 1.5-3 (documented over-prediction)

**Status:** open discrepancy, resolution-limited. **Code:** `master` at 0cdfdf02a. **Example:** `examples/2D_ibm_cylinder_standoff`.
**Last run:** 2026-11-22.

## Purpose

The distance of the detached bow shock from a bluff body is a sensitive check of an immersed boundary in supersonic flow: the shock
stands only a few cells from the wall at the resolutions a 3D run can afford. This measures the stand-off of a sphere at Mach 1.5, 2
and 3 and of a cylinder at Mach 2, against Billig's correlation \cite Billig1967 and against an independent body-fitted solution,
and shows how the error falls with the grid.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | sphere (`patch_ib` geometry 8) or cylinder (geometry 2), D = 1, centred on a cell-centre row of the stagnation line |
| Flow | ideal gas, \f$\gamma\f$ = 1.4, sound speed 1, Mach 1.5, 2, 3, impulsive start from the free stream, inviscid |
| Grid | sphere 20, 30, 40 and 60 cells/D at Mach 2, 20 and 40 at Mach 1.5 and 3; domain [-1.5 D, 2.5 D] x [-2 D, 2 D] in each lateral direction; cylinder 20, 40 and 80 cells/D, [-2.5 D, 2.5 D] x [-2 D, 2 D] |
| Boundaries | extrapolation on every face |
| Wall | fixed, no slip (inviscid flow, so the wall acts as a slip wall for the normal flux), inverse-distance image points |
| Numerics | WENO5-M, HLLC, SSP-RK3, acoustic CFL 0.4, fixed step |
| Length | t u/R = 40 (sphere), 60 (cylinder), 8 and 12 saves |
| Hardware | NVIDIA A100 80GB (Georgia Tech), NVHPC 25.11, double precision (sphere); GNU 13 CPU (cylinder and body-fitted solver) |

## Reference data

- Billig \cite Billig1967: \f$\Delta/R = 0.143\exp(3.24/M^2)\f$ (sphere), \f$0.386\exp(4.67/M^2)\f$ (cylinder). An inviscid, high-Reynolds-number fit; the paper's equation number was not re-read.
- The normal-shock pressure ratio of the Rankine-Hugoniot jump (p2/p1 = 4.5 at Mach 2) and the Rayleigh pitot pressure (1.253 p2 at Mach 2) for the nose-pressure check.
- Body-fitted reference: `examples/2D_ibm_cylinder_standoff/bodyfitted_sphere.py`, an axisymmetric Euler solver on a polar mesh around the sphere
  (finite volumes weighted by the cylindrical radius, MUSCL-minmod, HLLC on the rotated face normal, SSP-RK3; no immersed boundary).
  Outer radius 2.2 R, 80 and 160 cells in each direction; steady to 1e-5 R over the last 8 R/u at Mach 2.

## Quantities

- **Stand-off** \f$\Delta/R\f$: the distance from the body surface to the first point, scanning inward on the stagnation line from the inflow, where the pressure reaches the mean of p1 and the normal-shock p2; linear interpolation between cells.
  Mean of the last three saves; the range over them is listed.
- **p_nose/p2**: the pressure of the last cell before the body over p2 (theory: the pitot ratio, 1.253 at Mach 2, 1.388 at Mach 1.5, 1.167 at Mach 3).
- Script: `examples/2D_ibm_cylinder_standoff/standoff.py`.

## Results

Sphere, A100:

| M | cells/D | Delta/R | range over the last three saves | Billig | error vs Billig | error vs body-fitted |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.5 | 20 | 0.6988 | 0.001 | 0.6036 | +15.8 % | +14.0 % |
| 1.5 | 40 | 0.6650 | 0.002 | 0.6036 | +10.2 % | +8.5 % |
| 2 | 20 | 0.4167 | 0.008 | 0.3215 | +29.6 % | +19.1 % |
| 2 | 30 | 0.4025 | 0.003 | 0.3215 | +25.2 % | +15.1 % |
| 2 | 40 | 0.3939 | 0.008 | 0.3215 | +22.5 % | +12.6 % |
| 2 | 60 | 0.3776 | 0.0005 | 0.3215 | +17.5 % | +7.9 % |
| 3 | 20 | 0.2482 | 0.036 | 0.2050 | +21.1 % | +13.9 % |
| 3 | 40 | 0.2469 | 0.007 | 0.2050 | +20.4 % | +13.3 % |

Body-fitted solver (\f$\Delta/R\f$ on 80 and 160 cells): Mach 1.5 0.6130 and 0.6131 (1.6 % above Billig); Mach 2 0.3489 and 0.3498 (Billig 8.1 % low of the body-fitted value); Mach 3 0.2153 and 0.2179 (Billig 5.9 % low; changing by 1.2 % between the grids).

The nose pressure at Mach 2 is 1.15-1.26 p2 against the pitot value 1.253 (it oscillates with the shock position); at Mach 1.5, 1.29-1.41 against 1.388.

## Assessment

- The stand-off error against Billig falls with the grid at Mach 1.5 (+15.8 to +10.2 %) and 2 (+29.6, +25.2, +22.5 %), as in the earlier measurements of this configuration (+27.5, +22.0 and +17.5 % at 20, 40 and 60 cells/D in a viscous Re 5000 run; the 60 cells/D value agrees).
- Against the body-fitted solution the error at Mach 2 is +19.1, +15.1, +12.6 and +7.9 % at 20, 30, 40 and 60 cells/D, falling at an observed order of about 0.6 (20 to 40 cells/D) and 1.1 (40 to 60). Billig's correlation is 8.1 % low at Mach 2 against the converged body-fitted value, so errors quoted against it overstate the immersed boundary's by about that much.
- The sphere is resolution-limited, not a defect of the measurement: the stand-off is 4-12 cells at these grids and three definitions of the shock position agree in the earlier runs.
- Status is an open discrepancy: no grid here is in the asymptotic range, and Mach 3 is unsettled.

## Known limitations

- **Mach 3 is unsteady.** The stand-off oscillates between saves by 0.007 R (40 cells/D) to 0.036 R (20 cells/D) and the nose pressure is 1.04-1.20 p2 against 1.167. The body-fitted solver at Mach 3 shows a stand-off that still changes by 1.2 % between 80 and 160 cells.
  Whether the MFC oscillation is a carbuncle of the HLLC solver was not tested (an earlier Python prototype showed one, steady with HLL); no sensor-switched Riemann solver is proposed here.
- **Cylinder.** The stand-off of the cylinder does not settle in 60 R/u: Mach 2 at 20 cells/D 1.414 R (+14.0 % against Billig's 1.241 R, range 0.016), at 40 cells/D 1.464 R and still rising (1.425, 1.458, 1.508 over the last three saves), at 80 cells/D 1.377 R (+11.0 %, 1.354, 1.368, 1.408 over the last three saves, still rising). At Mach 1.5 (40 cells/D, domain from -3.5 D) it is still moving toward the inflow at the end (3.19 to 3.37 R over the last saves, Billig 3.08 R), and at Mach 3 it is 0.771 R (+18.9 %, range 0.026). The cylinder numbers are not a validation; the CI-sized case is a regression test.
- Viscous runs (Re 5000) were not repeated; in earlier runs the 2D stand-off changed by under 0.5 % between inviscid and Re 5000.
- First-order no-slip wall (the ghost velocity is the wall value, [#2389](https://github.com/MFlowCode/MFC/issues/2389)): it biases the stand-off the other way from the discrepancy; not quantified here.
- Run only to t u/R = 40, which settles the sphere at Mach 1.5 and 2 (60 cells/D: 0.3776 R in each of the last five saves to 0.0005); the Mach 2 sequence is not extrapolated.
- The body-fitted solver is axisymmetric inviscid Euler with an outer boundary at 2.2 R; at Mach 2 a 3 R boundary changed it by 0.0008 R in earlier runs.

## Reproduce

```
./mfc.sh run examples/2D_ibm_cylinder_standoff/case.py -t pre_process simulation -n 1 --gpu -- \
    --dim 3 --cells_per_D 40 --steps 0 --tstar 40 --frames 8
python3 examples/2D_ibm_cylinder_standoff/standoff.py <run directory>
python3 examples/2D_ibm_cylinder_standoff/bodyfitted_sphere.py --mach 2 --nr 160 --nphi 160
```

`--mach 1.5` and `--mach 3` for the other rows, `--dim 2 --xu 2.5 --steps 0 --tstar 60` for the cylinder (a build without MPI adds `--serial_io`).
The suite runs the CI-sized cylinder as `2D -> IBM -> Bow Shock Standoff`.
