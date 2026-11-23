@page validation_moving_body_budget Validation: Moving-Body Mass Budget and Galilean Twin

# Moving immersed boundary: gas mass created at the surface, and the heat and force of a body against its Galilean twin

**Status:** partially validated. **Code:** `master` at 0cdfdf02a. **Example:** `examples/2D_ibm_moving_mass_budget`.
**Last run:** 2026-11-22.

## Purpose

A ghost-cell immersed boundary is not conservative at a moving wall. Gas flows through the faces between the fluid and the body into
ghost cells whose values are then overwritten, and cells change side as the body moves, so the gas total drifts with the distance the
body travels. This measures that drift, checks the opt-in correction `ib_mass_relax` (the solver default for a moving body, 30), and
compares the conducted heat and the force of a translating body with the same body held in the same relative flow. Both matter for
particle-resolved runs, where a particle's drag and heating are read from these quantities.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | one circle, D = 1, centre off the grid lines by (dx/3, dx/5), in a doubly periodic square box |
| Gas | ideal gas, \f$\gamma\f$ = 1.4, \f$\rho\f$ = 1, p = 1/1.4 (sound speed 1), no chemistry |
| Flow | moving twin: the circle moves at Mach 0.3 through gas at rest (`moving_ibm = 1`); fixed twin: the circle is held and the gas moves at Mach -0.3 |
| Mass runs | inviscid, box 5 D, 16, 24 and 32 cells/D, 6 D travelled, 24 saves, no-slip wall |
| Twin runs | Re = U D/\f$\nu\f$ = 40, Pr 0.7, isothermal wall (`thermal_bc = 1`), wall/gas temperature 0.9 and 0.2, box 8 D, 16 and 32 cells/D, 6 D travelled, `ib_ip_interp = 2` |
| Numerics | WENO5-M, HLLC, SSP-RK3, fixed step from an acoustic CFL of 0.4, `fd_order` 4 and `weno_Re_flux` when viscous |
| Correction | `ib_mass_relax` = 0 (off) and 30 |
| Hardware | GNU 13 on an x86 CPU, double precision |

## Reference data

- Exact: the gas mass of a closed periodic box with an inert body is constant.
- Exact: a body moving through gas at rest and the same body held in gas moving at the opposite velocity feel the same relative flow, so
  their conducted heat and force agree (ratio 1). The twin is a self-consistency check, not a comparison with an external value.
- Qualitative background on the force noise and the mass error of fresh and dead cells: \cite Lee2011 and \cite SeoMittal2011. Neither
  is used as a number here.

## Quantities

- **Gas mass** M: \f$\sum\rho\,dV\f$ over the cells whose centre is outside the body at its prescribed position (periodic wrap), from the saves.
  **Steady rate**: the least-squares slope of 100 (M - M_0)/M_0 against the distance travelled after the first diameter (the impulsive
  start is excluded), in % of the gas per diameter.
- **Conducted heat** Q: \f$k\,R\!\int\partial T/\partial n\,d\theta\f$ at the wall, with \f$\partial T/\partial n\f$ from a cubic fit of T along 720 rays at 1.5, 2, 2.5 and 3 cell widths (bilinear samples). The same estimator
  is used for both twins, so only the ratio is meaningful. Means over the last 3 D travelled.
- **Force**: \f$F_x\f$ of `D/ib_forces.dat` (`ib_force_wrt`, one row per step), means over the same window; the RMS of the difference of the
  two histories (the fixed history interpolated to the moving times, mean removed) is given in percent of the mean fixed force.
- Script: `examples/2D_ibm_moving_mass_budget/check_budget.py`.

## Results

Steady rate of the gas mass, % of the gas per diameter travelled, and the total change after 6 D (inviscid, inverse-distance image
points):

| cells/D | `ib_mass_relax` = 0: rate | total | `ib_mass_relax` = 30: rate | total |
| ---: | ---: | ---: | ---: | ---: |
| 16 | +0.0503 | +0.357 % | +0.0022 | +0.016 % |
| 24 | +0.0237 | +0.192 % | -0.0008 | +0.002 % |
| 32 | +0.0130 | +0.134 % | -0.0023 | +0.003 % |

With TSC image points (`ib_ip_interp = 2`) the uncorrected rate is +0.0544 and +0.0166 %/D at 16 and 32 cells/D, and +0.0024 and
-0.0022 with the correction. The stair-step geometry makes the mass oscillate with a period of half a diameter at these grids, so the
slope uses all saves after the first diameter.

Galilean twin, viscous, moving over fixed (ratio 1 is exact agreement):

| wall T / gas T | cells/D | Q, relax 0 | Q, relax 30 | F_x, relax 0 | F_x, relax 30 | RMS of the F difference, relax 30 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.9 | 16 | 0.940 | 0.993 | 0.966 | 0.963 | 2.7 % |
| 0.9 | 32 | 0.966 | 1.002 | 0.969 | 0.968 | 1.6 % |
| 0.2 | 16 | 0.760 | 0.817 | 1.248 | 1.748 | 112 % |
| 0.2 | 32 | 0.882 | 0.917 | 0.690 | 0.682 | 7.5 % |

## Assessment

- A moving no-slip circle creates gas mass at 0.013-0.05 % of the gas per diameter travelled (inverse-distance image points), falling at about first order with refinement. Over the 6 D of the runs
  this is 0.13-0.36 %, and it would grow with the distance a particle travels.
- `ib_mass_relax = 30` holds the total to 0.016 % or less over 6 D at every grid and the steady rate to 0.002 %/D, with the start-up offset bounded.
- With the correction, the conducted heat of a wall at 0.9 of the gas temperature agrees with its twin to 0.7 % at 16 cells/D and 0.2 % at 32; the correction moves the ratio from
  0.94 and 0.97 and does not change the force (3 % low at both grids).
- A cold wall (0.2 of the gas temperature) is not reproduced: the heat is 18 % low at 16 cells/D and 8 % low at 32 with the correction,
  and the force of the moving body is not close to its twin's.
- Status is partially validated: the mass budget and the warm-wall twin close; the cold-wall twin and the force do not.

## Known limitations

- **Cold moving wall.** The conducted heat is low by 8-18 % at 16-32 cells/D even with the correction, converging about first order. The
  cause is a ghost-state error at the advancing nose that the correction does not remove (found in a separate Python prototype; not
  re-verified in MFC). The solver documentation lists measured twin ratios for the same case.
- **Force of a moving body.** `ib_forces.dat` is the lab-frame face momentum flux (`ib_face_ft`); it has no \f$-V f_m\f$ term and the cells a
  body covers or uncovers book no momentum. The twin ratio is 0.96-0.97 for a wall near the gas temperature and not meaningful for the cold wall
  (1.25 and 1.75 at 16 cells/D, 0.69 and 0.68 at 32, RMS difference 7-112 % of the fixed force). The force noise of a moving body (the cell-crossing frequency) is
  visible in the RMS column and was not separated here. Not changed in this example.
- Drift with the box: 2D, one body, Mach 0.3, inviscid; 3D bodies cross cells in three directions and a particle array was not run here
  (the fixed-body array is covered by `ib_flux_wall`, see the `2D_ibm_wall_mass_budget` example).
- The reacting moving body (gasification plus swept cells) is covered by the suite's `Moving Heat Output` tests (regression only), not by a validation case.
- The conducted-heat estimator is a cubic fit; its absolute value is not an independent measurement.

## Reproduce

```
./mfc.sh run examples/2D_ibm_moving_mass_budget/case.py -t pre_process simulation -n 1 --gpu -- \
    --cells_per_D 24 --mach 0.3 --box 5 --travel 6 --nsave 24 --relax 0
python3 examples/2D_ibm_moving_mass_budget/check_budget.py mass <run directory>
```

Twin: run the same case with `--cells_per_D 16 --mach 0.3 --box 8 --Re 40 --Twall 0.2 --ip_interp 2 --travel 6 --nsave 12 --stride 1`
once with `--relax 0`, once with `--relax 30`, and once with `--frame fixed`, then
`check_budget.py twin --moving <moving run> --fixed <fixed run>`. A build without MPI adds `--serial_io`.
The suite runs the CI-sized case as `2D -> IBM -> Moving Circle -> Mass Budget` and `... -> Cold Wall`.
