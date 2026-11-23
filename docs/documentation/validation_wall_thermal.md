@page validation_wall_thermal Validation: Solid-Wall Thermal Models

# Lumped cooling of a cylinder, and a cold wall in hot gas

**Status:** partially validated. **Code:** `master` at 0cdfdf02a. **Examples:** `examples/2D_ibm_cylinder_drag_nu` (lumped cooling,
`--thermal_bc 3`) and `examples/2D_ibm_cold_wall_transient`.
**Last run:** 2026-11-23.

## Purpose

A particle heats or cools while it sits in the gas, and the heat it absorbs in the first microseconds sets the ignition delay of a reacting
particle. Two wall models carry this in an immersed boundary: the lumped solid temperature (`thermal_bc = 3`, one temperature per body advanced
by the heat the faces carry) and the sub-grid conduction column (`patch_ib(1)%%wall_column`) that carries the flux of a wall layer thinner than a cell. This checks the first against the decay time
its own fixed-wall Nusselt number implies, and the second against a body-fitted 1D solution.

## Setup

| Item | Lumped cooling | Cold wall |
| :--- | :--- | :--- |
| Geometry | 2D cylinder, fixed | 2D cylinder R = 50 um, fixed |
| Gas | dry air, `h2o2.yaml` thermodynamics and transport, reactions off, 300 K, 1 atm, Mach 0.2, Re 20 and 40 | the same air frozen at 1500 K, 1 atm, at rest |
| Wall | `thermal_bc = 3`, T(0) = 330 K, heat capacity \f$\rho_s c_s\f$ chosen for a decay time \f$\tau\f$ of about 12 (Re 40) and 6 and 24 (Re 20) \f$D/U\f$ | `thermal_bc = 1`, 300 K, initial layer 0.004 R (a near-step), with and without `patch_ib(1)%%wall_column` |
| Grid | 12 and 16 cells per D (the 24 \f$D/U\f$ run: 12), box 16 D by 10 D, characteristic inflow and outflow | 10, 20 and 40 cells per R, box 8 R, zero-gradient faces |
| Numerics | WENO5-M, HLLC, SSP-RK3, `chem_params%%transport_dt`, `weno_Re_flux`, `fd_order = 4`, adaptive CFL | WENO5-M, HLLC, SSP-RK3, adaptive CFL 0.5, `chem_params%%transport_dt` |
| Hardware | GNU 13, CPU, one rank, double precision | GNU 13, CPU, one rank, double precision |

## Reference data

- **Lumped model.** A body of mass per area \f$\rho_s c_s D/4\f$ (2D cylinder) in a stream of heat-transfer coefficient \f$h = \mathrm{Nu}\,k_f/D\f$ relaxes to the adiabatic-wall
  temperature with \f$\tau = \rho_s c_s D^2/(4\,\mathrm{Nu}\,k_f)\f$. Nu and \f$T_{aw}\f$ come from two fixed-wall `thermal_bc = 1` runs of the same case and grid (wall at 1.0 and 1.1 times the free-stream
  temperature): MFC's own fixed-wall result, so this tests the lumped update against the face heat it is fed, not the Nusselt number against the literature (for that see @ref validation_cyl_drag_nu).
- **Cold wall.** `examples/2D_ibm_cold_wall_transient/ref1d.py` solves conduction into the cylinder's surrounding gas at uniform pressure in the mass
  coordinate \f$\psi = \int \rho r\,dr\f$ (the wall at \f$\psi = 0\f$, no gas crosses it), \f$c_p\,\partial T/\partial t = \partial_\psi(\rho r^2 \lambda\,\partial T/\partial\psi)\f$ with
  \f$r^2 = R^2 + 2\int d\psi/\rho\f$, by finite volumes on a grid stretched from the wall and BDF in time, with the Cantera
  \f$\lambda(T)\f$ and \f$c_p(T)\f$ of the same air; the heat into the body \f$Q(t) = \int 2\pi R\,q\,dt\f$ is carried as one more unknown. Q changes by under 2e-4 between 800 and 1600 cells and between outer radii of 8 R and 16 R, and the solver shares no code with MFC.
  The initial gas profile is MFC's own, \f$T = T_\infty + (T_w - T_\infty)(R/r)^{1/2}\mathrm{erfc}((r-R)/\delta)\f$.

## Quantities

- **Lumped.** \f$T_w(t)\f$ from column 6 of `D/ib_heat.dat`; \f$\tau\f$ from the formula above with the twins' Nu (heat of column 5 of the hot and cold twin, the difference removing the viscous heating) and \f$k_f\f$ at the film temperature;
  the fit \f$T_a' + A\exp(-t/\tau')\f$ over \f$tU/D\f$ above the fit start with a free asymptote and \f$\mathrm{Nu}_\mathrm{fit} = \mathrm{Nu}\,\tau/\tau'\f$ (`check_lumped_cooling.py`).
- **Cold wall.** \f$Q(t)\f$, the running sum of column 5 of `D/ib_heat.dat` (energy into the body per unit depth, the face heat including the enthalpy of exchanged mass), against the reference at times in units of \f$R^2/\alpha_\infty\f$ of the hot gas (`check_heat.py`).

## Results

**Lumped cooling** (fit window and run lengths in the table; the fixed-wall twins are averaged over their last 20 %):

| Re | cells per D | \f$\rho_s c_s\f$ [J/m^3/K] | run length \f$tU/D\f$ | model \f$\tau U/D\f$ | fit \f$\tau' U/D\f$ | Nu (twins) | \f$\mathrm{Nu}_\mathrm{fit}\f$ | difference | heat rate / twins' law, first and last of window |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 40 | 12 | 6.9e3 | 36 | 12.1 | 13.3 | 3.296 | 3.002 | -8.9 % | 0.936, 0.846 |
| 40 | 16 | 6.9e3 | 36 | 12.2 | 13.4 | 3.278 | 2.976 | -9.2 % | 0.932, 0.836 |
| 20 | 12 | 5.1e3 | 20 | 5.9 | 7.0 | 2.501 | 2.104 | -15.9 % | 0.887, 0.760 |
| 20 | 16 | 5.1e3 | 20 | 6.0 | 7.1 | 2.480 | 2.080 | -16.1 % | 0.882, 0.760 |
| 20 | 12 | 2.05e4 | 60 | 23.7 | 25.1 | 2.501 | 2.365 | -5.5 % | 0.955, 0.920 |

The fit window starts at \f$tU/D\f$ = 6 (Re 40, \f$\tau\f$ about 12), 3 (Re 20, \f$\tau\f$ about 6) and 8 (Re 20, \f$\tau\f$ about 24); the twins ran to \f$tU/D\f$ = 30 (Re 40) and 20 (Re 20). The last heat-rate entry of a window is sensitive: \f$T_w - T_{aw}\f$ has fallen to 2 K there, so 0.2 K in \f$T_{aw}\f$ moves the ratio by 10 %.

**Cold cylinder in 1500 K gas, heat absorbed against the reference** (7.4 us = 1 \f$R^2/\alpha_\infty\f$; error \f$Q_\mathrm{MFC}/Q_\mathrm{ref} - 1\f$):

| cells per R | wall | 0.003 | 0.01 | 0.03 | 0.1 |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 10 | ghost mirror (default) | -51.9 % | -33.7 % | -19.7 % | -9.6 % |
| 20 | ghost mirror (default) | -37.5 % | -23.1 % | -12.9 % | -6.1 % |
| 40 | ghost mirror (default) | -20.8 % | -10.3 % | -4.2 % | -1.4 % |
| 10 | wall column | +12.7 % | +6.4 % | +2.5 % | +0.7 % |
| 20 | wall column | +7.1 % | +0.8 % | -1.0 % | -0.4 % |
| 40 | wall column | +6.0 % | -0.3 % | -1.3 % | +0.0 % |

## Assessment

- The lumped cylinder decays as an exponential but 5.5 to 16 % more slowly than its fixed-wall Nusselt number implies. The gap does not depend on the grid (-8.9 % at 12 cells per D and -9.2 % at 16 for Re 40; -15.9 % and -16.1 % for Re 20) and falls as the decay time grows (-16 % at \f$\tau U/D\f$ = 6 to -5.5 % at 24 for Re 20): a short \f$\tau\f$ makes the thermal boundary layer lag the falling wall temperature, so the heat rate sits below the steady law. A gap of about 5 % remained at the longest decay time run (the Re 40 run at \f$\tau U/D\f$ about 48 was not completed). The heat the faces report is 0.89 to 0.96 of the steady law at the same wall temperature early in the window.
- Status of the lumped model: partially validated: a 5 % under-prediction of the heat rate against MFC's own fixed-wall result, cause not found (not verified: the face energy flux of the first-order wall, which also carries a small spurious mass flux, see the discussion of `ib_heat_wrt` in the case documentation).
- The default wall under-predicts the heat a cold wall absorbs from hot gas by 21 to 52 % at the earliest time and by 1.4 to 9.6 % at 0.1 \f$R^2/\alpha_\infty\f$, falling with refinement at roughly first order (21 % at 40 cells per R at 0.003 against 52 % at 10). The wall layer is thinner than a cell for most of this transient, and the ghost mirror carries at most the flux of a layer one cell thick.
- With `patch_ib(1)%%wall_column` the heat is within 1.3 % of the reference from 0.01 \f$R^2/\alpha_\infty\f$ on at 20 and 40 cells per R, and within 2.5 % at 0.03 on all three grids. It over-predicts at the earliest time (+6.0 % at 40 cells per R, +12.7 % at 10), a bias that does not fall with refinement in the range run.
- Partially validated: the lumped cooling is two Reynolds numbers on two coarse grids against MFC's own fixed-wall Nusselt number; the cold wall has a low-Mach reference only; a moving or reacting wall and a sphere were not run.

## Known limitations

- **Low-Mach column.** The conduction column is a low-Mach sub-grid model; its early over-prediction (+6 % at 0.003 \f$R^2/\alpha_\infty\f$, 22 ns for a 50 um particle in 1500 K air) is probably the thermoacoustic reduction of the heat of a compressible gas that a constant-pressure reference does not have (the pressure of the gas next to the wall dips while it cools). A compressible reference was not run here, so the size of this effect is not verified. It decays with time and is at or below 1.3 % by 0.03 \f$R^2/\alpha_\infty\f$ (0.22 us), a small part of the gas-side heating time of a particle.
- The heat booked on the body at the first step of a run differs from what the columns drain by several percent (observed on a reacting wall, not traced); the earliest entries of the cold-wall table may include it.
- The column applies to a fixed isothermal circle or sphere with chemistry, `chem_params%%diffusion`, `chem_params%%transport_model = 1` and an adaptive or constant CFL step, and was run on circles here.
- The lumped model is one temperature per body (small Biot number), explicit in the update and has no recession.
- The lumped reference is MFC's own fixed-wall Nusselt number: an error common to the fixed and lumped runs (for example the first-order wall) does not appear in it.

## Reproduce

```
./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation -n 1 -- --Re 40 --cells_per_D 16 --xu 5 --xd 11 --hy 5 --steps 0 --tstar 30 --tw_ratio 1.0
./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation -n 1 -- --Re 40 --cells_per_D 16 --xu 5 --xd 11 --hy 5 --steps 0 --tstar 30 --tw_ratio 1.1
./mfc.sh run examples/2D_ibm_cylinder_drag_nu/case.py -t pre_process simulation -n 1 -- --Re 40 --cells_per_D 16 --xu 5 --xd 11 --hy 5 --steps 0 --tstar 36 --tw_ratio 1.1 --thermal_bc 3 --rho_cs 6.9e3
python3 examples/2D_ibm_cylinder_drag_nu/check_lumped_cooling.py <lumped run> --cold <cold run> --hot <hot run> --tmin 6

./mfc.sh run examples/2D_ibm_cold_wall_transient/case.py -t pre_process simulation -n 1 -- --cells_per_R 40 --half 4 --tend 0.1 --wall_column 1
python3 examples/2D_ibm_cold_wall_transient/check_heat.py examples/2D_ibm_cold_wall_transient
```

Each run writes into the example directory, so copy it between runs. The 40 cells per R cold-wall run takes about 20 minutes on one core. The suite runs the CI-sized cold wall (12 cells per R, 0.03 \f$R^2/\alpha_\infty\f$) as `2D -> Chemistry -> IBM -> Cold Wall Transient -> Wall Column`.
