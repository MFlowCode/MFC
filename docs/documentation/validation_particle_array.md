@page validation_particle_array Validation: Periodic Particle Arrays

# Periodic arrays of fixed particles: drag, heat transfer and the reacting-array budget

**Status:** partially validated. **Code:** `master` at 0cdfdf02a. **Example:** `examples/2D_ibm_periodic_array`.
**Last run:** 2026-11-23.

## Purpose

Particle-resolved training and closure data for gas-solid flow come from periodic boxes holding many particles, driven at a fixed
bulk flow, with a force and a heat flux recorded for every particle. This checks that configuration: the global momentum balance of the
periodic box, the mean drag of random sphere arrays against the Tenneti, Garg and Subramaniam correlation \cite Tenneti2011, the
spread of the per-particle force, and the Nusselt number against Gunn \cite Gunn1978. It also records what is known of the burning-array
mass budget.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | `N` fixed no-slip isothermal bodies in a periodic box, placed by seeded random sequential addition (minimum surface gap 0.1 D, bodies may straddle the periodic faces); 2D cylinders or (`--dim 3`) spheres |
| Flow | dry air (`h2o2.yaml`, chemistry on for thermodynamics and mixture-averaged transport, reactions off), 300 K, 1 atm; the fluid-averaged mass flux is held at \f$Q = \mathrm{Re}_m \mu/((1-\phi) D)\f$ by `mean_mflux` with `mean_dpdx_gain = 0.05`, so \f$\mathrm{Re}_m = \rho(1-\phi) W D/\mu\f$ with \f$W\f$ the fluid-phase mean velocity; Mach 0.3 (3D), 0.2 (2D) |
| Grid | 16 and 24 cells/D, uniform; 3D: 20 spheres at \f$\phi\f$ = 0.2 in a cube of 3.74 D (0.26 and 0.88 M cells) |
| Boundaries | triply (2D: doubly) periodic |
| Wall | fixed, no slip, `thermal_bc = 1` at the free-stream temperature (they absorb the heat the driving force dissipates); `mass = 0` keeps a body-force term out of the force record |
| Numerics | WENO5-M, HLLC, SSP-RK3, `weno_Re_flux`, `fd_order = 4`, fixed step at acoustic CFL 0.2 (3D) or 0.3 (2D) because the controller needs a constant `dt`; `ib_ip_interp = 1`; these runs used `ib_wall_order = 1` and `ib_visc_compact = F`, the defaults then (the Reproduce commands pass `--wall_order 1 --visc_compact 0`; the example now follows the solver defaults 2 and T) |
| Length | \f$t W/D\f$ = 10 (3D), 43 (2D), averages over the last 40 % |
| Hardware | NVIDIA A100 80GB, NVHPC 25.11, double precision |

## Reference data

- Tenneti, Garg and Subramaniam \cite Tenneti2011, Eq. 43: the total force per particle over \f$F_\mathrm{St} = 3\pi\mu D(1-\phi)W\f$, fit range \f$0.1 \le \phi \le 0.5\f$, \f$\mathrm{Re}_m \le 300\f$. The immersed-boundary force has no mean-pressure-gradient part, so the total force is the recorded force over \f$(1-\phi)\f$ (their Eqs. 33 and 38). The authors report an average deviation of 2.5 % of the fit from their data.
- Gunn \cite Gunn1978 for the particle Nusselt number of a fixed bed, with \f$\epsilon = 1-\phi\f$ and \f$\mathrm{Re} = \mathrm{Re}_m\f$; and Sun, Tenneti and Subramaniam \cite Sun2015 for the thermally developed state (`mean_temp_wall`).
- Exact: the momentum balance of a periodic box, \f$\sum_i F_{x,i} = -\partial_x \bar p\, V_f\f$ at steady state.
- No reference exists for the 2D cylinder array: its checks are the balance and the spread.

## Quantities

All from `check_array.py` on the run directory (`D/ib_forces.dat`, `D/ib_heat.dat`, `mean_dpdx.dat`; no field output).

- **Q / Q\***: the held mass flux over its target.
- **balance**: \f$\sum_i F_{x,i}/(-\partial_x\bar p\,V_f)\f$, with \f$V_f\f$ the fluid cells of the staircase mask (the controller's definition).
- **F / F_St** (3D), **Cd** (2D, total force per length over \f$\tfrac12\rho((1-\phi)W)^2 D\f$).
- **spread**: standard deviation over the mean of the time-averaged streamwise force of the bodies.
- **Nu** (`--nu`, runs with `--thermostat`): the heat the bodies give the gas over \f$N \pi D^2 k_f (T_s - T_m)\f$, \f$k_f\f$ at the film temperature, \f$T_m\f$ = `mean_temp`.

## Results

2D, 12 cylinders at \f$\phi\f$ = 0.2, \f$\mathrm{Re}_m\f$ = 20, 14 cells/D (the CI-sized case run to \f$tW/D\f$ = 43): Q/Q\* = 1.000000, balance 1.0000, Cd 9.58, spread 0.44 (bodies between 0.40 and 1.72 of the mean).

3D, 20 spheres at \f$\phi\f$ = 0.2, mean drag (error against Tenneti Eq. 43 at the realized \f$\mathrm{Re}_m\f$):

| Re_m | seed | cells/D | F/F_St | Tenneti | error | balance | spread |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 1 | 16 | 7.075 | 7.366 | -4.0 % | 1.0001 | 0.20 |
| 20 | 2 | 16 | 7.354 | 7.355 | -0.0 % | 1.0001 | 0.18 |
| 100 | 1 | 16 | 13.973 | 12.621 | +10.7 % | 1.0001 | 0.25 |
| 100 | 2 | 16 | 14.337 | 12.605 | +13.7 % | 1.0001 | 0.24 |
| 100 | 1 | 24 | 13.894 | 12.606 | +10.2 % | 1.0000 | 0.25 |

### Nusselt number

Thermally developed state (`--thermostat`: `mean_temp_wall` at the body temperature, mixing-cup temperature held at 300 K), 20 spheres at \f$\phi\f$ = 0.2, \f$\mathrm{Re}_m\f$ = 20, 16 cells/D, Mach 0.1, \f$tW/D\f$ = 8, last 30 %, bodies at 290 K and 310 K (seed 1). Gunn gives 5.995 and the Sun et al. correlation 3.39 on the same (surface-area) basis, from the equation as fitted to their data (a mixing-cup Nusselt number per cross-section perimeter, times \f$\pi/4\f$).

The thermostat removes the viscous dissipation from the energy equation (`mean_temp_wall` with `viscous`, \cite Patankar1977); the first columns are with it removed (A100, NVHPC 25.11, one rank, last 40 %), the last column is before (last 30 %).

| bodies | Nu | Gunn | difference | F/F_St | Nu before the removal |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 290 K (cold) | 6.113 | 5.995 | +2.0 % | 7.361 | 6.647 |
| 310 K (hot) | 5.711 | 5.998 | -4.8 % | 7.432 | 5.039 |
| mean of the pair | 5.91 | 5.996 | -1.4 % | | 5.84 |

The hot and cold runs disagree by 7 % (28 % before the removal, which was mostly viscous heating), so the pair mean is the quotable number: 1.4 % below Gunn and 1.7 times the Sun et al. value. The drag does not depend on the temperature difference (7.36 and 7.43).

### Burning array (not repeated on `master`)

A reacting array (carbon spheres in hot air, the Bradley-Chelliah surface mechanism of `examples/2D_ibm_reacting_surface`, lumped solid wall `thermal_bc = 3`, `mean_temp` thermostat) was run in earlier development builds, not on `master`, and is not part of the example or the suite (11-species chemistry is too heavy for the suite). The measurements are quoted with their settings: 64 spheres at \f$\phi\f$ = 0.10 in a cube of 6.96 D, 2000 K and 10 atm air, Re_m 100 held by `mean_mflux`, 24 cells/D, `ib_ip_interp = 2`, one realization, two flow-throughs.

- The per-particle gasification rate agreed with a single-film (kinetic plus diffusion) estimate to +2 % at 24 and 32 cells/D and was 8 % low at 12 cells/D.
- The gas of the closed box gained 9.5 % more mass than the surface chemistry gasified at 24 cells/D, falling with resolution (28, 23, 9.5 and 5.2 % at 12, 16, 24 and 32 cells/D): a resolution error of the ghost-cell wall, which the opt-in `ib_flux_wall` closes (the face carbon equals the face mass flux to round-off).
- Drag was +7.6 % above Tenneti and the Nusselt number (conductive, mixing-cup temperature) +6 % above Gunn, within the scatter of the non-reacting arrays above.

The single-body closed-box budget is reproducible on `master` with `examples/2D_ibm_wall_mass_budget`.

## Assessment

- The momentum balance closes to 1e-4 in every run and the controller holds the flux exactly: the periodic box, the force record and the controller agree.
- At \f$\mathrm{Re}_m\f$ = 20 the mean drag of two 20-sphere realizations is 2 % below Tenneti (the realizations differ by 4 %, the size of the scatter between realizations of such small boxes). At \f$\mathrm{Re}_m\f$ = 100 it is 11 to 14 % above at 16 cells/D.
- Status is partially validated: boxes of 20 spheres, two realizations and at most two grids are not enough to separate the correlation's own 2.5 % deviation, the realization scatter and the grid error.

## Known limitations

- The Re 100 excess is a grid effect in part: it changes little from 16 to 24 cells/D (+10.7 % to +10.2 % for seed 1), so the 20-sphere boxes here do not show the 7 to 8 % grid error of the reference. Tenneti et al. report that their own drag converges from below with resolution and that 20 cells/D is about 7 to 8 % low at \f$\phi\f$ = 0.2 (their Fig. 3, read from the plot), so the correlation itself may sit below a converged solution; this is an argument about the reference, not a measurement of MFC.
- The boxes hold 20 spheres in a 3.74 D cube; the literature arrays hold 80 or more. Realization scatter is large (4 % at Re 20 here).
- Not run: 3D arrays above 1 M cells and the 80-particle boxes (0.8 to 5.8 M cells); the pre-merge drag, Mach and per-particle force tables of those runs are not repeated here.
- **Nusselt number, open.** The hot and cold bodies of the thermostat pair differ by 7 % (6.11 and 5.71 against 5.995 for Gunn). Against the Sun et al. correlation the pair mean is 1.7 times higher (earlier runs: 1.7 to 1.9), which is not explained (the methods, constant against variable properties and the Mach number differ). Treat the Nusselt numbers as plus or minus 15 %. At Mach 0.3 (Brinkman number 0.77) the hot run used to cool its core to 1 K and abort at \f$tW/D\f$ = 8.3 ([#2418](https://github.com/MFlowCode/MFC/issues/2418)); with the dissipation removed it completes \f$tW/D\f$ = 12 with the fluid between 253 K and 322 K. The cell-centre dissipation takes out about 85 % of the pump work at 16 cells/D, so the rest still biases high-Mach runs: use Mach 0.1 for a Nusselt number.
- A closed periodic box with isothermal bodies gains mass slowly (documented in @ref ib_accuracy_limits and the mass-budget example `examples/2D_ibm_wall_mass_budget`).
- A body force on moving particles is not a mean pressure gradient: with `mean_dpdx` the force is not applied inside the bodies, so use the controller for fixed or moving arrays alike.

## Reproduce

```
./mfc.sh run examples/2D_ibm_periodic_array/case.py -n 1 --gpu -- --tstar 40 --stride 20 --wall_order 1 --visc_compact 0
python3 examples/2D_ibm_periodic_array/check_array.py examples/2D_ibm_periodic_array
./mfc.sh run examples/2D_ibm_periodic_array/case.py -n 1 --gpu -- --dim 3 --n_bodies 20 --phi 0.2 --mach 0.3 \
    --Re 100 --cells_per_D 16 --tstar 10 --stride 50 --cfl 0.2 --seed 1 --wall_order 1 --visc_compact 0
```

The suite runs the CI-sized case (12 cylinders, 96 x 96 cells, 200 steps, 2 ranks) as `2D -> Chemistry -> IBM -> Periodic Array`.
