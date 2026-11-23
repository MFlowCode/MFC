@page validation Validation

# Validation

This page records what MFC has been checked against: exact solutions, reference data and conservation laws.
Each section is one family of checks. It gives the status, the code version, the setup, the reference, the measured values with signed errors, and the limits.
Every check has an example in `examples/` with a checker script and a README, and a small version of it runs in the test suite (`./mfc.sh test`).

Status words: **validated** (agrees with the reference within the stated error), **partially validated** (agrees where measured; a regime or quantity is open), **open discrepancy** (disagrees, cause known or not), **comparison only** (reference data are qualitative).
Errors are (MFC − reference)/reference.

Other pages: @ref validation_gas_chemistry "Gas-phase chemistry" (ignition, flame speed, reactive shock tube).

## Moving-body kinematics {#sec-validation-moving-kinematics}

**Status:** validated. **Code:** `master` at 067f23074, with the stage-consistent position update of [#2040](https://github.com/MFlowCode/MFC/pull/2040). **Examples:** `examples/2D_mibm_constant_acceleration`, `examples/3D_mibm_wall_collision`.

Two-way coupled bodies (`moving_ibm = 2`) are advanced from the force on them. Particle-resolved runs depend on the integration being exact where an exact answer exists, and on collisions returning the prescribed restitution.

### Constant acceleration {#sec-validation-moving-kinematics-accel}

| Item | Value |
| :--- | :--- |
| Body | circle, R = 50 µm, 1000 kg/m, `moving_ibm = 2` |
| Gas | air at rest, 1e5 Pa, 1.2 kg/m³, inviscid; the gas force is negligible against the body force |
| Drive | `bf_x = T`, `g_x = 1e5` m/s², `k_x = w_x = p_x = 0` |
| Numerics | first-order reconstruction, HLLC, SSP-RK3, fixed \f$\Delta t\f$ at acoustic CFL 0.4, 200 steps, 6R half-width, extrapolation boundaries |
| Reference | exact kinematics: \f$x = g t^2/2\f$, \f$v = g t\f$ |
| Measured | final centroid and velocity from `restart_data/ib_state_*.dat`, checked by `check_kinematics.py` |

| Cells per radius | Position error | Velocity error |
| ---: | ---: | ---: |
| 5 | −3.6e-13 | −2.3e-13 |
| 10 | −2.9e-13 | −1.9e-13 |

Both are round-off. A build with the stage weights of the position update from before #2040 gave a position error of 1.2e-3 with an exact velocity (measured earlier on a build without #2040; not repeated). The check cannot see a wrong stage velocity in the velocity, only in the position, which is why it uses a position.

Limits: only constant acceleration is checked. For a prescribed body (`moving_ibm = 1`) with an analytic velocity such as \f$U_0\sin\omega t\f$, the velocity is evaluated at the step time rather than the stage time, and the centroid is off by about 1e-3 relative (measured earlier; the cause is not verified and is not addressed here).

### Sphere-wall collision {#sec-validation-moving-kinematics-collision}

| Item | Value |
| :--- | :--- |
| Body | sphere, R = 5 m, mass 1e6 (gas force negligible), impact speed 3.9 m/s |
| Wall | `y = 0`, `bc_y%%beg = -15`; other boundaries extrapolation |
| Model | `collision_model = 1`, `collision_time` = N\f$\Delta t\f$, \f$\Delta t\f$ = 5e-6 s, `coefficient_of_restitution` e, `ib_coefficient_of_friction` \f$\mu\f$ |
| Numerics | first-order reconstruction, HLLC, SSP-RK3, fixed \f$\Delta t\f$, 16³ cells, air at 1e5 Pa (the gas plays no role) |
| Reference | rigid-body impact of a sphere on a fixed wall with Coulomb friction and no tangential elasticity (\cite Maw1976): \f$v_n' = -e v_n\f$; sliding throughout (\f$\tan\theta \ge 3.5\mu(1+e)\f$): \f$v_t' = v_t - \mu(1+e)|v_n|\f$; sticking (\f$\tan\theta < 3.5\mu(1+e)\f$): \f$v_t' = \frac{5}{7}v_t\f$. Setup after \cite DiMaio2004 |
| Measured | \f$e_n = -v_y^\text{out}/v_y^\text{in}\f$ and \f$v_x^\text{out}\f$ from the last `ib_state` record, by `check_collision.py` |

| Case | Prescribed | Measured | Error |
| :--- | ---: | ---: | ---: |
| normal, e = 0.5 | 0.5000 | 0.4765 | −4.70 % |
| normal, e = 0.8 | 0.8000 | 0.7877 | −1.54 % |
| normal, e = 0.98 | 0.9800 | 0.97948 | −0.05 % |
| 30°, \f$\mu\f$ = 0, e = 0.98: \f$v_x^\text{out}\f$ | 1.9500 | 1.94998 | 0.00 % |
| 30°, \f$\mu\f$ = 0.092, e = 0.98: \f$v_x^\text{out}\f$ (sticking) | 1.3929 | 1.3958 | +0.21 % |
| 50°, \f$\mu\f$ = 0.092, e = 0.98: \f$v_x^\text{out}\f$ (sliding) | 2.5309 | 2.5310 | +0.00 % |

Error in \f$e_n\f$ against the steps per collision:

| e | 20 | 40 | 80 |
| ---: | ---: | ---: | ---: |
| 0.5 | −4.70 % | −2.33 % | −1.16 % |
| 0.8 | −1.54 % | not run | −0.37 % |

The restitution error is first order in the steps per collision. At the default of 20 steps, expect a few percent in e at low e. Its size and sign depend on where in a step the contact begins: an earlier run with the contact beginning on a step boundary gave +2.17 % and +0.69 % for e = 0.5 and 0.8 at 20 steps. Without friction the tangential velocity is conserved exactly. The sticking and sliding regimes both match the rigid-body result.

Limits: one sphere against a wall, three angles. Body-body collisions, spin and other speeds are not covered. The first-order behavior is in the collision force; whether it is the force evaluation or the contact detection is not separated.

### Reproduce {#sec-validation-moving-kinematics-run}

```
./mfc.sh run examples/2D_mibm_constant_acceleration/case.py -t pre_process simulation -n 1
python3 examples/2D_mibm_constant_acceleration/check_kinematics.py examples/2D_mibm_constant_acceleration
./mfc.sh run examples/3D_mibm_wall_collision/case.py -t pre_process simulation -n 1 -- --angle 30 --e 0.98 --mu 0.092
python3 examples/3D_mibm_wall_collision/check_collision.py examples/3D_mibm_wall_collision --angle 30 --e 0.98 --mu 0.092
```

The small cases run in the test suite as `2D -> IBM -> Moving Circle -> Constant Acceleration` and `3D -> IBM -> Wall Collision -> Normal` and `-> Sliding Friction`.
The tables above were measured with GNU 13.3 (CPU, one rank). The golden tolerance is 1e-6 for both. Against goldens generated on the current master, the GNU, Intel, NVHPC and macOS CI lanes agreed to 1e-13 absolute on the fields and 5e-9 relative on the contact force; the checker scripts validate the physics, and the golden is a regression guard. Both tests are skipped in the single-precision lane. The sticking case is checked by `check_collision.py` and is not a golden, because the friction direction there depends on round-off.

## Thermally developed channel {#sec-validation-mean-forcing-channel}

**Status:** validated (plane channel; the same source in particle arrays is not). **Code:** `master` at 1e3f1be0b, with the thermally developed source `mean_temp_wall` of #2200. **Example:** `examples/2D_thermal_channel`.

Periodic flow with a heated or cooled wall never reaches a uniform temperature, so the mean-temperature source (`mean_temp`, `mean_temp_wall`) removes the streamwise decay: it solves for the periodic part \f$\psi\f$ of the thermally developed state \f$T = T_s + \psi\,e^{-\beta x}\f$ \cite Sun2015 and holds the mixing-cup temperature at its target. A plane channel with isothermal walls has an exact developed Nusselt number, so it tests the source, its decay rate \f$\beta\f$ and the wall heat transfer together.

| Item | Value |
| :--- | :--- |
| Geometry | 2D plane channel, gap \f$H\f$ = 34.28 µm (\f$D_h = 2H\f$), periodic length 25 cells |
| Gas | air (O2 0.233, N2 0.767 by mass, `h2o2.yaml`, reactions off), 1 atm, mixing-cup \f$T_m\f$ held at 300 K, bulk Mach 0.1 (\f$U_b\f$ = 34.8 m/s), \f$\mathrm{Re}_{D_h}\f$ = 150, \f$\mathrm{Pe}_H\f$ = 53.3 |
| Boundaries | `x` periodic; `y` no-slip isothermal walls (`bc_y%%beg = bc_y%%end = -16`) at \f$T_s\f$ = 330 K (hot) or 270 K (cold) |
| Numerics | WENO5 (mapped), HLLC, SSP-RK3, `fd_order = 4`, mixture-averaged transport, fixed \f$\Delta t\f$ = 0.4 \f$\Delta y/(U_b + a)\f$ |
| Drive and source | `mean_mflux` and `mean_dpdx` (gain 0.05); `mean_temp` = 300 K (gain 0.05) with `mean_temp_wall` = \f$T_s\f$ |
| Reference | constant-property laminar flow between parallel plates, constant wall temperature, fully developed: Nu = 7.5407 \cite ShahLondon1978. With axial conduction at the run's \f$\mathrm{Pe}_H\f$, the eigenproblem \f$\psi'' + (6\,\mathrm{Pe}_H\,\eta(1-\eta)\,b + b^2)\psi = 0\f$, \f$\psi(0) = \psi(1) = 0\f$, \f$b = \beta H\f$, gives Nu = 7.5417 and \f$\beta H\f$ = 0.14122; `check_nu.py` solves it |
| Measured | `mean_temp.dat` (\f$T_m\f$, \f$\dot q\f$, \f$V_f\f$, \f$\beta\f$) and the last save, by `check_nu.py` |

\f$\mathrm{Nu}_\mathrm{src} = q_w D_h/(k(T_s - T_m))\f$ takes the wall heat from the source total, \f$q_w = -\dot q V_f/(2L_x)\f$, which equals it at steady state. \f$\mathrm{Nu}_\mathrm{grad}\f$ takes it from the temperature gradient at the walls (second-order one-sided) and does not use the source. The pair value uses half-differences between the hot and the cold run, which cancels viscous heating and the first-order property asymmetry. Averages are over the last 10 % of the steps; the state is developed by about 0.1 \f$H^2/\alpha\f$.

| Cells across \f$H\f$ | Run [\f$H^2/\alpha\f$] | \f$T_m\f$ hold (max error) | \f$\beta H\f$ hot / cold (exact 0.1412) | \f$\mathrm{Nu}_\mathrm{src}\f$ pair | \f$\mathrm{Nu}_\mathrm{grad}\f$ pair |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 0.3 | 2e-8 | 0.1457 / 0.1371 | 7.552 (+0.13 %) | 7.554 (+0.16 %) |
| 64 | 0.15 | 7e-9 | 0.1456 / 0.1371 | 7.546 (+0.05 %) | 7.546 (+0.05 %) |

- A single run is biased by viscous heating and by the property change across the gap, with opposite sign for the hot and the cold wall (32 cells: \f$\mathrm{Nu}_\mathrm{src}\f$ 7.507 (−0.46 %) hot, 7.600 (+0.77 %) cold; \f$\mathrm{Nu}_\mathrm{grad}\f$ 7.141 (−5.3 %) and 7.998 (+6.1 %), which includes the dissipation that the source removes). The pair removes both.
- \f$\beta H\f$ of the pair mean is 0.1414 at 32 cells (+0.15 % against 0.1412) and 0.1413 at 64 (+0.06 %); the hot and the cold value differ from it by the property effect on the Peclet number.
- The source total and the wall gradient agree within 0.1 % after pairing: the source removes the heat that crosses the walls, which closes the fluid energy budget.
- The mixing-cup temperature is held to 2e-8 (32 cells) and 7e-9 (64 cells) of its target.

Limits:
- **Plane channel only.** In periodic arrays of spheres \f$\beta\f$ drifted by more than 1 % and the hot-wall and cold-wall runs disagreed (12 cells per diameter, cause not understood); do not take an array Nusselt number from this source (note under `mean_temp_wall` in the case-file documentation). Viscous heating also biases a single-wall-temperature run at high Mach number.
- **Near-wall velocity sawtooth (open).** With the no-slip isothermal wall (`bc_y%%beg = -16`) the streamwise velocity in the first four cells next to each wall alternates about the Poiseuille profile by up to 12 % of \f$U_b\f$ (32 cells, hot wall: \f$u/U_b\f$ = 0.195, 0.147, 0.459, 0.543 against 0.092, 0.268, 0.432, 0.584). It is steady, uniform in \f$x\f$ and symmetric, and the temperature is smooth; the Nusselt numbers above are not affected at the 0.2 % level. The cause is not verified.
- **Round-off sensitivity.** With the default mapped WENO5 and `weno_eps = 1e-16`, the transverse momentum in the wall rows differs by 2 to 3 % between a GNU 13 `-O3` build and a debug build after 20 to 200 steps. The registered test therefore uses non-mapped WENO3 with `weno_eps = 1e-6`, where the builds agree to 5e-8 (relative), and an AMD MI250X (Cray Fortran, OpenMP offload, one GCD) run agrees with the GNU CPU golden to 8e-8.

### Reproduce {#sec-validation-mean-forcing-channel-run}

```
./mfc.sh run examples/2D_thermal_channel/case.py -t pre_process simulation -n 1 -- --ny 32 --nx 25 --tw 330 --tend 0.3
python3 examples/2D_thermal_channel/check_nu.py <hot run directory> <cold run directory>
```

The cold run uses `--tw 270` in a copy of the directory (the output is written next to `case.py`). The tables were measured with GNU 13.3 (CPU, one rank, `-O3`). The suite runs the first 100 steps of the hot case as `2D -> Chemistry -> Thermal Channel -> Wall Temperature`; its golden holds the conservative fields, and the tolerance is 1e-6 (GNU 13 double -O3, reldebug and debug agree to 5e-8, AMD MI250X and an AMD CPU build to 8e-8 and 6e-8; single precision differs from double by 4e-4, so the test is skipped in that lane).

## Circular cylinder drag and Nusselt number {#sec-validation-cyl-drag-nu}

**Status:** partially validated. **Example:** `examples/2D_ibm_cylinder_drag_nu`. The face-flux immersed-boundary force and the face heat of a no-slip isothermal cylinder at Re 20, 40 and 100 against Dennis and Chang, Williamson and the Nusselt-number correlations, with a momentum balance as an independent drag. Results and limits: @ref validation_cyl_drag_nu.

## Moving-body gas-mass budget and Galilean twin {#sec-validation-moving-body-budget}

**Status:** partially validated. **Example:** `examples/2D_ibm_moving_mass_budget`. The gas mass a translating circle creates in a closed periodic box, with and without `ib_mass_relax`, and the conducted heat and force of the moving body against the same body held in the same relative flow. Results and limits: @ref validation_moving_body_budget.

## Sphere drag and Nusselt number {#sec-validation-sphere-drag-nu}

**Status:** partially validated. **Example:** `examples/3D_ibm_sphere_drag_nu`. The face-flux drag and face heat of a no-slip isothermal sphere at Mach 0.2 and Re 20-300 against Schiller-Naumann, Ranz-Marshall and Whitaker, with a momentum balance as an independent drag. Results and limits: @ref validation_sphere_drag_nu.

## Bow-shock stand-off {#sec-validation-standoff}

**Status:** open discrepancy (resolution-limited). **Example:** `examples/2D_ibm_cylinder_standoff`. The detached-shock distance of a cylinder and a sphere at Mach 1.5-3 against Billig's correlation and a body-fitted Euler solution. Results and limits: @ref validation_standoff.

## Other families

- @ref validation_tsc_annulus "Isothermal conduction annulus": image-point interpolation (inverse distance and TSC) against the exact radial conduction solution.

## Adding a section

One `##` section per family, with an anchor `{#sec-validation-<family>}`: status and code version, one setup table per check, the reference and how it was obtained, a results table with signed errors, and the limits.
Put each check's example, checker script and README under `examples/`, and register a small version in `toolchain/mfc/test/cases.py`.
Cite references with `\cite` keys from `docs/references.bib`.

## Keeping the goldens current

The registered tests are regression guards: their goldens hold the output of the code at the time they were generated. When a default or an algorithm in the numerics a test exercises changes (for moving bodies, for example the `ib_mass_relax` default), regenerate the goldens with `./mfc.sh test --generate --only <UUID>` in the same change, and re-run the checker scripts to confirm the physics results above. A case whose tolerance must be looser than the default takes `override_tol` chosen from the spread between compiler and device lanes ([Testing](@ref testing)).
