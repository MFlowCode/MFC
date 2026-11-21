@page validation_gas_chemistry Gas-phase chemistry

# Gas-phase chemistry validation

**Status:** 0D ignition and equilibrium: validated. 1D H2/air flame speed: validated at the resolutions run. Reactive shock tube: partially validated (shock structure agrees, reactive features differ by the amount a 0D comparison of the two mechanisms predicts). **Code:** `master` at dd9778b26. **Examples:** `examples/1D_h2_air_ignition`, `examples/1D_h2_air_flame_speed`, `examples/1D_reactive_shocktube`. See also @ref validation.

Gas-phase kinetics, thermodynamics and mixture-averaged transport are checked against Cantera on the same mechanism file (`h2o2.yaml`, 10 species, 29 reactions), so that later results can be attributed to the flow and boundary models and not to the chemistry. Both reaction integrators are covered: the source inside the Runge-Kutta stages (`chem_params%%reaction_substeps = 0`) and the operator-split alpha-QSS integrator (`reaction_substeps > 0`).

All tables were measured with GNU 13.3 (CPU, one rank, double precision) on `master` at dd9778b26.

## Constant-volume ignition {#sec-validation-gas-chemistry-ignition}

| Item | Value |
| :--- | :--- |
| Geometry | 1D periodic box of 26 cells, 1 m long: a uniform state, so no fluxes |
| Gas | stoichiometric H2/air (X: H2 2, O2 1, N2 3.76), 1 atm, at rest, \f$T_0\f$ = 1000, 1100, 1200 K, `h2o2.yaml` |
| Numerics | `model_eqns = 2`, SSP-RK3, fixed \f$\Delta t\f$ = 1e-8 s (split, `reaction_substeps = 2`) or 5e-10 s (in-stage); 0.8 ms, 2000 saves |
| Reference | Cantera `IdealGasReactor` (constant volume, adiabatic), rtol 1e-12, atol 1e-22; `equilibrate("UV")` of the initial state |
| Measured | ignition delay \f$\tau\f$ (time of max \f$dT/dt\f$, parabola through the peak); \f$T(t)\f$; the end state against the UV equilibrium; checked by `compare_cantera.py` |

The error in \f$\tau\f$ compares MFC with Cantera sampled at the same times with the same estimator, which isolates the time-integration error.

| \f$T_0\f$ [K] | \f$\tau\f$ Cantera [µs] | \f$\tau\f$ error, split | max \f$|\Delta T|\f$, split [K] | \f$T_\text{end} - T_\text{eq}\f$ [K] |
| ---: | ---: | ---: | ---: | ---: |
| 1000 | 305.37 | +2.8e-6 | 0.29 | +0.060 |
| 1100 | 86.166 | +1.1e-5 | 0.29 | +0.051 |
| 1200 | 44.227 | +2.0e-5 | 0.28 | +0.043 |

Time-step study, \f$T_0\f$ = 1100 K, split integrator (200 saves, so \f$\tau\f$ errors are against the coarser sampling):

| \f$\Delta t\f$ [s] | \f$\tau\f$ error | max \f$|\Delta T|\f$ [K] | \f$T_\text{end} - T_\text{eq}\f$ [K] |
| ---: | ---: | ---: | ---: |
| 5e-8 | +4.5e-5 | 6.46 | +0.92 |
| 2e-8 | +6.9e-6 | 1.13 | +0.18 |
| 1e-8 | +1.7e-6 | 0.29 | +0.051 |

The split integrator is second order in \f$\Delta t\f$ (max \f$|\Delta T|\f$ falls by 4 per halving). With the in-stage integrator (`reaction_substeps = 0`), \f$T_0\f$ = 1200 K and \f$\Delta t\f$ = 5e-10 s, \f$\tau\f$ agrees to 2e-11 and max \f$|\Delta T|\f$ is below 1e-4 K over 120 µs. Its stability limit is \f$\Delta t < 2.51/|\lambda|_\text{max}\f$ of the reactor Jacobian, about 6e-10 s for burnt H2/air at 1 atm (measured earlier on a pre-merge build: NaN aborts at 1e-9 s and above; not repeated here), so a flow time step above that needs `reaction_substeps > 0`.

Limits: the split integrator renormalizes the mass fractions in the sub-step, so it does not conserve elements exactly (the inert N2 drifts by about \f$-4\times10^{-5}\f$ at \f$\Delta t\f$ = 1e-8 s, measured earlier, not repeated); it shifts the final equilibrium by \f$O(\Delta t^2)\f$ (+0.05 K here). A second mechanism (11 species, CO/moist air) was validated on a pre-merge build only; it is not an example because it would add a chemistry build.

## One-dimensional premixed flame speed {#sec-validation-gas-chemistry-flame}

| Item | Value |
| :--- | :--- |
| Geometry | 1D, \f$x \in [0, 8]\f$ mm; flame started at 6.5 mm and moving toward \f$x = 0\f$ into gas at rest |
| Gas | stoichiometric H2/air, \f$T_u\f$ = 300 K, 1 atm, `h2o2.yaml`, `chem_params%%transport_model = 1` (mixture-averaged), no Soret, `viscous = F` |
| Boundaries | \f$x = 0\f$: reflective wall; \f$x = 8\f$ mm: characteristic subsonic outflow (`bc_x%%end = -8`) with the pressure target `bc_x%%grcbc_out`, `bc_x%%pres_out` = 101325 Pa |
| Numerics | mapped WENO5, HLLC, SSP-RK3, \f$\Delta t = 0.25\,\Delta x/1100\f$ m/s, `reaction_substeps = 2`, run to 1.5 ms |
| Initial condition | Cantera `FreeFlame` profile (`hcid = 170`), shifted to the frame of the unburnt gas |
| Reference | Cantera `FreeFlame`, mixture-averaged, Soret off, 2 cm wide, refine slope = curve = 0.02: \f$S_L\f$ = 2.3317 m/s (2.3316 m/s with 0.01); thermal thickness 329 µm |
| Measured | `compare_flame_speed.py`: \f$S_c = \int\omega_{H_2}dx/(\rho_u(Y_u - Y_\text{out}))\f$ from the H2 mass balance of the domain, and the front speed \f$S_f\f$ from the \f$c = 0.5\f$ isoline, averaged over \f$t \ge 0.75\f$ ms |

The outlet value \f$Y_\text{out}\f$ is used because the domain ends 3 to 5 mm behind the front, inside the slow recombination zone; the Cantera burnt value 2 cm downstream would read 3 % low.

| \f$\Delta x\f$ [µm] | cells per \f$\delta_L\f$ | \f$S_c\f$ [m/s] | error | \f$S_f\f$ [m/s] | error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | 8.2 | 2.3140 ± 0.0010 | -0.76 % | 2.313 ± 0.047 | -0.79 % |
| 20 | 16.5 | 2.3268 ± 0.0008 | -0.21 % | 2.327 ± 0.042 | -0.22 % |
| 10 (earlier\*) | 32.9 | 2.3308 ± 0.0007 | -0.03 % | 2.330 ± 0.023 | -0.08 % |
| 5 (earlier\*) | 65.7 | 2.3318 ± 0.0007 | +0.01 % | 2.331 ± 0.023 | -0.04 % |

The ± is the standard deviation in time. \*The 10 and 5 µm rows were measured earlier on a pre-merge build (same settings and the same Cantera reference within 1e-4) and are not repeated; they are quoted to show the convergence (observed order of \f$S_c\f$ 2.0 between 20, 10 and 5 µm). The 40 and 20 µm rows were run on `master` and agree with the earlier values (-0.74 % and -0.20 %). With `reaction_substeps = 0` at 20 µm the earlier build gave the same \f$S_c\f$ to 0.004 %; not repeated.

The pressure target at the outflow is needed: with extrapolation at both ends the whole domain accelerates, and with a wall and an outflow without a target it pressurizes (unburnt density +13 % by 1 ms, measured earlier on a pre-merge build, trial runs not repeated).

Limits: the domain ends inside the recombination zone (burnt \f$T_\text{max}\f$ 2294 K against 2382 K at 2 cm), so a thermal thickness measured on the MFC profile does not compare with the Cantera value. Not covered: Soret, multicomponent transport, pressures other than 1 atm, the unity-Lewis model, GPU builds.

## Reactive shock tube {#sec-validation-gas-chemistry-shocktube}

A shock reflects from a closed end and the reflected shock ignites the gas; the reaction wave catches the shock and merges with it. Case of Martinez Ferrer et al. \cite MartinezFerrer2014 (section 4.7), `examples/1D_reactive_shocktube`.

| Item | Value |
| :--- | :--- |
| Geometry | 1D, \f$x \in [0, 12]\f$ cm, initial discontinuity at 6 cm |
| Gas | H2/O2/Ar 2/1/7 molar. Left: 0.072 kg/m³, at rest, 7173 Pa. Right: 0.18075 kg/m³, \f$u\f$ = -487.34 m/s, 35594 Pa. Euler equations |
| Boundaries | \f$x = 0\f$: reflective wall; \f$x = 12\f$ cm: extrapolation |
| Numerics | `model_eqns = 2`, mapped WENO5 with `mp_weno`, HLLC, SSP-RK3, \f$\Delta t\f$ = 1.2285e-8 s at 401 cells (acoustic CFL 0.07), chemistry `h2o2.yaml`, `reaction_substeps = 4` as shipped (also 2) |
| Reference | temperature curves at 170, 190 and 230 µs of Fig. 4.9 of the thesis \cite MartinezFerrer2013, which treats the same case (the 2014 paper is paywalled and was not read). The polyline vertices were read from the vector figure with the axes fitted to the tick labels (quantization 0.035 mm in \f$x\f$, 0.9 K in \f$T\f$) and are in `ref_T.txt`. The thesis mechanism is Kee et al. (1993), 9 species and 18 reactions, not `h2o2.yaml`; the curves match the markers of a later replot of the 2014 data by eye (not a fit) |
| Measured | `compare_reference.py`: leading-shock position (rightmost crossing of 975 K), reaction front (1600 K), L1 difference of \f$T\f$ in percent of a 2400 K range, wall temperature, mean offset in the burnt gas (1 to 5 cm at 190 µs, 2 to 8 cm at 230 µs) |

Positions in mm (reference: \f$x_s\f$ = 43.25, 59.81, 99.89; \f$x_r\f$(170 µs) = 36.81; 230 µs front 99.60), `reaction_substeps = 4`:

| Cells | \f$x_s\f$ 170 µs | \f$x_r\f$ 170 µs | \f$x_s\f$ 190 µs | \f$x_s\f$ 230 µs | L1 of \f$T\f$, 170 / 190 / 230 µs [%] | \f$T\f$ wall 170 / 190 / 230 µs [K] (ref. 2459 / 2622 / 2839) | burnt \f$\Delta T\f$ 190 / 230 µs [K] |
| ---: | ---: | ---: | ---: | ---: | :--- | :--- | :--- |
| 401 | 43.16 | 34.17 | 58.76 | 100.62 | 1.01 / 2.50 / 4.16 | 2396 / 2625 / 2837 | +102 / +107 |
| 801 | 43.25 | 34.42 | 58.32 | 100.27 | 1.00 / 2.47 / 3.92 | 2405 / 2625 / 2849 | +86 / +100 |

With `reaction_substeps = 2` the 401-cell run differs by at most 0.02 mm in position and 0.03 in L1 (1.01 / 2.49 / 4.13 %).

The shock position at 170 µs and the induction-zone state agree with the reference; the reaction front trails by 2.5 mm at 170 µs and the burnt gas is 90 to 110 K hotter. A 0D comparison of the two mechanisms from the wall state (`h2o2.yaml` against Kee et al.), made on a pre-merge build with the Kee set, gave +85 K at 20 µs and +96 K at 50 µs after ignition and an ignition time later by 6.8 %, the same sign and size; MFC was not run with the reference mechanism, so the attribution to the mechanism is verified in 0D only. The 230 µs front position moves by -0.35 mm from 401 to 801 cells and is not converged. The reference curves have fewer vertices (107 to 327) than grid points, so narrow peaks may be clipped; only temperature is shipped, and the velocity, \f$\gamma\f$ and \f$Y_H\f$ comparisons of the write-up on a pre-merge build are not repeated. On a pre-merge build, in-stage reactions (`reaction_substeps = 0`) at this \f$\Delta t\f$ ended in NaN at 165 to 175 µs once burnt gas formed, because the step is 4 to 17 times the in-stage stability limit; `master` ships the example with `reaction_substeps = 4`, which runs to 230 µs here (not repeated with 0).

### Reproduce {#sec-validation-gas-chemistry-run}

```
./mfc.sh run examples/1D_h2_air_ignition/case.py -t pre_process simulation post_process -n 1 -- --T0 1100 --nsaves 2000
python3 examples/1D_h2_air_ignition/compare_cantera.py examples/1D_h2_air_ignition --case-args "--T0 1100 --nsaves 2000"
./mfc.sh run examples/1D_h2_air_flame_speed/case.py -t pre_process simulation post_process -n 1 -- --dx_um 20
python3 examples/1D_h2_air_flame_speed/compare_flame_speed.py examples/1D_h2_air_flame_speed --case-args "--dx_um 20"
./mfc.sh run examples/1D_reactive_shocktube/case.py -t pre_process simulation post_process -n 1 -- --t_save_us 5
python3 examples/1D_reactive_shocktube/compare_reference.py examples/1D_reactive_shocktube
```

The checkers need Cantera and the Silo output. The suite runs the first test as `1D -> Chemistry -> H2-Air Ignition -> Split Reactions` (\f$T_0\f$ = 1200 K, \f$\Delta t\f$ = 2e-8 s, 100 µs); the flame and the shock-tube checks are not registered because of their length.

The golden of the registered test was generated with GNU 13 (CPU, one rank, double). An NVIDIA A100 with NVHPC 25.11 agrees with it to 2e-12 relative on every value above 1e-9. NVHPC 25.5 on the same A100 differs by up to 1e-3 relative in the minor species, growing through the induction period (first visible at 14 to 30 µs) as the radical pool grows, while \f$\rho\f$ and the total energy agree. A CPU build with NVHPC 25.5 and a GPU build with NVHPC 25.11 (CUDA 13.0 or 12.9) reproduce the GNU result bit for bit in the printed digits, and the 25.5 GPU difference does not change with `-Mfprelaxed` removed, `-Mnofma`/`-gpu=nofma` or `-O1`, so it comes from the 25.5 GPU compilation or runtime and not from those flags; the cause is not identified. The test is therefore skipped on NVHPC < 25.11 GPU builds (issue #2181); CPU builds run it. This is not chaotic ignition: the box is uniform and the run is smooth.
