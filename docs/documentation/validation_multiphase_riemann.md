@page validation_multiphase_riemann Validation: Two-Fluid Shock Tube, Interface Advection and Phase-Change Equilibrium

# Two-fluid flows: shock tube, interface advection, phase-change equilibrium

**Status:** water-air shock tube: validated. Interface advection: validated (pressure equilibrium to round-off at every grid). Phase-change equilibrium: partially validated (equilibrium of the relaxed state only). **Code:** `master` at 5838b641c (interface advection: 3e4c53258). **Examples:** `examples/1D_kapilashocktube`, `examples/2D_advection`, `examples/1D_exp_tube_phasechange`. **Last run:** 2026-11-24. See also @ref validation.

Except for the interface-advection table, all tables were measured with GNU 13.3 (CPU, one rank, double precision) and the A100 rows with NVHPC 25.11 (OpenACC, one GPU); the interface-advection table states its own compilers.

## Water-air shock tube {#sec-validation-mp-riemann}

| Item | Value |
| :--- | :--- |
| Geometry | 1D, x in [0, 1] m, interface at 0.7 m, transmissive ends |
| Fluids | left: water, stiffened gas gamma 4.4, pi_inf 6e8 Pa, 1000 kg/m^3, 1e9 Pa; right: air, gamma 1.4, 50 kg/m^3, 1e5 Pa; both at rest, pure (volume fractions 0 and 1) |
| Numerics | 5-equation model \cite Kapila01 with `alt_soundspeed`, WENO5-M, HLLC, SSP-RK3, dt = 4e-8 s at N = 1000 (CFL about 0.2), scaled with 1/N; t = 241 us |
| Reference | exact Riemann solution of two stiffened gases (Toro's iterative solver \cite Toro09 with p + pi_inf in each material), sampled at the cell centres |
| Measured | L1 of rho, u, p over the tube (normalized by the left-right jump of each), the alpha_1 = 1/2 crossing against the contact, the mean over the middle fifth of each star region; `check_riemann.py` |

Exact star state: p* = 1.4190e7 Pa, u* = 482.61 m/s, rho* = 804.44 kg/m^3 (water) and 288.17 kg/m^3 (air).

| N | L1 rho | L1 u | L1 p | order (rho, u, p) | interface [cells] | air plateau rho, p |
| ---: | ---: | ---: | ---: | :--- | ---: | :--- |
| 250 | 3.69e-3 | 7.98e-3 | 2.08e-3 | | +0.02 | +4.8e-2, +1.2e-3 |
| 500 | 2.32e-3 | 4.52e-3 | 1.06e-3 | 0.67, 0.82, 0.98 | +0.01 | -5.5e-3, -4.0e-3 |
| 1000 | 1.11e-3 | 1.61e-3 | 5.28e-4 | 1.06, 1.48, 1.00 | +0.01 | -1.8e-2, +9.1e-4 |
| 2000 | 7.18e-4 | 1.16e-3 | 2.69e-4 | 0.63, 0.48, 0.97 | +0.01 | -1.8e-3, -2.3e-4 |
| 1000, A100 | 1.11e-3 | 1.61e-3 | 5.28e-4 | | +0.01 | -1.8e-2, +9.1e-4 |

- The solution converges to the exact one: the pressure error is first order, as expected with a shock, and the material
  interface sits within 0.02 cells of the exact contact on every grid. Density and velocity converge irregularly (order 0.5 to 1.5
  between grids) because their error is concentrated at the shock and the contact.
- The plateau errors are the post-shock oscillation of the air state (up to 5 % in density at N = 250, 0.2 % at N = 2000); in the
  water plateau p is within 0.5 % and u within 1e-4 on all grids.

## Interface advection {#sec-validation-mp-advection}

| Item | Value |
| :--- | :--- |
| Geometry | 2D periodic unit box, N x N cells; water disk of radius 0.2 centred at (0.4627, 0.5381), off the cell faces and the box centre |
| Fluids | water (stiffened gas, gamma 2.35, pi_inf 1e9 Pa, 1000 kg/m^3) in air (gamma 1.4, 1 kg/m^3); 1e-8 of the other fluid in each pure region; initial interface smoothed (`smooth_coeff` 0.5, about 9 cells wide at every N) |
| State | p = 1e5 Pa, u = v = 100 m/s everywhere |
| Numerics | 6-equation model, `mpp_lim = F`, WENO5-M, HLLC, SSP-RK3, acoustic CFL 0.4; one flow-through along the diagonal (t = 0.01 s, the disk returns to its initial position) |
| Reference | exact: p and u stay uniform, the volume fraction translates; after one flow-through it equals the initial field |
| Measured | max abs(p/p0 - 1) and max abs(u - U)/U over all saves and cells; L1 of alpha_water against the initial field per unit interface length, in cells; the width of the mixed zone (0.01 < alpha < 0.99) per unit interface length, in cells; `check_advection.py` |

| N | max dp/p0 | max du/U | L1 [cells] | L1 [m] | width, initial [cells] | width, final [cells] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 3.7e-11 | 1.1e-11 | 0.066 | 1.3e-3 | 9.17 | 9.25 |
| 64 | 2.5e-11 | 1.4e-12 | 0.075 | 1.2e-3 | 9.21 | 9.26 |
| 80 | 2.4e-11 | 2.1e-12 | 0.084 | 1.0e-3 | 9.26 | 9.24 |
| 100 | 2.5e-11 | 5.6e-12 | 0.095 | 9.5e-4 | 9.18 | 9.26 |
| 200 | 1.5e-10 | 1.9e-10 | 0.142 | 7.1e-4 | 9.16 | 9.54 |
| 50, `model_eqns = 2` | 1.9e-11 | 9.4e-14 | 0.066 | 1.3e-3 | 9.17 | 9.25 |
| 100, `model_eqns = 2` | 2.9e-11 | 3.9e-13 | 0.095 | 9.5e-4 | 9.18 | 9.26 |
| 50, A100 | 3.8e-11 | 1.1e-11 | 0.066 | 1.3e-3 | 9.17 | 9.25 |
| 100, A100 | 3.4e-11 | 2.1e-12 | 0.095 | 9.5e-4 | 9.18 | 9.26 |

This table was measured at 3e4c53258 on the CPU with NVHPC 25.11 (four ranks) and on the A100 with NVHPC 25.5 (OpenACC, one GPU).

- The uniform pressure and velocity hold to round-off (1e-11 to 1e-10) at every grid, across a density ratio of 1000 and a
  pi_inf contrast of 1e4 p0, over the whole flow-through and with either model: the scheme keeps the interface equilibrium.
- The interface is not diffused beyond its initial smoothing: the mixed zone widens by 0.1 cells or less per flow-through up to
  100 cells (0.4 cells at 200, where the disk travels 200 cells in each direction). The L1 error per unit interface length is
  about a tenth of a cell and falls as N^-0.4 in physical units; the initial profile is a fixed number of cells wide, so this is
  not the convergence of one solution.
- **`mpp_lim` must stay off here** ([#2476](https://github.com/MFlowCode/MFC/issues/2476)). With `mpp_lim = T`, which the deck
  used before, the pressure error is 3e-6 at 100 and 200 cells. At 50 and 64 cells it jumps after 0.4 flow-throughs to 2.5e-3 p0
  (250 Pa; 2.7e-2 p0 with `model_eqns = 2`) in a patch of water cells that hold traces of air. `mpp_lim` clips each volume
  fraction to [0, 1] and renormalizes the set before it forms the mixture \f$\Gamma\f$ and \f$\Pi_\infty\f$, in the
  conversion to primitives and on the Riemann-solver face states. The energy is not changed with them, so where alpha_air
  undershoots by \f$\delta\f$ the clip moves the water pressure by about \f$\delta\,\gamma\pi_\infty\f$ = 2.35e9
  \f$\delta\f$ Pa: \f$\delta\f$ = 1e-7 gives the 250 Pa. Without the clip the volume fractions do undershoot at 50 cells
  (min alpha_air = −2.7e-5, min alpha_water = −3.4e-7, partial densities down to −3.4e-4 kg/m³; none at 100 and 200 cells),
  but the equation of state stays consistent with the advected alpha and energy, and the run is stable. `face_bound = 1` (bounds
  on the face density and p + pi_inf, `model_eqns = 2`) does not change the result with or without `mpp_lim`.

## Phase-change equilibrium {#sec-validation-mp-phase-change}

| Item | Value |
| :--- | :--- |
| Geometry | 1D, x in [0, 1] m, 201 cells, transmissive ends |
| Fluids | liquid water, water vapor, air (stiffened gas with `cv`, `qv`, `qvp` of \cite LeMetayer04); volume fractions 0.850, 0.029, 0.121 |
| State | 1.175e5 Pa and 378.146 K, on the saturation curve of the parameters; u = -2 m/s for x < 0.5, +2 m/s for x > 0.5 |
| Numerics | 6-equation (`-me 3`) or 5-equation (`-me 2`) model, `relax_model = 6` (infinitely fast p-T-Gibbs relaxation \cite Saurel08 after each step), WENO3-M, HLLC, SSP-RK3; t = 511 us |
| Reference | equilibrium conditions: T of each phase from its own EOS equal, and T = T_sat(p) where liquid and vapor coexist, T_sat the root of g_liquid(p, T) = g_vapor(p, T) |
| Measured | in every cell of the save, `check_saturation.py` |

| Model | Save | max abs(T_k/T_vapor - 1) | max abs(T - T_sat(p)) [K] | max abs(sum alpha - 1) | vapor mass |
| :--- | :--- | ---: | ---: | ---: | ---: |
| 6-equation | initial | 7.2e-10 | 2.5e-7 | 5.7e-14 | 2.0413e-2 kg/m^2 |
| 6-equation | last | 9.4e-7 | 7.0e-4 | 8.0e-11 | 2.1792e-2 kg/m^2 (+6.76 %) |
| 5-equation | last | 9.4e-7 | 6.2e-4 | 8.0e-11 | 2.1792e-2 kg/m^2 (+6.76 %) |
| 6-equation, A100 | last | 9.4e-7 | 7.0e-4 | 8.0e-11 | 2.1792e-2 kg/m^2 |

- The relaxed state is in thermal and chemical equilibrium: every cell lies within 7e-4 K of the saturation curve after the
  expansion has evaporated 6.8 % more vapor, with all phases at one temperature to 1e-6.
- The parameters give T_sat(101325 Pa) = 373.79 K against 373.12 K measured (a property of the fit, not of the solver).

Limits: the expansion is weak (the pressure falls by 420 Pa, 0.04 %), and only the equilibrium of the state is checked, not the
wave speeds or the evaporated mass against a reference solution. Metastable and flashing cases are not covered.

## Reproduce {#sec-validation-mp-run}

```
./mfc.sh run examples/1D_kapilashocktube/case.py -t pre_process simulation post_process -n 1 -- -N 1000
python3 examples/1D_kapilashocktube/check_riemann.py examples/1D_kapilashocktube --case-args "-N 1000"
./mfc.sh run examples/2D_advection/case.py -t pre_process simulation post_process -n 1 -- -N 100
python3 examples/2D_advection/check_advection.py examples/2D_advection --case-args "-N 100"
./mfc.sh run examples/1D_exp_tube_phasechange/case.py -t pre_process simulation post_process -n 1 -- -me 3
python3 examples/1D_exp_tube_phasechange/check_saturation.py examples/1D_exp_tube_phasechange --case-args "-me 3"
```

The checkers take several run directories (and `--case-args` for each) to print the observed order. The shock tube and the expansion
tube run in the test suite at a reduced size, as every registered example does; `2D_advection` is on the skip list. No test was added.
