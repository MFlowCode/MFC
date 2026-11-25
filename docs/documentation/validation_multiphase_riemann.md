@page validation_multiphase_riemann Validation: Two-Fluid Shock Tube, Interface Advection and Phase-Change Equilibrium

# Two-fluid flows: shock tube, interface advection, phase-change equilibrium

**Status:** water-air shock tube: validated. Interface advection: partially validated (open: a localized pressure error on coarse grids). Phase-change equilibrium: partially validated (equilibrium of the relaxed state only). **Code:** `master` at 5838b641c. **Examples:** `examples/1D_kapilashocktube`, `examples/2D_advection`, `examples/1D_exp_tube_phasechange`. **Last run:** 2026-11-24. See also @ref validation.

All tables were measured with GNU 13.3 (CPU, one rank, double precision); the A100 rows used NVHPC 25.11 (OpenACC, one GPU).

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
| Numerics | 6-equation model with `mpp_lim`, WENO5-M, HLLC, SSP-RK3, acoustic CFL 0.4; one flow-through along the diagonal (t = 0.01 s, the disk returns to its initial position) |
| Reference | exact: p and u stay uniform, the volume fraction translates; after one flow-through it equals the initial field |
| Measured | max abs(p/p0 - 1) and max abs(u - U)/U over all saves and cells; L1 of alpha_water against the initial field per unit interface length, in cells; the width of the mixed zone (0.01 < alpha < 0.99) per unit interface length, in cells; `check_advection.py` |

| N | max dp/p0 | max du/U | L1 [cells] | L1 [m] | width, initial [cells] | width, final [cells] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 2.5e-3 | 1.7e-4 | 0.066 | 1.3e-3 | 9.17 | 9.25 |
| 64 | 2.5e-3 | 7.2e-7 | 0.075 | 1.2e-3 | 9.21 | 9.26 |
| 80 | 1.2e-4 | 5.0e-8 | 0.084 | 1.0e-3 | 9.26 | 9.24 |
| 100 | 3.1e-6 | 1.9e-9 | 0.095 | 9.5e-4 | 9.18 | 9.26 |
| 200 | 3.2e-6 | 2.2e-9 | 0.142 | 7.1e-4 | 9.16 | 9.54 |
| 50, A100 | 2.5e-3 | 1.7e-4 | 0.066 | 1.3e-3 | 9.17 | 9.25 |
| 100, A100 | 2.9e-6 | 1.8e-9 | 0.095 | 9.5e-4 | 9.18 | 9.26 |

- At 100 and 200 cells the uniform pressure and velocity hold to 3e-6 and 2e-9 across a density ratio of 1000 and a pi_inf
  contrast of 1e4 p0 over the whole flow-through: the scheme keeps the interface equilibrium.
- The interface is not diffused beyond its initial smoothing: the mixed zone widens by 0.1 cells or less per flow-through up to
  100 cells (0.4 cells at 200, where the disk travels 200 cells in each direction). The L1 error per unit interface length is
  about a tenth of a cell and falls as N^-0.4 in physical units; the initial profile is a fixed number of cells wide, so this is
  not the convergence of one solution.
- **Open: localized pressure error on coarse grids.** At 50 and 64 cells the pressure error stays below 3e-6 for 0.4
  flow-throughs and then jumps to 2.5e-3 p0 (250 Pa) in a patch of about 30 water cells that hold 1e-6 to 1e-4 of air, and
  travels with the disk. The onset coincides with a cell where a volume fraction undershoots to zero (6-equation) or below zero
  (the same case with the 5-equation model, where the error reaches 2.5e-2 p0 at 50 cells and stays at 2e-11 at 100). It does not
  change with the CFL number (0.2), the trace fraction (1e-6 to 1e-10), `mp_weno` (1.7e-3) or a disk path that avoids the box
  corner, and the A100 reproduces it to the digits shown. Relative to the water stiffness p0 + pi_inf the error is 2.5e-7. Cause
  not identified.

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
