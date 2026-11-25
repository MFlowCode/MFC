@page validation_surface_tension Validation: Surface Tension

# Surface tension: Laplace pressure jump of a static drop

**Status:** validated (static 2D drop). **Code:** `master` at 5838b641c. **Example:** `examples/2D_laplace_pressure_jump`. **Last run:** 2026-11-24. See also @ref validation.

## Purpose

Surface tension (`surface_tension`, the capillary model of \cite Schmidmayer17 with a color function) must balance the pressure
jump across a curved interface at rest. A 2D drop in equilibrium has the exact jump \f$\Delta p = \sigma/R\f$; the residual flow
(spurious currents) measures how well the discrete capillary force and the pressure gradient balance.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | quarter of a 2D drop centred on the corner of [0, 0.375]^2 m (symmetry planes at x = 0 and y = 0, extrapolation on the far faces); r0 = 0.15 m; N x N cells |
| Fluids | water (stiffened gas, gamma 2.1, pi_inf 1e6 Pa, 1000 kg/m^3) in air (gamma 1.4, 1 kg/m^3), 1e5 Pa, at rest; 1e-9 of the other fluid in each pure region; interface smoothed (`smooth_coeff` 0.95) |
| Surface tension | sigma = 8 N/m, so sigma/r0 = 53.33 Pa (5.3e-4 of p0) |
| Numerics | 6-equation model, WENO5-M with `mp_weno` and `weno_avg`, HLLC, SSP-RK3, dt = 5e-6 s x 50/N; t = 0.5 s (about 60 acoustic crossings of the drop) |
| Reference | \f$\Delta p = \sigma/R\f$, R the equivalent radius of the water area, \f$\sqrt{4A/\pi}\f$ with \f$A = \sum \alpha_\mathrm{water}\,\Delta x\Delta y\f$ |
| Measured | mean p where alpha_water > 0.999 minus mean p where alpha_water < 0.001 and r > 1.5 R, averaged over the last 20 saves; the largest speed in the window; `check_laplace.py` |

## Results

GNU 13.3 (CPU, one rank) unless noted; the capillary speed \f$\sqrt{\sigma/(\rho_l R)}\f$ is 0.23 m/s.

| N | R [m] | dp [Pa] | vs sigma/R | vs sigma/r0 | spread of dp | max abs(u) [m/s] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 25 | 0.15068 | 53.633 | +1.02 % | +0.56 % | 1.3e-3 | 0.053 |
| 50 | 0.15017 | 53.279 | +0.010 % | -0.10 % | 2.6e-4 | 0.038 |
| 50, A100 | 0.15017 | 53.280 | +0.013 % | -0.10 % | 3.5e-4 | 0.037 |
| 100, A100 | 0.15004 | 53.185 | -0.25 % | -0.28 % | 8.9e-3 | 0.042 |

## Assessment

- The pressure jump matches \f$\sigma/R\f$ to 1.0 % at 25 cells, 0.01 % at 50 and 0.25 % at 100 cells across the box (10, 20 and 40
  cells per drop radius). The error does not decrease monotonically past 50 cells: at 100 cells the jump still varies by 0.9 % over
  the averaging window, which is larger than its error, so the 50 and 100 cell values agree within that unsteadiness.
- R moves by 0.4 % from the initial radius at 25 cells and by 0.03 % at 100 (smoothing of the initial interface); measured against
  the initial radius the error is below 0.6 % on every grid.
- The A100 and CPU runs at 50 cells agree to 3e-5 in the jump.
- The spurious currents are 0.04 to 0.05 m/s, 16 to 23 % of the capillary speed, and do not fall with the grid.

## Known limitations

- Only the static drop is checked. The oscillation frequency of a perturbed drop (Rayleigh and Lamb) and capillary waves were not
  run: no example sets them up.
- The spurious currents (16 to 23 % of the capillary speed) do not fall with the grid; the jump is averaged over 20 saves to
  reduce their pressure signature, which at 100 cells is still 0.9 % of the jump.
- 2D only; a 3D drop (\f$\Delta p = 2\sigma/R\f$) was not run.

## Reproduce

```
./mfc.sh run examples/2D_laplace_pressure_jump/case.py -t pre_process simulation post_process -n 1 -- -N 50
python3 examples/2D_laplace_pressure_jump/check_laplace.py examples/2D_laplace_pressure_jump --case-args "-N 50"
```

`case.py` takes `-N` (default 50), `--tend` (0.5 s) and `--saves` (100). The 50-cell run takes 27 minutes on one CPU core and 22
minutes on an A100 with `--case-optimization` (the run is small enough that the GPU does not help); 100 cells were run on the GPU only.
The example is on the test-suite skip list; no test was added.
