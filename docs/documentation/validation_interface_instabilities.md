@page validation_interface_instabilities Validation: Interface Instabilities

# Interface instabilities: viscous Rayleigh-Taylor growth rate

**Status:** validated (growth rate within 1 % of linear theory; a -0.7 % offset remains at the finest grid). **Code:** `master` at 5838b641c. **Example:** `examples/2D_rayleigh_taylor`. **Last run:** 2026-11-24. See also @ref validation.

## Purpose

A heavy fluid above a light one under gravity is unstable; while the interface amplitude is small each Fourier mode grows
exponentially at a rate fixed by the density and viscosity profiles \cite Chandrasekhar61. This checks the coupling of the
two-fluid model, gravity (`bf_y`) and the viscous stresses against that rate.

## Setup

| Item | Value |
| :--- | :--- |
| Geometry | 2D, one wavelength lambda = 0.2 m wide (periodic in x), 1.2 m tall, no-slip walls at y = 0 and 1.2 m (`bc_y = -16`); N x 6N cells |
| Fluids | ideal gases (gamma 1.4), rho = 3 kg/m^3 above, 1 below (Atwood number 0.5); mu = 0.0219 and 0.0073 Pa s (equal kinematic viscosity 0.0073 m^2/s) |
| Initial state | hardcoded IC 204: interface y = 0.7 - a0 cos(k x), k a0 = 0.05, tanh profile of the volume fraction with thickness delta = 2.5e-3 m; hydrostatic pressure, 1e5 Pa at the top; at rest |
| Drive | g = 9.81 m/s^2 (Mach number of the growth about 1e-4) |
| Numerics | 6-equation model (the table was run with `mpp_lim = T`; the example now sets `mpp_lim = F` (#2476), and the 32-cell CPU case rerun with it off over the same window gives n = 6.4335 1/s, -1.12 % (-1.16 % and -1.30 % for the other two caps), Phoenix GNU 12.3), WENO5-M with `mp_weno`, `weno_Re_flux`, HLLC, SSP-RK3, dt = 0.5 dx/c; t = 0.6 s, 60 saves |
| Reference | growth rate of the linearized incompressible Navier-Stokes equations about the same rho(y) and mu(y) = alpha mu_H + (1 - alpha) mu_L: the quadratic eigenproblem \f$n^2 M w + n C w + G w = 0\f$ for the vertical velocity, \f$M w = k^2\rho w - (\rho w')'\f$, \f$C w = (\mu(w'' + k^2 w))'' + k^2\mu(w'' + k^2 w) - 4k^2(\mu w')'\f$, \f$G w = -g k^2\rho' w\f$, \f$w = w' = 0\f$ at 0.3 m from the interface, second-order finite differences on 3200 points |
| Measured | the interface height (alpha = 1/2 crossing in each column) and its cos(k x) amplitude a(t); n from a least-squares fit of \f$a = A e^{nt} + B e^{-mt}\f$ over the saves with k a <= 0.25 (the second term absorbs the start from rest); `check_rt.py` |

The eigen solver gives n = 6.5067 1/s (3200 points; 6.5070 with 800). For a sharp interface it gives 6.534 against 6.535 from the
exact two-fluid dispersion relation (a 4 x 4 determinant of the semi-infinite solutions), and without viscosity 12.03 for the
diffuse and 12.37 for a ten times thinner profile, against \f$\sqrt{Agk}\f$ = 12.41 for a sharp one. The diffuse interface lowers
the viscous rate by 0.4 %; viscosity halves it.

## Results

Linear theory n = 6.5067 1/s. Fit over k a <= 0.25 (t <= 0.29 to 0.33 s, 30 to 34 saves); the two other columns move the cap:

| N (cells per wavelength) | Platform | n [1/s] | error | error, k a <= 0.15 | error, k a <= 0.35 |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 32 | GNU 13.3 CPU | 6.4337 | -1.12 % | -1.16 % | -1.30 % |
| 32 | A100 | 6.4337 | -1.12 % | | |
| 64 | A100 | 6.4560 | -0.78 % | | |
| 100 | A100 | 6.4611 | -0.70 % | -1.05 % | -0.75 % |

The A100 runs used NVHPC 25.11 with `--case-optimization`. The 32-cell CPU and GPU runs agree to the digits shown.

## Assessment

- The growth rate is within 1.2 % of the viscous linear rate on every grid and approaches it with refinement (-1.12, -0.78,
  -0.70 %). Viscosity halves the rate against \f$\sqrt{Agk}\f$, so the agreement checks the viscous stresses of both fluids,
  not only the buoyancy.
- The remaining -0.7 % at 100 cells is twice the spread between fit windows (0.35 %). Part of it may come from the fit, which
  assumes one growing and one decaying exponential for a start from rest; it was not separated from a solver error.
  Compressibility (\f$n/(kc)\approx 6\times10^{-4}\f$) is too small to explain it.
- At 32 cells the initial amplitude is 0.04 instead of 0.05 in k a: the tanh profile (2.5 mm) is thinner than a cell (6.25 mm)
  and is sampled at the cell centres. The growth rate does not depend on the amplitude in the linear regime.

## Known limitations

- The Kelvin-Helmholtz example (`examples/2D_kelvin_helmholtz`) is not checked: its shear layers are sharp, so the linear growth
  rate of the perturbation depends on the numerical thickness of the layer, and its perturbation (0.1 of the velocity jump) is not
  small.
- Only one wavelength, Atwood number and viscosity; the late-time (bubble and spike) regime is not compared with anything.
- The 3D cases (`examples/3D_rayleigh_taylor`) are not run.

## Reproduce

```
./mfc.sh run examples/2D_rayleigh_taylor/case.py -t pre_process simulation post_process -n 1 --gpu --case-optimization -- -N 100 --tend 0.6 --saves 60
python3 examples/2D_rayleigh_taylor/check_rt.py examples/2D_rayleigh_taylor --case-args "-N 100 --tend 0.6 --saves 60"
```

`case.py` takes `-N` (cells per wavelength, default 200), `--tend` (default 2 s) and `--saves` (default 100); the defaults are the
original deck. On NVIDIA GPUs build with `--case-optimization`: without it the 6-equation pressure relaxation allocates on the device
heap, and the default 240k-cell deck aborts ([#2466](https://github.com/MFlowCode/MFC/issues/2466)). The example runs in the test
suite at a reduced size, as every registered example does; no test was added.
