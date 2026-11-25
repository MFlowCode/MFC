@page validation_vortex_order Validation: Isentropic Vortex Order of Accuracy and Taylor-Green Decay

# Smooth-flow accuracy: the translating isentropic vortex and the viscous Taylor-Green vortex

**Status:** open discrepancy (isentropic vortex: second order for every WENO order,
[#2464](https://github.com/MFlowCode/MFC/issues/2464)); validated (2D Taylor-Green kinetic-energy decay). **Code:** `master` at
5838b641c. **Examples:** `examples/2D_isentropicvortex` (`--convect`), `examples/2D_TaylorGreenVortex`. **Last run:** 2026-11-25.

## Purpose

Shock tubes test the discontinuities; a smooth exact solution tests the design order of the schemes. The isentropic vortex is an
exact solution of the Euler equations that translates without change, so the error after a translation measures the accuracy of
the whole scheme on a nonlinear smooth flow. The 2D Taylor-Green vortex decays at an exact rate set by the viscosity, which tests
the viscous terms.

## Isentropic vortex

| Item | Value |
| :--- | :--- |
| Flow | the vortex of `examples/2D_isentropicvortex` \cite Coralic14 (after \cite Yee1999), \f$\varepsilon\f$ = 5, on a uniform state (1, 1, 1), \f$\gamma\f$ = 1.4: \f$T = 1 - \frac{\varepsilon^2}{16(\gamma+1)\pi^2}e^{2(1-r^2)}\f$, \f$\rho = T^{\gamma}\f$, \f$p = T^{\gamma+1}\f$, \f$u_\theta = \frac{\varepsilon}{2\pi} r e^{1-r^2}\f$ (min \f$\rho\f$ = 0.39), translating at (1, 0.5). Its pressure balances the rotation exactly, so it is an exact solution of the Euler equations (its entropy varies with r) |
| Domain | [-5, 5]², periodic, N × N uniform cells |
| Initial condition | exact conservative cell averages (`hcid = 283`, 3 × 3 Gauss) |
| Numerics | WENO (orders 3, 5, 7; JS, mapped, Z) or TENO5, HLLC, SSP-RK3, fixed \f$\Delta t\f$ at CFL 0.4 on a bound of \|u\| + c, t = 2 |
| Reference | the initial vortex moved by (2, 1), averaged over each cell with 5 × 5 Gauss points |
| Measured | density error, L1 and L2 normalized by the box area, and the max norm, by `check_vortex.py` |

L1 density error (order between successive grids in parentheses):

| N | WENO3 | WENO5-M | WENO7-M | WENO5-JS | WENO5-Z | TENO5 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 5.14e-3 | 8.13e-4 | needs N \f$\ge\f$ 36 | 1.36e-3 | 7.80e-4 | 7.69e-4 |
| 64 | 1.20e-3 (2.10) | 2.15e-4 (1.92) | 2.16e-4 | 2.32e-4 (2.55) | 2.15e-4 (1.86) | 2.14e-4 (1.84) |
| 128 | 2.57e-4 (2.22) | 5.47e-5 (1.98) | 5.47e-5 (1.98) | 5.48e-5 (2.08) | 5.47e-5 (1.97) | 5.47e-5 (1.97) |
| 256 | 4.37e-5 (2.55) | 1.371e-5 (2.00) | 1.371e-5 (2.00) | 1.371e-5 (2.00) | 1.371e-5 (2.00) | 1.371e-5 (2.00) |

Max-norm density error of WENO5-M: 3.25e-2, 8.17e-3, 2.36e-3, 6.07e-4 (order 2.0, 1.8, 2.0).

- From N = 64 on, all fifth- and seventh-order schemes give the same error to three digits, and it falls at second order.
- It does not depend on the time step: at N = 128 with WENO7-M, CFL 0.8, 0.4, 0.2 and 0.1 give 5.462e-5, 5.468e-5, 5.468e-5 and
  5.469e-5; CFL 0.05 gives the values of CFL 0.4 for WENO5-M and WENO7-M at N = 32 to 128.
- On an A100 (NVHPC 25.11) WENO5-M at N = 128 and 256 and WENO7-M at N = 256 agree with the CPU to all printed digits (final densities
  within 2e-13 of each other).

The common second-order error is consistent with the conversion of cell averages to primitive variables before the reconstruction:
the primitive of the averaged state differs from the average of the primitive by O(h²). A numpy replica of this vortex gives
\f$\max|\overline{\rho u}/\bar\rho - \bar u|\f$ = 9.5e-3, 3.5e-3, 9.5e-4, 2.4e-4 at N = 32 to 256 (order 2), the size of the
measured max-norm error. The 1D and 2D advection convergence tests of the suite, in which u and p are constant so the primitives are
linear in the conservative variables, reach the design orders. The cause is not verified by a run that reconstructs conservative or
characteristic variables (that needs a code change).

## 2D Taylor-Green vortex

| Item | Value |
| :--- | :--- |
| Flow | \f$u = U_0 \sin(x/L_0)\cos(y/L_0)\f$, \f$v = -U_0\cos(x/L_0)\sin(y/L_0)\f$ \cite TaylorGreen1937 (`patch_icpp%%geometry = 20`), \f$U_0\f$ = 0.1, \f$L_0\f$ = 1, \f$\rho\f$ = 1.22, p = 1e5, \f$\mu\f$ = 1e4 (\f$U_0L_0/\nu\f$ = 1.2e-5, Mach 3e-4) |
| Domain | [-\f$\pi\f$, \f$\pi\f$]², periodic, N × N cells |
| Numerics | WENO5-M, HLLC, SSP-RK3, `weno_Re_flux`, `weno_avg`, \f$\Delta t\f$ = 1e-8, 10000 steps (\f$4\nu t/L_0^2\f$ = 3.28), 100 saves |
| Reference | exact incompressible solution, velocity decaying as \f$e^{-2\nu t/L_0^2}\f$, so the kinetic energy decays at the rate \f$4\nu/L_0^2\f$ = 3.27869e4 s⁻¹ |
| Measured | rate fitted to ln KE over all saves, and KE at the last save against the exact decay, by `check_tgv.py` |

| N | KE drop \f$4\nu t\f$ | Fitted rate | Rate error | KE(t) error |
| ---: | ---: | ---: | ---: | ---: |
| 32 | 3.28 | 3.27869e4 | -7.3e-7 | +2.4e-6 |
| 64 | 3.28 | 3.27869e4 | -9.7e-9 | +3.2e-8 |
| 128 | 3.28 | 3.27869e4 | -1.1e-10 | +3.6e-10 |
| 200 | 3.28 | 3.27869e4 | -8.3e-12 | +3.3e-11 |

The rate error falls at sixth order or faster with N (6.2, 6.5 and 5.7 between successive grids; at 200 cells it is near round-off): the single Fourier mode is resolved by
the high-order viscous fluxes, and the time error is negligible at this \f$\Delta t\f$. On an A100 (NVHPC 25.11, `parallel_io = T`) the 64-cell
run agrees with the CPU to 1e-12 in density and to all printed digits in the table, and the 200-cell run gives -7.9e-12 and +3.2e-11.

## Assessment

- **Isentropic vortex: open discrepancy.** Every WENO and TENO variant converges at second order on this nonlinear smooth flow,
  with the same error constant for fifth and seventh order from N = 64. The high-order schemes still have a 3 to 6 times smaller
  error than WENO3 at moderate N, but the design order is not reached. Filed as
  [#2464](https://github.com/MFlowCode/MFC/issues/2464).
- **Taylor-Green: validated.** The viscous decay rate is exact to 1e-6 at 32 cells.

## Known limitations

- The Taylor-Green check is at a Reynolds number of 1.2e-5 and a Mach number of 3e-4, where the flow is a single decaying mode. The
  3D Re = 1600 transition (`examples/3D_TaylorGreenVortex`) and the IGR Taylor-Green examples are not checked here.
- The initial pressure of the 2D Taylor-Green patch is a quarter of the exact one ([#2465](https://github.com/MFlowCode/MFC/issues/2465)).
  At this Mach number the start-up acoustic transient does not affect the decay rate measurably.
- The time-step order is not measured here: the spatial error dominates the vortex, and the suite covers the RK1-3 orders on 1D
  advection (`Convergence -> Temporal`).
- WENO7 needs at least 36 cells per direction (`num_stcls_min * weno_order`).

## Reproduce

```
./mfc.sh run examples/2D_isentropicvortex/case.py -t pre_process simulation -n 1 -- --convect -N 128 --order 5
python3 examples/2D_isentropicvortex/check_vortex.py <run N=64> <run N=128> <run N=256>
./mfc.sh run examples/2D_TaylorGreenVortex/case.py -t pre_process simulation -n 1 -- -N 64
python3 examples/2D_TaylorGreenVortex/check_tgv.py <run N=32> <run N=64>
```

`--convect` also takes `--scheme wenojs|wenom|wenoz|teno`, `--cfl`, `--t-end` and `--time-stepper`; without it the example runs its
original 10-step stationary vortex on a stretched grid. The N = 256 vortex takes about 70 s on one CPU core, the 128-cell
Taylor-Green run about 12 minutes. Both examples are registered in the suite at reduced size in their default mode; no test is
added here.
