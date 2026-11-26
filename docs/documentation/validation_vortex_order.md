@page validation_vortex_order Validation: Isentropic Vortex Order of Accuracy and Taylor-Green Decay

# Smooth-flow accuracy: the translating isentropic vortex and the viscous Taylor-Green vortex

**Status:** isentropic vortex: validated against the expected second order, a known limit of the method
([#2464](https://github.com/MFlowCode/MFC/issues/2464)); 2D Taylor-Green kinetic-energy decay: validated. **Code:** `master` at
5838b641c (vortex, Taylor-Green), 3e4c53258 (1D smooth waves). **Examples:** `examples/2D_isentropicvortex` (`--convect`), `examples/2D_TaylorGreenVortex`, `examples/1D_smooth_wave_order`. **Last run:** 2026-11-25.

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

### Why the vortex is second order

MFC reconstructs **primitive** variables. Each cell's conserved averages are converted to primitives cell by cell, and the WENO
or TENO reconstruction treats those values as cell averages of the primitives. This is deliberate: in the diffuse-interface
multiphase models, reconstructing \f$\alpha\f$, \f$u\f$ and \f$p\f$ keeps a material interface in pressure and velocity
equilibrium, whereas conservative reconstruction creates pressure oscillations there (Abgrall's condition; \cite JohnsenColonius2006,
\cite Coralic14). The price is on nonlinear smooth flows. The primitive of an average differs from the average of the primitive
by \f$O(h^2)\f$, for example \f$\overline{\rho u}/\bar\rho - \bar u = \frac{h^2}{12}\,\rho_x u_x/\rho + O(h^4)\f$. This is a
smooth \f$O(h^2)\f$ perturbation of the data the reconstruction receives, so the face states, and with them the flux
divergence, carry an \f$O(h^2)\f$ error whatever the WENO order. Wherever the primitives are linear in the conserved
variables the conversion is exact and the design order is kept.

In more than one dimension a second, independent limit applies. MFC's finite-volume scheme works dimension by dimension and
evaluates each face flux at one point, the face midpoint, from face-averaged states. For a nonlinear flux this is second-order
accurate: the face average of the flux is not the flux of the face average. Zhang, Zhang and Shu \cite ZhangZhangShu2011 show that
such schemes are second order on nonlinear multidimensional problems, while for linear problems they keep the full order of the
1D reconstruction. Removing both limits would need a Gaussian quadrature of the flux on each face and a conservative or
characteristic reconstruction. MFC does neither, on purpose.

### 1D smooth waves: which flows reach the design order

The two mechanisms can be separated in 1D, where there is no face quadrature. Two smooth periodic problems were run with the
initial condition set to the exact conservative cell averages (10-point Gauss, read through `hcid = 170`; `examples/1D_smooth_wave_order`) and the density
measured against exact cell averages at the end. HLLC and SSP-RK3 were used, with CFL 0.02 (0.005 for seventh order); halving
the CFL changes no digit shown.

- **Linear:** a density wave \f$\rho = 1 + 0.2\sin 2\pi x\f$ carried at u = 1 through p = 1 (\f$\gamma\f$ = 1.4) for one period.
  u and p are uniform, so the primitives are linear in the conserved variables.
- **Nonlinear:** the isentropic \f$\gamma = 3\f$ wave \f$\rho_0 = 1 + 0.2\sin 2\pi x\f$, \f$u_0 = 0\f$, \f$p = \rho^3/3\f$. For
  \f$\gamma = 3\f$ the Riemann invariants \f$u \pm c\f$ each obey Burgers' equation, so the exact solution follows from the
  characteristics; it stays smooth until t = 0.80 and is compared at t = 0.2.

L1 density error, CPU (NVHPC 25.11), one rank (order between successive grids in parentheses):

| problem | scheme | initial condition | N = 32 | 64 | 128 | 256 |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| linear | WENO5-M | cell averages | 3.88e-6 | 1.21e-7 (5.00) | 3.80e-9 (5.00) | 1.20e-10 (4.99) |
| linear | TENO5 | cell averages | 3.85e-6 | 1.21e-7 (4.99) | 3.80e-9 (5.00) | 1.20e-10 (4.99) |
| linear | WENO5-M | point values | 2.01e-4 | 5.10e-5 (1.98) | 1.28e-5 (2.00) | 3.20e-6 (2.00) |
| linear | WENO3 | cell averages | 9.67e-3 | 2.34e-3 (2.05) | 5.33e-4 (2.13) | 1.16e-4 (2.20) |
| \f$\gamma\f$ = 3 | WENO5-M | cell averages | 1.96e-5 | 4.27e-6 (2.20) | 1.05e-6 (2.03) | 2.61e-7 (2.00) |
| \f$\gamma\f$ = 3 | TENO5 | cell averages | 1.96e-5 | 4.27e-6 (2.20) | 1.05e-6 (2.03) | 2.61e-7 (2.00) |
| \f$\gamma\f$ = 3 | WENO5-M | point values | 1.09e-4 | 2.77e-5 (1.98) | 6.94e-6 (2.00) | 1.73e-6 (2.00) |
| \f$\gamma\f$ = 3 | WENO3 | cell averages | 5.82e-4 | 1.07e-4 (2.44) | 1.59e-5 (2.76) | 1.84e-6 (3.11) |

Seventh order needs N \f$\ge\f$ 36:

| problem | scheme | N = 48 | 64 | 96 | 128 |
| :--- | :--- | ---: | ---: | ---: | ---: |
| linear | WENO7-M | 1.87e-9 | 2.51e-10 (6.98) | 1.58e-11 (6.83) | 3.69e-12 (5.05) |
| linear | TENO7 | 1.87e-9 | 2.51e-10 (6.98) | 1.57e-11 (6.83) | 3.69e-12 (5.05) |
| \f$\gamma\f$ = 3 | WENO7-M | 7.44e-6 | 4.18e-6 (2.01) | 1.86e-6 (2.00) | 1.04e-6 (2.00) |
| \f$\gamma\f$ = 3 | TENO7 | 7.44e-6 | 4.18e-6 (2.01) | 1.86e-6 (2.00) | 1.04e-6 (2.00) |

(The linear WENO7 error at 128 cells, 4e-12, is at the level of the time and round-off error.)

- The linear wave reaches the design order (5 and 7) when the initial condition is the exact cell averages, and only second
  order from point values: a point value differs from the cell average by \f$O(h^2)\f$, and the scheme carries that difference
  along.
- The nonlinear wave is second order for every order and both WENO and TENO, starting from exact averages, in 1D. The fifth- and
  seventh-order errors agree to three digits (2.61e-7 at N = 256; 1.04e-6 at N = 128 for both). A numpy replica of the same
  scheme (WENO5, LLF flux) gives 2.61e-7 at N = 256 with primitive reconstruction and 3.3e-10 at fifth order with conservative
  reconstruction (posted on [#2464](https://github.com/MFlowCode/MFC/issues/2464)). So the primitive conversion alone sets the
  second-order limit of the vortex, and the face quadrature adds a second one in 2D.
- The WENO order still matters. It sets the error at moderate resolution: on the \f$\gamma = 3\f$ wave WENO3 is 30 times less
  accurate than WENO5 at N = 32 and 7 times at 256, and on the vortex 5.6 to 3.2 times (N = 64 to 256). It also sets the
  shock-capturing quality (@ref validation_shock_tubes). Above about 64 cells per wavelength, fifth and seventh order give the
  same error on nonlinear smooth flows.

The convergence tests of the suite (`Convergence -> 1D/2D/3D -> WENO5, WENO7, TENO5, TENO7` in `toolchain/mfc/test/cases.py`,
run by `run_h_sweep` in `toolchain/mfc/test/convergence.py` on `examples/1D_euler_convergence`, `2D_advection_convergence` and
`3D_advection_convergence`) assert the design order and pass. They are the linear case above: a density sine carried at uniform
u and p, where the average-to-primitive conversion is exact and the flux is linear in the face state, so neither limit applies.
They also start from point values but compare the solution with the initial array shifted by a whole number of cells (or after
a whole period), not with exact averages, so the \f$O(h^2)\f$ initial-condition difference cancels. They test the
reconstruction stencils and weights, not the accuracy on a nonlinear smooth flow.

The 1D Alfven wave of @ref validation_mhd reaches fifth order once the time step is small: its density is uniform, so the
conversion is exact there too. Its order of 3.0 at CFL 0.4 is the SSP-RK3 time error.

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

- **Isentropic vortex: second order, as expected.** Every WENO and TENO variant converges at second order on this nonlinear
  smooth flow, with the same error constant for fifth and seventh order from N = 64. This is the expected order of
  primitive-variable reconstruction from cell averages (also in 1D, shown above) and of dimension-by-dimension finite volumes
  with midpoint face fluxes \cite ZhangZhangShu2011. The high-order schemes still have a 3 to 6 times smaller error than WENO3
  at moderate N.
- **Taylor-Green: validated.** The viscous decay rate is exact to 1e-6 at 32 cells.

## Known limitations

- On nonlinear smooth flows the spatial order is two for every `weno_order` (see above). The design orders hold where the
  primitives are linear in the conserved variables and the flux is linear in the face state (uniform-velocity, uniform-pressure
  advection), and the reconstruction order still sets the error at moderate resolution and the shock capturing.

- The Taylor-Green check is at a Reynolds number of 1.2e-5 and a Mach number of 3e-4, where the flow is a single decaying mode. The
  3D Re = 1600 transition (`examples/3D_TaylorGreenVortex`) and the IGR Taylor-Green examples are not checked here.
- The 2D Taylor-Green patch sets the exact initial pressure since [#2465](https://github.com/MFlowCode/MFC/issues/2465) (it was a
  quarter of it); the table above, measured before the fix, is unchanged at the printed digits (32 cells, rerun after the fix).
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

The 1D smooth waves are `examples/1D_smooth_wave_order` (`--problem advect|gamma3`, `--ic avg|point`, `--order`, `--teno`,
`--cfl`); its README gives the sweep commands. It is not in the test suite.

`--convect` also takes `--scheme wenojs|wenom|wenoz|teno`, `--cfl`, `--t-end` and `--time-stepper`; without it the example runs its
original 10-step stationary vortex on a stretched grid. The N = 256 vortex takes about 70 s on one CPU core, the 128-cell
Taylor-Green run about 12 minutes. Both examples are registered in the suite at reduced size in their default mode; no test is
added here.
