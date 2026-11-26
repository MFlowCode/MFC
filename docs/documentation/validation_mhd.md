@page validation_mhd Validation: Ideal MHD

# Ideal MHD: Alfven wave, Brio-Wu shock tube, Orszag-Tang divergence error

**Status:** Alfven wave: validated (exact solution, design order). Brio-Wu: partially validated (grid convergence and the compound-wave states against published values; there is no exact solution to compare with). Orszag-Tang: partially validated (divergence error and mass conservation measured; no reference flow data). **Code:** `master` at 5838b641c. **Examples:** `examples/1D_mhd_smooth_alfven_wave`, `examples/1D_brio_wu_hlld`, `examples/1D_brio_wu`, `examples/2D_orszag_tang_hyper_cleaning`. See also @ref validation.

All tables were measured with GNU 13 (CPU, one rank, double precision) on `master` at 5838b641c unless a row says otherwise.

## Circularly polarized Alfven wave {#sec-validation-mhd-alfven}

| Item | Value |
| :--- | :--- |
| Case | Tóth \cite Toth2000 section 6.3.1: \f$\rho = 1\f$, \f$p = 0.1\f$, \f$\gamma = 5/3\f$, \f$B_x = 1\f$, \f$v_\perp = B_\perp = 0.1(\sin 2\pi x, \cos 2\pi x)\f$ (`hcid = 150`) |
| Domain | \f$x \in [0, 1]\f$, periodic, one period (t = 1) |
| Numerics | HLLD, SSP-RK3; WENO3 (the example) or WENO5, `weno_eps` = 1e-7; sweep at fast-wave CFL 0.4 |
| Reference | exact: the wave moves at \f$B_x/\sqrt{\rho} = 1\f$ in \f$-x\f$ with \f$\rho\f$, \f$p\f$, \f$|B_\perp|\f$ constant, so the state at t = 1 is the initial state |
| Measured | `check_alfven.py`: L1 error (mean \f$|\cdot|\f$) of \f$B_y\f$ at t = 1; \f$v_y\f$, \f$v_z\f$, \f$B_z\f$ agree with it to 1 %; max \f$|\rho - 1|\f$ |

| cells | WENO3 L1(\f$B_y\f$) | order | WENO5 L1(\f$B_y\f$) | order | WENO5 max \f$|\rho - 1|\f$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 5.69e-3 | | 1.99e-5 | | 2.1e-6 |
| 64 | 1.39e-3 | 2.03 | 1.21e-6 | 4.04 | 6.5e-8 |
| 128 | 2.34e-4 | 2.57 | 1.11e-7 | 3.44 | 1.8e-9 |
| 256 | 2.05e-5 | 3.51 | 1.27e-8 | 3.13 | 6.2e-11 |
| 512 | 1.16e-6 | 4.14 | 1.56e-9 | 3.03 | 1.9e-12 |

The example deck itself (1024 cells, WENO3, \f$\Delta t\f$ = 1e-5, CFL 0.011) gives L1(\f$B_y\f$) = 6.3e-8 and \f$|\rho - 1| \le 2.5\f$e-8.
WENO5 at fixed CFL approaches third order, the order of SSP-RK3: at CFL 0.4 the time error dominates from 128 cells on. With the
CFL number reduced, WENO5 reaches its design order (CPU, NVHPC 25.11; L1(\f$B_y\f$) at 32, 64, 128 and 256 cells):

| fast-wave CFL | 32 | 64 | 128 | 256 | orders |
| ---: | ---: | ---: | ---: | ---: | :--- |
| 0.4 | 1.99e-5 | 1.21e-6 | 1.11e-7 | 1.27e-8 | 4.04, 3.44, 3.13 |
| 0.1 | 1.37e-5 | 4.38e-7 | 1.48e-8 | 6.03e-10 | 4.97, 4.89, 4.62 |
| 0.025 | 1.36e-5 | 4.26e-7 | 1.33e-8 | 4.13e-10 | 5.00, 5.00, 5.01 |

The primitive-variable reconstruction that limits nonlinear smooth flows to second order (@ref validation_vortex_order) does not
act here: \f$\rho\f$ is uniform, so the velocity is linear in the momentum, and the pressure recovered from the cell averages
differs from p by a constant (\f$|v_\perp|\f$ and \f$|B_\perp|\f$ are uniform), which has no gradient.
WENO3 is second order on coarse grids, where its weights are
far from optimal at the extrema of the sine, and steepens toward its optimal-weight order as the grid is refined (the 4.1 between
256 and 512 cells is pre-asymptotic). The wave amplitude 0.1 is not small, so this checks the nonlinear HLLD fluxes as well as the
linear Alfven speed.

Limits: one dimension, so \f$\nabla\cdot B = 0\f$ holds trivially (the divergence error is measured on the Orszag-Tang vortex below); the 2D and 3D versions of the wave are not run.

## Brio-Wu shock tube {#sec-validation-mhd-briowu}

| Item | Value |
| :--- | :--- |
| Case | Brio and Wu \cite BrioWu1988: left \f$(\rho, p, B_y) = (1, 1, 1)\f$, right \f$(0.125, 0.1, -1)\f$, \f$B_x = 0.75\f$, at rest; \f$\gamma = 2\f$, and \f$\gamma = 5/3\f$ as in Ryu and Jones \cite RyuJones1995 |
| Domain | \f$x \in [0, 1]\f$, jump at 0.5, extrapolation ends, t = 0.2 (the right fast rarefaction reaches the end at about t = 0.13; all runs share this) |
| Numerics | WENO5-M, HLLD (`1D_brio_wu_hlld`) or HLL (`1D_brio_wu`), SSP-RK3, \f$\Delta t = 2\times10^{-4} \cdot 801/N\f$ |
| Reference | no exact solution: with the compound wave the Riemann problem is not unique, and dissipative schemes converge to the solution with the compound wave. (a) A 25632-cell HLLD run of the same deck, block-averaged onto each grid. (b) For \f$\gamma = 5/3\f$, the states on both sides of the intermediate shock of the compound wave that Ryu and Jones read from their 8192-cell calculation (section 5 of \cite RyuJones1995): preshock \f$\rho\f$ = 0.647, p = 0.484, \f$B_y\f$ = 0.536; postshock \f$\rho\f$ = 0.827, p = 0.758 |
| Measured | `check_brio_wu.py`: L1 errors against (a); preshock state (median over cells 12 to 4 before the sign change of \f$B_y\f$) and postshock state (the density peak within 10 cells after it) for (b) |

Grid convergence, HLLD, L1 errors against the 25632-cell run:

| cells | \f$\gamma = 2\f$: L1(\f$\rho\f$) | order | L1(\f$B_y\f$) | \f$\gamma = 5/3\f$: L1(\f$\rho\f$) | order | L1(\f$B_y\f$) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 801 | 2.20e-3 | | 3.47e-3 | 2.19e-3 | | 3.78e-3 |
| 1602 | 1.25e-3 | 0.81 | 2.04e-3 | 1.15e-3 | 0.93 | 1.87e-3 |
| 3204 | 7.64e-4 | 0.71 | 1.19e-3 | 6.87e-4 | 0.74 | 1.13e-3 |
| 6408 | 4.90e-4 | 0.64 | 7.07e-4 | 4.22e-4 | 0.70 | 6.76e-4 |

HLL (\f$\gamma = 2\f$) against the same HLLD reference: L1(\f$\rho\f$) = 2.87e-3 at 801 cells and 1.23e-3 at 3204 (order 0.61), 1.3 and 1.6 times the HLLD error. An L1 order below one is usual for a solution with contact discontinuities and a compound wave; the last ratio is also lowered by the reference being only four times finer. On the A100, the 1602-cell \f$\gamma = 2\f$ run against its own 25632-cell reference gives the same L1 errors to three digits.

Compound-wave states, \f$\gamma = 5/3\f$ (intermediate shock at x = 0.4435 at t = 0.2):

| cells | preshock \f$\rho\f$ | p | \f$B_y\f$ | postshock \f$\rho\f$ | p |
| ---: | ---: | ---: | ---: | ---: | ---: |
| Ryu and Jones | 0.647 | 0.484 | 0.536 | 0.827 | 0.758 |
| 801 | 0.6457 (−0.2 %) | 0.4823 (−0.3 %) | 0.5336 (−0.5 %) | 0.8035 (−2.8 %) | 0.7190 (−5.2 %) |
| 6408 | 0.6461 (−0.1 %) | 0.4829 (−0.2 %) | 0.5344 (−0.3 %) | 0.8019 (−3.0 %) | 0.7151 (−5.7 %) |
| 25632 | 0.6460 (−0.2 %) | 0.4827 (−0.3 %) | 0.5339 (−0.4 %) | 0.8298 (+0.3 %) | 0.7489 (−1.2 %) |

The preshock state agrees with Ryu and Jones to 0.5 % at every resolution, which is within the precision of values read from a plot. The postshock state is the tip of a narrow peak (the intermediate shock is followed directly by the slow rarefaction), so a cell-averaged maximum underestimates it until the peak is resolved: it is within 0.3 % (\f$\rho\f$) and 1.2 % (p) only at 25632 cells, and does not converge monotonically below that. MFC, like the TVD code of Ryu and Jones, produces the compound wave rather than the rotational discontinuity plus slow shock of the other exact solution.

Limits: no exact solution is compared; the reference values for (b) were read from a figure; only the coplanar problem (\f$B_z = 0\f$) is run.

## Orszag-Tang vortex: divergence error {#sec-validation-mhd-orszag-tang}

| Item | Value |
| :--- | :--- |
| Case | \f$\rho = 25/(36\pi)\f$, \f$p = 5/(12\pi)\f$, \f$\gamma = 5/3\f$, \f$v = (-\sin 2\pi y, \sin 2\pi x)\f$, \f$B = (-\sin 2\pi y, \sin 4\pi x)/\sqrt{4\pi}\f$ (`hcid = 250`), periodic unit square |
| Numerics | 513 x 513 cells, WENO5-M, HLL, SSP-RK3, \f$\Delta t\f$ = 4e-4, 2500 steps to t = 1; `hyper_cleaning` with speed 2.5 and damping time 0.004 (the example), and the control without it (`--no-clean`) |
| Reference | \f$\nabla\cdot B = 0\f$ exactly (the initial field is divergence free); total mass conserved |
| Measured | `check_divb.py`: \f$\epsilon = \overline{|\nabla\cdot B|}\,\Delta x/\overline{|B|}\f$ with \f$\nabla\cdot B\f$ by second-order central differences, its maximum \f$\max|\nabla\cdot B|\,\Delta x/\max|B|\f$, the relative change of the total mass |
| Hardware | NVIDIA A100 (NVHPC 25.11, OpenACC); the first 200 steps also with GNU on the CPU |

| t | \f$\epsilon\f$, cleaning | max, cleaning | \f$\epsilon\f$, no cleaning | ratio |
| ---: | ---: | ---: | ---: | ---: |
| 0.04 | 7.0e-7 | 3.8e-5 | 2.5e-7 | 0.35 |
| 0.12 | 3.5e-6 | 1.9e-4 | 3.2e-6 | 0.93 |
| 0.16 | 2.0e-5 | 3.1e-3 | 1.0e-4 | 5.1 |
| 0.24 | 1.7e-4 | 1.4e-2 | 7.0e-4 | 4.2 |
| 0.32 | 2.0e-4 | 1.9e-2 | 2.7e-3 | 13.3 |
| 0.36 | 2.6e-4 | 1.8e-2 | NaN at step 900 | |
| 0.6 | 4.2e-4 | 1.2e-2 | | |
| 1.0 | 9.1e-4 | 3.0e-2 | | |

The total mass changes by −1.3e-13 over the cleaned run (round-off). The first 200 steps on the CPU (GNU 13.3) give the same \f$\epsilon\f$ and maximum to four digits. From t = 0.16 on the cleaning keeps \f$\epsilon\f$ 4 to 13 times below the uncleaned run, and the uncleaned run stops with NaN at step 900 (t = 0.36), so this deck needs the cleaning. Before that the two runs differ little: the L1 density difference at t = 0.32 is 4.7e-3 of the mean density.
The divergence error with cleaning grows to 9e-4 at t = 1 and does not saturate over this time. The maximum, 1 to 3 % from t = 0.24 on, is a cell-centred central difference and includes the jump of B across captured discontinuities, so it bounds the divergence the scheme carries rather than measuring it; where it occurs was not mapped.

Limits: one resolution and one cleaning speed and damping time; the flow field is not compared with published Orszag-Tang profiles; the WENO1 deck `2D_orszag_tang` is not measured.


## Reproduce {#sec-validation-mhd-run}

```
./mfc.sh run examples/1D_mhd_smooth_alfven_wave/case.py -t pre_process simulation post_process -n 1 -- -N 128 --cfl 0.4 --weno 5
python3 examples/1D_mhd_smooth_alfven_wave/check_alfven.py examples/1D_mhd_smooth_alfven_wave
```
