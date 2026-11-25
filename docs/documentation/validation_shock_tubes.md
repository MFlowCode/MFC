@page validation_shock_tubes Validation: Shock Tubes and the Shu-Osher Problem

# 1D Riemann problems against the exact solution, and the Shu-Osher problem

**Status:** validated (Sod, Lax, 123 with WENO, TENO and MUSCL; Shu-Osher against a converged reference); partially validated for
IGR (Sod and Lax only). **Code:** `master` at 5838b641c. **Examples:** `examples/1D_sod_convergence`, `examples/1D_sodshocktube`,
`examples/1D_laxshocktube`, `examples/1D_shuosher_*`. **Last run:** 2026-11-25.

## Purpose

Shock tubes are the basic check of a shock-capturing scheme: the exact solution of the Riemann problem is known, so the error of the
shock, contact and rarefaction can be measured, not only seen. This checks the single-fluid ideal-gas path (`model_eqns = 2`, one
fluid) with each reconstruction and Riemann solver against the exact solution, and the shock-entropy-wave interaction of Shu and
Osher against a converged reference.

## Setup

| Item | Value |
| :--- | :--- |
| Problems | Sod \cite Sod1978: (1, 0, 1) / (0.125, 0, 0.1), t = 0.2. Lax \cite Lax1954: (0.445, 0.698, 3.528) / (0.5, 0, 0.571), t = 0.14. 123 (Toro test 2 \cite Toro09): (1, -2, 0.4) / (1, 2, 0.4), t = 0.15. States are (\f$\rho\f$, u, p), \f$\gamma\f$ = 1.4, jump at x = 0.5 on [0, 1] |
| Boundaries | extrapolation (`bc_x = -3`) at both ends |
| Numerics | default WENO5-M (`mapped_weno`, `weno_eps` = 1e-40), HLLC, SSP-RK3, fixed \f$\Delta t\f$ at CFL 0.3 on a bound of the signal speed (`examples/1D_sod_convergence/case.py`) |
| Reference | exact Riemann solution (Toro chapter 4, `riemann_exact.py`), averaged over each cell with 64 sub-points |
| Measured | L1 norm \f$\sum_i |q_i - q_i^\mathrm{exact}|\,\Delta x\f$ of density, velocity and pressure at the final time, by `check_riemann.py` |
| Hardware | GNU 13.3, CPU, one rank; the A100 rows below are NVHPC 25.11, OpenACC |

`riemann_exact.py` reproduces the star states of Toro's table 4.2 (Sod: p* = 0.30313, u* = 0.92745; 123: p* = 0.00189, u* = 0).

## Results

### Grid convergence (WENO5-M, HLLC)

L1 density error, observed order between successive grids in parentheses:

| Cells | Sod | Lax | 123 |
| ---: | ---: | ---: | ---: |
| 128 | 2.70e-3 | 9.28e-3 | 4.21e-3 |
| 256 | 1.22e-3 (1.14) | 5.80e-3 (0.68) | 2.25e-3 (0.91) |
| 512 | 7.80e-4 (0.65) | 2.80e-3 (1.05) | 1.54e-3 (0.55) |
| 1024 | 4.19e-4 (0.90) | 1.82e-3 (0.63) | 1.06e-3 (0.53) |
| 2048 | 2.18e-4 (0.94) | 9.23e-4 (0.98) | 6.82e-4 (0.64) |

The order of the L1 error of a solution with discontinuities is at most one: the shock and the contact are each smeared over a few
cells. The 123 problem converges more slowly (0.5 to 0.65): its error sits at the centre, where the two rarefactions leave a
near-vacuum (\f$\rho\f$ = 0.022, p = 0.0019) and the density shows the usual dip of shock-capturing schemes there; its velocity error is
the largest of the three (L1 4.0e-3 at 2048 cells against 3.0e-4 for Sod).

### Reconstruction and Riemann solver (Sod)

L1 density error at 512 and 2048 cells, and the order from 1024 to 2048:

| Scheme | 512 | 2048 | order |
| :--- | ---: | ---: | ---: |
| WENO1 (first order) | 7.06e-3 | 2.80e-3 | 0.67 |
| WENO3 | 1.39e-3 | 4.03e-4 | 0.92 |
| WENO5-M | 7.80e-4 | 2.18e-4 | 0.94 |
| WENO7-M | 7.71e-4 | 2.07e-4 | 0.99 |
| TENO5 (`teno_CT` = 1e-6) | 7.15e-4 | 2.41e-4 | 0.78 |
| MUSCL, minmod | 1.85e-3 | 5.99e-4 | 0.82 |
| MUSCL, MC | 1.02e-3 | 3.18e-4 | 0.86 |
| WENO5-M, HLL | 8.19e-4 | 2.28e-4 | 0.93 |
| WENO5-M, Lax-Friedrichs | 9.93e-4 | 2.60e-4 | 1.00 |
| IGR, order 5 (2D strip, see below) | 5.21e-3 | 1.39e-3 | 0.97 |
| IGR, order 3 | 5.02e-3 | 1.37e-3 | 0.95 |

At a fixed grid the error falls with the order of the reconstruction up to fifth order; seventh order gains little because the
error is set by the discontinuities. HLL and Lax-Friedrichs, which do not resolve the contact, are 4 % and 19 % above HLLC at 2048
cells. Lax with IGR (order 5): 6.45e-3 at 512 cells and 1.95e-3 at 2048 (2.1 times the WENO5-M error).

### The examples as shipped

| Example | Cells | t | L1 \f$\rho\f$ | L1 u | L1 p |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `1D_sodshocktube` (WENO5-M, HLLC) | 400 | 0.1 | 9.53e-4 | 1.94e-3 | 7.31e-4 |
| `1D_laxshocktube` (WENO5-M, HLLC) | 201 | 0.14 | 9.08e-3 | 1.26e-2 | 1.34e-2 |

### CPU and GPU

On an A100 (NVHPC 25.11, OpenACC, `parallel_io = T`) the L1 errors of Sod, Lax and 123 at 1024 cells, of IGR Sod at 1024, and of both
examples above agree with the GNU CPU values to all printed digits; the final densities differ by at most 1.4e-9 of the maximum
(4e-10 for Sod and Lax, 1e-12 and below for the examples and IGR).

## Shu-Osher problem

| Item | Value |
| :--- | :--- |
| Problem | a Mach 3 shock running into the density wave \f$\rho = 1 + 0.2\sin 5x\f$ \cite ShuOsher1989: (3.857, 2.629, 10.333) for x < 1, (\f$1 + 0.2\sin 5x\f$, 0, 1) for x > 1 on [0, 10], t = 1.8, extrapolation boundaries |
| Numerics | as in `examples/1D_shuosher_*` (WENO5-JS, -M, -Z, TENO5, TENO7; HLLC, SSP-RK3, 2000 steps at 1001 cells, CFL about 0.4) |
| Reference | `examples/1D_shuosher_wenoz5` run with `-N 16016` (WENO5-Z, 32032 steps), averaged onto each run's cells |
| Measured | L1 density error over the domain, and over 5.5 \f$\le\f$ x \f$\le\f$ 7.5 (the steepened entropy waves behind the shock and the shock), by `check_shuosher.py` |

| Example (1001 cells) | L1 domain | L1 waves |
| :--- | ---: | ---: |
| `1D_shuosher_wenojs5` | 8.02e-2 | 5.99e-2 |
| `1D_shuosher_wenom5` | 8.41e-2 | 5.85e-2 |
| `1D_shuosher_wenoz5` | 6.53e-2 | 4.57e-2 |
| `1D_shuosher_teno5` | 8.87e-2 | 6.41e-2 |
| `1D_shuosher_teno7` | 1.23e-1 | 9.81e-2 |

WENO5-Z, grid sweep: 7.19e-1 (143 cells), 3.47e-1 (286), 9.69e-2 (572), 5.92e-2 (1144); the error drops quickly once the
post-shock waves are resolved (about 500 cells) and then at order 0.7, where the shock and the remaining wave error compete.
The reference itself: WENO5-Z at 2288, 4004 and 8008 cells differs from it by 3.83e-2, 1.70e-2 and 6.23e-3 (order 1.45), so its
own error is about 2.5e-3, under 4 % of the errors at 1001 cells.

At 572 and 1144 cells WENO5-Z puts the steepest density jump (the shock, x = 7.39) in the same cell as the reference. WENO5-Z is
the most accurate at 1001 cells; TENO7 with `teno_CT = 1e-9` is the least accurate and the least robust to round-off: the CPU and
the A100 final densities differ by up to 0.10 (in the cells next to the shock and in the waves; 5e-4 after 200 steps), and the L1 error is 1.23e-1 on the CPU and 1.29e-1 on the A100. TENO5 behaves the same way at a
smaller scale (densities differ by 1.1e-2; L1 8.871e-2 and 8.864e-2). WENO5-JS and WENO5-Z agree between the two to 5e-12, and the
16016-cell reference to 3e-9.

## Assessment

- The exact-solution errors of all three Riemann problems decrease at the order expected of discontinuous solutions, for every
  reconstruction and Riemann solver, and are the same on the CPU and the A100.
- The ordering of the schemes is as expected: first-order WENO1 is 13 times less accurate than WENO5 at 2048 cells, MUSCL-minmod
  2.7 times, and HLLC is the most accurate of the three Riemann solvers.
- The Shu-Osher examples agree with a converged reference to 6.5e-2 to 1.2e-1 (L1 over a domain of length 10) at their default
  1001 cells.

## Known limitations

- **IGR has no 1D path** ([#2463](https://github.com/MFlowCode/MFC/issues/2463)): a 1D case with `igr = T` passes validation and
  segfaults. The IGR rows use a 2D strip of 8 cells in y (periodic) with `elliptic_smoothing` (10 iterations); without the smoothing
  the Sod strip gives NaN within 10 steps. With `alf_factor = 10` its error is 6 to 7 times that of WENO5-M (the regularization
  widens the shock over several cells), and it converges at the same order. The 123 problem fails with IGR (NaN, also with
  `igr_pres_lim`).
- **Shu-Osher examples differ from the textbook problem.** `1D_shuosher_wenojs5`, `-wenom5`, `-wenoz5`, `-teno5` and `-teno7` put the
  shock at x = 1 on [0, 10] and use \f$\sin 5x\f$ (the textbook problem is on [-5, 5] with the shock at x = -4, which shifts the
  phase of the wave by 0.13 rad), and round the left state to three decimals. `1D_shuosher_old` and `1D_shuosher_analytical` use the
  textbook domain but a left density of 3.957143 instead of 3.857143. The reference above is for the shipped problem, so the
  comparison is consistent; it is not a comparison with published figures.
- TENO5 and TENO7 on the Shu-Osher problem are sensitive to round-off at the shock (the CPU and A100 runs differ by 1e-2 to 1e-1 in
  density near the shock, see above), presumably through the sharp stencil cut-off (`teno_CT`); not investigated further.
- Only ideal gas and one fluid are checked here; multi-fluid shock tubes and stiffened gases are not.

## Reproduce

```
./mfc.sh run examples/1D_sod_convergence/case.py -t pre_process simulation -n 1 -- -N 512 --problem lax
python3 examples/1D_sod_convergence/check_riemann.py <run 128> <run 256> ... <run 2048>
./mfc.sh run examples/1D_shuosher_wenoz5/case.py -t pre_process simulation -n 1 -- -N 16016
python3 examples/1D_shuosher_wenoz5/check_shuosher.py <reference run> <run> [<run> ...]
```

Each run writes next to its `case.py`, so copy the example directory per run. Other options of `1D_sod_convergence/case.py`:
`--order`, `--teno`, `--muscl --muscl-lim`, `--riemann hll|hllc|lax_friedrichs`, `--igr 3|5`. The 2048-cell Sod run takes about
10 s on one CPU core; the 16016-cell Shu-Osher reference about 9 minutes.
The suite already runs the L1 self-convergence of `1D_sod_convergence` for ten reconstructions (`Convergence -> Sod`) and the
shipped examples at reduced size; no test is added here.
