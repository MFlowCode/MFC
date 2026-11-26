@page validation_acoustics Validation: Acoustic Sources

# Acoustic sources: planar pulses and a focused transducer

**Status:** validated for 1D planar sources. The emitted pulses arrive on time and converge to linear acoustics at third order, down to a floor of about 2e-4. **Code:** `master` at 3e4c53258 with the stage-time fix of [#2462](https://github.com/MFlowCode/MFC/issues/2462). **Examples:** `examples/1D_acoustic_gauss_sigmatime`, `examples/1D_acoustic_gauss_sigmadist`, `examples/1D_acoustic_sine_frequency`, `examples/1D_acoustic_dipole`. See also @ref validation and @ref validation_coordinates (the same pulse on a stretched grid).

The tables were measured on the CPU with NVHPC 25.11 (one rank, double precision), except the stage-time rows of the time-step refinement, which ran on the A100. On an NVIDIA A100 (NVHPC 25.5, OpenACC) the 400-cell Gaussian run and the dipole deck give the CPU values to the digits shown.

## Source model

The source adds \f$S = \mathrm{mag}\,g(t)\,f(\mathbf{x})\f$ to the momentum along the emission direction, \f$S/c\f$ to the mass and \f$S c/(\gamma - 1)\f$ to the energy, so \f$\partial_t p = c S\f$ and only the invariant travelling along the emission direction, \f$p + \rho c u\f$, is driven.
With `dipole = T` the mass and energy terms are doubled and the momentum term is dropped, so both invariants are driven and two pulses of the same sign leave the source, one each way.
In 1D (`support = 1`) the support \f$f\f$ is a Gaussian of standard deviation \f$1.25\Delta x\f$ normalized to unit integral, and linear acoustics gives, away from the support,
\f[ p'(x, t) = \mathrm{mag} \sum_j f_j \Delta x_j\, g\big(t - t_d - |x - x_j|/c\big), \f]
which `check_acoustic.py` evaluates with the discrete weights \f$f_j\f$ of the run. The support width therefore enters the reference, and the comparison isolates the propagation and the time integration of the source.

## Planar Gaussian and sine pulses {#sec-validation-acoustics-planar}

| Item | Value |
| :--- | :--- |
| Medium | water, stiffened gas \f$\gamma\f$ = 4.4, \f$\pi_\infty\f$ = 5.57e8 Pa, 1100 kg/m³, 1e5 Pa at rest; c = 1492.8 m/s |
| Domain | \f$x \in [0, 1]\f$ mm, non-reflecting ends (`bc_x = -6`), source at 0.1 mm (0.5 mm for the dipole) |
| Pulses | Gaussian: \f$\sigma_t\f$ = 20 ns (`gauss_sigma_time`, or `gauss_sigma_dist` = 1493 \f$\sigma_t\f$), delay 100 ns; sine: 20 MHz, 3 cycles (15 cells per wavelength); `mag` = 1 Pa |
| Numerics | WENO5 with TENO (`teno_CT` = 1e-8), HLLC, SSP-RK3, \f$c\Delta t/\Delta x\f$ = 0.6 |
| Time stepping | every RHS term that depends explicitly on time, the source among them, is evaluated at its own stage time: t, t + \f$\Delta t\f$ and t + \f$\Delta t/2\f$ for SSP-RK3 (t and t + \f$\Delta t\f$ for RK2) |
| Reference | the expression above at the last save with the whole signal inside the domain (t = 500 ns; 324 ns for the dipole) |
| Measured | `check_acoustic.py`: relative L2 error of \f$p'\f$ outside 5 support widths; amplitude A and shift \f$\tau\f$ of the least-squares fit \f$p'_\text{MFC} = A\,p'_\text{lin}(t + \tau)\f$; the L2 error left after the fit |

Gaussian pulse (`1D_acoustic_gauss_sigmatime`), grid sweep at fixed \f$c\Delta t/\Delta x\f$:

| cells | \f$\Delta t\f$ [ns] | A − 1 | \f$\tau/\Delta t\f$ | rel. L2 error | order | rel. L2 after the fit | `mag` = 1e3 Pa: rel. L2 error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 2.0 | −6.5e-4 | +0.002 | 2.33e-3 | | 2.24e-3 | 2.38e-3 |
| 400 | 1.0 | −1.2e-4 | +0.001 | 3.08e-4 | 2.92 | 2.80e-4 | 2.99e-4 |
| 800 | 0.5 | −1.1e-4 | +0.002 | 2.53e-4 | | 2.27e-4 | 2.28e-4 |
| 1600 | 0.25 | −1.8e-4 | +0.002 | 4.67e-4 | | 4.30e-4 | 9.5e-5 |

The other decks, at their own resolution (200 cells):

| example | A − 1 | \f$\tau/\Delta t\f$ | rel. L2 error | rel. L2 after the fit |
| :--- | ---: | ---: | ---: | ---: |
| `1D_acoustic_gauss_sigmadist` | −6.5e-4 | +0.002 | 2.33e-3 | 2.23e-3 |
| `1D_acoustic_dipole` (sine, both directions) | −3.8e-4 | 0.000 | 1.35e-3 | 1.30e-3 |
| `1D_acoustic_sine_frequency` (15 cells per wavelength) | −4.6e-2 | +0.021 | 6.25e-2 | 4.21e-2 |

The pulses arrive on time (\f$|\tau| \le 0.002\,\Delta t\f$, against \f$-0.50\,\Delta t\f$ before the fix), so the plain L2 error now equals the error left after the fit and falls at third order (2.92 from 200 to 400 cells; it was first order, 3.45e-2, 1.75e-2, 8.8e-3 and 4.4e-3, while the pulse lagged).
From 800 cells on the error stops falling. At `mag` = 1 Pa the floor rises with the number of steps (2.3e-4 at 1000 steps, 4.3e-4 at 2000). This is round-off: a 1 Pa signal is 1.4e-9 of the energy \f$(p + \gamma p_\infty)/(\gamma - 1)\f$ = 7.2e8 Pa it is computed from, and the `master` code shows the same floor after its shift is removed (2.2e-4 and 4.1e-4). With `mag` = 1e3 Pa, still linear (\f$p'/\rho c^2\f$ = 4e-7), the error is 2.3e-4 at 800 cells and 9.5e-5 at 1600 (last column), the level the unfixed code reaches after its shift is removed (1.9e-4 and 1.7e-4).
The \f$-4.6\f$ % amplitude of the sine deck is the dissipation of WENO5 at 15 cells per wavelength over 9 wavelengths of travel; it was not refined (the deck has no resolution argument).
The dipole sends one pulse each way with the amplitude `mag`, as the doubled mass and energy terms predict.

Time-step refinement, 800 cells, `mag` = 1e3 Pa, \f$\Delta t\f$ = 0.5 ns / q (`master`: CPU; stage times: A100): the relative L2 difference of p' between successive \f$\Delta t\f$.

| scheme | q = 1 to 2 | 2 to 4 | 4 to 8 |
| :--- | ---: | ---: | ---: |
| `master` (step-start source time), WENO5 with TENO | 4.41e-3 | 2.21e-3 (order 1.00) | 1.10e-3 (1.00) |
| stage times, WENO5 with TENO (the deck) | 9.1e-5 | 3.1e-5 | 1.3e-4 |
| stage times, WENO5-M | 2.7e-4 | 6.4e-5 (2.06) | 5.2e-6 (3.63) |

The step-start source time made the whole solution first order in \f$\Delta t\f$. With the stage times the time error of WENO5-M falls at second to third order. With TENO (`teno_CT` = 1e-8) it stays at the 1e-4 level and does not fall regularly, consistent with the stencil selection switching between runs of different \f$\Delta t\f$; that was not traced further.

Limits: 1D planar sources only on this page; two linear amplitudes (1 Pa and 1e3 Pa on 1e5 Pa); the square and broadband pulses (`pulse = 3, 4`) and the 2D and 3D planar supports are not checked. Multifluid sources are the subject of [#1921](https://github.com/MFlowCode/MFC/issues/1921).

## Focused transducer

The axisymmetric focused transducer (`support = 6`, `examples/2D_acoustic_support6_axisym`) is on @ref validation_coordinates: its focal peak pressure converges to the exact on-axis solution of a spherical cap (−0.4 % at the finest grid) and arrives on time (the 1.5 \f$\Delta t\f$ lag measured before the fix was \f$\Delta t/2\f$ of source time and \f$\Delta t\f$ of probe time axis in the checker).

## Reproduce {#sec-validation-acoustics-run}

```
./mfc.sh run examples/1D_acoustic_gauss_sigmatime/case.py -t pre_process simulation post_process -n 1 -- -N 400
python3 examples/1D_acoustic_gauss_sigmatime/check_acoustic.py examples/1D_acoustic_gauss_sigmatime --case-args "-N 400"
```

`case.py` also takes `--mag` (the `mag` = 1e3 Pa column), `--wenom` (WENO5-M instead of TENO) and `--dt-div q` (dt/q at the
same end and save times). For the time-step refinement, run copies of the directory with `-N 800 --mag 1000 --dt-div q` for
q = 1, 2, 4, 8 and pass them in that order to `check_acoustic.py --dt-sweep`.
