@page validation_acoustics Validation: Acoustic Sources

# Acoustic sources: planar pulses and a focused transducer

**Status:** partially validated. The amplitude and shape of the emitted pulses converge to linear acoustics, but every pulse lags the prescribed waveform by \f$\Delta t/2\f$ (open discrepancy, [#2462](https://github.com/MFlowCode/MFC/issues/2462)). **Code:** `master` at 5838b641c. **Examples:** `examples/1D_acoustic_gauss_sigmatime`, `examples/1D_acoustic_gauss_sigmadist`, `examples/1D_acoustic_sine_frequency`, `examples/1D_acoustic_dipole`. See also @ref validation and @ref validation_coordinates (the same pulse on a stretched grid).

All tables were measured with GNU 13.3 (CPU, one rank, double precision) on `master` at 5838b641c. On an NVIDIA A100 (NVHPC 25.11, OpenACC) the 400-cell Gaussian run gives \f$\tau/\Delta t\f$ = −0.499, the same L2 errors to three digits and A − 1 = −1.7e-4 (GNU −1.3e-4), and the dipole deck A − 1 = −5.8e-4 (GNU −5.6e-4) with the other columns unchanged. The CPU-GPU spread of A − 1 is about 5e-5.

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
| Reference | the expression above at the last save with the whole signal inside the domain (t = 500 ns; 324 ns for the dipole) |
| Measured | `check_acoustic.py`: relative L2 error of \f$p'\f$ outside 5 support widths; amplitude A and shift \f$\tau\f$ of the least-squares fit \f$p'_\text{MFC} = A\,p'_\text{lin}(t + \tau)\f$; the L2 error left after the fit |

Gaussian pulse (`1D_acoustic_gauss_sigmatime`), grid sweep at fixed \f$c\Delta t/\Delta x\f$:

| cells | \f$\Delta t\f$ [ns] | A − 1 | \f$\tau\f$ [ns] | \f$\tau/\Delta t\f$ | rel. L2 error | order | rel. L2 after the fit |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 2.0 | −8.2e-4 | −0.992 | −0.496 | 3.45e-2 | | 2.5e-3 |
| 400 | 1.0 | −1.3e-4 | −0.499 | −0.499 | 1.75e-2 | 0.98 | 3.3e-4 |
| 800 | 0.5 | −2.7e-5 | −0.249 | −0.498 | 8.8e-3 | 1.00 | 9.6e-5 |
| 1600 | 0.25 | −6.9e-6 | −0.125 | −0.499 | 4.4e-3 | 1.00 | 8.8e-5 |

The other decks, at their own resolution (200 cells):

| example | A − 1 | \f$\tau/\Delta t\f$ | rel. L2 error | rel. L2 after the fit |
| :--- | ---: | ---: | ---: | ---: |
| `1D_acoustic_gauss_sigmadist` | −8.2e-4 | −0.496 | 3.45e-2 | 2.5e-3 |
| `1D_acoustic_dipole` (sine, both directions) | −5.6e-4 | −0.498 | 3.46e-2 | 1.5e-3 |
| `1D_acoustic_sine_frequency` (15 cells per wavelength) | −4.8e-2 | −0.443 | 1.29e-1 | 4.6e-2 |

A − 1 falls at second order or faster (2.7, 2.2, 2.0) and the L2 error left after removing the shift at about third order down to a floor near 1e-4.
The shift does not converge in units of the time step: every pulse arrives \f$\Delta t/2\f$ late, which makes the plain L2 error first order.
The source is evaluated at the time of the start of the step in all three Runge-Kutta stages (`sim_time = mytime`, and `mytime` advances after `s_tvd_rk`), so its effective time is \f$t^n\f$ instead of the stage average \f$t^n + \Delta t/2\f$ (filed as [#2462](https://github.com/MFlowCode/MFC/issues/2462)).
The \f$-4.8\f$ % amplitude of the sine deck is the dissipation of WENO5 at 15 cells per wavelength over 9 wavelengths of travel; it was not refined (the deck has no resolution argument).
The dipole sends one pulse each way with the amplitude `mag`, as the doubled mass and energy terms predict.

Limits: 1D planar sources only on this page; a single linear amplitude (1 Pa on 1e5 Pa); the square and broadband pulses (`pulse = 3, 4`) and the 2D and 3D planar supports are not checked. Multifluid sources are the subject of [#1921](https://github.com/MFlowCode/MFC/issues/1921).

## Focused transducer

The axisymmetric focused transducer (`support = 6`, `examples/2D_acoustic_support6_axisym`) is on @ref validation_coordinates: its focal peak pressure converges to the exact on-axis solution of a spherical cap (−0.4 % at the finest grid), with a 1.5 \f$\Delta t\f$ lag.

## Reproduce {#sec-validation-acoustics-run}

```
./mfc.sh run examples/1D_acoustic_gauss_sigmatime/case.py -t pre_process simulation post_process -n 1 -- -N 400
python3 examples/1D_acoustic_gauss_sigmatime/check_acoustic.py examples/1D_acoustic_gauss_sigmatime --case-args "-N 400"
```
