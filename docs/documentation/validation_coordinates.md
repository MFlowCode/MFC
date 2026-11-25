@page validation_coordinates Validation: Axisymmetric Focusing and Stretched Grids

# Coordinates and grids: axisymmetric focusing and a stretched 1D grid

**Status:** partially validated. The focal pressure of an axisymmetric transducer converges to the exact on-axis solution (−0.4 % at the finest grid), but it arrives 1.5 \f$\Delta t\f$ late and 4 % of the waveform is left after removing that shift. On a stretched grid the 1D pulse is as accurate as on the uniform grid, with the same \f$\Delta t/2\f$ lag ([#2462](https://github.com/MFlowCode/MFC/issues/2462)). **Code:** `master` at 5838b641c. **Examples:** `examples/2D_acoustic_support6_axisym`, `examples/1D_acoustic_gauss_sigmatime`. See also @ref validation and @ref validation_acoustics.

## Axisymmetric focused transducer {#sec-validation-coordinates-focus}

| Item | Value |
| :--- | :--- |
| Medium | water, stiffened gas \f$\gamma\f$ = 4.4, \f$\pi_\infty\f$ = 5.57e8 Pa, 1100 kg/m³, 1e5 Pa; c = 1492.8 m/s |
| Geometry | `cyl_coord = T`, \f$x \in [0, 1]\f$ mm, \f$r \in [0, 0.5]\f$ mm, axis at r = 0 (`bc_y%%beg = -2`), non-reflecting elsewhere |
| Source | `support = 6`: spherical cap, radius of curvature F = 0.54 mm, aperture 2a = 0.8 mm (depth h = 0.177 mm), apex at x = 0.06 mm; Gaussian pulse \f$\sigma_t\f$ = 20 ns, delay 100 ns, `mag` = 1 Pa |
| Numerics | WENO5 with TENO, HLLC, SSP-RK3; 100 x 50 cells (the example) times k = 1, 2, 4, 8, \f$\Delta t\f$ = 2 ns / k |
| Reference | exact on-axis linear field of a cap with uniform normal velocity (O'Neil \cite ONeil1949): \f$p(z,t) = \frac{F}{F - z}[p_0(t - z/c) - p_0(t - R_e/c)]\f$, \f$R_e = \sqrt{a^2 + (z - h)^2}\f$, z from the apex; at the focus \f$(h/c)\,dp_0/dt\f$, peaks \f$\pm\,\mathrm{mag}\,h/(c\sigma_t\sqrt{e})\f$ = ±3.60 Pa |
| Measured | `check_focus.py`: pressure of a probe in the cell just past the focus, every step (`--probe`); positive and negative peak values and times, relative L2 error of the signal, least-squares time shift and the L2 error left after it |

| k | cells per \f$c\sigma_t\f$ | positive peak [Pa] | error | negative peak [Pa] | error | shift [ns] | shift / \f$\Delta t\f$ | rel. L2 | rel. L2 after the shift |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.0 | 2.947 | −17.5 % | −2.986 | −16.4 % | −2.59 | −1.30 | 0.252 | 0.178 |
| 2 | 6.0 | 3.425 | −4.5 % | −3.428 | −4.4 % | −1.47 | −1.47 | 0.110 | 0.059 |
| 4 | 11.9 | 3.554 | −1.1 % | −3.557 | −1.0 % | −0.76 | −1.52 | 0.062 | 0.041 |
| 8 | 23.9 | 3.584 | −0.37 % | −3.585 | −0.35 % | −0.38 | −1.51 | 0.046 | 0.039 |

The exact peak at the probe (0.62 to 5 µm past the focus) is 3.572 to 3.597 Pa; the errors are against that value.
k = 1, 2 and 4 were run with GNU 13.3 on the CPU (one rank) and k = 1, 2, 4 and 8 on an NVIDIA A100 (NVHPC 25.11, OpenACC); where both ran, the peaks agree to four digits and the L2 errors to three.
The peak pressure converges at second order (4.5 %, 1.1 %, 0.37 %), so the axisymmetric source with its curvature correction, the geometric source terms and the axis reproduce the focal gain of the cap.
Two differences remain. The focal signal lags by 1.5 \f$\Delta t\f$ at every resolution: \f$\Delta t/2\f$ of it is the step-start source time of [#2462](https://github.com/MFlowCode/MFC/issues/2462); the other \f$\Delta t\f$ was not traced (the probe row is taken at \f$(k+1)\Delta t\f$, the time of the state it is written from). After removing the shift, a 4 % L2 difference stays as the grid is refined; it sits in the tails of the signal, and whether it comes from the free (unbaffled) edge of the MFC source against the baffled cap of the O'Neil solution was not separated.

Limits: one transducer geometry and one probe point on the axis; the off-axis field, the element arrays (`support = 9, 10`) and the 3D transducers are not compared.

## Stretched grid {#sec-validation-coordinates-stretched}

The 1D Gaussian pulse of @ref validation_acoustics run on a stretched grid (`stretch_x`, `a_x` = 4, clustered on 0.3 to 0.6 mm, one loop): the same number of cells cover 0 to 1.68 mm with \f$\Delta x\f$ from 0.81 to 3.6 times the uniform spacing of 1 mm/N, and the largest change between neighbouring cells 0.8 %. The reference is the same linear solution, evaluated with the run's own cell widths, and \f$\Delta t\f$ is that of the uniform grid (c \f$\Delta t/\Delta x_\text{min}\f$ = 0.73).

| cells | \f$\Delta x\f$ range [µm] | A − 1 | \f$\tau/\Delta t\f$ | rel. L2 error | rel. L2 after the fit | uniform grid, same cells: A − 1 | uniform: rel. L2 after the fit |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 4.07-17.9 | −1.65e-3 | −0.488 | 3.39e-2 | 4.0e-3 | −8.2e-4 | 2.5e-3 |
| 400 | 2.03-8.96 | −3.0e-4 | −0.491 | 1.73e-2 | 5.7e-4 | −1.3e-4 | 3.3e-4 |
| 800 | 1.02-4.49 | −8.1e-5 | −0.494 | 8.7e-3 | 1.9e-4 | −2.7e-5 | 9.6e-5 |

The stretched grid converges like the uniform one (A − 1 at order 2.5 and 1.9, the L2 error after the fit at 2.8 and 1.6), with errors about twice the uniform-grid values, with cells up to 3.6 times coarser than the uniform spacing. The lag is \f$\Delta t/2\f$ as on the uniform grid. The 400-cell run on the A100 gives A − 1 = −3.5e-4 and the same lag and L2 errors (5.6e-4 after the fit).

Limits: one stretching and one direction; a pulse that crosses the clustered region once; no shocks.

## Reproduce {#sec-validation-coordinates-run}

```
./mfc.sh run examples/2D_acoustic_support6_axisym/case.py -t pre_process simulation -n 1 -- --refine 4 --probe
python3 examples/2D_acoustic_support6_axisym/check_focus.py examples/2D_acoustic_support6_axisym --case-args "--refine 4 --probe"
./mfc.sh run examples/1D_acoustic_gauss_sigmatime/case.py -t pre_process simulation post_process -n 1 -- -N 400 --stretch 4
python3 examples/1D_acoustic_gauss_sigmatime/check_acoustic.py examples/1D_acoustic_gauss_sigmatime --case-args "-N 400 --stretch 4"
```
