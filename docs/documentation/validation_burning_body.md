@page validation_burning_body Burning carbon body

# Burning carbon body validation

**Status:** cylinder burn rate against a 1D radial reference: validated with `ib_ip_interp = 2` at 24 cells per diameter (+0.3 % at 10 us, +0.2 % at 20 us), +0.8 to +1.4 % at 12; with the default image points +0.6 to +0.9 % at 24 and +1.3 to +1.9 % at 12; sphere: partially validated with `ib_ip_interp = 2` (+2.5 to +4 % at 12 cells per diameter), +6 to +9 % with the default. **Code:** `master` at 0cdfdf02a. **Example:** `examples/nD_ibm_carbon_burn`. See also @ref validation.

A carbon body at fixed temperature gasifies in hot, high-pressure moist air with a CO flame off the surface. The reference is an independent 1D radially symmetric transient solution of the same model (`ref1d.py`, shipped), so the comparison checks the immersed-boundary surface chemistry (`ib_flux_wall`, the surface solve's species fluxes through the wall faces) and the reacting transport around a staircase body, with the mechanisms (`2D_ibm_reacting_surface`) shared. Tables marked "earlier" were measured on a pre-merge build and were not repeated.

## Setup {#sec-validation-burning-body-setup}

| Item | Value |
| :--- | :--- |
| Body | fixed cylinder (2D) or sphere (3D), D = 100 um (R = 50 um), at the centre of a box of half-width 6 R; the body appears at t = 0 |
| Gas | quiescent moist air (mass fractions N2 0.7606, O2 0.233, H2O 0.0058, CO2 0.00064), 10 atm, \f$T_\infty\f$ = 2000 K; reduced 11-species GRI gas mechanism, reactions on, mixture-averaged diffusion and conduction, viscous |
| Wall | isothermal (`thermal_bc = 1`) at 2000 K, surface reactions on (Bradley/Chelliah carbon mechanism, `surface_reaction = 1`), `ib_flux_wall = T`, no slip |
| Boundaries | extrapolation on all faces (the diffusion front reaches about 2 R by 60 us) |
| Numerics | mapped WENO5, HLLC, SSP-RK3, `weno_Re_flux`, adaptive CFL 0.5 with `transport_dt`, split alpha-QSS reactions (2 substeps, adaptive to 16) |
| Reference | `ref1d.py`: finite volumes on r in [R, 200 R] with geometric stretching (N = 400, first cell 0.002 R), cell unknowns T and Y; Cantera thermodynamics, mixture-averaged transport and kinetics; the wall receives the surface-solve fluxes; BDF time integration. Tabulated in `ref_cylinder.txt`, `ref_sphere.txt` |
| Measured | `compare_burn_rate.py`: `D/ib_heat.dat` column 4 (face mass outflow of the flux wall) and column 3 (surface-solve rate) against the reference \f$2\pi R\sum_k W_k s_k\f$ (cylinder) or \f$4\pi R^2\sum_k W_k s_k\f$ (sphere), averaged over \f$\pm\f$0.5 us; flame peak temperature and radius from the temperature binned in r |

### Cylinder

Burn rate \f$\dot m'\f$ [kg/(m s)] against the reference (default settings of `master`) (column 4 of `ib_heat.dat`, error = MFC/reference - 1; column 3 in parentheses). Measured on `master` at 0cdfdf02a with NVHPC 25.11 on an NVIDIA A100 (12 cells per diameter to 60 us, 24 cells per diameter to 30 us); the 12-cell run on the CPU (NVHPC 25.11) agrees with the A100 to 1e-5 relative and the 24-cell run at 10 us to 1e-4. The right-hand columns are the errors of the earlier study (pre-merge build whose defaults differed from `master`: triangular-shaped-cloud image points and the compact viscous face gradient).

| t [us] | reference | 12 cells/D | 24 cells/D | earlier, 12 / 24 |
| ---: | ---: | ---: | ---: | ---: |
| 10 | 1.4312e-4 | +1.87 % (+2.02 %) | +0.86 % (+0.92 %) | +0.92 % / +0.16 % |
| 20 | 1.1045e-4 | +1.36 % (+1.48 %) | +0.64 % (+0.69 %) | +0.08 % / -0.11 % |
| 30 | 9.680e-5 | +1.34 % (+1.46 %) | +0.93 % (+0.98 %) | -0.02 % / -0.14 % |
| 60 | 7.876e-5 | +1.39 % (+1.50 %) | not run | +0.03 % / -0.03 % |

Flame: the peak of the temperature binned in r is 2404.5 K at r = 2.09 R at 60 us on 12 cells per diameter (reference 2399.1 K at 2.04 R) and 2308.0 K at 1.80 R at 30 us on 24 cells per diameter (reference 2305.6 K at 1.79 R).

The cause of the gap to the earlier study was traced by rerunning the 12-cell cylinder with each default set back one at a time (CPU, 10 and 20 us; a run is compared at times at least 0.5 us before its end, since the averaging window of `compare_burn_rate.py` is cut off at the last row):

| Build and setting (12 cells/D) | 10 us | 20 us |
| :--- | ---: | ---: |
| earlier study | +0.92 % | +0.08 % |
| the earlier study's build (pre-merge, `ib_ip_interp = 2` by default), this case, GNU 13 CPU | +0.92 % | +0.08 % |
| `master`, defaults | +1.87 % | +1.36 % |
| `master`, `ib_ip_interp = 2` (TSC image points) | +1.39 % | +0.80 % |
| `master`, `ib_smooth_contact = T` | +1.87 % | not run |
| `master`, `ib_mass_relax = 30` or `1` | +1.86 % | not run |
| `master`, `weno_Re_flux_compact = T` | +1.86 % | +1.85 % (end of run) |
| `master` before #2248 (wall conduction through the unbounded mirror temperature) | burn rate equal to `master` to 5 digits through 5 us | |

The earlier study's build reproduces its numbers exactly with this case file, so the case and the reference are not the difference. The default that moves the burn rate is the image-point interpolation: the earlier build used triangular-shaped-cloud weights (`ib_ip_interp = 2`), `master` uses inverse-distance weights (`ib_ip_interp = 1`). Setting it recovers about half of the gap at 12 cells per diameter and nearly all of it at 24 (below). The other options of the earlier build do not matter here (`ib_smooth_contact`, `ib_mass_relax`, which the flux wall makes redundant, `weno_Re_flux_compact`). What remains at 12 cells per diameter (+0.5 and +0.7 percentage points at 10 and 20 us) comes from changes of the earlier build that are not in `master` (for example the mirrored temperature of the advected ghost of an isothermal or reacting wall); this was not isolated. The recommended setting for these cases is `--set ib_ip_interp=2`.

With `ib_ip_interp = 2` (A100, NVHPC 25.11):

| t [us] | 12 cells/D | 24 cells/D | earlier, 12 / 24 |
| ---: | ---: | ---: | ---: |
| 10 | +1.39 % | +0.32 % | +0.92 % / +0.16 % |
| 20 | +0.80 % | +0.16 % | +0.08 % / -0.11 % |

### Sphere

Same setup in 3D (`--ndim 3`), 12 cells per diameter (72 x 72 x 72 cells) to 30 us, `master` at 0cdfdf02a, NVHPC 25.11 on an NVIDIA A100, reference in kg/s:

| t [us] | reference | MFC, column 4 | error | column 3 | error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 1.8906e-8 | 2.0296e-8 | +7.35 % | 2.0378e-8 | +7.78 % |
| 20 | 1.5750e-8 | 1.6730e-8 | +6.22 % | 1.6785e-8 | +6.57 % |
| 30 | 1.4490e-8 | 1.5416e-8 | +6.39 % | 1.5464e-8 | +6.72 % |

With TSC image points (`ib_ip_interp = 2`, the recommended setting) the error is +4.05 % at 5 us, +3.17 % at 10 us and +2.49 % at 20 us (default: +9.14 %, +7.35 %, +6.22 %), against +2.07 %, +1.14 % (10, 20 us) of the earlier study (pre-merge build with `ib_ip_interp = 2` as the default and an older flux-wall face rule). The flame peak is 2247.5 K at r = 1.59 R at 30 us with the default settings (reference 2236.1 K at 1.61 R). The remaining 1 to 1.3 percentage points are not isolated (see the cylinder). Finer grids were not run (a 24-cell sphere is 144^3 cells). Status: partially validated with `ib_ip_interp = 2` (within 3.2 % at 12 cells per diameter), open discrepancy with the default image points.

## Known limitations {#sec-validation-burning-body-limits}

- The 1D reference converges at order 1.3 in its own grid (stretched mesh, first-order wall closure); its N = 400 cylinder value is within 0.08 % of its Richardson limit and the box (6 R and 8.5 R) changes it by 0.1 %, which bounds how far MFC can be checked (earlier).
- Element conservation of the flux wall (closed-box O, H and N inventories to round-off) is checked by `2D_ibm_wall_mass_budget` (see @ref validation); the 64-particle array budget of the earlier study is not repeated here.
- The film-theory series and the stagnation-flow comparisons (Makino and co-workers) of the earlier study are not part of this page: they need the stagnation reference of a separate study and were not re-measured.
