@page validation_burning_body Burning carbon body

# Burning carbon body validation

**Status:** cylinder burn rate against a 1D radial reference: partially validated (within 2 % at 12 cells per diameter and 1 % at 24, converging with the grid; the earlier study's 0.2 % at 24 cells per diameter is not reproduced on `master`); sphere: open discrepancy on `master` at 12 cells per diameter (+6 to +7 %). **Code:** `master` at 0cdfdf02a. **Example:** `examples/nD_ibm_carbon_burn`. See also @ref validation.

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

Burn rate \f$\dot m'\f$ [kg/(m s)] against the reference (column 4 of `ib_heat.dat`, error = MFC/reference - 1; column 3 in parentheses). Measured on `master` at 0cdfdf02a with NVHPC 25.11 on an NVIDIA A100 (12 cells per diameter to 60 us, 24 cells per diameter to 30 us); the 12-cell run on the CPU (NVHPC 25.11) agrees with the A100 to 1e-5 relative and the 24-cell run at 10 us to 1e-4. The right-hand columns are the errors of the earlier study (pre-merge build whose defaults differed from `master`: triangular-shaped-cloud image points and the compact viscous face gradient).

| t [us] | reference | 12 cells/D | 24 cells/D | earlier, 12 / 24 |
| ---: | ---: | ---: | ---: | ---: |
| 10 | 1.4312e-4 | +1.87 % (+2.02 %) | +0.86 % (+0.92 %) | +0.92 % / +0.16 % |
| 20 | 1.1045e-4 | +1.36 % (+1.48 %) | +0.64 % (+0.69 %) | +0.08 % / -0.11 % |
| 30 | 9.680e-5 | +1.34 % (+1.46 %) | +0.93 % (+0.98 %) | -0.02 % / -0.14 % |
| 60 | 7.876e-5 | +1.39 % (+1.50 %) | not run | +0.03 % / -0.03 % |

Flame: the peak of the temperature binned in r is 2404.5 K at r = 2.09 R at 60 us on 12 cells per diameter (reference 2399.1 K at 2.04 R) and 2308.0 K at 1.80 R at 30 us on 24 cells per diameter (reference 2305.6 K at 1.79 R).

On `master` the burn rate sits 1.3 to 1.9 % above the reference at 12 cells per diameter and 0.6 to 0.9 % above it at 24, falling by about 2 between the grids at 10 us (+1.87 to +0.86 %), where the earlier study measured -0.1 to +0.9 %. Two options of the earlier build were tried at 12 cells per diameter (CPU, to 20 us): triangular-shaped-cloud image points (`ib_ip_interp = 2`) give +1.39 % and +1.25 % at 10 and 20 us, and the compact viscous face gradient (`weno_Re_flux_compact = T`) +1.86 % and +1.85 %; neither closes the gap to the earlier numbers. The cause of the remaining difference is not identified.

### Sphere

Same setup in 3D (`--ndim 3`), 12 cells per diameter (72 x 72 x 72 cells) to 30 us, `master` at 0cdfdf02a, NVHPC 25.11 on an NVIDIA A100, reference in kg/s:

| t [us] | reference | MFC, column 4 | error | column 3 | error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 1.8906e-8 | 2.0296e-8 | +7.35 % | 2.0378e-8 | +7.78 % |
| 20 | 1.5750e-8 | 1.6730e-8 | +6.22 % | 1.6785e-8 | +6.57 % |
| 30 | 1.4490e-8 | 1.5416e-8 | +6.39 % | 1.5464e-8 | +6.72 % |

The flame peak is 2247.5 K at r = 1.59 R at 30 us (reference 2236.1 K at 1.61 R). The earlier study (pre-merge build, one flux-wall share rule not yet merged) measured the burn rate within 0.9 to 2 % of the reference at 12 to 16 cells per diameter and 0.3 % at 24; this is not reproduced on `master`: the 3D burn rate is 6 to 7 % above the reference at 12 cells per diameter, with columns 3 and 4 within 0.6 % of each other. Finer grids and the options of the earlier build were not run (a 24-cell sphere is 144^3 cells). Status: open discrepancy on `master` at this resolution; the cylinder shows the same sign with a smaller error.



## Known limitations {#sec-validation-burning-body-limits}

- The 1D reference converges at order 1.3 in its own grid (stretched mesh, first-order wall closure); its N = 400 cylinder value is within 0.08 % of its Richardson limit and the box (6 R and 8.5 R) changes it by 0.1 %, which bounds how far MFC can be checked (earlier).
- Element conservation of the flux wall (closed-box O, H and N inventories to round-off) is checked by `2D_ibm_wall_mass_budget` (see @ref validation); the 64-particle array budget of the earlier study is not repeated here.
- The film-theory series and the stagnation-flow comparisons (Makino and co-workers) of the earlier study are not part of this page: they need the stagnation reference of a separate study and were not re-measured.
