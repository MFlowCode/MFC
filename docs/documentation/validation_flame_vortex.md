@page validation_flame_vortex Flame-vortex interaction

# Flame-vortex interaction validation

**Status:** 2D vortex pair: partially validated (self-convergence with the characteristic outflow reproduces the earlier study to 6e-4; flat-flame speed +0.33 % against Cantera; the extrapolation outflow is first order). 3D vortex ring: verified by self-convergence and the flat-flame speed in the earlier pre-merge study; on `master` only a 32^3 run was made, to check the case runs. **Code:** `master` at 0cdfdf02a. **Examples:** `examples/2D_flame_vortex_pair`, `examples/3D_flame_vortex_ring`. See also @ref validation and @ref validation_gas_chemistry.

A laminar premixed H2/air flame wrapped and pierced by a vortex (the flame-vortex interaction of \cite Poinsot1991) couples transport, chemistry and flow with no immersed boundary. No external reference exists for the interaction, so the check is self-convergence of integral metrics (heat release, H2O mass, flame length or area) at a fixed time step, plus the flat-flame speed against Cantera. Tables marked "earlier" were measured on a pre-merge build and were not repeated.

## Setup {#sec-validation-flame-vortex-setup}

| Item | Value |
| :--- | :--- |
| Geometry | square (2D) or cube (3D) of side 3.2 mm; transverse directions periodic; slip wall at x = 0 |
| Gas | stoichiometric H2/air, 1 atm, \f$T_u\f$ = 300 K, `h2o2.yaml` (10 species, 29 reactions), reactions on, mixture-averaged diffusion and conduction, viscous. Cantera \f$S_L\f$ = 2.3315 m/s, thermal thickness 0.33 mm |
| Flame | Cantera `FreeFlame` profile (`flame_profile.txt`) in the frame of the unburnt gas (reactants at rest), T = 1200 K at x = 1.8 mm |
| Vortex | 2D pair: two Lamb-Oseen vortices of circulation \f$\pm\f$0.05 m\f$^2\f$/s, core radius 0.25 mm, centres 0.8 mm apart, y-periodic image sums, cyclostrophic pressure, starting at x = 0.9 mm. 3D ring: radius 0.4 mm, Gaussian core 0.15 mm, circulation 0.022 m\f$^2\f$/s, velocity from the exact axisymmetric Stokes stream function of the core |
| Boundaries | x = 0 slip wall (`bc_x%%beg = -2`); x = L characteristic subsonic outflow (`bc_x%%end = -8`, `bc_x%%grcbc_out`, `bc_x%%pres_out` = 1 atm) |
| Numerics | mapped WENO5, HLLC, SSP-RK3, `weno_Re_flux`, split alpha-QSS reactions (2 substeps, adaptive to 16); fixed time step on every grid (2.632e-9 s in 2D, 5e-9 s in 3D), so the first-order split-chemistry error does not change with the grid |
| Measured | `analyze_fv.py`, `analyze_fv3d.py`: Q = \f$\int \dot q\,dA\f$ (Cantera rates at each cell state), H2O mass, length (2D) or area (3D) of the T = 1200 K isoline; observed order \f$p = \log_2((f_2 - f_1)/(f_3 - f_2))\f$ and Richardson value from three grids |

## 2D vortex pair {#sec-validation-flame-vortex-2d}

Characteristic outflow (`bc_x%%end = -8`), fixed time step. Measured on `master` at 0cdfdf02a; "earlier" is the pre-merge study with the same case.

| Quantity | t [us] | 64 x 64 | 128 x 128 | 256 x 256 | earlier (64, 128, 256, 512) |
| :--- | ---: | ---: | ---: | ---: | :--- |
| Q [W/m] | 40 | 18674.5 | 18660.7 | 18632.6 | 18675.1, 18662.2, 18632.8, 18618.6 |
| | 70 | 23317.2 | 23367.4 | 23401.6 | 23332.1, 23367.6, 23403.3, 23411.9 |
| H2O mass [1e-7 kg/m] | 40 | 1.52564 | 1.52500 | 1.52421 | 1.525664, 1.525013, 1.524224, 1.523798 |
| | 70 | 1.61155 | 1.61153 | 1.61110 | 1.611683, 1.611551, 1.611127, 1.610895 |
| flame length [mm] | 40 | 3.5284 | 3.5239 | 3.5243 | 3.52840, 3.52388, 3.52434, 3.52454 |
| | 70 | 5.0647 | 5.0801 | 5.0845 | 5.06566, 5.08060, 5.08497, 5.08562 |

The 64 x 64 run to 100 us gives Q = 32008.3 W/m, H2O mass 1.73517e-7 kg/m and flame length 7.14806 mm (earlier 32013.4, 1.735351e-7, 7.14902). Every value agrees with the earlier study to 6e-4 or better. The 64, 128 and 256 runs are the 64 x 64 and 128 x 128 runs on the CPU (NVHPC 25.11) and the 256 x 256 run on an NVIDIA A100 (NVHPC 25.11), all in double precision.

Observed orders from 64, 128 and 256 are not in the asymptotic range (the heat release at 40 us is non-monotone in the grid: -1.0 at 40 us, +0.6 at 70 us; the flame length at 70 us gives 1.8), as in the earlier study, where the orders 1.0 (40 us) and 2.1 (70 us) came from the three finest grids (128, 256, 512). The 512 x 512 grid was not run here.

Flat flame (`--N 256 --ny 32 --gamma 0`, a 256 x 32 strip at the 256 x 256 spacing, 26.6 thousand steps to 70 us, NVIDIA A100, NVHPC 25.11): the displacement speed over 30 to 70 us (41 saves one microsecond apart, x_f fitted linearly, u_a and rho_a at T_up + 10 K) is S_d = 2.3393 m/s against the Cantera S_L = 2.3315 m/s, +0.33 % (earlier: +0.07 % over 30 to 90 us, +0.05 % in 3D).

## 3D vortex ring {#sec-validation-flame-vortex-3d}

The ring case (`examples/3D_flame_vortex_ring`) was run here only at 32^3 cells to 20 us on an NVIDIA A100 (NVHPC 25.11, MPI build, `pre_process`, `ic3d.py`, `simulation`, `post_process`) to check that it runs end to end: at 20 us the H2O mass is 4.73474e-10 kg, the flame area 1.02521e-5 m^2 and the heat release 59.59 W. The 64^3 and 128^3 convergence runs were not completed: on the shared development machine they would have taken longer than the 100 us of physical time justifies (about 0.8 s per step at 32^3 on the A100 under contention, 20000 steps). The table is from the earlier pre-merge study (fixed dt = 5e-9 s, 64^3, 128^3 on one machine and 256^3 on another; the 256^3 run agrees with a run of the same case on a second GPU architecture to 2e-13 at 52 us, so the result does not depend on the architecture), and its numbers were not re-measured.

| Quantity | t [us] | 64^3 | 128^3 | 256^3 | order | Richardson | 256^3 vs R |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| H2O mass [1e-10 kg] | 40 | 4.88629 | 4.87827 | 4.87602 | 1.84 | 4.87515 | +1.8e-4 |
| | 70 | 5.11127 | 5.10359 | 5.10228 | 2.55 | 5.10200 | +5.3e-5 |
| | 100 | 5.37166 | 5.36390 | 5.36204 | 2.07 | 5.36146 | +1.1e-4 |
| flame area [1e-5 m^2] | 40 | 1.03664 | 1.03586 | 1.03585 | (6.1) | 1.03585 | +1.7e-7 |
| | 70 | 1.14723 | 1.14843 | 1.14873 | 2.00 | 1.14883 | -8.7e-5 |
| | 100 | 1.31105 | 1.31249 | 1.31256 | (4.3) | 1.31257 | -3.1e-6 |
| heat release [W] | 40 | 57.809 | 57.326 | 57.339 | non-monotone | | |
| | 70 | 60.349 | 60.509 | 60.641 | 0.28 | 61.26 | -1.0 % |
| | 100 | 67.187 | 67.712 | 67.884 | 1.61 | 67.97 | -0.12 % |

The H2O mass converges at observed order 1.8 to 2.6 and the 256^3 value is within 2e-4 of its Richardson value at every time. The flame area converges at order 2.0 at 70 us, where the flame is most wrinkled (at 40 and 100 us it is converged to 1e-5 already on 128^3, so the order is not meaningful). The instantaneous heat release is not in the asymptotic range on 64^3 (6.6 cells per flame thickness). A flat-flame control (256 x 32 x 32 cells, no ring) gave a displacement speed within 0.05 % of Cantera (earlier). Status: verified by self-convergence in the earlier study, not re-measured on `master`.


## Known limitations {#sec-validation-flame-vortex-limits}

- The `bc_x%%end = -3` extrapolation outflow lets the pressure of this wall-closed domain drift upward by a grid-dependent, first-order amount (upstream excess pressure at 40 us of 1878, 1338, 1030 and 866 Pa at 64 to 512 cells, earlier), which turns the convergence of Q and the H2O mass first order. Use `bc_x%%end = -8` with `bc_x%%grcbc_out` for a closed domain with an outlet.
- The y-period (3.2 mm) is four pair separations: the pair's images strain the whole flame, so the vortex run's consumption speed is stretch-affected; the flat strip is the flame-speed check.
- The split chemistry is first order in the time step (Lie); its error is the same on all grids here (fixed step) and was 0.05 to 0.1 % in Q at \f$\Delta t\f$ = 1.2e-8 s (earlier, 2D).
- The window consumption speed of the flat strip reads 0.9 to 1.0 % below Cantera on every grid (earlier). It is converged in the grid, so it is the relaxation of the Cantera profile toward MFC's flame within 0.6 flame times, not discretization error (not verified by a longer run).
- The instantaneous heat release is not in the asymptotic range on 64 cells (6.6 cells per flame thickness): it is non-monotone in 3D at 40 us.
- MFC needs at least 25 cells per WENO direction, so the flat strips are 32 cells wide.
- The 3D case needs three commands (the ring is written into the restart file by `ic3d.py`) and an MPI build; it is not in the test suite.
