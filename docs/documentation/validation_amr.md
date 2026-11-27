@page validation_amr Adaptive mesh refinement

# AMR validation

AMR (@ref amr) must give the uniform fine-grid answer where the blocks cover the flow, keep the uniform grid's order of accuracy through coarse/fine seams, and must not create or destroy mass, momentum or energy at coarse/fine faces, in regrids or in subcycling. See also @ref validation.

The Sod and conservation checks read the fine data from the serial restart files (level 0 in `q_cons_vf*.dat`, the blocks in `amr_fine.dat`, via `mfc.viz.read_amr_restart`).

## Sod shock tube against the exact Riemann solution {#sec-validation-amr-sod}

**Status:** validated. **Example:** `examples/1D_amr_sod`.

| Item | Value |
| :--- | :--- |
| Problem | Sod: \f$(\rho, u, p)\f$ = (1, 0, 1) for x < 0.5, (0.125, 0, 0.1) for x > 0.5, \f$\gamma\f$ = 1.4, t = 0.2, outflow ends |
| Numerics | 5-equation model, one fluid, mapped WENO5, HLLC, SSP-RK3, fixed \f$\Delta t\f$ at CFL 0.4 on the cells that set it (finest level lock-step, base level subcycled) |
| AMR | 1 or 2 levels at 2:1, `amr_regrid_int` 2 (subcycled) or \f$2 \cdot 2^L\f$ (lock-step), `amr_tag_eps = 0.05` (shock and contact tagged, rarefaction not), `amr_buf = 8` (derived `amr_snap` 2 with one level, 0 with two), `amr_max_grid_size = N`, linear ghost fill |
| Compared | AMR on base N; uniform N; uniform at the finest AMR spacing (2N, 4N); N = 100, 200, 400 |
| Reference | exact Riemann solution (pressure function and Newton iteration, \cite Toro09 ch. 4), computed in `check_sod.py`: \f$p^*\f$ = 0.30313, \f$u^*\f$ = 0.92745, contact at 0.68549, shock at 0.85043 |
| Measured | on the composite grid (finest level covering each point): \f$L_1 = \sum_i \vert\rho_i - \bar\rho_i^{exact}\vert\,\Delta x_i\f$ with the exact cell average, and the contact and shock positions (where \f$\rho\f$ crosses the middle of its exact jump), by `check_sod.py` |

L1 density error, and in brackets (AMR − uniform fine)/uniform fine:

| N | uniform N | 1 level | 1 level, subcycled | uniform 2N | 2 levels | 2 levels, subcycled | uniform 4N |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 3.42e-3 | 1.89e-3 (+0.0 %) | 1.90e-3 (+0.1 %) | 1.89e-3 | 9.82e-4 (+4.9 %) | 9.87e-4 (+5.3 %) | 9.37e-4 |
| 200 | 1.89e-3 | 9.35e-4 (−0.2 %) | 9.36e-4 (−0.1 %) | 9.37e-4 | 4.87e-4 (+8.2 %) | 4.90e-4 (+9.0 %) | 4.50e-4 |
| 400 | 9.37e-4 | 4.43e-4 (−1.6 %) | 4.43e-4 (−1.6 %) | 4.50e-4 | 2.70e-4 (+2.9 %) | 2.71e-4 (+3.4 %) | 2.62e-4 |

Contact and shock position errors (x − x_exact), lock-step; the subcycled runs agree to 2e-6 except the N = 100 two-level contact (−4.03e-4):

| N | uniform N: contact, shock | 1 level | uniform 2N | 2 levels | uniform 4N |
| ---: | :--- | :--- | :--- | :--- | :--- |
| 100 | −2.17e-3, +1.25e-3 | −9.49e-4, +6.87e-4 | −9.49e-4, +6.87e-4 | −3.42e-4, +3.55e-4 | −4.31e-4, +3.55e-4 |
| 200 | −9.49e-4, +6.87e-4 | −4.29e-4, +3.55e-4 | −4.31e-4, +3.55e-4 | −2.84e-4, +2.33e-4 | −2.83e-4, +2.33e-4 |
| 400 | −4.31e-4, +3.55e-4 | −2.85e-4, +2.34e-4 | −2.83e-4, +2.33e-4 | −1.76e-4, +1.24e-4 | −1.75e-4, +1.26e-4 |

With one level the error is the uniform 2N error within 1.6 %, with two levels within 9 % of the uniform 4N error, on 40 to 48 % of its cells; the AMR error halves per doubling of N, as the uniform one does. The positions match the uniform fine grid to a fraction of a finest cell. Subcycling changes the L1 error by at most 0.7 %.

Found by this check, two levels with `amr_buf = 4`: the toolchain used to derive `amr_snap = 2` (regrid hysteresis) for every case, and with two levels the shock then ran ahead of the exact one by +3.8e-3, +3.3e-3 and +3.4e-3 at N = 100, 200, 400. That is a shock-speed error near 1 % that does not converge: L1 was 3.7 times the uniform 4N error at N = 200 and 5.9 times at N = 400. A level-2 box is clipped to its parent inset by `amr_cpat_mar` = 3 coarse cells, so at `amr_buf = 4` a box held by the snap let the shock reach the level-2 edge. With `amr_snap = 0` the same run matches the uniform fine grid (L1 4.87e-4, shock +2.32e-4 at N = 200). `amr_snap` is now derived only for one-level runs (#2474). Reproduce: `case.py -N 200 --levels 2 --buf 4`.

Limits: 1D, one problem, the linear ghost fill only (`amr_prolong = 2` is not run), one rank, `amr_tag_eps` and `amr_regrid_int` fixed. The rarefaction is not refined, so the two-level error contains the base grid's rarefaction error; that is part of the 8-9 % at N = 200. Only density is compared.

## Conservation in a periodic box {#sec-validation-amr-conservation}

**Status:** validated. **Example:** `examples/2D_amr_conservation`.

| Item | Value |
| :--- | :--- |
| Problem | periodic unit square; a square of \f$\rho = p = 2\f$ (edges on base faces, side 1/8) in \f$\rho = p = 1\f$, uniform velocity (0.5, 0.25), \f$\gamma\f$ = 1.4, t = 0.4: the contact advects, the pressure jump sends shocks across block edges and the periodic boundaries |
| Numerics | 5-equation model, mapped WENO5, HLLC, SSP-RK3, CFL 0.3, 64² base |
| AMR | 2 levels at 2:1, lock-step (684 steps, regrid every 8) and subcycled (172 steps, regrid every 2), `amr_tag_eps = 0.02`, `amr_buf = 4`, up to 52544 refined cells |
| Reference | totals at t = 0 (pre_process's level-0 data; the square's edges on base faces make it the blocks' starting average too) |
| Measured | \f$(S(t) - S(0))/S(0)\f$ for mass, both momenta and energy at four saves, summed exactly (`math.fsum`) over level 0; and the gap to the totals over the leaf cells of every level, by `check_conservation.py` |

Level 0 equals the composite grid because at the end of every step `s_amr_step_fold` (`m_amr_stage.fpp`, called at the end of `s_tvd_rk` in `m_time_steppers.fpp`) restricts the fine solution bottom-up into the covered cells; the checker confirms it at every save.

| Mode | mass | x-momentum | y-momentum | energy | level 0 − composite |
| :--- | ---: | ---: | ---: | ---: | ---: |
| lock-step | −6.6e-16 | 0 | 0 | −1.6e-16 | ≤ 4.4e-16 |
| subcycled | 0 | 0 | +2.2e-16 | +1.6e-16 | ≤ 2.2e-16 |

The largest drift over the four saves is given; all are one to three units in the last place of the totals.

Limits: 2D, one rank, one fluid, inviscid, ratio 2, two levels, linear ghost fill, periodic boundaries only; tagging skips the global boundary cells, so no block spans the periodic seam. Viscous, chemistry, multi-fluid, multi-rank and immersed-boundary conservation are not covered here.

### Reproduce {#sec-validation-amr-run}

```
./mfc.sh run examples/1D_amr_sod/case.py -t pre_process simulation -n 1 -- -N 200 --levels 2
python3 examples/1D_amr_sod/check_sod.py examples/1D_amr_sod
./mfc.sh run examples/2D_amr_conservation/case.py -t pre_process simulation -n 1 -- --subcycle
python3 examples/2D_amr_conservation/check_conservation.py examples/2D_amr_conservation
```

The checkers import `mfc.viz`, so run them in MFC's Python environment (`source build/venv/bin/activate` after any `./mfc.sh` command). `case.py -N <N> --levels 0` and `--fine` give the uniform grids; `--subcycle` the subcycled AMR runs. Each run needs its own directory (copy `case.py`), since the checker reads the restart files of the case directory. The suite runs `AMR -> 1D -> Example -> Sod two levels` (995DE2F6) and `AMR -> 2D -> Example -> conservation subcycle` (2352E14C).

## Smooth flow across block seams {#sec-validation-amr-vortex}

**Status:** validated (second order in L1 for every block layout and fill, after #2482). **Example:** `examples/2D_amr_isentropic_vortex`.

A smooth flow carried through coarse/fine seams has an exact solution only when nothing at the seams disturbs it, so the error at the end shows what the ghost fill, refluxing and restriction leave behind.

| Item | Value |
| :--- | :--- |
| Flow | isentropic vortex \cite Shu1998, \f$\epsilon = 5\f$, \f$\gamma = 1.4\f$, \f$\rho_\infty = p_\infty = 1\f$, carried by \f$u = 1\f$ from \f$x = -3\f$ to \f$x = 3\f$ (\f$t = 6\f$) on a periodic \f$[-10, 10]^2\f$ |
| Grids | base N = 128, 256, 512, 1024; one 2:1 level. Static: one block over \f$[-1.25, 1.25] \times [-2.5, 2.5]\f$ (3.1 % of the domain), which the vortex enters, crosses and leaves. Regrid: density-jump tags (\f$\epsilon_\text{tag} = 0.006\,\Delta x\f$, down to \f$r \approx 2.5\f$), `amr_buf = 3`, regrid every 10 steps (6.5-9 % of the domain refined) |
| Variants | lock-step and subcycled, `amr_prolong = 1` (linear) and `2` (quartic); uniform runs at N and 2N |
| Numerics | WENO5 (mapped, \f$\epsilon_W = 10^{-16}\f$), HLLC, SSP-RK3; one fine-level \f$\Delta t\f$ for all runs at a given N (CFL 0.4 on the 2N grid), so they differ only in space |
| Initial data | exact cell averages of the conserved variables on every level (3x3-point Gauss in the patch expressions) |
| Reference | the initial vortex translated by \f$ut\f$; exact cell averages of \f$\rho\f$ (4x4-point Gauss) |
| Measured | L1 (area mean) and Linf of \f$\rho\f$ on the composite grid (fine cells where a block covers, level 0 elsewhere) from `restart_data/lustre_<step>.dat` and `lustre_amr_<step>.dat`, by `check_vortex.py` |

L1 density error and observed order (512 → 1024):

| Run | N = 128 | 256 | 512 | 1024 | order |
| :--- | ---: | ---: | ---: | ---: | ---: |
| uniform N | 4.37e-5 | 1.11e-5 | 2.78e-6 | 6.95e-7 | 2.00 |
| uniform 2N | 1.11e-5 | 2.78e-6 | 6.95e-7 | 1.74e-7 | 2.00 |
| static, lock-step, linear | 1.08e-4 | 2.29e-5 | 4.95e-6 | 1.13e-6 | 2.13 |
| static, lock-step, quartic | 5.91e-5 | 1.56e-5 | 4.01e-6 | 1.03e-6 | 1.97 |
| static, subcycled, linear | 1.13e-4 | 2.35e-5 | 5.04e-6 | 1.14e-6 | 2.14 |
| static, subcycled, quartic | 5.91e-5 | 1.56e-5 | 4.01e-6 | 1.03e-6 | 1.97 |
| regrid, lock-step, linear | 1.40e-5 | 3.71e-6 | 9.84e-7 | 2.49e-7 | 1.98 |
| regrid, lock-step, quartic | 1.32e-5 | 3.39e-6 | 8.90e-7 | 2.32e-7 | 1.94 |
| regrid, subcycled, linear | 1.41e-5 | 3.87e-6 | 1.05e-6 | 2.70e-7 | 1.96 |
| regrid, subcycled, quartic | 1.32e-5 | 3.38e-6 | 8.89e-7 | 2.31e-7 | 1.94 |

Linf density error:

| Run | N = 256 | 512 | 1024 | order 512 → 1024 | on a seam at 1024 |
| :--- | ---: | ---: | ---: | ---: | :--- |
| uniform N | 8.26e-4 | 2.08e-4 | 5.21e-5 | 2.00 | – |
| uniform 2N | 2.08e-4 | 5.21e-5 | 1.30e-5 | 2.00 | – |
| static, lock-step, linear | 1.16e-3 | 2.34e-4 | 7.24e-5 | 1.69 | yes |
| static, lock-step, quartic | 6.17e-4 | 1.75e-4 | 5.57e-5 | 1.65 | yes |
| static, subcycled, linear | 1.13e-3 | 2.36e-4 | 7.27e-5 | 1.70 | yes |
| static, subcycled, quartic | 6.17e-4 | 1.75e-4 | 5.58e-5 | 1.65 | yes |
| regrid, lock-step, linear | 1.58e-4 | 8.49e-5 | 3.73e-5 | 1.19 | yes |
| regrid, lock-step, quartic | 1.60e-4 | 4.01e-5 | 1.82e-5 | 1.14 | yes |
| regrid, subcycled, linear | 1.64e-4 | 8.46e-5 | 2.77e-5 | 1.61 | yes |
| regrid, subcycled, quartic | 1.60e-4 | 4.02e-5 | 1.70e-5 | 1.24 | yes |

- **The uniform scheme is second order here**, in L1 and Linf, also from exact cell averages. The likely cause (not verified) is that MFC reconstructs primitive variables from cell averages of the conserved ones, which is second order when density, velocity and pressure all vary; the 1D tests of the fifth-order rate advect a density wave at constant velocity and pressure. The AMR rates are to be read against 2, not 5.
- **Blocks that follow the vortex:** L1 converges at 1.94-1.98 with either fill, 1.3 (quartic) to 1.6 (linear, subcycled) times the uniform 2N error. Linf converges more slowly at N = 1024 (1.1-1.6) and its largest value sits on a seam, 1.3-2.9 times the uniform 2N value.
- **Static block:** at t = 6 the vortex has crossed the block and sits on level 0, so the error is a level-0 error: L1 1.5-1.6 times and Linf 1.07-1.4 times the uniform N value, both converging (L1 at 2.0). Lock-step and subcycled agree to two digits.
- **Found by this check (#2473, fixed in #2482):** the volume fractions were refluxed on their flux alone, without the \f$\alpha\,\nabla\cdot u\f$ source that pairs with it, so a coarse cell beside a block face lost its uniform \f$\alpha\f$ where the velocity varies along the face. The static block then left a density error in the vortex core that did not converge (Linf 1.9e-3 at N = 1024, quartic; L1 at first order). Mass, momentum and energy were conserved throughout.

Limits: one vortex strength, one flow direction (aligned with the grid), one refinement ratio and one level. Linf on the composite grid compares cells of two sizes. The tag threshold scales with \f$\Delta x\f$ so the refined region is about the same at every N, but its edges move by a cell between resolutions, which adds noise to the seam Linf rates.

### Reproduce {#sec-validation-amr-vortex-run}

```
./mfc.sh run examples/2D_amr_isentropic_vortex/case.py -t pre_process simulation -n 1 -- --N 512 --mode static --prolong 2
./mfc.sh run examples/2D_amr_isentropic_vortex/case.py -t pre_process simulation -n 1 -- --N 512 --mode regrid --prolong 2 --subcycle
./mfc.sh run examples/2D_amr_isentropic_vortex/case.py -t pre_process simulation -n 1 -- --N 512 --mode uniform --fine
python3 examples/2D_amr_isentropic_vortex/check_vortex.py <run directories, increasing N>
```

Each run writes next to `case.py`, so copy the directory per run. Measured on one MI250X GCD per run (amdflang, OpenMP offload, double); the 40 runs take 3.6 minutes of one node (8 at a time), the largest (uniform 2048², 4608 steps) 140 s. The suite runs 40 steps at N = 64 as `AMR -> 2D -> Example -> isentropic vortex static block` (C409AA2D) and `... regrid subcycle` (D209FD0C, linear fill) as regression guards only: they start from the hardcoded vortex `hcid` 283 (\f$\epsilon = 3.5\f$, a different profile from Shu's) so that the default build runs them, and its fine blocks start from interpolation rather than their own cell averages. Goldens: MI250X GPU build (amdflang, double, Frontier); the Cray CPU build passes them at the default tolerance.

## Sedov-Taylor blast {#sec-validation-amr-sedov}

**Status:** validated (shock radius; AMR equals the uniform fine grid). **Example:** `examples/3D_amr_sedov`. The AMR runs are with #2482; it moves the radius by at most 0.05 percentage points (192³ base) and leaves the 96³ base unchanged.

The blocks follow a strong spherical shock through regrids. The self-similar solution gives the shock radius exactly, and an AMR run should match a uniform run at its fine resolution.

| Item | Value |
| :--- | :--- |
| Blast | energy \f$E = 1\f$ in a sphere of radius 0.04 at the center of \f$[-0.6, 0.6]^3\f$; \f$\rho_0 = 1\f$, \f$p_0 = 10^{-3}\f$, \f$\gamma = 1.4\f$, extrapolation boundaries; to \f$t = 0.16\f$ (\f$R = 0.50\f$, shock Mach number about 33) |
| Grids | base 96³ and 192³; one 2:1 level, absolute density-jump tags (`amr_tag_abs`, 0.1), `amr_buf = 4`, regrid every 4 steps, lock-step and subcycled; uniform 96³, 192³ and 384³ |
| Numerics | MUSCL-2 (minmod), HLLC, SSP-RK3, CFL 0.3 (adaptive) |
| Reference | \f$R = \xi_0 (E t^2/\rho_0)^{1/5}\f$ \cite TaylorBlast1950 \cite Sedov1959 with \f$\xi_0 = 1.0328\f$ for \f$\gamma = 1.4\f$, from the energy integral of the similarity equations (`check_sedov.py` integrates them; \f$\gamma = 5/3\f$ gives 1.1517); post-shock density \f$(\gamma + 1)/(\gamma - 1) = 6\f$ |
| Measured | per save, on the composite grid: \f$E\f$ from the data (total energy above the ambient), the shell-averaged density in shells one finest cell thick, \f$R\f$ at the outermost half-jump between \f$\rho_0\f$ and the peak, by `check_sedov.py` |

Shock radius error (R − R_exact)/R_exact, with \f$E\f$ measured in each run (0.9908 at 96³ and 192³ fine cells, 0.9954 at 384³):

| t | R exact | uniform 96³ | AMR 96³ lock-step | AMR 96³ subcycled | uniform 192³ | AMR 192³ lock-step | AMR 192³ subcycled | uniform 384³ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.02 | 0.216 | +5.60 % | +2.57 % | +2.57 % | +2.58 % | +1.36 % | +1.35 % | +1.38 % |
| 0.06 | 0.335 | +3.54 % | +1.92 % | +1.86 % | +1.85 % | +0.94 % | +0.94 % | +0.95 % |
| 0.10 | 0.411 | +2.96 % | +1.58 % | +1.58 % | +1.58 % | +0.78 % | +0.79 % | +0.79 % |
| 0.14 | 0.470 | +2.55 % | +1.35 % | +1.33 % | +1.35 % | +0.72 % | +0.74 % | +0.72 % |
| 0.16 | 0.495 | +2.50 % | +1.27 % | +1.27 % | +1.26 % | +0.74 % | +0.74 % | +0.73 % |

Peak shell-averaged density over \f$\rho_0\f$ at t = 0.16: 3.25 (uniform 96³), 3.97 (AMR 96³, both; uniform 192³), 4.56 (AMR 192³, both; uniform 384³), against 6.

- AMR matches the uniform grid at its fine resolution: on a 96³ base the radius is within 0.07 percentage points of uniform 192³ and the peak density within 0.5 % at every save (lock-step and subcycled); on a 192³ base within 0.03 percentage points of uniform 384³ (lock-step and subcycled; 0.01 at t = 0.16) and the peak density within 0.2 %.
- The radius error falls with the fine cell size, from +2.50 % (96³) to +1.26 % (192³) and +0.73 % (384³) at t = 0.16 (first order, as expected for a captured shock). It is positive: the half-jump of a shock smeared over 2-3 cells lies ahead of the discontinuity.
- The post-shock density is resolution-limited: the peak of the shell average is 66-76 % of the strong-shock value at these grids.

Limits:
- **Early times.** The energy sits in a finite sphere (0.04, 3-13 cells) with a staircase surface, which shifts the virtual origin; the first save (t = 0.02, R = 5.4 r0) carries that and the startup transient. The patch stair-step loses 0.9 % of the energy at 96³/192³ and 0.5 % at 384³; the reference uses the measured energy.
- **Counter-pressure.** \f$p_0 = 10^{-3}\f$ leaves the shock at Mach 33 at the end (\f$\dot R = 0.4 R/t\f$ against \f$c_0 = \sqrt{\gamma p_0/\rho_0}\f$); the strong-shock solution ignores \f$p_0\f$ (an effect of order \f$1/M^2\f$, about 0.1 %, not separated from the discretization error here).
- **Scheme.** On a uniform 256³ octant (no AMR), WENO5 with `mp_weno = T` goes NaN at t = 0.0024, at the collapsing centre in the reflective corner: WENO5 has no positivity safeguard there (#2502, as in #2469). With `face_bound = 1` or `2` (see [Positivity](validation_positivity.md)) the same run reaches t = 0.02. The check uses MUSCL with minmod. An octant with symmetry walls cannot be refined at the walls (a block must lie `buff_size` cells inside the domain), so the whole sphere is simulated.
- Save times come from `run_time.inf`, which has four significant digits (0.04 % in R).
- Only the radius and the peak density are checked; the interior profiles are not compared.

### Reproduce {#sec-validation-amr-sedov-run}

```
./mfc.sh run examples/3D_amr_sedov/case.py -t pre_process simulation -n 1 -- --N 96 --mode amr --subcycle
./mfc.sh run examples/3D_amr_sedov/case.py -t pre_process simulation -n 1 -- --N 96 --mode uniform --fine
python3 examples/3D_amr_sedov/check_sedov.py <run directories>
```

Measured on MI250X GCDs (amdflang, OpenMP offload, double): one GCD per run, 3-4 minutes for the 96³ AMR runs, 21 (subcycled) and 25 (lock-step) minutes for the 192³ AMR runs, 2.4 minutes for uniform 192³; the uniform 384³ run on 8 GCDs in 6.4 minutes. The full check took 32 minutes of one node (the 384³ run first, then the rest side by side). The suite runs the 32³ base, subcycled, to t = 0.005 as `AMR -> 3D -> Example -> Sedov regrid subcycle` (163DBF54), golden from the AMD CPU build (amdflang); the MI250X GPU build passes it at the default tolerance.

## Cylinder wake {#sec-validation-amr-wake}

**Status:** validated (Re 100; AMR against the uniform grid it replaces, and both against published values). **Example:** `examples/2D_amr_cylinder_wake`. Measured on `master` at 5838b641c; with #2482 the lock-step and subcycled runs on the 24 D x 16 D domain give the same Cd and St to four digits.

A static AMR block over a body and its near wake must give the drag and shedding frequency of a uniform grid at the block's resolution. Under lock-step AMR the immersed-boundary force sums each face on the finest level covering it, so Cd is a fine-grid quantity; the lift history gives St.

| Item | Value |
| :--- | :--- |
| Body | 2D cylinder, D = 1, no slip, fixed (`patch_ib` geometry 2) |
| Flow | ideal gas, \f$\gamma\f$ = 1.4, U = 0.2, c = 1 (M = 0.2), Re = 100; a transverse-velocity patch (0.2 U, 2 D x 1 D) behind the body seeds shedding |
| Domain | [-12 D, 24 D] x [-16 D, 16 D] (3.1 % blockage); Dirichlet inflow (`-17`), extrapolation elsewhere (`-3`) |
| Grids | uniform 16 and 32 cells/D; AMR: 16 cells/D base and one static 2:1 block over -1.5 D <= x <= 6 D, \|y\| <= 2.5 D (3.3 % of the 32 cells/D grid), lock-step and `amr_subcycle` |
| Numerics | WENO5-M with `mp_weno`, HLLC, SSP-RK3, `weno_Re_flux`, `fd_order = 2`, fixed \f$\Delta t\f$ at CFL 0.38 on the finest cells (level-0 cells when subcycled), to t = 650 (130 D/U) |
| Reference | Henderson \cite Henderson1995 (Cd 1.35), Park et al. \cite Park1998 (Cd 1.33, St 0.165), Williamson's fit \cite Williamson1989 (St 0.1643); all incompressible and unconfined |
| Measured | `D/ib_forces.dat` by `check_wake.py`: mean Cd and St = D/(U T), T the mean spacing of the upward zero crossings of Cl, over the whole periods in t = 400-650 (8 periods); an FFT of the window agrees to 0.2 % |

| Run | Cd | vs Henderson | vs Park | vs uniform 32 | St | vs Williamson | vs Park | vs uniform 32 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform 32 cells/D | 1.3687 | +1.38 % | +2.91 % | | 0.1674 | +1.86 % | +1.45 % | |
| uniform 16 cells/D | 1.3672 | +1.27 % | +2.79 % | −0.11 % | 0.1647 | +0.23 % | −0.17 % | −1.60 % |
| AMR, lock-step | 1.3682 | +1.35 % | +2.87 % | −0.03 % | 0.1673 | +1.83 % | +1.42 % | −0.03 % |
| AMR, subcycled | (1.6970) | | | (+24 %) | 0.1673 | +1.83 % | +1.42 % | −0.03 % |

- Lock-step AMR reproduces the uniform fine grid to 0.03 % in Cd and St, and its Cl amplitude to 0.7 % (0.292 against 0.294).
- The subcycled run sheds at the same frequency, but its force is summed on coarse faces only (as documented), and its Cd is 24 % high: not a drag measurement.
- On a 24 D x 16 D domain (6.25 % blockage, inflow 6 D upstream) every grid gave Cd and St 3-4 % higher (32 cells/D: 1.4187, 0.1730); AMR still matched the uniform fine grid (−0.31 %, −0.17 %). The domain, not the resolution, sets most of the remaining error.
- The remaining +1.4 % in Cd and +1.9 % in St against the incompressible references is of the size expected from the 3.1 % blockage and M = 0.2 (Prandtl-Glauert factor 1.02 on the pressure loads); neither was separated.
- Cost on one MI250X GCD: @ref amr "When AMR pays".

Limits:
- One Reynolds number and one block; the block's downstream extent (6 D) and the 16 cells/D base were not varied.
- The mean Cd rises by about 0.6 % through the window (uniform 32: 1.365 at t = 400-450, 1.373 at 600-650) while the Cl amplitude is steady; probably the far-field pressure creep of the extrapolated outflow seen in @ref validation_cyl_drag_nu, not verified.
- Cl amplitude 0.29 against 0.33 (Park), as in the uniform-grid cylinder check.
- Not in the test suite: a run is 130 D/U long (the uniform fine run took 42 min on one GCD). `2D_amr_cylinder_wake` is on the example skip list; the AMR immersed-boundary path is covered by `amr_golden_tests`.

### Reproduce {#sec-validation-amr-wake-run}

```
./mfc.sh run examples/2D_amr_cylinder_wake/case.py -t pre_process simulation -n 1 --gpu
python3 examples/2D_amr_cylinder_wake/check_wake.py examples/2D_amr_cylinder_wake [<uniform 32 cells/D run>]
```

The uniform and subcycled runs are edits at the top of `case.py` in a copy of the directory (`amr = False` with `N = 32` or `N = 16`; `subcycle = True`). Measured on one MI250X GCD, amdflang with OpenMP offload, double precision.
