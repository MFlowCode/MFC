@page validation_blast Validation: Circular Blast Wave

# Circular blast wave against a radial reference, uniform and adaptive grids

**Status:** validated (density and shock radius, uniform grids and two-level AMR on the CPU); the AMR run fails on an NVIDIA GPU
([#2467](https://github.com/MFlowCode/MFC/issues/2467)). **Code:** `master` at 5838b641c.
**Example:** `examples/2D_amr_blast_multilevel` (default: AMR; `--no-amr -N`: uniform grid). **Last run:** 2026-11-25.

## Purpose

A cylindrical explosion has a shock that weakens as it expands, a contact and an inward rarefaction, and its exact solution is
radially symmetric, so a fine 1D radial solution is an independent reference for the 2D Cartesian run. It checks the
multidimensional shock capturing on a grid not aligned with the flow, and the adaptive-mesh refinement of the example against
uniform grids.

## Setup

| Item | Value |
| :--- | :--- |
| Problem | disc of radius 0.1 at (0.5, 0.5) with (\f$\rho\f$, p) = (1, 1) in gas at (0.125, 0.1), at rest, \f$\gamma\f$ = 1.4 (the 2D Sod explosion, Toro section 17.1 \cite Toro09); t = 0.15625 |
| Domain | [0, 1]², extrapolation boundaries (the shock stays inside, at r = 0.34) |
| Numerics | WENO5-M, HLLC, SSP-RK3. AMR (default): 128² level 0, two 2:1 levels (`amr_max_level = 2`) regridded every 8 steps, lock-step \f$\Delta t\f$ = 0.1 \f$\Delta x_0\f$/4, 800 steps. Uniform (`--no-amr`): N², \f$\Delta t\f$ = 0.1 \f$\Delta x\f$ |
| Reference | radially symmetric Euler equations, \f$\partial_t(rU) + \partial_r(rF) = (0, p, 0)\f$, second-order finite volume (MUSCL-minmod, HLL, SSP-RK2) on 16000 cells over 0 \f$\le\f$ r \f$\le\f$ 0.75, in `check_blast.py`. Its L1 change from 4000 to 8000 cells is 1.1e-4 (a fifth of the smallest 2D error below) |
| Measured | L1 density error \f$\sum|\rho - \rho_\mathrm{ref}|\,\Delta x\,\Delta y\f$ against the reference averaged over each cell (4 × 4 samples), and the shock radius (outermost radius where the azimuthal-mean density exceeds the mean of the ambient and the peak density), measured the same way on both, by `check_blast.py` |

AMR output is the level-0 field with the finer levels restricted into it, so the AMR run is compared at 128².

## Results

| Grid | L1 \f$\rho\f$ | order | Shock radius | Reference |
| :--- | ---: | ---: | ---: | ---: |
| uniform 128² | 1.98e-3 | | 0.3400 | 0.3400 |
| uniform 256² | 1.08e-3 | 0.87 | 0.3419 | 0.3421 |
| uniform 512² | 5.90e-4 | 0.87 | 0.3431 | 0.3432 |
| AMR, 128² + 2 levels (512² effective at the fronts) | 5.47e-4 | | 0.3402 | 0.3400 |

On an A100 (NVHPC 25.11) the uniform 512² run agrees with the CPU to 7e-12 in density (the same L1 error and shock radius). The
AMR run aborts on the A100 at step 665 of 800 when its flux registers grow ([#2467](https://github.com/MFlowCode/MFC/issues/2467)),
so the AMR row is CPU only.

## Assessment

- The density converges at order 0.87 toward the radial reference, as expected with a shock and a contact; the shock radius agrees
  within a tenth of a cell at every resolution.
- The two-level AMR run, refined to 512² resolution around the fronts, is as accurate as the uniform 512² grid (5.5e-4 against
  5.9e-4) on its level-0 output.

## Known limitations

- The disc is a staircase on the Cartesian grid (`geometry = 2` assigns whole cells), which adds an O(\f$\Delta x\f$) error to the initial
  mass; it is part of the error above.
- Only the cylindrical explosion is checked. A Sedov point blast (self-similar, \f$r_s \propto t^{1/2}\f$ in 2D) has no example; the 3D
  spherical blast `examples/3D_amr_blast` is not checked here.
- The AMR comparison uses the restricted level-0 field; the fine-level solution is not compared directly.
- AMR on an NVIDIA GPU with NVHPC aborts in this example ([#2467](https://github.com/MFlowCode/MFC/issues/2467)).

## Reproduce

```
./mfc.sh run examples/2D_amr_blast_multilevel/case.py -t pre_process simulation -n 1
./mfc.sh run examples/2D_amr_blast_multilevel/case.py -t pre_process simulation -n 1 -- --no-amr -N 256
python3 examples/2D_amr_blast_multilevel/check_blast.py <AMR run> <uniform run> [...]
```

The 512² uniform run takes about 10 minutes on one CPU core; the reference about 1 minute. The example is not in the example suite (its AMR blocks do not fit
the capped test grid); AMR blast waves are covered by the `AMR -> 2D` golden tests. No test is added here.
