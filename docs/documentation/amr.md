@page amr Adaptive Mesh Refinement

# Adaptive Mesh Refinement

> **Experimental.** AMR is off by default (`amr = F`). The case validator rejects unsupported combinations by name.

## Overview {#amr-overview}

AMR adds refined rectangular blocks on top of the base grid (level 0) where the flow needs resolution: shocks,
interfaces, wakes, bodies. The base grid is advanced by the ordinary solver; refined blocks run the same solver on
their own grids, nested up to `amr_max_level` levels at `amr_ref_ratio`:1 per level.

With `amr_regrid_int > 0`, the blocks are rebuilt around the cells where the density gradient exceeds `amr_tag_eps`:
once per level before the first step of a fresh run, so they start on the initial features, and then every
`amr_regrid_int` steps. With `amr_regrid_int = 0` a single block stays at `amr_block_beg`..`amr_block_end`.

On a fresh run the blocks start from `patch_icpp` evaluated at their own cell centers, as on a uniform grid at
that resolution. If the patches alone do not reproduce the level-0 data (hard-coded or STL patches, a perturbed or
smoothed initial condition, old data), the simulation prints a note and the blocks start from data interpolated from
level 0. The step-0 files are pre_process's level-0 data: under a block the run starts from the average of the fine
data instead, which differs from them where a patch varies within a coarse cell. For a smooth analytic profile
sampled at cell centers that difference is second order in the cell width; an expression that is the exact cell
average (it may use `dx`) agrees on every level.

AMR lives in `simulation` (`src/simulation/m_amr*.fpp`); post_process can overlay the blocks. Internals:
@subpage amr_implementation.

## Algorithm {#amr-algorithm}

- **Tagging and clustering.** A coarse cell is tagged when its normalized density jump exceeds `amr_tag_eps`
  (`amr_tag_abs = T` uses the absolute jump instead). `amr_tag_var` tags on the volume fraction of fluid 1, the
  velocity magnitude, the pressure, the temperature or a species mass fraction (`amr_tag_species`) instead of the
  density, for interfaces, shear layers, pressure waves and flames with no density jump. The vorticity
  (`amr_tag_var = 5`) is not a jump: a cell is tagged when its vorticity magnitude times its cell size
  `dx` (central velocity differences; a velocity, halving per level for a resolved vortex) exceeds `amr_tag_eps`,
  whatever `amr_tag_abs`. Tags are clustered
  by Berger-Rigoutsos bisection into boxes that reach `amr_cluster_eff` tag efficiency, are at least `amr_blocking_factor` cells wide, are padded by `amr_buf`
  cells and are capped at `amr_max_grid_size` cells per dimension. Levels nest properly (a level-2 box lies inside one
  level-1 box).
- **Prolongation.** New blocks and fine ghost cells are filled from the parent level by conservative minmod-limited
  linear interpolation (`amr_prolong = 1`) or conservative quartic interpolation (`amr_prolong = 2`). The quartic fill
  falls back to limited linear in a cell whose children would leave the range of the neighbouring coarse averages (widened
  only where the second differences share a sign: at a smooth local extremum and next to one, so a smooth extremum keeps
  the quartic, a steep monotone front stays monotone, and an overshoot beside a front is not amplified), and wherever it would give
  non-positive density or internal energy. Volume fractions are
  renormalized to sum to one and species mass fractions to `sum(Y) = 1, Y >= 0`. A block that survives a regrid
  keeps its own fine solution.
- **Advance.** Lock-step (default): every level takes each step at `dt`, so `dt` must satisfy the finest level's CFL
  limit. With a CFL-based time step (`cfl_adap_dt`, `cfl_const_dt`) each coarse cell's limit is that of the finest cell
  covering it, so `dt` drops where blocks appear; the reported ICFL is the finest level's too. With a fixed `dt` the run
  aborts at start-up or after a regrid if `dt` exceeds the finest cells' inviscid limit (ICFL = 1), and prints the
  largest stable `dt`. With regrid on, start-up also aborts if `dt` would exceed that limit were the fastest cells
  refined to `amr_max_level`, so a run cannot fail on it later. Subcycled (`amr_subcycle = T`): each refined level takes two steps per parent step, with ghost
  cells interpolated in time; `dt` is the level-0 step. Needs `time_stepper = 2` or `3`.
- **Refluxing and restriction.** Fine face fluxes (advective, viscous, diffusive) are accumulated in flux registers and
  replace the coarse fluxes at coarse/fine faces, so mass, momentum and energy are conserved across the interface. Each
  covered coarse cell is then set to the average of its fine children.

## Conservation and accuracy {#amr-conservation}

Conservation across coarse/fine faces holds to round-off for single- and multi-fluid, viscous and chemistry cases, and
a uniform state stays uniform with regrid on. Interiors keep the WENO order; the coarse/fine ghost fill is second
order (linear) or fifth order (quartic). For a smooth wave crossing a static block, the measured global L2 rate is 2.5 with
the linear fill (its second-order ghost values leave a second-order error in the few cells at each interface) and 5 with
the quartic fill, as for WENO5 on a uniform grid; the `Convergence -> 1D -> AMR` tests track it. Results do not depend on
message timing, and agree across rank counts to round-off when `amr_max_grid_size` is set (with `amr_max_grid_size = 0` the
block size cap, and therefore the box set, depends on the decomposition).

## Parallelism {#amr-parallel}

Each block is owned whole by one rank. Owners come from a cost-weighted space-filling-curve cut over the blocks (cost:
one per coarse cell plus `K_ib` per immersed-boundary cell), recomputed at every regrid. `amr_place_owner = T` instead
places each block on the rank that owns most of the coarse cells under it, capped at 5% above the mean load, which
cuts communication. A block may straddle rank boundaries; its coarse-side coupling moves point-to-point between the
owner and the ranks it overlaps. All exchanges are posted in a fixed order, so runs are deterministic.

On GPUs every block array is device-resident. A rank advances its blocks of equal level in batches stacked along the
last dimension, one solver call per batch; this is why AMR needs a uniform Cartesian grid. The toolchain enables
`amr_device_pack` (fused exchange kernels) when `amr_max_grid_size <= 64` and sets `amr_snap = min(2, amr_buf - 2)`
(regrid hysteresis) under dynamic regrid, unless they are set explicitly.

`rank_time_wrt = T` prints a per-phase time budget of the AMR step, the block-balance report, and at the end of the run
the memory each group of AMR arrays holds per rank (`[amr-mem]`: block store, advance scratch, coarse patches, flux
registers, the fine IB ghost-point park, and the host-side per-slot IB marker copies).

## Supported physics {#amr-physics}

Required: `model_eqns = 2` or `3`, `time_stepper = 1`, `2` or `3`, WENO or MUSCL (`recon_type = 1` or `2`) or IGR, a uniform Cartesian grid, and
`mpp_lim = T` when `num_fluids > 1`.

| Feature | Status |
| :--- | :--- |
| Single- and multi-fluid Euler, 6-equation model | Supported |
| Viscous, hypoelastic (with continuum damage), body forces | Supported |
| Chemistry (reactions, diffusion) | Lock-step only; no reaction substepping |
| Immersed boundaries (fixed, prescribed or force-driven motion; not STL, `thermal_bc`, `surface_reaction`, `ib_surface_wrt`, `ib_heat_wrt`, `ib_mass_relax` or `ib_uncover_fill`) | Supported; blocks always contain each body; moving bodies are lock-step only; multi-level with IB is single-rank; forces and torques are summed on the finest level covering each face (the coarse grid only with `amr_subcycle`); the default `ib_mass_relax = -1` resolves to off, so a moving body differs from a uniform-grid run, where it is on |
| Lagrangian bubbles | Lock-step; the cloud stays on the base grid (blocks are kept clear of it) |
| Acoustic sources | On the base grid; blocks are kept clear of the source support |
| Body forces (`bf_x/y/z`, `bf_spatial_support`) | Supported; the source is evaluated on each block's own grid |
| MHD / RMHD | 1D only |
| IGR | Single-rank, lock-step; coarse/fine conservation is truncation-order |
| Synthetic turbulence forcing | Supported; evaluated on each block's own grid |
| Phase change (`relax`) | Supported; each block is relaxed at the end of the step, then restricted |
| Euler-Euler bubbles, QBMM, surface tension | Not supported |
| Riemann-extrapolation boundaries (`bc = -4`) | Supported, except with `null_weights` |
| Cylindrical coordinates, stretched grids | Not supported |

## Parameters {#amr-parameters}

Full descriptions: @ref case.

| Parameter | Default | Description |
| :--- | :---: | :--- |
| `amr` | F | Enable AMR |
| `amr_block_beg(i)`, `amr_block_end(i)` | 0 | Initial block, level-0 cell indices (inclusive) |
| `amr_regrid_int` | 0 | Steps between regrids; 0 keeps the initial block |
| `amr_tag_eps` | 0.1 | Tagging threshold on the normalized density jump (> 0 with regrid) |
| `amr_tag_abs` | F | Tag on the absolute jump; `amr_tag_eps` then has the units of the tagged variable |
| `amr_tag_var` | 1 | Tagged variable: 1 total density, 2 volume fraction of fluid 1, 3 velocity magnitude (needs `amr_tag_abs`), 4 pressure (five-equation model), 5 vorticity (magnitude times the cell size against `amr_tag_eps`; 2D/3D), 6 temperature (chemistry), 7 mass fraction of species `amr_tag_species` (chemistry, needs `amr_tag_abs`) |
| `amr_tag_species` | 1 | Species index tagged by `amr_tag_var = 7` |
| `amr_buf` | 3 | Padding around tagged cells, coarse cells (>= 1 with regrid) |
| `amr_cluster_eff` | 0.7 | Minimum tag efficiency of a clustered box |
| `amr_blocking_factor` | 4 | Minimum box width, coarse cells |
| `amr_max_grid_size` | 0 | Maximum box width, coarse cells; 0 derives it from the decomposition |
| `amr_max_blocks` | 1024 | Maximum number of blocks |
| `amr_max_level` | 1 | Number of refined levels |
| `amr_ref_ratio` | 2 | Refinement ratio, 2 or 4 (4: one level only) |
| `amr_subcycle` | F | Subcycle refined levels in time (ratio 2, up to two levels) |
| `amr_prolong` | 1 | 1: limited linear, 2: quartic (ratio 2) |
| `amr_place_owner` | F | Owner-aware block placement |
| `amr_snap` | 0 | Regrid hysteresis, coarse cells (<= `amr_buf - 2`) |
| `amr_device_pack` | F | Fused device pack/unpack of the coarse-patch exchange |

## Example {#amr-example}

```python
# 2D shock tracking: one refined level, regrid every 5 steps
"amr": "T",
"amr_block_beg(1)": 48, "amr_block_end(1)": 79,
"amr_block_beg(2)": 48, "amr_block_end(2)": 79,
"amr_regrid_int": 5,
"amr_tag_eps": 0.05,
"amr_buf": 3,
"amr_max_grid_size": 32,
```

Complete cases: `examples/2D_amr_droplet` (a tracked interface), `examples/2D_amr_blast_multilevel` (two refined
levels), `examples/2D_amr_ibm_cylinder` (an immersed body) and `examples/3D_amr_blast` (3D, subcycled).

## When AMR pays {#amr-payoff}

AMR costs more per cell than a uniform grid: extra kernels for ghost fills, flux registers and restriction, plus the
regrids. On a GPU it pays only when the grid it replaces is large enough that the cell count, and not that fixed cost,
sets the time per step.

Measured on one MI250X GCD (amdflang, OpenMP offload): `examples/2D_amr_droplet` with one refined level, against the
same case on a uniform grid at the fine resolution. Time is the solver's time to the same physical time; below 1 the
AMR run is faster.

| Base grid | Uniform fine grid | Fine cells held by AMR | Lock-step / uniform fine | Subcycled / uniform fine |
| :--- | :--- | ---: | ---: | ---: |
| 128^2 | 256^2 (0.07 M cells) | 11% | 2.98 | 2.19 |
| 512^2 | 1024^2 (1.0 M cells) | 7% | 2.05 | 1.54 |
| 2048^2 | 4096^2 (16.8 M cells) | 6% | 1.05 | 0.70 |

The last row is from the time per step over 400 steps. A 3D spherical blast (`examples/3D_amr_blast`, 64^3 base
against a uniform 128^3, the blocks holding 54% of the fine cells) took 1.86 (lock-step) and 1.40 (subcycled) times the uniform run.

- **Size.** Below a few million cells on the grid AMR would replace, run the uniform grid: a GPU advances those
  cells at nearly no extra cost.
- **Coverage.** The blocks should hold a small part of the fine grid. The droplet's 6% breaks even at 16.8 M cells; the
  blast's 54% would need a far larger grid.
- **Subcycling.** `amr_subcycle = T` took a quarter to a third less time than lock-step in every case above.
- **Accuracy.** In the droplet runs the AMR solution reproduces the uniform fine one: the L1 density difference from
  it is 6.1e-8 at the 512^2 base, against 5.0e-4 for the uniform 512^2 grid. The blocks cover the interface from the
  first step and start from the patches evaluated at the fine resolution.

## Limitations {#amr-limitations}

- **Block size.** A block spans at most about half a rank's subdomain per dimension unless `amr_max_grid_size` is set;
  larger features are covered by several blocks.
- **Restart.** `parallel_io` restarts work on any rank count when `amr_max_grid_size` is set; otherwise restart on no
  more ranks than wrote the checkpoint. The serial format needs the rank count that wrote it.
  The file header records `amr_ref_ratio`, `amr_max_level` and the storage precision, which must match (the level count may
  only grow); a missing fine-level file aborts the restart.
- **Nesting.** A level-2 or deeper block lies inside one block of the level above; regrid builds each block's children
  inside it, a few coarse cells from its edge. Where a feature crosses the seam between two parent blocks, a strip that
  wide on each side of the seam stays at the parent's resolution. A larger `amr_max_grid_size` gives fewer seams.
- **Regrid.** When no cell is tagged, regrid replaces the hierarchy with one small dormant level-1 block (one coarse cell
  at the domain center, padded by `amr_buf`), since the solver needs at least one block. At every level, the part of a
  new block that an old block of the same level covered keeps the old data; the rest is prolonged from the parent.
- **Output.** Field output is at level-0 resolution with the fine solution restricted into it; post_process with
  `amr = T` overlays the blocks as separate SILO domains. Each save also writes `D/amr_levels.<step>.dat`: one line per
  refined level with its fine cell count and the mean of `|q|` over its cells for each conserved variable.
