@page amr Adaptive Mesh Refinement

# Adaptive Mesh Refinement

> **Experimental.** AMR is off by default (`amr = F`). The case validator rejects unsupported combinations by name.

## Overview {#amr-overview}

AMR adds refined rectangular blocks on top of the base grid (level 0) where the flow needs resolution: shocks,
interfaces, wakes, bodies. The base grid is advanced by the ordinary solver; refined blocks run the same solver on
their own grids, nested up to `amr_max_level` levels at `amr_ref_ratio`:1 per level.

With `amr_regrid_int > 0`, every `amr_regrid_int` steps the blocks are rebuilt around the cells where the density
gradient exceeds `amr_tag_eps`. With `amr_regrid_int = 0` a single block stays at `amr_block_beg`..`amr_block_end`.

AMR lives in `simulation` (`src/simulation/m_amr*.fpp`); post_process can overlay the blocks. Internals:
@subpage amr_implementation.

## Algorithm {#amr-algorithm}

- **Tagging and clustering.** A coarse cell is tagged when its normalized density jump exceeds `amr_tag_eps`
  (`amr_tag_abs = T` uses the absolute jump instead). Tags are clustered by Berger-Rigoutsos bisection into boxes
  that reach `amr_cluster_eff` tag efficiency, are at least `amr_blocking_factor` cells wide, are padded by `amr_buf`
  cells and are capped at `amr_max_grid_size` cells per dimension. Levels nest properly (a level-2 box lies inside one
  level-1 box).
- **Prolongation.** New blocks and fine ghost cells are filled from the parent level by conservative minmod-limited
  linear interpolation (`amr_prolong = 1`) or conservative quartic interpolation (`amr_prolong = 2`). The quartic fill
  falls back to limited linear wherever it would give non-positive density or internal energy. Volume fractions are
  renormalized to sum to one and species mass fractions to `sum(Y) = 1, Y >= 0`. A block that survives a regrid
  keeps its own fine solution.
- **Advance.** Lock-step (default): every level takes each step at `dt`, so `dt` must satisfy the finest level's CFL
  limit. Subcycled (`amr_subcycle = T`): each refined level takes two steps per parent step, with ghost cells
  interpolated in time; `dt` is the level-0 step.
- **Refluxing and restriction.** Fine face fluxes (advective, viscous, diffusive) are accumulated in flux registers and
  replace the coarse fluxes at coarse/fine faces, so mass, momentum and energy are conserved across the interface. Each
  covered coarse cell is then set to the average of its fine children.

## Conservation and accuracy {#amr-conservation}

Conservation across coarse/fine faces holds to round-off for single- and multi-fluid, viscous and chemistry cases, and
a uniform state stays uniform with regrid on. Interiors keep the WENO order; the coarse/fine ghost fill is second
order (linear) or higher (quartic). Results do not depend on message timing, and agree across rank counts to round-off
when `amr_max_grid_size` is set (with `amr_max_grid_size = 0` the block size cap, and therefore the box set, depends
on the decomposition).

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

`rank_time_wrt = T` prints a per-phase time budget of the AMR step and the block-balance report.

## Supported physics {#amr-physics}

Required: `model_eqns = 2` or `3`, `time_stepper = 3`, WENO (`recon_type = 1`) or IGR, a uniform Cartesian grid, and
`mpp_lim = T` when `num_fluids > 1`.

| Feature | Status |
| :--- | :--- |
| Single- and multi-fluid Euler, 6-equation model | Supported |
| Viscous, hypoelastic (with continuum damage), body forces | Supported |
| Chemistry (reactions, diffusion) | Lock-step only; no reaction substepping |
| Immersed boundaries (fixed, prescribed or force-driven motion; not STL) | Supported; blocks always contain each body; moving bodies are lock-step only; multi-level with IB is single-rank |
| Lagrangian bubbles | Lock-step; the cloud stays on the base grid (blocks are kept clear of it) |
| Acoustic sources | On the base grid; blocks are kept clear of the source support |
| MHD / RMHD | 1D only |
| IGR | Single-rank, lock-step; coarse/fine conservation is truncation-order |
| Euler-Euler bubbles, QBMM, phase change (`relax`), surface tension, synthetic turbulence | Not supported |
| Cylindrical coordinates, stretched grids, `bc = -4` | Not supported |

## Parameters {#amr-parameters}

Full descriptions: @ref case.

| Parameter | Default | Description |
| :--- | :---: | :--- |
| `amr` | F | Enable AMR |
| `amr_block_beg(i)`, `amr_block_end(i)` | 0 | Initial block, level-0 cell indices (inclusive) |
| `amr_regrid_int` | 0 | Steps between regrids; 0 keeps the initial block |
| `amr_tag_eps` | 0.1 | Tagging threshold on the normalized density jump (> 0 with regrid) |
| `amr_tag_abs` | F | Tag on the absolute density jump; `amr_tag_eps` is then a density |
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

See `examples/2D_amr_droplet` for a complete case.

## Limitations {#amr-limitations}

- **Block size.** A block spans at most about half a rank's subdomain per dimension unless `amr_max_grid_size` is set;
  larger features are covered by several blocks.
- **Restart.** `parallel_io` restarts work on any rank count when `amr_max_grid_size` is set; otherwise restart on no
  more ranks than wrote the checkpoint. The serial format needs the rank count that wrote it.
  The file header records `amr_ref_ratio`, `amr_max_level` and the storage precision, which must match (the level count may
  only grow); a missing fine-level file aborts the restart.
- **Nesting.** A level-2 block must lie inside one level-1 block; the run aborts otherwise (reduce
  `amr_max_grid_size`).
- **Regrid.** When no cell is tagged, regrid keeps the current blocks, so refinement does not retract until something is
  tagged again. Level-2 and deeper blocks are re-prolonged from their parent at every regrid (only level-1 data is
  carried over), so frequent regrids cost accuracy at those levels.
- **Output.** Field output is at level-0 resolution with the fine solution restricted into it; post_process with
  `amr = T` overlays the blocks as separate SILO domains.
