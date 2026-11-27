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
  level-1 box); boxes that end on either side of a seam between two parent blocks with the same extents along it are
  joined across the seam.
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
  replace the coarse fluxes at coarse/fine faces, so mass, momentum and energy are conserved across the interface. A
  volume fraction's flux pairs with an \f$\alpha\,\nabla\cdot u\f$ source on the same faces, so the registers also hold
  the face velocity and the volume fractions' reflux corrects both: uniform \f$\alpha\f$ stays uniform, and a material
  interface on a block face keeps pressure and velocity equilibrium (HLLC, HLL with `hll_u_interface` and LF, with two or
  more fluids). The same correction covers the `alt_soundspeed` \f$K\,\nabla\cdot u\f$ term, the 6-equation internal
  energies' \f$\alpha_k p\,\nabla\cdot u\f$, and HLL Method 1, whose per-fluid interface volume fractions the registers
  also hold. Each covered coarse cell is then set to the average of its fine children.

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
`amr_device_pack` (fused exchange kernels) and sets `amr_snap = min(2, amr_buf - 2)`
(regrid hysteresis) under dynamic regrid with one refined level, unless they are set explicitly.

`rank_time_wrt = T` prints a per-phase time budget of the AMR step, the block-balance report, and at the end of the run
the memory each group of AMR arrays holds per rank (`[amr-mem]`: block store, advance scratch, coarse patches, flux
registers, the fine IB ghost-point park, and the host-side per-slot IB marker copies). Every block's store column is
sized to the block size cap (`amr_max_grid_size`, or the rank-derived cap) plus ghost layers, so it also prints that column's
size and, per level, the most blocks any rank held over the run, their store, and how much of its column (ghost layers
included, so below 1 even at the cap) a block's interior fills (min, mean, max). Each level's row also splits its store into the ghost-layer part and gives its flux registers (one register slot per block), and GPU builds add one line with the device and host totals: every group but the IB marker copies is device-resident and keeps a host mirror.

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
| Condensed-phase burn (`reactive_burn`) | Supported; with `rburn%%substeps > 0` each block takes the split burn at the end of the step, before the fold |
| Euler-Euler bubbles, QBMM, surface tension | Not supported |
| Riemann-extrapolation boundaries (`bc = -4`) | Supported, except with `null_weights` |
| Characteristic boundaries (`bc = -5` to `-12`) | Supported; the coarse grid applies them at the domain boundary (blocks stay `buff_size` cells inside it) |
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
| `amr_max_grid_size` | 0 | Maximum box width, coarse cells; 0 derives it from the decomposition (and, on several ranks, tiles regrid boxes smaller, not below 16, until there is one per rank) |
| `amr_max_blocks` | 1024 | Maximum number of blocks |
| `amr_max_level` | 1 | Number of refined levels |
| `amr_ref_ratio` | 2 | Refinement ratio, 2 or 4 (4: one level only) |
| `amr_subcycle` | F | Subcycle refined levels in time (ratio 2, up to two levels) |
| `amr_prolong` | 1 | 1: limited linear, 2: quartic (ratio 2) |
| `amr_place_owner` | F | Owner-aware block placement |
| `amr_snap` | 0 | Regrid hysteresis, coarse cells (<= `amr_buf - 2`; one refined level only under regrid) |
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
levels), `examples/2D_amr_ibm_cylinder` (an immersed body), `examples/2D_amr_cylinder_wake` (a viscous wake, drag and St), `examples/3D_amr_blast` (3D, subcycled) and
`examples/1D_amr_reactive_shocktube` (chemistry: the blocks follow a reacting front tagged on temperature).

## When AMR pays {#amr-payoff}

AMR costs more per cell than a uniform grid: extra kernels for ghost fills, flux registers and restriction, plus the
regrids. On a GPU it pays only when the grid it replaces is large enough that the cell count, and not that fixed cost,
sets the time per step.

Measured on one MI250X GCD (amdflang, OpenMP offload): `examples/2D_amr_droplet` with one refined level, against the
same case on a uniform grid at the fine resolution. Time is the wall time of the step loop to the same physical time;
below 1 the AMR run is faster.

| Base grid | Uniform fine grid | Fine cells held by AMR | Lock-step / uniform fine | Subcycled / uniform fine |
| :--- | :--- | ---: | ---: | ---: |
| 128^2 | 256^2 (0.07 M cells) | 11% | 2.49 | 1.77 |
| 512^2 | 1024^2 (1.0 M cells) | 7% | 1.49 | 1.04 |
| 2048^2 | 4096^2 (16.8 M cells) | 6% | 0.51 | 0.32 |

The last row runs 400 fine steps. A 3D spherical blast (`examples/3D_amr_blast`, 64^3 base against a uniform 128^3 to
t = 0.1875, the blocks holding 54% of the fine cells) took 1.17 (lock-step) and 0.87 (subcycled) times the uniform run.

On CPUs the cost per cell is higher and the fixed cost matters less, so AMR pays at every size measured. Same droplet
cases on 8 MPI ranks of one Frontier node (amdflang, CPU), relative to the uniform fine grid:

| Base grid | Uniform fine grid | Lock-step / uniform fine | Subcycled / uniform fine |
| :--- | :--- | ---: | ---: |
| 128^2 | 256^2 | 0.54 | 0.37 |
| 512^2 | 1024^2 | 0.41 | 0.27 |
| 2048^2 | 4096^2 | 0.36 | 0.22 |

A 2D viscous cylinder wake (`examples/2D_amr_cylinder_wake`: Re 100, Mach 0.2, 16 cells per diameter, a static block
over the body and near wake; the variants are edits at the top of `case.py`), run to 130 D/U on one MI250X GCD (whole-run wall time):

| Domain | Uniform fine grid | Fine cells held by AMR | Lock-step / uniform fine | Subcycled / uniform fine | Uniform coarse / uniform fine |
| :--- | :--- | ---: | ---: | ---: | ---: |
| 24 D x 16 D | 768 x 512 (0.39 M cells) | 9.8% | 1.94 | 1.37 | 0.38 |
| 36 D x 32 D | 1152 x 1024 (1.18 M cells) | 3.3% | 1.40 | 0.98 | 0.29 |

On the larger domain the lock-step run matches the uniform fine grid's mean drag coefficient and Strouhal number to
0.03%; the uniform coarse grid is off by -0.1% and -1.6%. The subcycled run's St matches as well; its drag is not
comparable (forces on coarse faces only).

- **Size.** On a GPU, AMR breaks even near a million cells on the grid it replaces (the subcycled droplet at 1.0 M
  cells, 1.04) and pays well above that (0.32 at 16.8 M). Below that, run the uniform grid: a GPU advances those cells
  at little extra cost. On CPUs AMR paid off from the smallest case above.
- **Coverage.** The smaller the part of the fine grid the blocks hold, the sooner AMR pays: the wake's blocks hold
  3.3% of a 1.2 M-cell grid and its subcycled run breaks even (0.98); the blast's hold 54% of a 2.1 M-cell grid and
  only its subcycled run is faster than the uniform one.
- **Subcycling.** `amr_subcycle = T` took 26% to 39% less time than lock-step in every case above.
- **Accuracy.** In the droplet runs the AMR solution reproduces the uniform fine one: the L1 density difference from
  it is 6.1e-8 at the 512^2 base, against 5.0e-4 for the uniform 512^2 grid. The blocks cover the interface from the
  first step and start from the patches evaluated at the fine resolution. For an integral quantity check the coarse grid
  first: in the wake it already has the drag to 0.1%; refinement buys the 1.6% in St.

## Limitations {#amr-limitations}

- **Block size.** A block spans at most about half a rank's subdomain per dimension unless `amr_max_grid_size` is set;
  larger features are covered by several blocks. On several ranks, regrid tiles the level-1 boxes smaller (not below 16
  coarse cells) until there is at least one per rank (or per block-pool slot, if `amr_max_blocks` is smaller), so a small
  refined region still spreads over the ranks; a size that would leave a tile narrower than the ghost layer or overflow the
  block pool is skipped.
- **Restart.** `parallel_io` restarts work on any rank count when `amr_max_grid_size` is set; otherwise restart on no
  more ranks than wrote the checkpoint. The serial format needs the rank count that wrote it.
  The file header records `amr_ref_ratio`, `amr_max_level` and the storage precision, which must match (the level count may
  only grow); a missing fine-level file aborts the restart.
- **Nesting.** A level-2 or deeper block lies inside one block of the level above; regrid builds each block's children
  inside it, a few coarse cells from its edge. Where a feature crosses the seam between two parent blocks, the children
  on the two sides are joined across it when their extents along it are equal (the tag buffer usually makes them so,
  also for a disc crossing the seam off its centre):
  each is extended to the seam, or, where that would exceed its level's size cap, the gap is filled by a thin block of
  its own. Otherwise a strip that wide on each side of the seam stays at the parent's resolution; so does every seam with
  immersed boundaries, or with `amr_prolong = 2` on a scheme with fewer
  than 4 ghost cells (WENO3 and below). A larger `amr_max_grid_size` gives fewer seams.
- **Regrid.** When no cell is tagged, regrid replaces the hierarchy with one small dormant level-1 block (one coarse cell
  at the domain center, padded by `amr_buf`), since the solver needs at least one block. At every level, the part of a
  new block that an old block of the same level covered keeps the old data; the rest is prolonged from the parent.
- **Output.** Field output is at level-0 resolution with the fine solution restricted into it; post_process with
  `amr = T` overlays the blocks as separate SILO domains. Each save also writes `D/amr_levels.<step>.dat`: one line per
  refined level with its fine cell count and the mean of `|q|` over its cells for each conserved variable.
  Step 0 is pre_process's output: level 0 only, with no blocks and no `amr_levels` file, since simulation regrids at
  step 0 but first saves at step `t_step_save`. A run started at step 0 is therefore always a fresh start: it redoes
  the initial regrid and fills the blocks from the patches.
  The overlay carries density, pressure, velocity and volume fractions only; probes and derived post_process fields read
  level 0, `run_time_info` reads the level-0 state but scales its CFL numbers to the finest level covering each cell
  (its Mie-Gruneisen compression check also reads the fine blocks), and `sim_data` is rejected.
