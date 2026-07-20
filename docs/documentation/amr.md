@page amr Adaptive Mesh Refinement

# Adaptive Mesh Refinement

> **Experimental.** AMR is off by default (`amr = F`) and requires explicit opt-in.
> The physics support matrix is enforced at run time by the checker; unsupported
> combinations abort with a clear message rather than producing silently wrong results.

## Overview {#amr-overview}

Block-structured adaptive mesh refinement (AMR) concentrates resolution where the flow
demands it — around shocks, interfaces, and bubble clouds — while leaving the rest of the
domain at the coarser base-grid resolution. MFC implements a multi-level block hierarchy: the
unmodified base (level-0) solve runs as usual, and one or more refined rectangular blocks advance
alongside it on a finer grid — nested recursively to `amr_max_level` levels when enabled.

The fine blocks are dynamically repositioned every `amr_regrid_int` coarse steps using a
gradient-based cell tagger and Berger–Rigoutsos block clustering, so they follow moving
features automatically. When `amr_regrid_int = 0` the block is fixed at the initial
`amr_block_beg`/`amr_block_end` position for the whole run.

AMR lives entirely in the `simulation` executable (`src/simulation/m_amr.fpp` and
`m_amr_registers.fpp`) and is the only part of MFC that modifies the solver's grid
mid-run.

---

## The Block-Structured Model {#amr-model}

The hierarchy spans levels `0` through `amr_max_level` (default `1`, i.e. two levels):

- **Level 0** — the base grid with cell spacing `dx`, `dy`, `dz`. The ordinary MFC solver
  advances this level every step.
- **Level 1** — a list of up to `amr_max_blocks` rectangular refined blocks, each covering
  a sub-region of the level-0 domain at `ref_ratio`:1 refinement (`ref_ratio` = 2 or 4;
  the cell spacing shrinks by `ref_ratio` in every active direction).
- **Levels 2 … `amr_max_level`** — when `amr_max_level > 1`, blocks nest recursively: a
  level-`l` block refines a region of its parent level-(`l-1`) block by a further
  `ref_ratio`, tracking a moving feature to arbitrary depth. Multi-level nesting requires
  `ref_ratio = 2`. See @ref amr_multilevel for the nesting and reflux details.

Each block is described by its bounding box in level-0 cell-index space
(`amr_block_beg(1:num_dims)` to `amr_block_end(1:num_dims)`). The initial (level-1)
fine-block extents must satisfy:

```
ref_ratio*(amr_block_end(i) - amr_block_beg(i) + 1) - 1 <= N_i
```

where `N_i` is the global cell count in direction `i`. This ensures the fine scratch
(which is sized to the base grid) is never overflowed. A level-`l` block's fine extent
grows as `ref_ratio**l`, so the nested boxes are sized accordingly.

**Fixed-slot storage.** All block slots are pre-allocated at init to the maximum possible
block size (half the per-rank subdomain in each dimension). Setting `amr_max_blocks = N`
requires roughly `N` times the device memory of one max-size block. Regrid reuses the
existing allocations by adjusting the active geometry metadata; no re-allocation occurs
at run time.

---

## Algorithm {#amr-algorithm}

### Prolongation (coarse-to-fine interpolation) {#amr-prolongation}

When a block is first created or moved by regrid, the fine solution is initialized from
the coarse level by **conservative-linear prolongation**: a piecewise-linear fit to the
coarse cell averages is evaluated at each fine cell centre, ensuring the fine-cell averages
are consistent with the coarse parent to second order. Physics-specific closures
are applied after prolongation:

- **Multi-fluid**: volume fractions are renormalized so they sum to one on the fine level.
- **Euler-Euler bubbles**: the radius moment is floored to a small positive fraction of
  the parent so reconstructed radius and number density stay positive.
- **Chemistry**: species partial densities are rescaled so `sum(Y_k) = 1` and `Y_k >= 0`
  on the fine level by construction.

### Per-block advance {#amr-advance}

The fine block advances using the same WENO/Riemann/RK3 solver as the coarse level. The
solver globals (`m`, `n`, `p`, `dx`, ...) are swapped to the block geometry before the
advance and restored afterward. Ghost cells at the coarse/fine boundary are filled by
conservative-linear prolongation from the current coarse solution (or time-interpolated
between the coarse `t^n` and `t^{n+1}` states when subcycling is enabled).

Without subcycling, the fine block takes one step at the case `dt` (same as the coarse
step). With subcycling it takes two steps at `dt/2`.

### Flux registers and refluxing {#amr-reflux}

A **flux register** accumulates the fine-level face fluxes at every coarse/fine interface
face during the fine advance. After the fine advance completes, the **reflux correction**
replaces the corresponding coarse-level fluxes with the summed fine fluxes (Berger–Colella
refluxing). This correction is what makes the scheme exactly conservative: the coarse and
fine levels exchange mass, momentum, and energy at machine precision across the block
boundary, regardless of the solution inside the block.

Viscous stress and work enter as face-centred source fluxes of the same form as the
advective flux, so they travel through the same registers. The reflux matches the full
advective + viscous flux, and energy (including viscous work) is conserved.

### Restriction (fine-to-coarse averaging) {#amr-restriction}

After refluxing, the fine cell averages inside the block are **volume-averaged** back to
the coarse level (each coarse cell equals the average of its 2^d fine children). This
overwrites the coarse solution inside the block with the fine-level values. The coarse
level outside the block is unchanged. On uniform Cartesian grids the children share a cell
volume, so this is a plain arithmetic mean; under `cyl_coord` (2D axisymmetric) cell volume
scales with radius, so the fold-back is weighted by each fine child's cell-center radius and
the flux-register reflux area-weights the coarse/fine faces (radial faces by the outside
cell's `r_face/r_cell`, axial faces by the covering fine faces' radii), keeping the
radius-weighted conserved quantity exact to machine precision.

### Dynamic regrid {#amr-regrid}

Every `amr_regrid_int` coarse steps (when `amr_regrid_int > 0`):

1. **Tag** every coarse cell whose normalized density gradient exceeds `amr_tag_eps`.
2. **Cluster** the tagged cells using Berger–Rigoutsos recursive bisection into a list of
   rectangular boxes, each grown by `amr_buf` coarse cells of buffer padding on each side.
3. **Merge** boxes whose padded extents come within a ghost-cell buffer width of each
   other (guaranteeing no fine–fine adjacency, which the solver does not support).
4. **Split** each box further if its tag efficiency (tagged cells / total cells) has not
   yet reached `amr_cluster_eff`. Stop splitting when the box count reaches
   `amr_max_blocks`.
5. **Prolongate** the new block geometry from the current coarse solution, then continue.

Setting `amr_buf >= 1` and `amr_tag_eps > 0` is required when regridding is active.

### Subcycling {#amr-subcycle}

`amr_subcycle = T` enables Berger–Colella dt/2 subcycling:

- The coarse level takes one step at the case `dt`.
- The fine level takes two steps at `dt/2`.
- Ghost values at the coarse/fine boundary are time-interpolated between the saved
  `t^n` and `t^{n+1}` coarse states for each fine substep.
- Accumulated fine fluxes are refluxed back to the coarse level after each coarse step.

Subcycling improves temporal accuracy at the block boundary and is the recommended mode
for production runs. It requires a fixed time step (`cfl_dt = F`).

---

## Conservation and Accuracy {#amr-conservation}

**Exact conservation.** The flux-register reflux mechanism ensures per-fluid mass,
momentum, and total energy are conserved to machine precision across the coarse/fine
boundary. Defects are consistently ~1e-15 in validated runs (single-fluid, multi-fluid,
viscous, bubble, chemistry, and phase-change cases).

**Free-stream preservation.** A uniform-state run with AMR active (including
subcycling and regrid) preserves the free stream to machine precision — no spurious
velocities or pressure drift at the block boundary.

**Element-exact multi-rank.** An `np=1` run and an `np=2` run with the block spanning
the rank boundary produce bit-identical results (element-exact seam), confirming that
the mirror decomposition and fine halo exchange introduce no asymmetry.

**Accuracy posture.** Inside the block the fine-level solution converges at WENO order.
At the coarse/fine boundary the conservative-linear ghost fill is second order. The
viscous seam carries a bounded (~1e-6) np-dependent error only at the *prolongation
ghost* layer (the inherently-approximate coupling zone of any c/f boundary); bulk
accuracy is unaffected. For viscous cases with strong shear or boundary layers, a static
or generously buffered block is recommended (the density-gradient tagger does not sense
shear well; error-estimator taggers are future work).

---

## Parallelism {#amr-parallel}

### Multi-rank (MPI) {#amr-mpi}

The fine level uses **mirror decomposition**: each MPI rank holds the fine cells
that cover the intersection of the block with its own base-level subdomain. A block
may span rank boundaries and move across them freely under dynamic regrid. Fine ghost
cells at rank seams are exchanged using a dedicated halo exchange
(`s_mpi_sendrecv_amr_fine_halo`), which mirrors the base-level halo path but operates
over the fine geometry. Flux registers are distributed in the same pattern — each rank
owns the register faces that border its fine subdomain — and the reflux correction is
applied rank-locally after an MPI reduction of the distributed fluxes.

The fine block may cover at most about half of any rank's subdomain per dimension, since
the fine advance reuses the rank-local solver scratch (sized to the base subdomain).

Restart with `parallel_io` repartitions the fine blocks across any rank count; the
serial (per-rank-file) restart path requires the same rank count as the run that wrote
the fine-level file and aborts otherwise (see the Restart section below).

### Load-balance coupling {#amr-loadbalance}

When `load_balance = T` (see the load-balance parameters in @ref case), the weighted
Cartesian decomposition accounts for the extra work of the fine advance: cells inside
the block footprint are given a higher weight, biasing the rank boundaries toward a
balanced allocation of fine work. A deterministic feasibility clamp ensures the weighted
split never violates the half-subdomain constraint.

### GPU {#amr-gpu}

The fine-block field arrays (`q_cons`, `q_prim`, `rhs`, and the ghost-lerp sources) are
device-resident from allocation onward. The ghost fill, per-block RHS/RK advance, and
restriction are all performed on-device. GPU correctness has been validated on NVIDIA
V100 (OpenACC / nvhpc) for all supported physics rungs: single-fluid, multi-fluid,
viscous, bubbles, multi-block, phase-change, and chemistry.

The GPU build uses `src/simulation/` OpenACC/OpenMP macros; see @ref gpuParallelization
for the macro API.

---

## Supported Physics {#amr-physics}

The table below summarises what is and is not supported under AMR. The checker
(`src/simulation/m_checker.fpp`) enforces every restriction at run time and aborts with
a diagnostic message for unsupported combinations.

| Physics | Status | Notes |
| :--- | :---: | :--- |
| Single-fluid Euler (`num_fluids = 1`) | Supported | Base configuration |
| Multi-fluid Euler (`num_fluids > 1`) | Supported | Requires `mpp_lim = T`; volume fractions sum-preserved on prolongation |
| Viscous (`viscous = T`) | Supported | Viscous fluxes refluxed; bounded seam error at prolongation ghost layer |
| Euler-Euler bubbles (`bubbles_euler`; polytropic or non-polytropic; `nb >= 1`, polydisperse) | Supported | Moment realizability floor applied to all positive moments on prolongation |
| Phase change / relaxation (`relax = T`) | Supported | Per-cell relaxation runs on the fine block before restriction |
| Chemistry — reactions + advection + diffusion (`chemistry = T`) | Supported | Species sum/positivity closure on prolongation; temperature ghost exchanged at rank seams; diffusion fluxes refluxed like viscous |
| Surface tension (`surface_tension = T`) | **Not supported** | The capillary force depends on the interface-normal direction; the prolonged fine ghost color cannot reproduce the coarse normal across a 2:1 boundary, producing a growing spurious seam current. See @ref case section 7.1. |
| Hypoelasticity (`hypoelasticity = T`, incl. continuum damage) | Supported | Stress components prolong on the generic conservative path; the fine swap recomputes the spacing-dependent FD coefficients |
| Hyperelasticity | **Not supported** | Gated (no upstream test coverage to validate against) |
| MHD / RMHD, 1D | Supported | div(B) = d(Bx)/dx and 1D evolves only By/Bz (Bx is the uniform `Bx0` parameter), so div(B) = 0 by construction - the 2D/3D seam failure mode is structurally absent; By/Bz reflux and restrict as ordinary conserved scalars (HLL and HLLD; incl. relativistic) |
| MHD, 2D/3D | **Not supported** | Attempted and measured: per-component B prolongation/reflux is not divergence-preserving - on a magnetized 2D Brio-Wu the coarse/fine seam is a continuous O(1) monopole source that GLM cleaning spreads but cannot remove (max abs div(B) 0.53 block-interior / 0.36 far-field vs the no-AMR 1.4e-3 cleaning background; HLLD, which has no GLM coupling, NaNs outright). Needs constrained-transport-class B prolongation and reflux |
| QBMM bubbles (polytropic) | Supported | Bubble moments live in `q_cons`, injected piecewise-constant at prolongation to preserve CHyQMOM realizability |
| QBMM bubbles (non-polytropic) | Supported | Each block carries its own `pb`/`mv` quadrature side-state: prolonged piecewise-constant (realizability), advanced with the block's own rhs scratch, restricted back with the moments; dynamic regrid and `amr_subcycle` are both supported (the side-state bounces through the regrid and time-lerps its ghost shell) |
| Lagrangian bubbles (`bubbles_lagrange = T`) | Supported (cloud excluded from blocks) | Two-way coupling lives on the coarse grid: regrid suppresses tags and clips candidate boxes around the cloud's padded bbox (positions + `mapCells` smearing + stencil + drift margin, recomputed collectively each regrid), the fine advance skips the EL hooks, EL volume fractions prolong WITHOUT the sum-to-one closure (their sum is the local liquid fraction), and a per-stage guard aborts if the cloud reaches an active block |
| Immersed boundaries (`ib = T`; one or more non-STL bodies, static or prescribed-motion `moving_ibm=1`) | Supported | Per-block fine-grid IB markers/ghost points, rebuilt each fine substage at the body's sub-time position for a moving body; non-conservative ghost-cell forcing at the body; with dynamic regrid candidate boxes expand to fully contain each body at its live position plus margin, the fine IB state rebuilds after every regrid, and a per-substage guard aborts if a moving body reaches its block boundary between regrids; force-driven (`moving_ibm=2`)/STL gated; a body spanning a rank seam is rejected at startup |
| IGR solver (`igr = T`) | Supported (restriction-only coupling) | The fine block runs its own fixed-iteration sigma solve, seeded and Dirichlet-bounded by the converged coarse sigma (frozen ghost ring; the per-iteration BC populate is skipped); the coarse warm-start state is saved/restored across the fine advance. The Berger-Colella reflux is NOT captured from the fused IGR flux kernels, so seam conservation is truncation-order rather than exact; free-stream preservation is exact. `amr_subcycle` is gated |
| 2D axisymmetric (`cyl_coord = T`, `p = 0`) | Supported (single level) | Geometric sources read the live (swapped) grid; the axis half-width cell's per-cell WENO coefficients are recomputed for each block on swap/restore; blocks stay `buff_size` off the axis (domain-edge clamp); the axis-singularity viscous treatment runs on the coarse pass only. Cell volume scales with radius, so the fold-back is radius-weighted (fine `y_cc`) and the reflux area-weights radial faces (`r_face/r_cell`) and axial fine-flux averages (fine-face radius) - conservation is exact to machine precision. Restricted to `amr_max_level = 1` and not combined with non-polytropic QBMM (their parent-frame / `pb`-`mv` fold-backs are not yet radius-weighted; both are checker-gated) |
| 3D cylindrical (`cyl_coord = T`, `p > 0`) | **Not supported** | The per-stage azimuthal Fourier filter is a global operation incompatible with the block-local fine advance |
| Grid stretching (`stretch_x[y,z] = T`) | Supported | Fine ghost-shell coordinates extend by exact parent-cell bisection, and the spacing-dependent WENO coefficients are recomputed for the active grid on every block swap/restore (`amr_weno_coef_recompute`, armed automatically when the grid is nonuniform); prolongation stays conservative but its slope estimate is first-order on nonuniform parents. Stretched grids do NOT combine with Lagrangian bubbles or dynamic regrid with immersed bodies (their position-to-cell-index conversions assume uniform spacing; init abort) |
| Riemann-extrapolation BCs (`bc = -4`) | **Not supported** | Boundary-adjusted WENO coefficient rows cannot be inherited by interior blocks (checker gate) |
| `active_box` | Supported (single-rank, per active_box's own MPI gate) | Blocks must sit strictly inside the monotonically-growing active window (init abort + regrid clamp: the windowed coarse update would drop reflux corrections at faces outside it); the fine advance disables the coarse-indexed windowing and treats its whole block as active; the frozen exterior is valid ambient data for ghost prolongation |
| `acoustic_source` | Supported | The source acts on the coarse grid only: its support must not overlap the initial block (startup abort), and dynamic regrid keeps its boxes clear of the support (tags suppressed, candidate boxes clipped); emitted waves enter blocks through the coarse/fine coupling |

**Mandatory solver settings.**

AMR requires:
- WENO reconstruction (`recon_type = 1`, any order) or the IGR solver (`igr = T`)
- SSP-RK3 time-stepping: `time_stepper = 3`
- 5- or 6-equation model: `model_eqns = 2` or `3` (for 6-eq the per-stage pressure relaxation also runs on each fine block)

---

## Parameters {#amr-parameters}

The table below summarises the AMR parameters. For the full parameter descriptions,
default values, and cross-parameter constraints see @ref case section 7.1.

| Parameter | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| `amr` | Logical | F | Enable AMR (off by default) |
| `amr_block_beg(i)` | Integer | 0 | Initial block start cell index in direction `i` (level-0 index space) |
| `amr_block_end(i)` | Integer | 0 | Initial block end cell index in direction `i` (level-0 index space) |
| `amr_regrid_int` | Integer | 0 | Coarse steps between regrid events; 0 = static block |
| `amr_tag_eps` | Real | 0.1 | Normalized density-gradient threshold for refinement tagging; required `> 0` when `amr_regrid_int > 0` |
| `amr_buf` | Integer | 3 | Coarse-cell padding around tagged cells; required `>= 1` when `amr_regrid_int > 0` |
| `amr_subcycle` | Logical | F | Advance fine level at `dt/2` (two substeps per coarse step) with Berger-Colella refluxing |
| `amr_max_blocks` | Integer | 4 | Number of fixed refined-block slots preallocated; each slot is max-block sized (~N times device memory for N slots) |
| `amr_max_level` | Integer | 1 | Maximum refinement depth: `1` = single refined level, `> 1` = recursive multi-level nesting (needs `amr_max_blocks >= 2` and `ref_ratio = 2`). See @ref amr_multilevel |
| `ref_ratio` | Integer | 2 | Cell-refinement ratio between adjacent levels; must be 2 or 4. `ref_ratio = 4` is single-level only (no nesting, no subcycling) |
| `amr_cluster_eff` | Real | 0.7 | Berger-Rigoutsos min tag efficiency a clustered box reaches before splitting stops; must satisfy `0 < amr_cluster_eff <= 1` |

---

## Usage Example {#amr-example}

The minimal configuration for a 1D run with a static refined block spanning cells 16
through 47 (level-0 index space):

```python
# case.py excerpt — 1D single-fluid AMR, static block
{
    # ... base grid, patch, time-stepping settings ...
    'recon_type'     : 1,   # WENO (required)
    'time_stepper'   : 3,   # SSP-RK3 (required)
    'model_eqns'     : 2,   # 5- or 6-equation model (2 or 3)

    # AMR
    'amr'            : 'T',
    'amr_block_beg(1)': 16,
    'amr_block_end(1)': 47,
    'amr_regrid_int' : 0,   # static block; set > 0 for dynamic regrid
}
```

For a dynamic run that tracks shocks, add:

```python
    'amr_regrid_int' : 10,    # regrid every 10 coarse steps
    'amr_tag_eps'    : 0.1,   # normalized density-gradient threshold
    'amr_buf'        : 3,     # 3-cell buffer padding
    'amr_subcycle'   : 'T',   # fine level at dt/2
```

For multi-fluid (5-equation), additionally set:

```python
    'num_fluids'     : 2,
    'mpp_lim'        : 'T',   # required for num_fluids > 1 under AMR
```

---

## Limitations and Notes {#amr-limitations}

- **Fixed-slot memory.** All `amr_max_blocks` slots are allocated at init at maximum size.
  More blocks = more device memory. Compact per-block memory pools are future work.
- **Restart across rank counts.** `parallel_io` restart repartitions the fine blocks
  across any `num_procs` (each block is one contiguous region under whole-block ownership).
  The serial (per-rank-file) restart path requires the same `num_procs` and aborts with a
  clear message otherwise.
- **Half-subdomain limit.** Each block may span at most half of any rank's local subdomain
  per dimension, because the fine advance reuses the rank-local solver scratch.
- **Multi-level constraints.** Recursive multi-level nesting (`amr_max_level > 1`) requires
  `ref_ratio = 2`; with immersed boundaries it is single-rank only, and a moving body is
  not yet supported. `ref_ratio = 4` is single-level only.
- **Level-0 output only.** Standard visualization output (HDF5/SILO) is written at
  level-0 resolution; the restricted fine solution is already folded into the coarse
  fields over the block region, so existing visualization workflows are unchanged.
  Fine-resolution output is future work.
- **Unsupported physics.** See the [physics matrix](#amr-physics) above. The restrictions
  are enforced by the checker and reflect the validation state at the time of writing.

<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-07-04</div>

## Design notes {#amr-design-notes}

Internal design / implementation records (how the code works, not how to configure it):

- @subpage amr_multilevel — multi-level nesting design and reflux
- @subpage amr_fine_distribution — fine-block distribution across MPI ranks
