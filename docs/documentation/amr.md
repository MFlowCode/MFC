@page amr Adaptive Mesh Refinement

# Adaptive Mesh Refinement

> **Experimental.** AMR is off by default (`amr = F`) and requires explicit opt-in.
> The physics support matrix is enforced at run time by the checker; unsupported
> combinations abort with a clear message rather than producing silently wrong results.

## Overview {#amr-overview}

Block-structured adaptive mesh refinement (AMR) concentrates resolution where the flow
demands it — around shocks, interfaces, and bubble clouds — while leaving the rest of the
domain at the coarser base-grid resolution. MFC implements a two-level hierarchy: the
unmodified base (level-0) solve runs as usual, and one or more 2:1 refined rectangular
blocks advance alongside it on a finer grid.

The fine blocks are dynamically repositioned every `amr_regrid_int` coarse steps using a
gradient-based cell tagger and Berger–Rigoutsos block clustering, so they follow moving
features automatically. When `amr_regrid_int = 0` the block is fixed at the initial
`amr_block_beg`/`amr_block_end` position for the whole run.

AMR lives entirely in the `simulation` executable (`src/simulation/m_amr.fpp` and
`m_amr_registers.fpp`) and is the only part of MFC that modifies the solver's grid
mid-run.

---

## The Block-Structured Model {#amr-model}

The hierarchy has exactly two levels:

- **Level 0** — the base grid with cell spacing `dx`, `dy`, `dz`. The ordinary MFC solver
  advances this level every step.
- **Level 1** — a list of up to `amr_max_blocks` rectangular refined blocks, each covering
  a sub-region of the level-0 domain at 2:1 refinement (half the cell spacing in every
  active direction).

Each block is described by its bounding box in level-0 cell-index space
(`amr_block_beg(1:num_dims)` to `amr_block_end(1:num_dims)`). The fine-block extents must
satisfy:

```
2*(amr_block_end(i) - amr_block_beg(i) + 1) - 1 <= N_i
```

where `N_i` is the global cell count in direction `i`. This ensures the fine scratch
(which is sized to the base grid) is never overflowed.

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
level outside the block is unchanged.

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

Restart requires the same rank count as the run that wrote the fine-level file.

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
| Hypoelasticity / hyperelasticity / MHD | **Not supported** | Unvalidated with the fine-level advance |
| QBMM bubbles (polytropic) | Supported | Bubble moments live in `q_cons`, injected piecewise-constant at prolongation to preserve CHyQMOM realizability |
| QBMM bubbles (non-polytropic) | **Not supported** | `pb`/`mv` quadrature side-state is a global array the fine advance would corrupt through the swap |
| Lagrangian bubbles | **Not supported** | Lagrangian tracking not advanced on the fine level |
| Immersed boundaries (`ib = T`; single non-STL body, static or prescribed-motion `moving_ibm=1`; `amr_regrid_int = 0`) | Supported | Per-block fine-grid IB markers/ghost points, rebuilt each fine substage at the body's sub-time position for a moving body; non-conservative ghost-cell forcing at the body; force-driven (`moving_ibm=2`)/multi-body/STL/dynamic-regrid gated; a body spanning a rank seam is rejected at startup |
| IGR solver | **Not supported** | Unvalidated with the fine-level advance |
| Cylindrical / axisymmetric coordinates | **Not supported** | Unvalidated with the fine-level advance |
| `active_box` | **Not supported** | Unvalidated combination |
| `hybrid_weno` / `hybrid_riemann` | **Not supported** | Unvalidated combination |
| `acoustic_source` | **Not supported** | dt-dependent RHS source; unvalidated with the fine-level advance |

**Mandatory solver settings.**

AMR requires:
- WENO reconstruction: `recon_type = 1` (any order)
- SSP-RK3 time-stepping: `time_stepper = 3`
- 5-equation model: `model_eqns = 2`

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
    'model_eqns'     : 2,   # 5-equation model (required)

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
- **Same-rank-count restart.** The AMR restart file encodes the block geometry and fine
  solution per rank. Restarting with a different `num_procs` is not supported and aborts
  with a clear message; np-flexible restart is future work.
- **Half-subdomain limit.** Each block may span at most half of any rank's local subdomain
  per dimension, because the fine advance reuses the rank-local solver scratch.
- **Single refinement level.** Only one refined level is supported. Multi-level AMR
  (levels 2, 3, ...) is not implemented.
- **Level-0 output only.** Standard visualization output (HDF5/SILO) is written at
  level-0 resolution; the restricted fine solution is already folded into the coarse
  fields over the block region, so existing visualization workflows are unchanged.
  Fine-resolution output is future work.
- **Unsupported physics.** See the [physics matrix](#amr-physics) above. The restrictions
  are enforced by the checker and reflect the validation state at the time of writing.

<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-07-04</div>
