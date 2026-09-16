# AMR Implementation Reference

`amr.md` describes *what* MFC's AMR does. This document describes *how it is built*: the index
spaces, the data structures, the exact call sequences, and the invariants that hold between them.

It exists because the conceptual description is not sufficient to modify the code safely. Defects
in this code typically come from a wrong belief about one of the items below, not from a wrong
algorithm. Sections 2 and 11 are the two most load-bearing.

Routines are cited by name rather than line number so this document ages gracefully; every name
here is greppable in `src/simulation/m_amr*.fpp` (see the module map below),
`src/simulation/m_global_parameters.fpp`, or `src/simulation/m_time_steppers.fpp`.

---

## 1. Module map

| File | Responsibility |
|---|---|
| `src/simulation/m_global_parameters.fpp` | Slot pool metadata, working-slot mirrors, `s_amr_select_slot` |
| `src/simulation/m_amr.fpp` | AMR driver: module init/finalize; re-exports the eight per-block modules below so `use m_amr` reaches all of them |
| `src/simulation/m_amr_state.fpp` | All per-block module state (slot pool, flat store, wave buffers, `sw_*` swap mirrors) and leaf helpers |
| `src/simulation/m_amr_distribution.fpp` | Rank decomposition, box arithmetic, SFC cut, block ownership and owned-block lists |
| `src/simulation/m_amr_store.fpp` | Flat store slot machinery: dense index, reserve, alloc/free/reconcile, prim and bridge loaders |
| `src/simulation/m_amr_exchange.fpp` | Gather plans, pack/unpack, seams, ghost fills, coarse halo, and the stage/parent fill waves |
| `src/simulation/m_amr_frame.fpp` | The grid-state swap (`s_amr_swap_to_fine`/`s_amr_restore_coarse`) and pb/mv side-state services |
| `src/simulation/m_amr_transfer.fpp` | Prolongation, restriction, reflux and the restrict/reflux/freg waves |
| `src/simulation/m_amr_advance.fpp` | Fine stage advance (fused and batched), IB/Lagrange fine services |
| `src/simulation/m_amr_l0.fpp` | Level-0 tiling (`s_l0_*`) |
| `src/simulation/m_amr_regrid.fpp` | Tagging, clustering, nesting, box shaping, slot rebuild |
| `src/simulation/m_amr_registers.fpp` | Flux registers (capture and application) |
| `src/simulation/m_amr_restart.fpp` | Checkpoint of the hierarchy |
| `src/simulation/m_time_steppers.fpp` | The driver: where the AMR calls are sequenced within an RK stage |

The eight `m_amr_*` modules form a chain (each may `use` only the ones listed above it); `m_amr_exchange.fpp` is the
largest. Every symbol they export is re-exported by `m_amr`, so callers keep writing `use m_amr, only: ...`.

---

## 2. The four index spaces

**This is the section to read before touching anything.** Four distinct integer index spaces are
in play, they are all `integer`, and nothing in the type system distinguishes them. Confusing two
of them produces a silent wrong answer, never a compile error.

### 2.1 Global slot index: `k`, and the working value `amr_cur`

Range `1 : amr_max_blocks`. **Replicated on every rank**: every rank knows the level, owner, and
region of every block in the hierarchy, whether or not it owns it. This is what indexes:

```
amr_block_level(k)      amr_block_owner(k)      amr_slot_live(k)
amr_region_lo_all(:,k)  amr_region_hi_all(:,k)
amr_isect_lo_all(:,k)   amr_isect_hi_all(:,k)   amr_owns_all(k)
```

The pool is **tiles-prefix partitioned**:

```
[1 : l0_slot_off]                              level-0 L0 tiles
[l0_slot_off+1 : l0_slot_off+amr_max_fine]     regrid-managed fine blocks
```

With the L0-tile path disabled (`l0_ntile = 0`, the common configuration) `l0_slot_off = 0` and
`amr_max_fine = amr_max_blocks`.

`amr_cur` is the *working* slot. `s_amr_select_slot(islot)` sets `amr_cur = islot` and copies that
slot's region/intersection/ownership into the scalar working mirrors (`amr_region_lo`,
`amr_isect_lo`, `amr_rank_owns_block`, ...) that the per-block machinery reads. Most per-block
routines take no slot argument: **they read `amr_cur` implicitly**. Calling one without a
preceding `s_amr_select_slot` operates on whatever block was selected last.

`amr_num_blocks` is the number of currently active slots, and the loops are
`do islot = 1, amr_num_blocks` **on every rank**, with a `cycle` for blocks the rank does not own.
The loop is therefore collective in structure even where the work is not.

### 2.2 Local dense index: `loc = amr_loc_of(k)`

Range `1 : amr_loc_n`, **per-rank and private**. Defined only for blocks this rank owns;
`amr_loc_of(k) == 0` means "not live here". This is the only thing that indexes the flat field
store's slot dimension.

The two spaces are related solely through `amr_loc_of`. There is no inverse map. `loc` values are
handed out by a counter with a recycle stack, **not derived from the live set**; see §5 and
§11.3 for why that distinction matters.

> **The two per-block arrays disagree on index space.** The field store is indexed by the *local*
> dense index; the **flux registers are indexed by the *global* slot index**
> (`freg(d)%%lo(:,:,:,amr_cur)` in `m_amr_registers.fpp`) and sized `1:amr_reg_cap`, which grows
> toward `amr_max_blocks`. So `amr_cons_st(...,amr_loc_of(k))` and `freg(d)%%lo(...,k)` refer to the
> same block through different subscripts. Passing one where the other is expected compiles
> cleanly and silently reads another block's data.

### 2.3 Box index: `k` within a regrid, `ks = f_l0_slot(k)`

During `s_amr_regrid_*`, freshly clustered boxes are numbered `1 : nboxes`. This is a *third*
space, mapped to global slots by `f_l0_slot(k)`. The regrid code holds **old** and **new** box sets
simultaneously (`old_np`, `old_ilo`, `old_ext`, `old_level`, `old_owns` describe the pre-regrid
set), and both are indexed in this space. Reading `old_*` with a new box index, or vice versa, is
the classic regrid bug.

### 2.4 Cell index

Block regions are stored **in level-0 cell indices at every level**. A level-`l` block's fine
extent per dimension is `amr_ref_ratio**l * region_width - 1`. So `amr_region_lo_all` is *not* in
the block's own resolution, and comparing a region bound against a fine-grid loop bound requires
the scale factor.

---

## 3. The flat field store

```fortran
real(stp), allocatable, dimension(:,:,:,:,:) :: amr_cons_st, amr_stor_st, amr_gst_a, amr_gst_b
!                                               (x, y, z, var, LOCAL slot)
$:GPU_DECLARE(create='[amr_cons_st, amr_stor_st, amr_gst_a, amr_gst_b]')
integer :: amr_st_cap = 0
```

One contiguous device-resident array per role, replacing a per-slot vector of independently
allocated `scalar_field`s. This is the layout AMReX's `MultiFab` single-chunk path uses and the
prerequisite for batching one kernel over all blocks.

| Array | Role |
|---|---|
| `amr_cons_st` | Conserved state. **Authoritative** for `q_cons` on fine blocks. |
| `amr_stor_st` | Regrid stash: the pre-regrid state, read to seed new blocks |

Every slot carries the same buffered extents (`mbuf*_lo : mbuf*_hi`), so one array serves them all.
The cost of that uniformity is that every slot is sized for the *largest* block, not its own.

`amr_st_cap` is the number of local slots the store is sized for.

### 3.1 The copy bridge

`s_compute_rhs`, `s_ibm_correct_state`, `s_pressure_relaxation_procedure` and
`s_infinite_relaxation_k` all take `type(scalar_field), dimension(sys_size)` and also serve the
monolithic (non-AMR) path, so the flat store cannot be passed to them. A pointer view into the
store is **not attachable on the OpenMP-offload backend** (see `misc/amr_ledger/block_batching.md`).

A single block-shaped `scalar_field` array bridges instead:

```
amr_cons_br  --s_amr_br_load-->  call  --s_amr_br_store-->  amr_cons_st
```

All four dummies are `intent(inout)` (`s_compute_rhs` writes the buffer region through
`s_populate_variables_buffers`), so **both directions are required at every crossing**. This
round trip, once per block or batch per RK stage, is a first-order cost.

`amr_br_batch` is the number of blocks the bridge spans: `amr_bat_max` (8) under
`amr_batched_advance`, otherwise 1. In the batched advance (`s_amr_fine_stage_advance_batched`)
owned blocks of equal level and extent are stacked two ghost shells apart along the last active
dimension and served by one `s_compute_rhs` call; a smaller block joins a batch led by a larger
one when padding it to the leader wastes at most a fixed fraction (0.10) of its cells. The
per-block advance (`s_amr_fine_stage_advance`) remains the path for physics whose fine-advance
hooks are per-block.

---

## 4. Per-block geometry and the swap

Because the shared solver reads module-scope grid state (`x_cb`, `dx`, `idwint`, `idwbuff`, ...),
running it on a fine block requires **temporarily overwriting the global grid state** with that
block's:

```
s_amr_swap_to_fine     -> rebuild coords, extents, WENO coefficients; push to device
   ... solver call ...
s_amr_restore_coarse   -> put the coarse grid state back
```

`s_amr_swap_to_fine` performs several host array copies plus heap allocations and calls
`s_amr_sync_grid_state_to_device` (several `GPU_UPDATE` groups); `s_amr_restore_coarse` does the
reverse. The swap is **not** a leading cost, but it is why the fine advance cannot simply be
hoisted out of the per-block loop.

---

## 5. Slot lifecycle and store sizing

Six routines, in `m_amr_store.fpp`:

| Routine | Effect |
|---|---|
| `s_amr_alloc_slot(islot)` | Full slot: dense index + geometry (+ QBMM side-state; + per-slot `rhs`/gated `q_prim` for L0 **tile** slots only). Idempotent for full slots; **upgrades a live stash-only slot in place** (keeps its index and stor data, adds the arrays). |
| `s_amr_alloc_slot_stash(islot)` | **Stash-only slot**: dense index + `amr_slot_live` only, no geometry or field arrays. Used for migration replicas, which only ever touch their `amr_stor_st` half. |
| `s_amr_free_slot(islot)` | Idempotent, handles both flavors (each array family's teardown guarded on its own `allocated`; the full-vs-stash discriminator is `allocated(x_cb)`). Pushes the index onto `amr_loc_free`. |
| `s_amr_st_reserve(nloc)` | Grows the store, increments capped at 16 slots. Never shrinks the allocation. |
| `s_amr_scr_init()` | Allocates the **pooled `q_prim`/`rhs` advance scratch** (`amr_scr_prim`/`amr_scr_rhs`), once, on every rank. Called where `mbuf*` are final per mode: module init (pure AMR) or after `s_l0_tiles_init`'s mbuf union (any tiles). |
| `s_amr_compact_store()` | Re-densifies the local index space **in place on the device**, every reconcile. Does NOT realloc. |

**Fine blocks carry NO per-slot `q_prim`/`rhs`.**
The fused per-block advance (`s_amr_fine_stage_advance`: rhs then rk on one block) leaves no
cross-block lifetime, so all fine blocks share one slot-shaped scratch pair; the stage routines
take the target arrays as dummies and the caller chooses (fine: scratch; L0 tiles: per-slot,
because all owned tiles' rhs coexist across the MPI reflux point, and a tile's `q_prim` is read
in the RK pass after other tiles' RHS work; allocated per-slot exactly when the `m_rhs`
copy-out gate writes it: `run_time_info|probe_wrt|ib|bubbles_lagrange`). Sharing the scratch
keeps the per-slot footprint to its share of the flat store and avoids per-slot alloc/free
churn on the device allocator. Per-slot device cost is the slot's share of the flat
store (`amr_cons_st` + `amr_stor_st`) plus, under QBMM only, the per-slot side-state.

### 5.1 How the store grows and the host-coherence contract

`s_amr_st_reserve` early-returns when `nloc <= amr_st_cap`. Otherwise it grows by
`max(min(oldcap/4, 16), 8)` slots and reallocates, staging **through the host**:

```
GPU_UPDATE(host=...) -> tmp = store -> DEALLOCATE -> ALLOCATE(newcap) -> zero -> copy -> GPU_UPDATE(device=...)
```

The increment cap matters at scale: a proportional increment is itself store-sized, and its
transient (old plus new allocation) can exceed device memory.

**The host round trip is LOAD-BEARING, not incidental.** The rebuild's overlap carry-forward
(`s_amr_regrid_rebuild_slots`) reads `amr_stor_st` on the HOST without pulling first. It relies
on every store resize leaving host == device, so that the migration stash's host writes (pushed
to device at write time) survive a mid-rebuild grow. A device-only resize breaks this path
silently (wrong answers under regrid churn, not a crash). Until the carry-forward is converted
to a device kernel, any change to `s_amr_st_reserve` must preserve: *after a resize, the full
host copy equals the device copy*.

### 5.2 Re-densification of the local index space

`amr_loc_n` is `(live slots) + (recycle-stack depth)` mid-rebuild, and
`s_amr_regrid_rebuild_slots` allocates new blocks while old stashes are live. Four mechanisms
keep the store's capacity bounded:

1. **`s_amr_compact_store()` runs at EVERY reconcile**, renumbering `amr_loc_of` densely and
   moving slot data **in place on the device** (`s_amr_st_move_slot`, no staging array, no
   realloc). In-place is safe because moves are processed in ascending source order: each
   destination (the rank of its source among live indices) is <= its source, and every pending
   source lies above the current destination. Moved slots' host copies go stale; that is the
   store's normal state between rebuilds, and §5.1's contract applies to *resizes*, not moves.
2. So `amr_loc_n` equals the live count at every reconcile, and `amr_st_cap` plateaus at
   the **rebuild-transient high-water** (old + new generations coexisting) instead of growing
   for the whole run. Without re-densification, at np >= 2 the owned set is a shifting SFC
   window plus received replicas, so the union of indices ever handed out grows without bound.
3. The rebuild loop **early-frees** old-only slots the moment their last covering new box is
   built (`last_use` per old block, geometric region overlap as a conservative superset of every
   stash read), so freed indices recycle into the very next allocations and the transient union
   itself stays small.
4. Replicas received in migration are **stash-only** slots (see table above).

**Instrumentation:** `[amr-cap] rank r live n cap c` prints on stderr at every reconcile when
`rank_time_wrt` is on. Per-rank store capacity is the invariant quantity (device memory =
f(live local boxes)), and wall time cannot see it. `cap - live` is the transient envelope.

> **Comparison.** AMReX and Parthenon have no analogous quantity, because neither *allocates* a
> local index. AMReX's `localindex` is a binary search over the sorted owned-box list; Parthenon's
> `lid = n - nbs` is recomputed contiguous every regrid. In both, the index space *is* the live
> set by construction, so it cannot grow. Re-densification gets MFC the same invariant at
> reconcile granularity; deleting `amr_loc_free`/`amr_loc_nfree` outright ("derive the index")
> is a possible cleanliness follow-up. See §11.3.

---

## 6. The timestep

Within each RK stage `s`, in `s_tvd_rk` (`m_time_steppers.fpp`):

```
PH_COARSE   s_compute_rhs on the coarse (level-0) grid                          [1 per stage]

if (amr):
  PH_HALO   s_amr_exchange_coarse_cons_halo                                     [1 per stage]

  PH_GATHER s_amr_stage_fill_wave     ! ALL level-1 fills as ONE F1+F3 wave     [1 per stage]
            do ilev = 2, amr_num_levels
  PH_GATHER   s_amr_parent_fill_wave(ilev)  ! level-lev fills as one F2 wave    [1 per LEVEL]

  PH_SEAM   s_amr_fine_fine_halo(0)   ! all levels together                     [1 per stage]

            do islot = 1, amr_num_blocks           ! skip level 0
  PH_RHS      s_amr_fine_stage_advance                                          [1 per BLOCK]

            do islot = 1, amr_num_blocks           ! skip level 0 AND level >= 2
  PH_REFLUX   PH_RFP2P  s_amr_p2p_reflux_faces                                  [1 per L1 BLOCK]
              PH_RFAPP  s_amr_apply_reflux
```

Then after the stage loop, in reverse slot order:

```
PH_RESTR    do islot = amr_num_blocks, 1, -1
              s_restrict_fine_to_coarse
              s_amr_reflux_to_parent
```

**Multiplicity is the thing to notice.** `PH_COARSE`, `PH_HALO` and `PH_SEAM` are once per stage;
everything else is once per *block* (or per batch under `amr_batched_advance`). With hundreds of
blocks that is the difference between a handful of calls and hundreds of calls per step.

Three facts that are easy to get wrong:

- **`PH_COARSE` is not part of `PH_RHS`.** The coarse `s_compute_rhs` is bracketed separately.
  Charging it to `PH_RHS` overstates the fine advance's share of device time.
- **The reflux loop skips level >= 2**, so `PH_REFLUX` covers level-1 blocks only.
- **The loops run on every rank** and `cycle` on non-owned blocks, so loop *trip count* is global
  even though the work is local.

---

## 7. The level-1 / level-2 divergence

`s_amr_gather_coarse_patch`, the routine that fills a fine block's coarse-side data, branches at
the top:

```fortran
if (amr_block_level(amr_cur) >= 2) then
    call s_amr_gather_from_parent(pull_host)   ! entirely different path
    return
end if
```

**Everything after that branch is level-1-only code.** In a deep hierarchy most blocks are
level >= 2, so most blocks never execute the routine's main body. Any instrumentation or
optimization of coarse-patch gathering **must cover both branches** or it reports on a minority
of the work.

The level >= 2 path (`s_amr_gather_from_parent`) does a parent-to-child point-to-point exchange:
the parent's owner packs the patch on device and sends it; the child's owner receives and unpacks.
The send side uses a deferred `ISEND` pool (§8.2).

---

## 8. Communication

### 8.1 Inventory

The `m_amr_*` modules (mostly `m_amr_exchange.fpp`, `m_amr_transfer.fpp` and `m_amr_l0.fpp`) mix blocking `MPI_SEND`/`MPI_RECV` sites with `MPI_ISEND`/`MPI_IRECV`/`MPI_WAITALL`
waves and a few `MPI_ALLREDUCE`s; `m_amr_regrid.fpp` adds the regrid collectives
(`MPI_ALLREDUCE`, `MPI_ALLGATHER(V)`).

**Static call sites are not expanded counts.** Many sit inside Fypp `#:for` loops over
dimensions, so a `grep` undercounts what actually executes: `s_amr_p2p_reflux_faces` shows 2
`MPI_RECV` call sites but issues **6** in 3D (three dimensions x lo/hi). Read the Fypp, not the
grep, when sizing a communication path.

**Blocking calls in a globally ordered block loop with fixed tags** produce head-of-line
blocking: rank A cannot progress past block *i*'s rendezvous even when block *i+1*'s partner is
ready. This convoy effect is why MPI time in the AMR exchange is dominated by *wait* rather than
bandwidth.

### 8.2 The deferred send pool

`amr_gsnd_pool` / `s_amr_gsnd_reserve` / `s_amr_gather_send_flush`, capacity `amr_gsnd_max = 64`.

The parent-gather send site packs each box on device into a pooled buffer and posts an
`MPI_ISEND`, with a single drain, instead of `allocate / MPI_SEND / deallocate` per box.

The safety rule: any call site whose original semantics were "the send has completed when this
returns" must be followed by `s_amr_gather_send_flush()`. Two such sites exist in
`s_amr_subtree_stage_advance` and carry that call with a comment. **A deferred send with no
downstream drain is a deadlock.**

Hoisting the drain to once per step does not help: deferring sends pays only where a downstream
synchronization already absorbs the timing drift.

### 8.3 Reflux exchange

`s_amr_p2p_reflux_faces`: the block owner posts `2 * num_dims` `ISEND`s per participating rank
followed by one `WAITALL`; **each participating non-owner does `2 * num_dims` blocking `MPI_RECV`s**
(6 in 3D) on fixed tags `2*d` and `2*d+1`, i.e. tags 2-7, then a `GPU_UPDATE`. Both the sends and
the receives are generated by a Fypp `#:for` over dimensions guarded by `if (d <= num_dims)`.
Reflux is almost entirely communication, and its cost is per-call *wait*, not bytes.

Reflux is the **sink** for load skew, not its source: imbalance originates in `PH_RHS` and is
absorbed here.

---

## 9. Regrid

`s_amr_regrid` (`m_amr_regrid.fpp`), in order, with three early exits:

```
PH_RGHALO    s_amr_exchange_coarse_cons_halo
             s_amr_compute_lag_supp                    (bubbles_lagrange only)
PH_RGTAG     s_amr_regrid_tag_cells                    -> tag_grid
PH_RGCLUS    s_amr_regrid_cluster_tags                 -> boxes, nboxes
             if (nboxes == 0) return                   ! nothing tagged anywhere
PH_RGSHAPE   s_amr_regrid_shape_boxes
             if (nboxes == 0) return                   ! all boxes died in the domain margin
             s_amr_regrid_nest_children
             s_amr_check_box_caps                      ! invariant: no box exceeds its level cap
             s_amr_regrid_boxes_unchanged -> same
             if (same) return                          ! identical box set: keep live slots
PH_RGMIG     s_amr_regrid_stash_migrate                -> amr_stor_st, ownership
PH_RGBUILD   s_amr_regrid_rebuild_slots                -> amr_cons_st, new slots
```

The `same` early exit is what makes regrid cost bimodal: a regrid that changes nothing is nearly
free, so averaging cost over regrid *calls* understates the cost of a real one. `amr_snap` widens
this exit by letting a new box within a few cells of a live same-level block take that block's
box.

### 9.1 stash_migrate

Writes each old block's state into `amr_stor_st` at its **old** local index, exchanges blocks whose
owner changed, sets the new `amr_num_blocks`, and calls `s_amr_assign_block_owners`
(SFC/work-balanced, `s_amr_sfc_cut`).

### 9.2 rebuild_slots and the ordering constraint

```fortran
do k = 1, nboxes
    ks = f_l0_slot(k)
    ...
    if (amr_block_owner(ks) == proc_rank) call s_amr_alloc_slot(ks)
    ...
    amr_cons_st(fi,fj,fk,i, amr_loc_of(ks)) = amr_stor_st(ofi,ofj,ofk,i, amr_loc_of(kks))
    !                       ^^ NEW slot                                  ^^ OLD slot, kks = f_l0_slot(kk)
```

**This line is why old slots cannot be freed wholesale before the rebuild.** The stash is read at
*old* local indices for `kk = 1, old_np` while new slots are being allocated, so both index sets
must be simultaneously valid. A blanket free-before-allocate would produce silent wrong answers,
not a crash. Instead the rebuild loop precomputes `last_use(kk)` (the last new box whose region
overlaps old block `kk`, a conservative geometric superset of every stash read) and frees each
old-only slot at the top of iteration `last_use(kk) + 1`. The freed dense indices recycle into
the very next allocations, which is what keeps the transient union (and with it the store's
capacity plateau) small.

`s_amr_reconcile_slots` runs at the end (`PH_RBREC`), frees the remaining now-dead slots, and
calls `s_amr_compact_store` (which re-densifies the index space every time, §5.2).

---

## 10. Phase instrumentation

`m_phase_timing.fpp` defines `PH_N` phase ids with `s_phase_tic` / `s_phase_toc`. Both ends
issue `GPU_WAIT()`, so brackets measure completed device work, not launch return.

Reporting is gated on `rank_time_wrt`. **Benchmark cases usually set it `.false.`** to keep I/O out
of the timings, so the phase budget is invisible unless it is turned on deliberately. The
counters accumulate regardless; enabling reporting costs two `MPI_ALLREDUCE`s at finalize.

The report emits mean/max/imbalance per phase, a calls column and ms/call (from an `MPI_ALLREDUCE`
over an `ncall` array), and per-rank lines for `PH_RHS`, `PH_REFLUX`, `PH_GATHER`, `PH_SEAM`.

Five traps:

1. **Zero-time phases are dropped from the report** (`if (gsum(i) <= 0._wp) cycle`). A phase id
   that is declared but never `tic`'d is indistinguishable from one whose code is missing. When
   adding a phase, check that `PH_N` matches the name count, that every used id is imported, and
   that every declared id is actually `tic`/`toc`'d.
2. `PH_NAME` is `len=8`; longer names silently truncate and can collide.
3. `RESIDUAL` is meaningless once brackets nest: nested phases double-count against wall.
4. A phase must bracket **both** branches of §7 or it reports on a minority of blocks.
5. Per-call cost, not total, is the diagnostic quantity when block counts differ between runs.

---

## 11. Invariants and traps

### 11.1 Checkable invariants

| Invariant | Where it can break |
|---|---|
| `amr_loc_n == live + amr_loc_nfree` | Recycle-stack bookkeeping in `alloc`/`free` |
| `amr_loc_of(k) > 0` iff rank owns live block `k` | `alloc`/`free` ordering in rebuild |
| `amr_st_cap >= amr_loc_n` | `s_amr_st_reserve` post-condition |
| Every `@:ALLOCATE` has a matching `@:DEALLOCATE` | The store's realloc paths |
| No box exceeds its level's slot cap | `s_amr_check_box_caps` |
| `amr_reg_cap >= amr_num_blocks` | Flux registers are indexed by *global* slot (§2.2) |
| Block regions are in L0 cell indices at all levels | Any new geometry code |

### 11.2 Traps

- **Implicit `amr_cur`.** Per-block routines take no slot argument. A missing `s_amr_select_slot`
  silently operates on the previous block.
- **Two-space confusion in regrid.** `old_*` arrays are indexed by old box index; the live arrays
  by new. Both are `integer` and both are in scope.
- **The level >= 2 early return** (§7) hides the majority of blocks from anything measuring the
  level-1 body.
- **Deferred send with no drain** deadlocks (§8.2).
- **`amr_buf` does double duty**: it sets the ghost width *and* enters the box-merge threshold
  `thr = buff_size + 2*amr_buf`. Changing it to tune ghosts also changes the box topology.
- **`rocm-smi` GPU[N] enumeration need not match HIP device N.** Rank-to-device is
  `dev = mod(local_rank, devNum)`; do not infer the slow rank from the busiest GPU row.
- **Recycle-stack depth is tautological when read at a high-water print**: it is printed only when
  the high-water rises, which requires the stack to be empty. Instrument where the value can vary.

### 11.3 Known structural weaknesses

1. **The local index is allocated, not derived** (§2.2, §5.2). AMReX derives it by binary search
   over the sorted owned-box list; Parthenon recomputes `lid` contiguous every regrid. The index
   space is re-densified at every reconcile, which gets the same invariant at reconcile
   granularity; deriving it outright (deleting `amr_loc_free`/`amr_loc_nfree`) is a cleanliness
   item, not a memory one.
2. **Store GROWTH stages through the host** (§5.1), and this is a documented contract, not merely
   a weakness: the rebuild's overlap carry-forward host-reads depend on resize leaving
   host == device. The obvious AMReX-style fix (full device-side remake into a staging array)
   does not work as-is: its old+staging transient doubles the device footprint at the moment of
   growth, and dropping the host round trip breaks the carry-forward. The path that remains:
   convert the carry-forward to a device kernel first, then growth can go device-side.
3. **Slots are sized for the maximum block, not their own** (§3). AMReX's single-chunk arena sums
   actual per-box bytes. This is an independent multiplier equal to the max/mean box-volume ratio;
   measure that ratio before deciding it is worth the change to non-uniform striding.
4. **Blocking-dominated communication** (§8.1) with fixed tags over a globally ordered loop.

---

## 12. Related documents

| Document | Content |
|---|---|
| `docs/documentation/amr.md` | User-facing overview, parameters, usage |
| `misc/amr_ledger/README.md` | What the performance campaign learned, with pointers into its design memos and records |

## AMR levels (silent-index traps)

- **A level-`l` block's fine extent is `amr_ref_ratio**l * (coarse-region width) - 1`, NOT `amr_ref_ratio*width`.**
  The `amr_ref_ratio*width` form is correct only for the level-1 initial block; nested boxes compound by
  `amr_ref_ratio` per level (`amr_ref_ratio**level`). Every fine-extent computation uses
  `amr_ref_ratio**amr_block_level`: geometry (`s_set_amr_fine_geometry`), the restart-reader extent check,
  load-weight, `fmul`. Assuming `amr_ref_ratio*width` rejects level>=2 blocks as corrupt (for example in a
  restart reader's extent check).
- **"coarse" in the AMR coupling routines means the block's PARENT level (`l-1`), not the base grid (level 0).**
  For a level-1 block the parent IS L0; for level>=2 the block folds to/from its parent block's fine array.
  `s_amr_gather_coarse_patch`, `s_interpolate_coarse_to_fine`, and the restrict/reflux path all operate in the
  parent-fine frame; assuming L0 silently corrupts level>=2 coupling.
- **The fine advance SWAPS the coarse grid globals (`m/n/p`, `idwint/idwbuff`, coords, `acoustic_source`,
  `ab_active`) to a fine block and restores them after; see the SWAP CONTRACT block at the `sw_*` declarations
  in `m_amr_state.fpp`.** Any module-level variable DERIVED from the grid that a kernel reads during the fine advance
  must be swapped there or refreshed per fine call at its use site; if it is `GPU_DECLARE`'d, its DEVICE copy
  must be refreshed too. A stale device copy of coarse bounds reads out of range on the fine grid under **CCE
  OpenACC only** (NVHPC/CCE-omp evaluate bounds host-side). This is why `ab_int` is refreshed by an
  unconditional `GPU_UPDATE` in `s_compute_rhs`; `amr_rvw` (cyl_coord radius weights) is instead kept safe by a
  `m_checker.fpp` gate. A CPU-only or NVHPC-acc pass does not exercise this class; it is CCE-acc-specific.
