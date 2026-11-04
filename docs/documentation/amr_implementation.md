@page amr_implementation AMR Implementation Reference

# AMR Implementation Reference

@ref amr says what AMR does; this page says how it is built: modules, index spaces, data layout, the step, regrid and
communication, and the invariants that are easy to break. Routines are cited by name (greppable in
`src/simulation/m_amr*.fpp`, `m_global_parameters.fpp`, `m_time_steppers.fpp`).

## Modules

Each module may `use` only those above it; `m_amr` re-exports the chain, so callers write `use m_amr, only: ...`.

| Module | Responsibility |
| :--- | :--- |
| `m_box` | Box type, per-rank partition arithmetic and the Morton key (common) |
| `m_amr_restart_io` | Fine-level restart file layout shared by the writer, the reader and post_process (common) |
| `m_global_parameters` | Block metadata pools, the working-slot mirrors, `s_amr_select_slot` |
| `m_amr_state` | Per-block module state: flat store, wave buffers, swap mirrors (`sw_*`) |
| `m_amr_distribution` | Block cost, SFC cut and owner placement, owned-block lists, box arithmetic |
| `m_amr_store` | Flat store: local index, reserve/compact, slot alloc/free, solver bridge load/store |
| `m_amr_wave` | Nonblocking exchange engine: pooled sides, post/send/wait |
| `m_amr_ghost` | Coarse-to-fine ghost fill: slab table, linear or quartic fill, quartic guard, volume-fraction closure |
| `m_amr_seam` | Fine-fine seams between same-level blocks, coarse conservative halo |
| `m_amr_exchange` | Coarse-patch gathers and the fill waves |
| `m_amr_frame` | Grid-state swap to a block or batch and back; 6-equation fine relaxation |
| `m_amr_registers` | Flux registers: capture during the RHS and reflux apply |
| `m_amr_transfer` | Prolongation, restriction, reflux and restriction waves |
| `m_amr_advance` | Batched fine RHS + RK advance, IB and Lagrangian fine services |
| `m_amr_cluster` | Berger-Rigoutsos clustering |
| `m_amr_regrid` | Tagging, box shaping, migration and slot rebuild |
| `m_amr_nest` | Child nesting for levels above 1, and the shared tag sensor |
| `m_amr_restart` | Hierarchy checkpoint (serial and MPI-IO) |
| `m_amr_stage` | Hooks called by `s_tvd_rk`: stage begin, fine stage, subcycled fine step, end-of-step fold |
| `m_phase_timing` | `rank_time_wrt` phase brackets and report |

## Index spaces

Four integer index spaces are in play and the compiler cannot tell them apart; mixing two gives a silent wrong answer.

1. **Global slot `k`** (`1:amr_num_blocks`), replicated on every rank: `amr_block_level(k)`, `amr_block_owner(k)`,
   `amr_region_lo_all(:,k)`, `amr_region_hi_all(:,k)`, `amr_owns_all(k)`. Per-block routines read the *working* slot
   `amr_cur` set by `s_amr_select_slot`; they take no slot argument. Block loops run on every rank and `cycle` on blocks
   the rank does not own. Flux registers are indexed by the global slot (through the dense map `amr_reg_of`).
2. **Local slot `loc = amr_loc_of(k)`** (`1:amr_loc_n`), per rank, defined only for owned blocks (`0` = not here). It
   indexes the last dimension of the flat store. Indices are recycled and re-densified at every reconcile
   (`s_amr_compact_store`).
3. **Box index** during regrid: new boxes `1:nboxes` and the old set (`old_*`) are both in scope; reading one with the
   other's index is the classic regrid bug.
4. **Cells.** Block regions are stored in level-0 cell indices at every level. A level-`l` block's fine extent is
   `amr_ref_ratio**l * (region width) - 1`.

"Coarse" in the coupling routines means the block's parent level, which is level 0 only for level-1 blocks.

## Flat store and the solver bridge

`amr_cons_st(x, y, z, var, loc)` holds the conserved state of every owned block (`amr_stor_st` is the regrid stash);
all slots share the buffered extents `mbuf_lo:mbuf_hi`, sized for the largest block. The shared solver
(`s_compute_rhs` and the physics hooks) takes `scalar_field` arrays, so a batch is copied into the bridge `amr_cons_br`
(`s_amr_br_load_batch`), advanced, and copied back. A batch stacks up to `amr_bat_lim` blocks of one level along the last
dimension (offsets `amr_bat_off`); a smaller block joins a larger leader when padding wastes at most `amr_bat_pad` of
its cells. Fine blocks share one pooled `q_prim`/`rhs` scratch.

The store grows by `s_amr_st_reserve` (never shrinks), on the device up to `amr_grow_dev_bytes` and through the host
above it.

## Grid-state swap

The solver reads module-level grid state (`m`, `n`, `p`, `idwint`, `idwbuff`, `x_cc`, `dx`, ...).
`s_amr_swap_to_fine` installs a block's (or batch's) grid, recomputes the spacing-dependent coefficients and pushes them
to the device; `s_amr_restore_coarse` puts the base grid back. Any module variable derived from the grid that a kernel
reads during the fine advance must be swapped there (see the swap contract at the `sw_*` declarations in
`m_amr_state`) or refreshed at its use site, including its device copy; CCE OpenACC reads a stale device copy.

AMR fails closed on physics switches: with `amr = T` the case validator rejects any logical switch set true that is not in
`AMR_VETTED_FLAGS` (`toolchain/mfc/case_validator.py`). New physics joins that list only after it meets the swap contract
above and has an AMR test.

## The step

`s_tvd_rk` calls, per RK stage:

```
s_amr_stage_begin          coarse cons halo for the coarse RHS
s_compute_rhs              level 0                                           [PH_COARSE]
s_amr_stage_fine           (lock-step only)
  fill waves               level-1 patches from level 0, then each level from its parent   [PH_HALO, PH_GATHER]
  seam halo                fine-fine faces between same-level blocks                       [PH_SEAM]
  batched advance          RHS + RK per batch, registers captured in the RHS                [PH_RHS, PH_RK, PH_SWAP]
  reflux                   level-1 faces as one wave, then one batched apply                [PH_REFLUX]
```

and after the stage loop:

```
s_amr_subcycle_fine        (subcycled only) level 1 takes two dt/2 steps with time-interpolated ghosts,
                           level 2 two steps per level-1 step; then the registers correct the coarse state
s_amr_step_fold            restrict each level into its parent, bottom-up, and reflux level>=2 into its parent  [PH_RESTR]
```

## Regrid

`s_amr_regrid`, every `amr_regrid_int` steps: halo, tag (`s_amr_regrid_tag_cells`), cluster
(`s_amr_regrid_cluster_tags`), shape (pad, clamp, cap, keep clear of bodies, sources and bubble clouds), nest, check
caps and disjointness, and return early if the box set is unchanged (`amr_snap` widens this). Otherwise
`s_amr_regrid_stash_migrate` stashes each old block at its old local index, assigns owners and ships stashes to new
owners, and `s_amr_regrid_rebuild_slots` builds the new blocks level by level: overlap with an old block is copied from
the stash, the rest is prolonged. Old slots are freed only after their last reader (`last_use`), then
`s_amr_reconcile_slots` compacts the store.

## Communication

Every exchange is a wave on `m_amr_wave`: all receives posted, sends packed and posted, one wait. Both sides derive the
same transfer list from replicated metadata, so no handshake is needed, and messages carry keyed tags, so matching does
not depend on posting order. Waves: coarse halo, level-1 and parent fill waves, seam halo, reflux faces, restriction and
level>=2 register waves, regrid migration. With `amr_device_pack` the fill waves pack and unpack in one fused kernel
per family. Regrid adds a few collectives (tag reductions, the box list, owners).

## Instrumentation

`rank_time_wrt = T` brackets the step phases (`PH_*` in `m_phase_timing`) with device-synchronizing timers and prints
mean, max and imbalance per phase, plus the block-balance report at each regrid. With it off the brackets return
immediately.

## Invariants and traps

- `amr_loc_of(k) > 0` exactly when this rank owns live block `k`; `amr_st_cap >= amr_loc_n`.
- No box exceeds its level's cap (`s_amr_check_box_caps`); same-level boxes are disjoint (`s_amr_check_box_disjoint`);
  a level>=2 box lies in exactly one parent.
- A missing `s_amr_select_slot` makes a per-block routine act on the previous block.
- `amr_buf` sets the tag padding and also enters the box-merge distance; changing it changes the box topology.
- The rank-to-GPU map is `mod(local_rank, ndevices)`; `rocm-smi` numbering need not match.
