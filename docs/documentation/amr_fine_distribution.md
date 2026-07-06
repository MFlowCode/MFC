# AMR fine-level distribution (design note)

**Status:** in progress, own branch (`amr-fine-dist`), separate from the AMR mega-PR.

## Problem

MFC's AMR uses a **mirror decomposition**: the fine level reuses the coarse rank
decomposition, so each rank owns exactly the fine cells over its coarse subdomain's
intersection with the blocks (`s_amr_compute_isect` → `amr_rank_owns_block`). This makes
coarse↔fine coupling entirely rank-local (zero communication) but leaves fine work
**unbalanced**: if refinement concentrates in a sub-region, only the ranks owning that
region do the fine advance while the rest idle. At scale this caps AMR speedup.

## Target

Give the fine level its **own distribution**, decoupled from the coarse decomposition, so
fine blocks spread across all ranks by measured work — the AMReX per-level
`DistributionMapping` idea. MFC already enforces an invariant that makes this far simpler
than the general case: **blocks are kept ≥ `buff_size` apart** (regrid merge), so with
**whole-block-per-rank** assignment there is *no fine–fine halo at all*
(`s_mpi_sendrecv_amr_fine_halo` disappears). The only coupling that becomes communication
is coarse↔fine, block-granular.

## Coupling changes (whole-block-per-rank)

| Step | Today (mirror, local) | Decoupled |
|---|---|---|
| Ownership | intersection, multi-owner (`amr_rank_owns_block`) | single owner (`amr_block_owner(k) == proc_rank`) |
| Ghost fill / prolong (coarse→fine) | read local coarse `q_cons` | **gather** coarse patch (block + `buff_size` halo) from coarse owner(s) to block owner |
| Reflux + restriction (fine→coarse) | write local coarse cells | **scatter** corrections from block owner back to coarse owner(s) |
| Fine–fine halo | coarse Cartesian neighbors | none (blocks separated) |
| Regrid | cluster (global) → mirror | cluster → **assign** (`amr_block_owner`) → **migrate** fine state |

The gather/scatter is the only new communication surface. SFC-ordered assignment keeps a
block's coarse patch on an SFC-nearby rank, bounding the cost. Conservation exactness is
preserved because the scatter applies the *same* reflux add / restrict overwrite the
mirror model applied locally — a pure data-movement change, not a numerics change.

## Reuse

`m_load_weight` (per-cell cost field) and `m_sfc_partition` (Morton order + chains-on-chains
balanced partition) already exist as an init-time diagnostic. They become the fine-level
distributor: feed the block list + per-block fine-work weight, return `amr_block_owner(:)`.

## Phasing (each independently mergeable)

1. **Distribution map, computed but not applied** (behavior-preserving). Add
   `amr_block_owner(:)` + SFC assignment + an imbalance diagnostic; mirror ownership and all
   coupling unchanged. Goldens bit-identical. *(this note's first increment)*
2. **Apply the map**: switch `amr_rank_owns_block` to single-owner, add the coarse↔fine
   gather/scatter, drop the fine–fine halo. Validate conservation + assignment-independence.
3. **Dynamic rebalance each regrid**; optionally drop the fixed max-size slot pool for
   right-sized boxes.

Phase 2 carries the correctness risk (gather/scatter exactness) and is the checkpoint gate.
