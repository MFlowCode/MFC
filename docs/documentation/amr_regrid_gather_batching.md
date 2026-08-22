# Batching the regrid coarse-patch gather

**Target:** `rb:gath` = 13.4% of wall. Regrid overall is 47.5% of wall at the matched operating
point, and `rg:build` is 28.6% of it.

> **CORRECTION (audit, 2026-08-20).** An earlier revision said "13.4% plus `rb:wait` 8.9% = 22.3%".
> **That double-counted.** `PH_RBWAIT` is bracketed *inside* `s_amr_gather_coarse_patch`
> (`m_amr.fpp:954`) while `PH_RBGATH` brackets the *call site* (`m_amr_regrid.fpp:1398`), so
> `rb:gath` (109.0 s) **already contains** `rb:wait` (72.9 s). The correct reading is: `rb:gath` is
> 13.4% of wall, **of which 67% is MPI wait**. Ten other phases nest inside it likewise
> (`PH_PGALL, PH_RBALLOC, PH_RBOWN, PH_RBPACK, PH_RBPOST, PH_RBRSV, PH_RBSEAM, PH_RBSEND,
> PH_RBUNPK, PH_RBUPD`). See
`docs/documentation/amr_action_plan.md` for the measurement.

## What happens today

`s_amr_regrid_rebuild_slots` loops over the new box set and, for **each box**, calls
`s_amr_gather_coarse_patch`, which is collective:

```
do k = 1, nboxes
    amr_cur = f_l0_slot(k)
    s_set_amr_fine_geometry(...)
    s_amr_gather_coarse_patch(q_cons_base, .false.)   ! <-- one full rendezvous PER BOX
    if (.not. amr_rank_owns_block) cycle
    s_interpolate_coarse_to_fine()                    ! consumes amr_cg
    ...overlap-copy from stashed old blocks...
end do
```

Inside that call, for a level-1 block:

- **owner**: unpacks its own overlapping coarse cells locally, then posts one `MPI_IRECV` per
  contributing rank (`amr_ovl_gather(:,amr_cur)`, tag = `amr_cur`), then `MPI_WAITALL`, then
  unpacks each slice into `amr_cg`.
- **every other overlapping rank**: packs its slice and sends it to the owner.

So the message count is `sum over boxes of (contributing ranks)` — roughly 224 boxes x 2-4
contributors = 450-900 messages per regrid — and there is a **separate `WAITALL` per box**. That
per-box rendezvous is the convoy: rank A cannot progress past box *i* even when box *i+1*'s partner
is ready.

SAMRAI fills an entire new level with **one `RefineSchedule`**; Chombo with **one `copyTo` +
`Copier`**. The gap is granularity, not algorithm.

## The constraint that shapes the design

The obvious fix — "gather everything first, then loop boxes" — **does not fit in memory**, and this
is the thing to understand before writing any code.

`amr_cg` is a *single* patch buffer sized for the largest block
(`amr_cg(i)%%sf(0:amr_cpat_hi(1), 0:amr_cpat_hi(2), 0:amr_cpat_hi(3))`). The gather fills it, and
`s_interpolate_coarse_to_fine` immediately consumes it. Holding all boxes' patches at once would
need `nboxes x patch`: a cap-64 block has a coarse patch of roughly `(64 + 2*margin)^3` cells, so with
`amr_cpat_mar = (buff_size + ref_ratio - 1)/ref_ratio + 1` and the matched case's **`sys_size = 6`**
(`model_eqns = 2`, `num_fluids = 1`) that is roughly `68^3 x 6 x 8 B` = **~15 MB per box, ~3.4 GB
for 224 boxes** (an earlier revision said 5.6 GB by assuming `sys_size ~ 10`). Against a working set
already near the 64 GiB device limit — and given that ungoverned store capacity is exactly what we
just spent a day removing — that is not acceptable.

## The design: chunked plan-then-execute

Process the box list in **chunks of `amr_gath_chunk` boxes** (default 32). Within a chunk:

1. **Plan** — for every box `k` in the chunk and every rank `r` in `amr_ovl_gather(:,k)`:
   - `r == me and owner(k) /= me` -> append `(k, bl, bh)` to `sendlist(owner(k))`
   - `owner(k) == me and r /= me` -> append `(k, bl, bh)` to `recvlist(r)`
   Both sides iterate boxes in the same global order and read the same cached overlap lists, so the
   two layouts agree by construction. **That agreement is the correctness invariant** and must be
   asserted, not assumed (see below).
2. **Exchange** — one packed `MPI_ISEND` per destination peer and one `MPI_IRECV` per source peer,
   then a **single `WAITALL` for the whole chunk**.
3. **Consume** — loop the chunk's boxes; for each, unpack its slices out of the peer buffers into
   `amr_cg`, prolong, and do the overlap-copy exactly as today.

Message count per regrid falls from ~450-900 to `ceil(nboxes/chunk) x peers` — about 7 chunks x <=7
peers = ~50. Staging memory is bounded at `chunk x patch` ~= 480 MB at chunk 32, and the chunk size
is the single knob trading memory against message count.

## Correctness invariants to assert, not assume

- **Layout agreement.** Sender and receiver must derive identical `(k, bl, bh)` sequences per peer.
  Deterministic iteration over the same cached lists gives this, but a mismatch would be a silent
  wrong answer, not a crash. Assert matching byte counts per peer before unpacking.
- **`amr_cur` is implicit.** Per-block routines read `amr_cur`; the plan phase must not leave it
  pointing at the wrong box for the consume phase.
- **Level >= 2 blocks take a different path** (`s_amr_gather_from_parent`) and are untouched by this
  work. They are the majority of blocks on the production case, so this change should be measured on
  the level-1 population specifically.
- **The `pull_host` path and the `num_procs == 1` shortcut** must keep their current semantics.
- **Deferred sends need a drain.** The existing `amr_gsnd_pool` already carries this rule; a batched
  send with no matching wait is a deadlock.

## Increments

1. **LANDED (2026-08-21 night).** Plan builder (`s_amr_build_gather_plan`, both families:
   level-1 contributor lists + sizes AND the level>=2 parent source/size) + always-on
   `@:ASSERT`s at all five derivation points in the per-box gathers, exchange still per-box
   (no behaviour change). Validated three ways: AMR subset 67/67 with asserts armed (zero
   trips = plan reproduces the live message set across the suite), and TWO seeded-bug
   tripwires that each ABORTED as required (level-1 size +1 -> "send entry mismatch" on a
   contributor rank; parent size +1 -> parent send-size assert on the parent owner; S0-style
   np=4, first rebuild). The plan may now be trusted by step 2.
2. Chunked exchange for BOTH families, `pull_host = .false.` only (see the S0 np=8 scope
   section above: pre-post the chunk's level-1 IRECVs and the level>=2 per-box IRECVs, sends
   stay pooled, consume in box order).
3. Extend to the pb/mv gather (`s_amr_gather_coarse_patch_pbmv`) if step 2 pays.

Step 1 is the safety net: if the plan does not reproduce the current message set box for box, the
batching is wrong and it is visible before any data moves.

## S0 np=8 REVERSES the level-1/level-2 priority (2026-08-21, post-P1 gate, logs/p1gate-0821_2144)

> **The section below this one concluded — correctly, for the MATCHED point — that the level-1
> WAITALL dominates and `pg:all` "cannot dominate." At the WEAK-SCALING point (S0, fixed
> 200^3/rank, amr_max_level=2) the ordering flips: the split is operating-point-dependent, and
> the increment must cover BOTH families.**

Measured sub-brackets of `rb:gath` (mean s):

| | np=4 | np=8 | mechanism |
|---|---|---|---|
| `pg:all` (level>=2 parent gather) | 15.6 | **108.0** | of which `pg:recv` 15.3 / **99.2** — the **blocking per-box `MPI_RECV`** (m_amr.fpp:1405), 252 / 208 ms/call, imb ~1.3-2.0 |
| `rb:wait` (level-1 WAITALL) | 1.2 | **58.9** | 85 calls x 693 ms at np=8 |
| everything else (pack/unpack/alloc/post/send) | ~0.15 | ~1.6 | dead — fixes aimed there stay dead |

At S0 the refined region is a deep blob: most boxes are level 2, and each one costs a pairwise
rendezvous — both ranks walk the same global box list, so the owner's blocking RECV for box k
absorbs whatever skew the parent rank accumulated (rb:slot/rb:ovl/rb:push differ by ownership).
R1 converted this site's SEND to the ISEND pool; **the RECV side is the unconverted half, and at
np=8 it is the single largest item inside `rg:build`.**

**Scope change:** the chunked plan-then-execute below covers BOTH families in one framework:
- level-1: peer send/recv lists exactly as designed below;
- level>=2: trivial plan (one known parent rank per box, sizes computable on both sides with no
  handshake) — pre-post one `MPI_IRECV` per owned level>=2 box in the chunk into a per-box
  buffer, keep the parent's pooled device-pack ISEND as-is, keep the device unpack;
- consume in box order as today; a box is exactly one level, so tag = box index stays unambiguous
  across the two families for pre-posted receives.
Chunk memory bound now covers both: level-1 staging (~480 MB at chunk 32) + per-box level-2
buffers (~15 MB each) — ~1 GB/rank at chunk 32, affordable post-P1 (>=12 GiB headroom).
**Validation tripwire:** increment 1's plan must reproduce today's message set exactly — the I1a
`XA_F1_*`/`XA_F2_*` conservation counters (m_amr_xchg_audit.fpp) must be IDENTICAL before/after.

Expected payoff at S0 np=8, stated as a bound: rb:wait + pg:recv = 158 s = 14.3% of wall; the
realistic target is the rendezvous share of it (not the bytes) — differenced runs required, and
the downstream wait shadows (reflux/gather/seam, which absorb regrid-side skew) may move too.

## The level-1 / level-2 split — RESOLVED without a run (MATCHED point; superseded for S0 above)

`PH_PGALL` (the level >= 2 parent-gather path) is nested inside `rb:gath`, so `rb:gath` covers both
the level-1 body this design batches and the level >= 2 path it does not. An earlier revision called
this a blocking unknown needing a fresh run. **It is not — the call counts settle it.**

| quantity | value | meaning |
|---|---|---|
| `rb:gath` calls / rebuild | 3615/16 = **226** | every box, every level |
| `rb:own` calls / rebuild | 128/16 = **8** | level-1 **owner branch only** |
| boxes owned per rank per rebuild | 226/8 = **28** | |
| => level-1 share of owned boxes | **28%** | level >= 2 is ~72% |

So ~72% of boxes *are* level >= 2. But the time splits the other way:

| | s | % of `rb:gath` | % wall |
|---|---|---|---|
| level-1 branch (bracketed phases, all after the line-875 early return) | **73.6** | **67.5%** | 9.0% |
| of which `rb:wait` | 72.9 | 66.9% | **8.9%** |
| remainder = `pg:all` + `rb:post` + `rb:send` + unbracketed | <= 35.4 | <= 32.5% | <= 4.3% |

**The level-1 path is the minority of boxes and the majority of the time**, because each level-1
gather performs a **593 ms owner-side `WAITALL`** (123 calls). `pg:all` is bounded above by 4.3% of
wall and cannot dominate.

`PH_RBWAIT` is at `m_amr.fpp:954`, after the level >= 2 early return at 871-875, so it is
unambiguously level-1 — verified from the source, not inferred.

**Conclusion: the design is aimed correctly.** Batching the level-1 gather targets the 8.9%-of-wall
wait term.

## Expected payoff, stated as a bound

Corrected for the nesting above:

| if we eliminated | wall | tax |
|---|---|---|
| all of `rb:gath` | x0.866 | 11.03x -> **9.56x** |
| only its wait term | x0.911 | 11.03x -> 10.04x |
| all of regrid | x0.525 | 11.03x -> 5.79x |

Batching removes the per-box *rendezvous*, not the bytes, so the realistic target is a fraction of
the wait term — call it 4-7% of wall, **which is close enough to the 4.96% noise floor that the
experiment must be differenced, not single-run.**

## One more correction from the audit

`regrid` shows 40 calls but `rg:build` and `rg:mig` show **16**: only 16 of 40 regrids actually
rebuild, the rest hit the `boxes_unchanged` early-out. So "9.7 s per regrid" (387.5/40) understates
the real thing — **a rebuilding regrid costs ~21.4 s** (`rg:build` 14.5 + `rg:mig` 6.9). Quote the
per-rebuild figure, not the per-call average.

Cross-check that validates the call accounting: `rb:gath` 3615 calls / 16 rebuilds = **225.9 boxes
per rebuild**, against a converged box count of ~224.
