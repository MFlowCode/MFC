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
2. **LANDED 01cc4318 + assert re-guard 3de4724e (2026-08-22). Verdict below.** Chunked
   exchange for BOTH families, `pull_host = .false.` only (pre-post the chunk's level-1
   IRECVs and the level>=2 per-box IRECVs, sends stay pooled, consume in box order).
3. Extend to the pb/mv gather (`s_amr_coarse_patch_pbmv`) — **DEPRIORITIZED** after the
   step-2 verdict: F3 is absent from the S0 operating point (QBMM gated), and the tag
   split (`k + amr_max_blocks`) it requires buys nothing until an F3-heavy case matters.

## STEP-2 VERDICT (2026-08-22, node k004-001, on-node differenced arms)

**Correctness: proven at every np.** AMR subset 67/67 goldens; XA exchange report
line-for-line identical in FOUR comparisons (np=4 tripwire, np=4 arm, np=8 both arms);
output BIT-IDENTICAL across binaries at np=4 (1.5 GB state + 14.8 GB hierarchy file) and
np=8 (3.1 + 31.3 GB) — which also demonstrates the pipeline is bitwise deterministic
cross-node, legitimizing bit-diff as the np>2 verification method. VRAM peaks identical
card-for-card (the chunk pool is host-side).

**Wall (same node, differenced): np=4 421.1 -> 413.0 s (-1.9%); np=8 1246.0 -> 1235.6 s
(-0.8%).** Both inside the 4.96% noise floor, both the right sign. Mechanism brackets at
np=8: `rb:wait` 63.9 -> 45.2 (-29%), `rb:gath` 180.1 -> 159.6, `rg:build` -13.4 s; ~11 s
reabsorbed at `rb:tail`/`rb:xchg` (skew moved one fence later, as the design's own caveat
predicted).

**THE ATTRIBUTION CORRECTION — read before proposing more exchange work.** `pg:recv` did
NOT move (np=4: 15.3 -> 16.2 s; np=8: 106.2 -> 101.7 s) even though most level-2 parents
sit in earlier chunks and DID get the early phase-B send. The banked "158 s rendezvous
bound" was wrong: the F2 payload is CREATED BY THE REBUILD ITSELF — the parent's store
exists only after the parent owner's own consume + prolong + device push — so the child's
wait measures the build-backlog difference between the two owners plus pack/D2H/copy.
Pre-posting removes none of that. `rb:wait` was the true rendezvous share (~20 s at np=8);
level-1 sends have no data dependency (`q_cons_base` is host-current before the loop), so
they overlap fully once posted early. **Further MPI batching of the parent-gather family is
a dead end; the residual ~100 s yields only to owner-local rebuild (removing the lockstep
all-ranks box walk) or rebuild-work rebalancing.**

**Node confound recorded:** k004-001's GCD 6 runs hot/slow (rank-6 rhs 237/216 s vs ~175
mean on BOTH binaries; +8% total GPU compute vs k004-004 on byte-identical work; the seven
healthy ranks absorb it as 90-129 s of reflux wait). Never compare walls across
k004-001/k004-004 — the +11.7% cross-node scare that triggered the same-node control was
entirely this.

**Follow-ups from the expert review round (MPI + AMR auditors, 2026-08-22):** the
`amr_gcr_*` chunk machinery is scaffolding — absorb and delete it when the v2 plan-based
exchange (cached per-peer schedules) lands; add a forced-drain counter to
`s_amr_gsnd_reserve` (cap-64 + width-regrowth triggers) when next instrumenting; the
program's next constraint is DEAD BYTES (the per-step fill family ships full patches whose
interior is never read — see the action plan's re-aimed increment list).

## Step-2 design — REVIEWED 2026-08-22: two independent adversarial reviewers (MPI/deadlock
## lens + state/lifetime lens) CONVERGED on one fatal defect and four bindings. The design
## below is AMENDED accordingly; implement only the amended form.

**FATAL (found by BOTH reviewers, D1): the phase-B level>=2 pack reads an UNBUILT parent.**
The parent's new-generation store content is produced only in the parent's own phase-C
iteration (alloc m_amr_regrid.fpp:1478, fill 1491-1521, device push 1524); the box list is
parents-first but not chunk-aligned, so a same-chunk parent/child pair makes phase B pack
stale or unallocated store memory (amr_loc_of = 0 -> out-of-bounds device read). The XA
identity invariant is PROVABLY BLIND to it (same counts/sizes, wrong payload).
**AMENDMENT: level>=2 sends split by parent position — phase B ONLY when parent index <
c_lo (parent consumed in an earlier chunk = the early-send win); otherwise the pack+ISEND
stays at the child's box position in phase C (parents-first guarantees the parent's consume,
index < k, has completed). Both reviewers confirm the deadlock induction survives this.**

Further bindings from the review (all mandatory):
- **Ownership predicate = `amr_block_owner(ks)` everywhere in phases A/B.** The mirrors
  (`amr_owns_all`, `amr_rank_owns_block`) hold the PREVIOUS generation until phase C's
  geometry call (writers: m_amr.fpp:2496/2519 only).
- **`amr_cpat_off` is per-level and per-phase**: level-1 = region - margin; level>=2 =
  parent-foot - margin (parent-fine frame). Recompute in every phase that calls a kernel or
  host loop reading it; NEVER inherit it across phases (phase B's last write is arbitrary).
- **Tag collision with QBMM pbmv (F3) is real**: pbmv reuses (src, tag=k) against F1/F2.
  Correctness rests on MPI NON-OVERTAKING plus fixed posting order (F1 recv in A before
  pbmv recv in C; F1 send in B before pbmv send in C) — now a STATED invariant. The pbmv
  per-box call stays in phase C, all ranks, box order. Step 3 MUST disambiguate the tag
  space (k + amr_max_blocks) before chunking pbmv.
- **Keep the per-box ``GPU_UPDATE(device='[amr_cg]')``** (m_amr.fpp:1125-1131) in phase C —
  device coherence of amr_cg is not proven dead; dropping coherence updates on a hunch is
  the restart-NaN bug class.
- **Recorded progress assumption**: any blocking MPI call progresses ALL traffic (true for
  Open MPI 4.1 ob1/vader; the induction survives even per-request progress). No collective
  may enter the chunk loop (the only rebuild collective is s_amr_reduce_xchg_flag at
  rb:tail — keep it there).
- Request-array reuse across chunks is safe ONLY because phase C consumes every box
  unconditionally (the owner-cycle comes after the WAITALL position) — preserve that.

## Step-2 design (original text; read with the amendments above)

Restructure `s_amr_regrid_rebuild_slots`'s box loop into chunks of `CHUNK` boxes (start 32):

```
call s_amr_build_gather_plan(nboxes)
do c_lo = 1, nboxes, CHUNK
    c_hi = min(c_lo + CHUNK - 1, nboxes)
    A: post   - for each box k in chunk I OWN: IRECVs per plan entry (level-1: one per
                contributor; level>=2 split: one from the parent owner), tag k, into the
                chunk recv pool; requests appended IN BOX ORDER so each box's requests are
                a contiguous run.
    B: send   - for each box k in chunk where I contribute (level-1: plan lists me;
                level>=2: I own the parent, split): pack + pooled ISEND exactly as today
                (s_amr_gsnd_reserve/amr_gsnd_pool unchanged, forced drains included).
    C: consume - per box k in chunk, IN ORDER: the existing per-box body (early-free,
                alloc, geometry) but the gather call replaced by: set amr_cpat_off for k;
                own-slice host copy; MPI_WAITALL on k's contiguous request run; unpack
                (host for level-1 slices, device for the level-2 buffer); co-located
                level>=2 parent copy happens HERE (device copy, as today); then
                interpolate/carry-forward/push unchanged.
end do
flush (existing rb:tail)  ! sends may still be in flight across chunks - unchanged
```

Load-bearing details, each verified against source during the step-1 audit:
- **Geometry without the swap.** The post/send phases need only plan data + `bl/bh`
  (recomputable from replicated caches, as the builder does) + `start_idx` (rank-constant).
  `amr_cpat_off` is consumed host-side by the pack/unpack kernels (captured into locals
  before the kernel), so setting it per box in the phase that calls the kernel suffices;
  the full `s_set_amr_fine_geometry` swap stays in phase C where it is today.
- **Buffers.** One flat `real(wp)` chunk recv pool + an offset table, sized from the plan
  (sum of the chunk's owned-box message sizes; bound ~CHUNK x patch ~ 0.5-1 GB at 32);
  one request array; a per-box (first-request, count) index. Reused across chunks, grown
  monotonically, freed at module finalize.
- **No (src, tag) ambiguity.** Tag = box index; a box is exactly one level; level-1 sources
  are distinct ranks; the level>=2 rebuild path sends only the `_cons` message. So every
  (source, tag) pair in a chunk is unique - pre-posting cannot mismatch.
- **Deadlock-freedom argument.** Recv-posting (A) contains no MPI waits. The only blocking
  points are per-box WAITALLs (C) and the send pool's forced drain (B, cap 64). Consider
  the rank at the globally minimal (chunk, phase) position: its phase A completes
  unconditionally; its drains wait on sends whose receivers are at-or-ahead and whose
  matching recvs are posted in THEIR phase A, which they reach without waiting on anything
  from the minimal rank's current chunk. So the minimum always advances - no cycle. (The
  refuted per-STEP drain hoist does not apply: that experiment deferred sends past a sync
  that could not absorb the drift; here recvs are pre-posted and the WAITALL that pays the
  skew is per-box but no longer serializes the WHOLE exchange behind it.)
- **`num_procs = 1` and co-located parents** degenerate naturally: empty plan entries, no
  posts/sends; phase C's own-copy and co-located device copy reproduce today's path.
- **Validation:** goldens + subset must stay green; XA_F1/F2 message and word totals must be
  IDENTICAL to the per-box path (same plan, same messages, different timing); S0 np=4/np=8
  differenced arms judge the win (target: pg:recv 99 s + rb:wait 59 s at np=8 collapse
  toward one skew absorption); watch the reflux/gather/seam shadows for the downstream
  effect. Step-1's per-box asserts are bypassed on the chunked path by construction - the
  XA identity is the replacement invariant.

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
