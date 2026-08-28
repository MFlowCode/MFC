# T1/S4: plan-based exchange and distributed block metadata — v2, post-review

v1 of this design was independently audited by four reviewers (code-truth vs source, MPI transport,
GPU/compiler portability, adversarial correctness). Every blocker they found is folded in below,
with the finding named where it changed the design. **v1's payoff mechanism was partly wrong and its
increment scoping was wrong; the architecture survived.** This version is the implementation
contract.

## What is wrong today, stated structurally

Every AMR data exchange is driven **per box**, inside loops every rank walks over **all global
boxes**. Two consequences:

1. **O(boxes) rendezvous and posting skew.** Measured: `rb:wait` 583 ms x 123, `mg:wait` 3227 ms x
   16, plus two barriers (`rb:xchg`, `rb:flush`) absorbing the skew.
2. **O(boxes) metadata per rank.** 21+ arrays at `amr_max_blocks` extent on every rank, plus
   O(boxes) stack automatics in the regrid path (`old_ilo`/`old_ext`/..., `t_box` lists) and
   O(boxes x ranks) request arrays (`rq(old_np*num_procs)`). At 10^6 boxes the tags alone break:
   box-id tags exceed Cray MPICH's `MPI_TAG_UB` (~2^21). Fatal at scale independent of the tax.

AMReX fills a level with one `FillPatch`; Chombo with one `copyTo`; SAMRAI with one
`RefineSchedule`. The gap is granularity. The same refactor removes both consequences.

## The exchange-family inventory — complete this time

v1 claimed six families / 26 call sites; the audit counted 17 call sites in the named families and
found the rest in families v1 never named. The full inventory:

| # | family | per | notes |
|---|---|---|---|
| F1 | level-1 coarse gather (`s_amr_gather_coarse_patch`) | box | regrid AND per-stage fill |
| F2 | level>=2 parent gather (`s_amr_gather_from_parent`) | box | source is the freshly built parent |
| F3 | **qbmm pb/mv twin** (`s_amr_gather_coarse_patch_pbmv`) | box | LIVE at np>1 (single-level qbmm passes the checker); blocking `MPI_SEND` (never got the R1 fix); shares tag `amr_cur` with F1 under a non-overtaking lockstep contract |
| F4 | migration (`s_amr_regrid_stash_migrate`) | old block | already posts-then-waits (see mechanism) |
| F5 | reflux faces (`s_amr_p2p_reflux_faces`) + `s_amr_p2p_freg_to_parent` | block | pure `wp`, no precision crossing |
| F6 | seam halo (`s_amr_fine_fine_halo`) | pair | **also called from the L0 tile advance** — v1's "does not touch L0" was false |
| F7 | **restrict** (`s_restrict_fine_to_coarse`, `s_amr_restrict_to_parent`, `s_amr_scatter_pbmv`) | box, per STEP | reverse slot order: finest level first; v1 named it and then implemented it in no increment |

The subcycle path (`amr_subcycle`) reaches F1/F2/F3/F5/F7 through its own call shapes (two parent
snapshots per child on the same tag, ordered by non-overtaking; forced `s_amr_gather_send_flush`
sites). **Scope decision, now explicit: the subcycle call sites are NOT converted in v2.** They keep
the per-box path behind their existing gates; each converted site asserts its `amr_subcycle`
handling. Conversion is a follow-on increment (I8) after the lockstep path is proven.

### Call-site inventory, verified against source 2026-08-21 (pre-I1 sweep)

A full sweep of every AMR p2p MPI call found **42 sites** and nine facts the family table above
misses. The ones that change I1/I2:

1. **Ten call sites in five `s_l0_*` tile-routing routines sit outside the seven families**
   (`s_l0_fill_tiles_from_coarse`, `s_l0_scatter_tiles_to_coarse`, `s_l0_add_reflux_to_tiles`,
   `s_l0_restrict_to_tiles`, `s_l0_migrate_tile` — all blocking pairs, tags `k`/`4300`/`4400+k`).
   **Out of scope for conversion pending the D-l0 deletion decision** (`amr_endstate.md` sec. 8):
   the machinery measures 28-35% cost for 0.4% recovery, so converting it first would be waste.
   The I1 validator instruments them read-only (they are inert at the `l0_ntile=0` default).
2. **F1 and F2 share ONE request pool and ONE drain** (`amr_gsnd_req`/`amr_gsnd_pool`, cap 64,
   `s_amr_gather_send_flush`); a reserve for either family can force-drain the other's sends.
   The two families must be converted or fenced TOGETHER at each wave boundary — a per-family
   `amr_gsnd_n == 0` assert is meaningless while the pool is shared.
3. **F2's receive is a blocking `MPI_RECV`** (m_amr.fpp:1388) against pooled nonblocking sends,
   and its send site is fypp-instantiated twice (`_cons`/`_stor` — the subcycle calls both per
   child). The doc's per-box count undercounts, and F3's blocking-send defect has a twin on F2's
   recv side.
4. **F5 uses literal tags 2-7 (reflux faces) and 42-47 (freg), which numerically collide with the
   `amr_cur` tag space used by F1/F2/F3/F7** at small block counts. Phase separation is the only
   thing preventing mispairing today — one more reason the runtime tag bases (above) must land
   before any family is converted.
5. **F7 is three routines with three different blocking disciplines** (ISEND+WAITALL/blocking-RECV;
   fully blocking pair; ISEND+WAITALL/blocking-RECV): the "restrict" row is not one shape.
6. **F6 is blocking `MPI_SENDRECV` only** (tags 4200/4201): the validator cannot count
   posts-vs-drains there and must instrument the SENDRECV calls directly.
7. F5's send-side `reqs` allocation reuses `nreq` with two meanings (rank count at allocation,
   request count at the drain) — correct today, but a naive posts-vs-drains check will trip on it.
8. Non-AMR p2p the validator's instrumentation must NOT capture: `m_ibm.fpp` force halo
   (MPI_PACKED), `m_mpi_proxy.fpp` Lagrangian particle exchange, `m_start_up.fpp` IB
   neighbor-table build.

## The mechanism, corrected (MPI review B1/B2; code-truth M2; convergent)

v1 claimed the waits were late-sender under rendezvous with no async progress. **Wrong for this
transport, and refuted by our own code:**

- Open MPI 4.1.8 vader uses **CMA single-copy** (confirmed on the production node:
  `single_copy_mechanism=cma`, `ptrace_scope=0`). The RGET rendezvous is **receiver-driven**: once
  an ISEND is posted, the owner in `WAITALL` pulls the payload itself via `process_vm_readv`. A
  sender computing elsewhere stalls nothing that is already posted.
- **The in-tree control:** migration already posts all IRECVs, packs, posts all ISENDs, then one
  `WAITALL` — the exact v1-proposed structure — and still measures `mg:wait` = 3.2 s/rebuild.

The recoverable cost decomposes as: (a) **posting-order / head-of-line skew** — a send not yet
POSTED because its rank is still walking earlier boxes (batching removes this); (b) **pack/arrival
skew** — serial host packs, per-slot `GPU_UPDATE` staging, uneven ownership work (batching does NOT
remove this; pack parallelisation and staging restructure do); (c) **node-aggregate bandwidth** —
each byte crosses host DRAM ~5x (D2H stage, wp pack, CMA copy, stp unpack, H2D): 8 ranks x ~1.1 GiB
x 5 against ~100-200 GB/s shared is a floor of **0.3-0.6 s per rebuild per node**, not v1's 130 ms.

Pending measurement: a CMA-off control run (two-copy vader, where sender progress genuinely gates)
is in flight; its result calibrates the split between (a) and (b).

**Payoff, re-derived with a floor:** ceiling ~20% of wall (regrid pool with `mg:wait` mostly
excluded, plus the per-step F1/F6 families at 13.7%, at partial recovery); **floor 8-12%** if
pack/arrival skew dominates everywhere as it provably does in migration. Tax 11.03x -> **9.0-9.9x
floor, ~8x ceiling** for T1 alone. The program's case does not rest on the tax number: the scale
argument (tags, O(boxes) metadata, O(boxes x ranks) requests) is unconditional.

## The abstraction (revised per GPU review B1/B2/M3)

**No derived types.** The host plan is SoA flat arrays of intrinsics (avoids the CCE module-scope
derived-type descriptor bug class already worked around at `m_amr.fpp:641`); the device form is the
same arrays:

```fortran
! per (family, level): module-scope, GPU_DECLARE(create=...), allocate-max-once, 1.25x growth,
! refilled at plan build via contiguous-prefix GPU_UPDATE(device='[pl_lo(:,1:n)]') updates.
integer              :: pl_npeer, pl_nxs, pl_nxr
integer, allocatable :: pl_peer(:), pl_soff(:), pl_scnt(:), pl_roff(:), pl_rcnt(:)
integer, allocatable :: pl_loc(:)          ! LOCAL dense slot, resolved at build (amr_loc_of is
                                           ! host-only; kernel-side blk->loc translation is
                                           ! impossible). Couples plan validity to slot
                                           ! reconciliation - see the epoch rule.
integer, allocatable :: pl_lo(:,:), pl_hi(:,:), pl_off(:,:)   ! (3, nx)
integer, allocatable :: pl_coff(:)         ! exclusive prefix of CELL counts per transfer,
                                           ! per-peer-sliced: the unit the pack kernel's flat
                                           ! index runs over (sys_size factor NOT included)
integer(8)           :: pl_epoch = -1
```

**The pack/unpack kernel form is mandated, not suggested** (GPU review B1): a gang loop over
transfers with an inner vector loop is **silently serial on CCE and AMD flang** (`OMP_LOOP` expands
empty there). The portable form — already shipping in this codebase at `m_amr.fpp:3403` and `4167`
— is one flattened index over the peer's concatenated cells, decoded by **binary search over
`pl_coff`**, `collapse=2` with the `sys_size` loop. Per-family Fypp instantiation (the `GSFX`
idiom) supplies the array and the precision conversion. All of it lives in a new `m_amr_plan.fpp`
(compile-time isolation from the 7k-line `m_amr.fpp`).

**Wire buffers are persistent module arrays** (GPU review M1): `GPU_DECLARE`d, allocated once at
high-water, drained by contiguous-prefix `GPU_UPDATE` per peer. Never per-launch `copyout` of pool
slices (re-imports the 2.00-copies-per-entity map tax), never strided-section updates (AMD flang
corrupts non-contiguous `target update` — documented three times in `m_amr.fpp`). Unpack is a
device kernel wherever the destination is device-resident; per-family residency:

| family | source | destination |
|---|---|---|
| per-step fill (F1/F2/F3) | device | device |
| regrid gather (F1/F2/F3) | **host** (`q_cons_base` host-current) | patch storage |
| migration (F4) | host stash | host stash, then device push (see the fixed bug) |
| reflux/freg (F5) | device registers (`wp`) | device registers |
| restrict (F7) | device | device |

Precision: the wire is always `real(wp)`/`mpi_p`; stp<->wp conversion is a per-family template
parameter. **F5 is `wp` end-to-end with no conversion — the one family a shared hard-wired
conversion would corrupt under `--mixed`, invisibly to default-precision goldens.**

## Execution: per-(family, level) waves — the central correctness rule

Three reviewers independently converged on v1's worst flaw: "exchange everything, then prolong
everything" is **impossible**, because a level-l block's exchange source is *produced* by the
level-(l-1) prolong in the same rebuild ("re-prolongs from its (freshly-built, parents-first)
parent", `m_amr_regrid.fpp:1428`). The rule:

> A level-l exchange wave may start only after the level-(l-1) prolong + overlap-copy **and its
> device push** have completed. (The F2 pack is a device kernel; prolong writes the parent on the
> host; batching the `PH_RBPUSH` pushes to the end would feed the pack stale device data even with
> correct level staging.)

So: `for lev = 1 .. num_levels: build/execute plan(family, lev); prolong(lev); push(lev)`. Restrict
(F7) is the mirror image: **finest level first**, reverse waves. The per-step stage *fill* is the
one place cross-level batching IS legal (children are inset by `amr_cpat_mar`, so fill gathers read
only parent interior stage-entry cells) — that asymmetry is deliberate and must not be "unified".

Two-pass state rule (adversarial M6): `amr_cpat_off` and the working mirrors are written by the
gather and consumed by prolong. Pass 2 **recomputes both per box** — never inherits pass 1's frame.
Invisible on single-block cases, where the frames coincide.

Order of operations within a wave: post all IRECVs FIRST, then pack, then ISENDs, then one WAITALL
(with real statuses + `MPI_Get_count` under `MFC_DEBUG`, not `MPI_STATUSES_IGNORE`), then unpack.
Self-transfers (`peer == proc_rank`) are device copy kernels, present in the no-MPI build.

## Staleness: the epoch, not the dirty flag (adversarial B2; GPU B2)

v1 folded `amr_seam_pairs_dirty` into the stamp. **That flag is consumed** — cleared by whichever of
five lazy cache-rebuild triggers fires first — and ownership changes with NO regrid exist
(`s_l0_rebalance` migrates tile ownership mid-run; restart reassigns owners). A cached plan can see
"clean" and execute with the old owner map: a hang under clean tags, silent corruption under
colliding ones.

Rule: a monotone `amr_mesh_epoch` (integer(8), module scope), incremented at every site that sets
the dirty flag (`m_amr_regrid.fpp:1271`, `m_amr_restart.fpp:450`, `m_amr.fpp:6478`), at every real
regrid (after the `same` early-out), and at every slot reconciliation (plans bake `pl_loc`, so
recycled local indices invalidate them). Plans compare epochs; the boolean remains for the seam
caches only.

## Tags (MPI M1/M5; adversarial M5)

Live tag spaces today: `amr_cur` in [1..amr_max_blocks] for SIX logical transfers, migration
[1..old_np], reflux 2-7, freg 42-47, seam 4200/4201, L0 move 4300 — already numerically colliding,
safe only via phase separation and non-overtaking. Rules:

- Family tag bases derived at runtime: `base_f = amr_max_blocks + 100*f` (never literals; asserted
  below `MPI_TAG_UB` at init — Cray MPICH's is ~2^21, not INT_MAX).
- The epoch folded into the tag: `tag = base_f + mod(pl_epoch, K)` — a skipped epoch then mismatches
  loudly instead of pairing silently.
- `@:ASSERT(amr_gsnd_n == 0)` at every plan-exchange entry (the deferred pool legitimately holds
  level>=2 sends tagged `amr_cur` until the rebuild flush).
- Chunk any per-peer message above ~256 MB (buffer footprint + completion granularity; the int32
  count limit is 16 GiB and is not the reason).

## The validator (I1) — hardened (MPI M3; adversarial M4)

The v1 validator (set comparison + per-peer byte counts) cannot see: same-size transposition (most
blocks are exactly slot-sized, so swapped xfers have equal bytes), cross-rank builder asymmetry,
order/multiplicity semantics, dropped self-transfers, or which family a message belongs to (three
families share tag `amr_cur`, disambiguated only by call-site order). Requirements:

1. **Instrument the real MPI call sites** — log what is actually sent, partitioned by call site,
   never re-derived from the same metadata the builder reads (that validates the builder against
   itself).
2. Compare **ordered multisets per (peer, family)**, not sets.
3. `MFC_DEBUG` per-xfer identity: header words `(blk, lo, hi, family)` prepended to each transfer,
   verified at unpack.
4. **Destination-coverage tiling assert**: local copies plus received transfers exactly tile each
   destination patch, no gaps, no overlaps — the only check that sees self-transfer bugs, and it
   runs at np=1 for serial CI.
5. Cross-rank plan-hash `MPI_Allreduce` (debug builds): sender-side and receiver-side plans agree.
6. Builders read ONLY the replicated `*_all` arrays — never `amr_isect_lo/hi` or `amr_cpat_off`
   (empty on non-owners; the exact trap the parent-gather comments warn about). Assert it.
7. Explicit per-increment family scope, so the F3 twin cannot silently keep running per-box inside
   a converted loop.

Verified property worth one assert: old blocks ARE pairwise disjoint (cluster partition + merge
threshold + IB overlap-merge), so per-peer unpack reordering is safe for migration/overlap
destinations. Enforce with `@:ASSERT` after `shape_boxes`, don't inherit it as folklore.

## STATUS (verified against commits and source, 2026-08-27)

**Landed:** I0, I1a, I1b, I2a, I3, I4a, I4b, I5.
**Outstanding:** I2b, I5b (~250 LOC), I6 (~200), I7 (~600), I8 (unpriced).

Verified against the code, not inferred: **19 of 41 AMR p2p call sites still tag per box** (22 use
plan tags `tq`). F1 retains an unconverted path that passes the block index `amr_cur` as the MPI tag,
and migration (F4) passes the column index. The `.not. amr_subcycle` assert in the stage-fill wave
confirms the subcycle deferral recorded below is still in force.

**CONTRADICTION TO RESOLVE.** `m_amr.fpp` states the `amr_max_blocks` term leaves the tag base "with
the last per-box family (increment **I7**)". That cannot be right as written: I7's own boundary below
says "any family left per-box (subcycle) keeps its tables", and subcycle conversion is **I8**. So the
tag space -- and with it the ~28k-rank W5 wall -- does not clear until **I8**, not I7. The source
comment has been corrected; this note records why.

**Relation to S3 (W4).** S3.1 deleted the level-1 tag ALLGATHERV. The clustering tag union is
explicitly OUT of scope for I7 ("tag-union/clustering stay global (that is S3)") and that boundary
still holds -- S3.2/S3.3 own it.

## Increments, re-staged and re-priced

| # | content | LOC | gate |
|---|---|---|---|
| **I0** | prep, no plan code: `amr_mesh_epoch`; tag-base constants + init assert; `amr_gsnd_n==0` asserts; **the migration-stash device-push fix** (found live by this review: the receive-unpack path lacked the push its sibling has, so a mid-rebuild store grow clobbered migrated data) + a ppn=2 churn+growth regression case; box-disjointness assert | ~120 | AMR subset + the new case |
| I1 | validator: call-site instrumentation across ALL families incl. F3, F5-freg, F7, plus the six checks above. **Split 2026-08-21: I1a = `m_amr_xchg_audit` site registry (30 ids over the 35 physical sites; fypp twins share ids), always-cheap per-site msg/word/tag-range counters at every AMR p2p call, per-family global send==recv conservation asserts at finalize, the stash-only-replica reconcile assert, and the `[amr-xa]` report. I1b = MFC_DEBUG per-xfer identity headers (blk, lo, hi, family verified at unpack) + the destination-coverage tiling assert, gated by a SEEDED-BUG tripwire test.** | ~500 | I1a: goldens + subset green with `[amr-xa]` totals sane at ppn=2 AND ppn=4, no behaviour change. I1b: the seeded bug is CAUGHT. |
| I2 | F1+F3 level-1 gathers via plans (both twins together — converting one breaks their tag-order contract), per-owned-box patch storage + `amr_cpat_off` threading through the prolong/fill chain | ~600 | message count; per-xfer identity; bitwise goldens |
| I3 | F2 as per-level waves with the device-push rule | ~250 | dynamic-multilevel ppn=2 golden, bitwise |
| I4 | migration: parallelise the serial host pack; per-peer aggregation; right-size `spack`/`rq` (kills two O(boxes) allocations). **Decision made now, not during: whole-block sends preserved in I4 so the bytes gate stays exact; overlap clipping is I4b with a recomputed expected-bytes gate** | ~300 | bytes exact vs expectation; pack time |
| I5 | F5 reflux+freg (wp, no conversion) + F6 seam including the L0-advance call site | ~300 | bitwise; seam validated under L0 coexist |
| I5b | F7 restrict as finest-first waves | ~250 | bitwise; reverse-order assert |
| I6 | plan caching on `amr_mesh_epoch`; per-step F1/F2/F3 fills through cached plans | ~200 | plan-build count == epoch increments |
| I7 | S4: distributed builder; shrink global arrays **whose consumers are all converted**; convert the O(boxes) stack automatics. Boundary stated: tag-union/clustering stay global (that is S3); any family left per-box (subcycle) keeps its tables | ~600 | `[amr-scale]` per-rank bytes vs size |
| I8 | subcycle call-site conversion (two-snapshot ordering via distinct tag bases) | later | subcycle goldens |

~3100 LOC total (v1 said 1900 — the delta is the pbmv twin, restrict, the patch-storage
restructure, and the validator hardening; better to know now).

### I1b implementation binding (2026-08-23, line numbers at commit 0b36c148)

Scope: **I1b-gather** — headers on the trio I2 converts (F1/F2/F3), per the validator's
own explicit-family-scope principle; remaining families get headers with their conversion
increments (F5/F6 with I5, F7 with I5b). Priced against the int=20 ladder: the gather
family is 23% of np=8 steady wall scaling 3.84x/doubling — I2 is the program's
highest-value increment and this is its gate.

Mechanics (the trick that keeps it ~100 LOC): `m_amr_xchg_audit` exports
`XA_NH` (= 8 under `MFC_DEBUG`, else 0) plus `s_xa_hdr_pack(buf, fam, blk, bl, bh)` /
`s_xa_hdr_check(buf, fam, blk, bl, bh)` (integers encoded as `real(wp)`, exact ≤ 2^53).
Every wire count/offset gains `+ XA_NH` UNCONDITIONALLY (zero when the header is off, so
production arithmetic is untouched and no call site needs an `#ifdef`); pack/verify calls
are `if (XA_NH > 0)` — dead-code-eliminated. Anchors are the existing `s_xa_rec` calls
(1:1 with wire ops by I1a's construction):
- plan sizes: `amr_gpl_sz`/`amr_gpl_psz` in `s_amr_build_gather_plan` gain +XA_NH per
  message (recv posts at m_amr.fpp:1063/1073 then need no size edits);
- chunked F1: pack/send 1140-1156 (header before the pack loop, `boxsz + XA_NH` on the
  wire), pool unpack 1266 (verify then offset); chunked F2: send via
  `s_amr_gather_from_parent_field_cons` (1876), device-unpack 1219 (host-verify the first
  XA_NH words before `s_amr_unpack_parent_patch_device`, then pass the offset slice);
- per-step F1: 1482 (recv)/1570 (send) + twins at 3403/3416 and 4307/4320 (fypp
  instantiations share XA ids — enumerate by grepping `XA_F1_`); F3: 1677/1766;
  per-step F2 blocking recv: 1911 (xbuf);
- pool reserves: `s_amr_gsnd_reserve(maxsz + XA_NH)` at 1136 and the `need` sum in
  `s_amr_gather_chunk_post` (~1040).
Header check failure → `call s_mpi_abort` with family/blk/expected-vs-got. Tiling assert
(the other I1b half): at each owner's unpack completion, assert the union of contributor
slabs plus the own-slice tiles the patch exactly (count cells, compare to patch volume) —
lives beside the existing gather asserts. Gate: the seeded-bug counterfactual (swap two
plan source entries locally → headers must abort; revert seed) + 75-golden subset +
`[amr-xa]` totals unchanged (headers are size-invisible to the F-family word counters:
count payload words only, i.e. record `cnt` not `cnt + XA_NH`).

### I2a implementation binding (2026-08-23)

I2 is staged. **I2a (landed)**: `s_amr_stage_fill_wave` (m_amr.fpp) converts the
non-subcycle per-stage LEVEL-1 fill (F1 + the F3 twin together) to one wave per RK stage —
SoA transfer lists built per wave from the replicated caches (both sides enumerate boxes
ascending with per-rank running offsets, so the per-(peer, family) wire layout — the
ascending-box concatenation of [XA_NH header | slab] — agrees with no metadata exchange),
recvs-then-packs-then-sends-then-one-WAITALL per the order-of-operations rule, tags
`amr_tag_base(1|3) + mod(amr_mesh_epoch, 100)`, wave audit sites XA_F1W_*/XA_F3W_* folding
into families F1/F3 (payload-words-only, so `[amr-xa]` family words stay exactly comparable
to the per-box baseline — the landed gate). Consume is box-major through the single
`amr_cg`/`amr_cpat_off`, reusing the per-box device kernels on contiguous pool slices:
**per-owned-box patch storage turned out unnecessary for the wave** — it only buys
cross-box batched unpack kernels, deferred to **I2b (contingent)** on the post-I2a budget
showing launch/map overhead rather than wait left in the gather share. The mandated
binary-search pack/unpack kernel form applies to I2b's kernels when/if they exist; I2a
adds zero device kernels. Plan caching on the epoch remains I6; the regrid-path F1
(chunked) and init/static/restart/subcycle sites keep the per-box path unchanged.

### I3 implementation binding (2026-08-23)

`s_amr_parent_fill_wave(lev)`: the per-step F2 gather as one wave per level per stage,
levels ascending from the driver (`do ilev = 2, amr_num_levels`). Each split child is
exactly one (parent-owner -> child-owner) transfer; both sides derive the pair list from
`f_amr_parent_block` + `s_amr_parent_foot` + `amr_block_owner` ONLY (the per-owner
mirrors lag a generation — the map's asymmetry finding). Tags `amr_tag_base(2)` + epoch
fold; audit sites XA_F2W_* (family F2, payload-words-only). Reuses the I2a wave's
scratch (q-side arrays only; the waves never overlap in time). The pack kernel reads
module `amr_cpat_off`, so the pack loop sets the CHILD's frame per transfer; consume
recomputes per box. `s_amr_fine_stage_fill` deleted with the conversion (no caller
remained). Restart gotcha fixed with it: `amr_num_levels` was regrid-only; the per-level
driver needs it truthful after restart too (set in `s_read_amr_restart`). Remaining
per-box F2: regrid chunked (converts with I6/I7 work if ever), subcycle (I8),
init/static.

### I5-F6 implementation binding (2026-08-23)

The seam wave lives INSIDE `s_amr_fine_fine_halo` (all four call sites covered,
subcycle shape-preserved). Plan = the replicated `amr_seam_pairs` list itself, walked
twice (sends, then recvs); each cross-rank pair is one transfer each way per owner;
per-peer aggregation on tag `amr_tag_base(6)` + epoch fold; audit sites XA_F6W_*
(payload-words-only keeps family F6 words exact vs the SENDRECV baseline — the landed
gate). Reuses the fill waves' fw scratch, with `amr_fw_spo/rpo` repurposed to carry
per-transfer `cnt` (the F6 payload is not derivable from bl/bh alone at consume). The
shared `amr_seambuf_x/y` and their tile-grow reconciliation are deleted. Header
convention: [site, sending slot, (d, pack dlo, pack dhi), (cnt, 0, 0)] — the receiver
derives the peer's pack bounds from the same replicated metadata. F5 remains per-box;
its conversion notes: faces recv directly into the mapped `freg` host mirror (zero-copy)
so headers must be DEBUG-ONLY COMPANION MESSAGES, not prefixes; `s_amr_reg_reserve`
must hoist ahead of any wave posting into `freg` (the apply can reallocate the
registers); the owner-side multicast membership is `cand ∩ f_amr_reflux_participates`
vs the receiver's bare predicate — a wave must reproduce the conjunction exactly.

- **No existing test runs np>2.** Every plan degenerates to <=1 remote peer: multi-peer slicing,
  peer ordering, and multi-contributor assembly are structurally unexercised. One ppn=4 dynamic
  regrid case is mandatory before I2 lands.
- A dynamic-regrid pbmv (qbmm non-polytropic) case at np=2 (the rebuild-path twin is uncovered).
- Bitwise golden diff mode: several AMR goldens carry `override_tol` up to 1e-5, so "tolerance
  zero" must be an explicit bitwise comparison, or drift hides inside existing tolerances.
- An `MFC_DEBUG` artificially-low chunk-threshold test (the >256 MB path never triggers at golden
  sizes).
- A periodic-seam np>1 case (a builder that forgets wrap pairs deadlocks only there).

## Merge gates

Four CI compilers + AMD flang; `GPU_*` macros only; `@:ALLOCATE` pairing; full 706-test suite at
merge (AMR subset per iteration); one increment = one commit, each independently green; bitwise
comparison per increment as above. Every wall-time claim measured against the 5.3% run-to-run floor
with >=3 repeats, or judged on exact byte/count gates instead.

### I5-F5 implementation binding (2026-08-23)

`s_amr_reflux_faces_wave` (level-1 walk, ascending slots) and `s_amr_freg_wave`
(level>=2 cowner->powner pairs) replace the per-box F5 exchanges in the lock-step
driver; the subcycle path keeps its per-box form (`do_xchg` on
`s_amr_reflux_to_parent`). Plan = the replicated owner/region tables: face-wave
participants are cand from `s_amr_ranks_overlapping` on region+-1 INTERSECTED with
`f_amr_reflux_participates(r)` (the conjunction the per-box form applied — a wave must
reproduce it exactly or a rank posts a recv no one sends). Receives are ZERO-COPY into
the `freg(d)%%lo/hi(:,:,:,slot)` host mirrors (each box owns a register slot — there is
no pool and MUST be none, or the zero-copy design is lost); owner D2H before the ISENDs,
receivers H2D after the WAITALL. Tags: `amr_tag_base(5) + mod(epoch, 50)`, freg wave at
`+50`. Identity headers are DEBUG-ONLY COMPANION 8-word messages carried in
`amr_fw_sq/rq` (never `s_xa_rec`'d, so [amr-xa] family words stay comparable to the
per-box baseline); `amr_fw_rblk` tracks the expected box order at consume.
`s_amr_reg_reserve(amr_num_blocks)` runs FIRST in both waves — the apply can reallocate
`freg`, and a recv posted into a stale mirror is a silent wrong-memory write.

### Ring-clip-on-waves implementation binding (2026-08-24)

The stepfill clip (`amr_stepfill_ring_clip.md` — the dead-byte proof survived the
revert) is applied inside the stage-fill wave's TWO plan walks: after each pair's
`s_amr_box_isect`, the slab is fed through `s_amr_shell_clip` against the box's shell
(`s_amr_shell_slabs` of the padded patch minus the open core [region_lo+1,
region_hi-1]; collapsed dims pass 0 so they never cut). Up to 6 sub-slabs replace the
one transfer; both walks derive the identical list from replicated metadata, so the
no-metadata-exchange property is untouched, and pack/unpack/consume — already generic
over (bl, bh) — need no change. Messages stay at the per-peer count; only payload
words drop (F1 -61% on the S0 probe). The consume's own-box copy becomes shell-only
(`s_amr_gather_own_shell_device`), which is what arms the debug NaN-poison gate
(`s_amr_poison_patch_device` floods the patch first; any consumer read of an unshipped
cell aborts within a step). The pbmv gather (F3) keeps its full-box wire contract, so
`do_pbmv` runs take the unclipped single-slab path — extending the clip to qbmm needs
its own dead-byte analysis. All four primitives are lifted verbatim from the reverted
implementation (archived: amr-bench/notes/ringclip_original_m_amr.fpp.txt).
