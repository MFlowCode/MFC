# T1/S4: plan-based exchange and distributed block metadata

This is the design for the only change that can move MFC's AMR tax materially. The incremental
items (T0) were measured and closed: each is a 4-8%% candidate against a **5.3%% run-to-run noise
floor**, so none is defensible. See `amr_action_plan.md` for that evidence.

## What is wrong today, stated structurally

Every AMR data exchange is driven **per box**, inside a loop that every rank walks over **all global
boxes**:

```
do k = 1, nboxes                       ! on EVERY rank, for EVERY box
    amr_cur = f_l0_slot(k)
    ... one collective rendezvous for box k ...
    if (.not. amr_rank_owns_block) cycle
end do
```

Two consequences, and they are the whole problem:

1. **Rendezvous count is O(boxes).** Rank A cannot pass box *i* even when box *i+1*'s partner is
   ready. Measured: `rb:wait` 583 ms x 123 calls, `mg:wait` 3227 ms x 16, and two downstream
   barriers (`rb:xchg`, `rb:flush`) that exist only to absorb the skew this creates.
2. **Metadata is O(boxes) per rank.** **21 distinct arrays** are allocated at `amr_max_blocks`
   extent on every rank (`amr_block_owner`, `amr_block_level`, the region/isect tables,
   `amr_owns_all`, `amr_slot_live`, `amr_slots`, `amr_loc_of`, the cached peer lists
   `amr_ovl_gather`/`amr_ovl_scatter`, `amr_seam_pairs`, ...), and **21 driver loops** walk the
   global box list. At 10^6 boxes this is fatal independently of the tax. (An earlier revision said
   "eight arrays"; the grep missed the peer caches and second allocation sites.)

AMReX fills an entire level with one `FillPatch`; Chombo with one `copyTo` + `Copier`; SAMRAI with
one `RefineSchedule`. **The gap is granularity, not algorithm** — and the same refactor fixes both
consequences, which is why T1 and S4 are one program rather than two.

## The abstraction

One plan type serves all six exchange families (`s_amr_gather_coarse_patch`,
`s_amr_gather_from_parent`, `s_amr_regrid_stash_migrate`, `s_amr_p2p_reflux_faces`,
`s_amr_fine_fine_halo`, `s_amr_restrict_to_parent` — 26 call sites in total):

```fortran
type :: t_amr_xfer                 !< one rectangular-subarray transfer (strided in memory, not contiguous)
    integer :: blk                 !< block this transfer belongs to
    integer :: lo(3), hi(3)        !< region in the SOURCE frame
    integer :: off(3)              !< offset into the DESTINATION frame
end type

type :: t_amr_plan
    integer                       :: npeer
    integer,          allocatable :: peer(:)              !< rank ids, ascending
    integer,          allocatable :: soff(:), scnt(:)     !< per-peer slice into xs
    type(t_amr_xfer), allocatable :: xs(:)                !< sends,   grouped by peer
    integer,          allocatable :: roff(:), rcnt(:)     !< per-peer slice into xr
    type(t_amr_xfer), allocatable :: xr(:)                !< receives, grouped by peer
    integer(8)                    :: stamp = -1           !< mesh generation this was built for
end type
```

Execution is generic; **receives are posted before anything else** (large messages use the
rendezvous protocol, and an unposted receive stalls the matching send):

```
s_amr_plan_exchange(plan, src, dst)
  1. size one buffer per peer from scnt/rcnt
  2. post ONE MPI_IRECV per peer                      <- FIRST, before pack
  3. ONE pack kernel/loop per peer over its xs slice
  4. post ONE MPI_ISEND per peer
  5. ONE MPI_WAITALL for the whole level
  6. ONE unpack kernel/loop per peer over its xr slice
```

**Pack location follows data residency, and residency differs per family** (`m_amr.fpp:830` states
this explicitly for the gather alone). The design's earlier "one device pack kernel" was wrong as a
universal statement:

| family | source residency | today's pack |
|---|---|---|
| per-step stage fill | device-current | device kernel |
| regrid coarse gather | **host-current** ("valid ghosts from the exchange at the top of s_amr_regrid") | host |
| migration (`rg:move`) | stash, host path | **serial host loop** — parallelising it is a free extra win |

MPI buffers are host memory throughout (`amr_gsnd_pool` is a plain host allocatable), i.e. the
transport is host-staged, not GPU-aware — so per-family the "bytes" cost includes the PCIe crossing,
and that stays true under the plan.

Message count per level falls from `sum over boxes of (peers per box)` to `<= num_procs`. On the
matched case that is **~450-900 messages per regrid -> <= 7**.

## Why this is also S4

A plan references only the blocks this rank actually sends or receives. Once every consumer reads
its plan instead of the global tables, the global arrays are needed **only by the plan builder** —
and the builder is then the single place to make distributed. That ordering matters: distributing
the metadata first would break 21 loops with nothing to replace them.

**The boundary, stated honestly:** I7 removes the global-table *consumers*. The tag union and
clustering (`s_amr_union_gtag`, `s_amr_cluster`) remain global — every rank still builds the whole
box list — and making *those* local is S3 (SAMRAI-style local clustering + boundary
reconciliation), a separate program. I7 makes S3 possible; it does not include it. One more piece of
S-track evidence found in this review: migration staging is allocated `spack(maxcnt, old_np)` —
**~11 GiB virtual per rank at the matched point, O(global box count) by construction**. Untouched
pages keep it harmless today, but at 10^6 boxes the allocate itself fails. Plan buffers are sized by
actual transfers and fix this for free.

## Correctness strategy — the part that makes this landable

**Increment 1 builds the plan and validates it against today's behaviour without using it.** The
validator asserts that the plan's derived `(peer, blk, lo, hi)` set is exactly the message set the
current per-box path produces. A mismatch is a silent wrong answer, not a crash, so this net comes
before any data moves through the new path.

Two invariants to assert, never assume:

- **Layout agreement.** Sender and receiver must derive identical per-peer sequences. Deterministic
  iteration over the same cached lists gives this; assert matching byte counts per peer before
  unpacking.
- **`amr_cur` is implicit.** Per-block routines read module state, not arguments. Plan execution
  must not leave `amr_cur` pointing at the wrong block for a later consumer.
- **Frames differ per family** (`amr_implementation.md` sec. 2): level-1 gathers are in L0 cell
  indices, level>=2 in the PARENT-fine frame, migration in old-block-local indices. `t_amr_xfer`
  carries lo/hi/off as bare integers with nothing marking the frame — each builder must document
  and assert its frame, or a cross-family copy-paste becomes a silent wrong index.

**Every increment is judged on bytes first, wall second.** Migration bytes were measured identical
to the digit across runs (`10748501376`), so aggregation can be proven exactly — message counts and
byte totals from the `[amr-scale]`/`[amr-mig]` counters — against a 5.3%% wall-noise floor that would
otherwise swamp the result.

## Increments

| # | content | LOC | verified by |
|---|---|---|---|
| I1 | plan type, builder for the level-1 coarse gather, **validator only** | ~400 | validator + AMR goldens; no behaviour change |
| I2 | route the level-1 gather through plan execution | ~200 | PRIMARY: message count and bytes (exact); wall only via >=3 repeats against the 5.3%% floor |
|    | *memory note: "one WAITALL per level" needs per-owned-box destination patches alive at once: ~28 owned x ~15 MiB = ~420 MiB/rank at the matched point. Acceptable there, but I2 must carry the chunk fallback from `amr_regrid_gather_batching.md` for cases where it is not.* | | |
| I3 | level >= 2 parent gather (`pg:all`) | ~150 | goldens; `pg:all` delta |
| I4 | migration (`rg:move`) | ~200 | **bytes exact**; `mg:wait` delta |
| I5 | reflux + seam | ~250 | goldens; `rf:wait`, `seam` deltas |
| I6 | cache plans across steps, invalidate at regrid | ~100 | plan-build count == regrid count |
| I7 | **S4**: distributed builder; shrink the global arrays | ~600 | `[amr-scale]` per-rank bytes vs problem size |

~1900 LOC total. I1-I2 alone are a complete, landable, independently valuable change.

## The MPI mechanism — why the waits are as large as they are, and what that predicts

`rb:wait` is 583 ms per call for a ~5-15 MB-per-contributor exchange. That is not bandwidth and not
queue depth; it is **late-sender under the rendezvous protocol with no asynchronous progress**:

- Slices this size use rendezvous, so the transfer starts only when the *sender* re-enters MPI.
- Since fix R1, contributors post `ISEND` into the pool and return to compute — and we measured the
  host spending ~40% of its CPU in Open MPI's progress engine (`amr-idle-is-mpi-progress`). The
  owner's per-box `WAITALL` therefore stalls until each contributor happens to make an MPI call.
- Plan exchange fixes this **because every rank posts everything and then sits in one `WAITALL`**,
  which is itself the progress engine. Nobody is late because nobody is elsewhere. The benefit is
  the post-then-wait-together structure, not the message count per se.

Sanity bound, not a promise: ~1.1 GiB moved per rank per rebuild (migration + gather) at a
host-staged ~8 GiB/s is **~130 ms**, against measured waits of ~4.5 s + 3.2 s per rebuild. An order
of magnitude of headroom exists **if** late-sender is the mechanism; if it is not, the floor is the
bytes and the win is small. I2's gate must distinguish these outcomes, which is why wall time alone
is not the gate.

Implementation guards from this analysis: per-peer messages can reach hundreds of MB — chunk any
message above ~256 MB (also guards the 2^31 element count); fan-out 1.048 means per-peer buffers
duplicate ~5% of bytes vs today's shared-buffer sends — negligible, accepted; MPI's per-(peer, tag,
comm) ordering makes the one-message-per-peer scheme robust; the plan `stamp` must fold in
`amr_seam_pairs_dirty`, not just the regrid generation, because the peer caches it consumes are
invalidated on that flag.

**Residual after batching, so the payoff is not oversold:** the skew that `rb:xchg`/`rb:flush`
absorb comes partly from uneven *ownership work* (prolong `rb:ovl` 22 s, `rb:slot` 27 s per-rank
spread), which batching does not touch. The sinks will shrink, not vanish.

## Expected payoff — and the honest limit

Plan execution removes the per-box **rendezvous**, not the **bytes**. Two pools, from the measured
budget (caveat: `mg:wait` is from a different run than the others — the percentages mix two
denominators 3.2%% apart, so the sums are indicative, not a single-run measurement):

**Regrid-path waits and their skew sinks** (increments I1-I4):

| term | %% wall |
|---|---|
| `rb:wait` | 9.0 |
| `mg:wait` | 6.3 |
| `rf:wait` | 3.3 |
| `rb:xchg` (skew sink) | 5.8 |
| `rb:flush` (skew sink) | 1.2 |
| subtotal | **25.6** |

**Per-STEP per-box families** (I5-I6 — these are in the six-family scope and an earlier revision of
this table omitted them):

| term | %% wall | evidence it is per-box |
|---|---|---|
| `gather` (stage fill) | 10.1 | 52614 calls/rank = 219 boxes x 240 stage-calls |
| `seam` | 3.6 | serial cross-rank `MPI_SENDRECV` loop |
| subtotal | **13.7** |

At an honest 80%% recovery of both pools, wall x ~0.69: tax **11.03x -> ~7.6x**. Counting only the
regrid pool (if I5-I6 disappoint): **-> ~8.8x**.

**That is not AMReX parity, and this design should not be sold as reaching it.** Parity at 3.13x
needs infrastructure at 6.6%% of wall against today's 78.2%% — an 11.9x reduction — and it additionally
requires the physics term: our fine advance is **2.41x less efficient per cell** than our own uniform
arm (ghost recompute, small-block GPU utilisation, the `scalar_field` copy bridge). T1 does not touch
that; T2 (fused advance) does, and T1's descriptor arrays are the prerequisite for it.

On T2's pool, a correction to the 2.41x figure: **34%% of the physics term is the coarse advance
(59.9 s), which re-advances the whole domain including under the 70%%-refined region.** That
redundancy is algorithmic in non-subcycled AMR — a fused fine advance does not touch it. T2's
reachable pool is the fine-advance inefficiency only (~115 s against ~60 ideal).

Realistic staging, arithmetic shown rather than asserted: **T1 full scope -> ~7.6x; T1+T2 ->
~6.9x** (T2 closing half its reachable pool). An earlier revision said "T1+T2 -> ~5-6x"; that
number is not supported unless T2 recovers nearly all of its pool AND the residual regrid work
(`rb:slot`, `rb:ovl`, `rg:clus`, pack/unpack volume) also shrinks. The honest claim for this
program is "closes most of the gap"; parity additionally needs the migration-volume and
residual-infrastructure work.

## Every corner — the implementation constraints that decide success

These are the traps that survive a correct-looking design review. Each is grounded in something
already measured or already worked around in this codebase; none is hypothetical.

### C1 — the plan must be POD-flat on the device (the trap we measured ourselves)

The `t_amr_plan` sketch above has allocatable components. **Handing that to a device pack kernel is
the exact descriptor trap this campaign measured**: every mapped array entity costs exactly 2.00
copy operations on the OMP backend (`amr-mapped-entity-law`), pointer views into the store are not
attachable (measured, `m_amr.fpp` copy-bridge comment), and CCE leaves module-scope derived-type
allocatable descriptors uninitialized (already worked around once at `m_amr.fpp` block-allocate).

**Therefore: `t_amr_plan` is HOST-side bookkeeping only.** The device-facing form is flat integer
arrays — `plan_blk(:)`, `plan_lo(3,:)`, `plan_hi(3,:)`, `plan_off(3,:)`, `plan_buf_off(:)` — mapped
once at plan build. The pack kernel takes bare arrays and scalars, nothing derived-typed. Parthenon
reaches the same conclusion from the other direction (POD `BndInfo` descriptors).

### C2 — precision crossings are per-family, not generic

The store is `stp`; every existing MPI buffer is `wp` (`spack`, `amr_gsnd_pool`), with the
conversion done inside the pack loop (`real(amr_stor_st(...), wp)`). Under `--mixed` (`wp` double,
`stp` half) this conversion is load-bearing, not cosmetic. The generic `s_amr_plan_exchange` must
let each family supply its pack/unpack element conversion; a shared kernel that assumes one
precision is a silent `--mixed` corruption that CI's default-precision goldens cannot see.

### C3 — the mesh-generation stamp does not exist yet

Nothing in the code counts regrids. I6's cache-invalidation needs `amr_mesh_gen`, incremented
exactly where the box set actually changes (after the `same` early-out, not per `s_amr_regrid`
call — 24 of 40 regrids change nothing). Getting this wrong silently reuses a stale plan: wrong
answers, no crash.

### C4 — tags currently carry meaning; per-peer messages need a tag discipline

Today's tags encode identity (box id `amr_cur`, old-block id `kk`) and that is what keeps concurrent
per-box messages from cross-matching. Plan exchange sends ONE message per (peer, family), so each
family needs a reserved tag base, and two families' exchanges must never be in flight concurrently
with the same (peer, tag) — the deferred-send pool (`amr_gsnd_*`) can still have level-2 sends
pending when a plan exchange starts. Either drain it first (the existing flush rule) or separate the
tag spaces; assert, don't assume.

### C5 — np=1, no-MPI, and self-transfers

The current per-box code has explicit `num_procs == 1` shortcuts and local-copy paths. In a plan,
self-transfers (`peer == proc_rank`) must bypass MPI as direct device copies, and the whole exchange
must compile and run in the no-MPI build (`#ifdef MFC_MPI` discipline: the plan builder and local
copies exist everywhere; only the ISEND/IRECV/WAITALL core is gated).

### C6 — the rebuild loop becomes two passes, and the stash constraint still binds

Whole-level exchange means restructuring `s_amr_regrid_rebuild_slots` from
"per box: gather -> prolong -> overlap-copy" into "exchange all patches, THEN per box: prolong ->
overlap-copy". The stash-read constraint (`amr_implementation.md` sec. 9.2: old local indices are
read throughout the rebuild) is unchanged — old slots still cannot be freed until the loop ends.
The per-box `s_set_amr_fine_geometry`/`amr_cur` sequencing survives in pass 2 unchanged.

### C7 — merge requirements, per this repo's contract

Four CI compilers + AMD flang; GPU code via `GPU_*` macros only; every `@:ALLOCATE` paired;
`wp`/`stp` discipline; format -> precheck -> build -> FULL test suite (not the AMR subset - that is
the iteration gate, not the merge gate); one logical change per commit, each increment
independently green. The bitwise expectation per increment: I1 changes no behaviour at all; I2-I6
deliver identical data in identical order, so goldens must pass at TOLERANCE ZERO drift - any diff
is a bug, never a regeneration.

## What this explicitly does not do

- It does not change numerics. Every increment is golden-verified.
- It does not touch the L0 tile path or the non-AMR solver.
- It does not require a new case parameter.
