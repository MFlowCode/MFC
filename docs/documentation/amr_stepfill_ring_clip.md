# Ring-clipping the runtime fill gathers

> **STATUS 2026-08-22: IMPLEMENTED, PROVEN CORRECT, and PARKED (reverted) on an amdflang
> whole-image codegen regression — NOT on any defect of this design.** The implementation
> (commits dc6d4129 + bd85c792 + a7970743, reverted immediately after) passed the full
> correctness bar: output bit-identity at np=4 and np=8, zero transport-assert trips,
> wire words −64 to −72%, gather family −33% at np=8. But adding its 7 target regions
> made amdflang's whole-image device link deterministically regenerate UNTOUCHED kernels
> with 2.4–4.5x worse ISA (weno 5x scratch, riemann VGPR→AGPR flip, LDS 2048→2560
> image-wide). Verdict, evidence, and the re-landing trigger: `amr_action_plan.md`
> "2026-08-22 (final)".

**Target:** the step-fill gather family = 14.3% of np=8 wall (185.8 s post-step-2 on
k004-001), whose words are **71.2% provably dead** (`[amr-cov]`, logs/p1gate-0822_1259:
48.2e9 of 67.7e9 words at np=8). The ghost-fill kernel reads only `floor(f/rr) +- 1`
(`s_amr_fill_fine_ghosts_*`), so of the full patch `[region-mar, region+mar]` only the
**hollow shell** — the patch minus the open core `[region_lo+1, region_hi-1]` — is ever
consumed. Ship the shell, not the patch.

## REVIEWED 2026-08-22: two adversarial reviewers (MPI/consumer lens + state/lifetime
## lens). The shell claim itself SURVIVED both line-by-line audits — every runtime
## consumer of amr_cg is a shell reader, the open core is exact with ZERO margin (reach
## is 1 on every branch incl. the multi-fluid closure; do NOT enlarge the shell), and the
## host/device asymmetry structurally insulates the rebuild's full-patch readers (runtime
## writes device-only; rebuild reads host truth after full-box host writes). But the
## original spec was defective in one place and validation-blind in four. This document
## is the AMENDED form; implement only this.

## Scope (amended per finding F1, both reviewers)

**ALL runtime gathers** — every `s_amr_gather_coarse_patch(..., .true.)` site and every
level>=2 `to_host = .false.` path — **subcycle included**. The original "lockstep only"
scope was unimplementable: the subcycle sites pass the identical flags (level-1
`m_amr.fpp:5527/5533`; level>=2 `5741-5751` vs the lockstep `1791/1796` — no argument
distinguishes them), and the expansion is proven safe because the subcycle consumers
(`fill_gsta/gstb`) are generated from the same Fypp body as `fill_cons` — the shell-read
proof covers them identically. Consequences: `s_amr_cov_note_fill` must also be called at
the subcycle fill sites (so accounting matches behavior), and the validation sweep must
include the subcycle goldens. Routing fact the original doc got wrong: the lockstep
level>=2 consumer is `fill_cons` (5387); `gsta/gstb` are subcycle-only.

Still out of scope, with named coherence walls (F11): the rebuild/init gathers
(`pull_host=.false.` / `to_host=.true.` — their consumers prolong the FULL patch), the
pbmv twin, and the four full-array `GPU_UPDATE`s of amr_cg (m_amr.fpp:1268, 1512, 1955,
1990) — they are the load-bearing host/device coherence walls; the clipped runtime path
must never add a host write of amr_cg, and no one may "optimize" those updates as part of
this work.

## Design

1. **One shared slab helper, used by every side** (sender pack, owner recv-size, owner
   unpack, own-box copy, np=1 shortcut, and the XA/cov shadow accounting):
   `(patch, core) -> <= 6 DISJOINT slabs`, fixed order (x-low/x-high full-transverse,
   then y-low/high restricted to the core's x-interval, then z-low/high restricted in x
   and y). Conventions (F6): collapsed dims define `core_d = [0,0]` (full interval, no
   transverse slabs); width<=2 regions have an empty core (shell = whole patch, legal);
   width-1 overlap resolved by clamped subtraction; always `@:ASSERT` slab disjointness
   AND `sum(slab words) == patch - core`. The transverse restriction is to the CORE
   interval, not the closed region (the closed-region variant double-covers face planes).
   The shell is the 3D box difference — never an independent per-face ring (transverse
   coordinates of shell reads go arbitrarily deep per-dim).
2. **Clipped exchange, same message set** (F8/Finding 2): iterate today's
   `amr_ovl_gather` contributor list UNCHANGED on both sides; a contributor whose overlap
   lies inside the core sends a ZERO-word message (legal MPI) — never a one-sided
   emptiness skip (that is an owner-side WAITALL hang; reachable at exascale operating
   points where rank slabs are narrower than region-2). One buffer, one message, today's
   tag per contributor; pack slices concatenated in slab order.
3. **amr_cg stays full-size**; only shell cells are written on the runtime path; the core
   is stale BY DESIGN — acceptable only because of the walls in Scope.
4. **Fused kernels** (F10): own-copy, pack, and unpack iterate the <= 6 slab
   intersections through the concatenated-flat-index single-kernel idiom the fill kernel
   already uses — one launch per message, never one per slab (the family's measured tax
   is launch-path serialization).
5. **Frame assert** (F7): the clip makes sender-frame (replicated region/foot) vs
   consumer-frame (`amr_isect_lo`) agreement load-bearing; `@:ASSERT(amr_isect_lo ==
   region_lo)` (level-1) / `== foot lo` (level>=2) on the owner in every clipped path.
6. **Precision/conversion points preserved verbatim** (F9): own-box/co-located copies
   stay direct stp:=stp; wire pack stays real(stp->wp); wire unpack stays real(wp->stp);
   keep the `(i, g3, g2, g1)` per-slab linear order and the fixed slab order both sides.

## Validation (amended per F2-F5 — the original plan was partially circular/blind)

- **Primary, non-circular: output BIT-IDENTITY** vs HEAD at np=1, 4, 8 (the step-2
  method) + goldens 67/67.
- **The MFC_DEBUG poison arm is MANDATORY and covers the WHOLE patch, not the core**
  (F3): a debug-gated device kernel writes quiet NaN (stp) over the box's entire patch
  extent at the top of every clipped gather, before any shell write — sites: the runtime
  branch of `s_amr_gather_coarse_patch` (covering the np=1 shortcut and the np>1 owner
  path), `s_amr_copy_parent_patch_*` (to_host=.false.), and
  `s_amr_unpack_parent_patch_device` (to_host=.false.). Any read of an unshipped cell —
  core OR a missed shell slab (the hole class the core-only poison is blind to) — NaNs
  the ghost fill and trips goldens/ICFL within a step. Sweep: the 67-case AMR subset
  under MFC_DEBUG at np=1 and the suite's ppn=2 cases (contributor-slab arithmetic exists
  only at np>1).
- **Transport asserts via MPI_GET_COUNT** (F4): the clipped recvs keep statuses (today
  they discard them) and assert `MPI_GET_COUNT == predicted live words` per message —
  sender-vs-receiver recomputation of the same replicated function is a tautology and is
  NOT the check. XA cannot see short sends (it records posted sizes).
- **Word accounting is shadow-counted at the send/recv sites** (F2): the `[amr-cov]`
  counter is the SAME arithmetic the slabs derive from (circular) and counts non-wire
  words (own-box, np=1) that XA never sees, and XA ids blend runtime/rebuild/subcycle.
  The word prediction therefore asserts TRANSPORT only; shell correctness rests on
  bit-identity + the poison arm.
- **Boundary-block coverage (F5): CLOSED BY EVIDENCE, no new test.** The attempted
  boundary golden (static block at 0, ppn=2) aborts in the RUNTIME checker: "amr block
  must lie at least buff_size cells inside the domain boundaries" — the reviewer's
  premise (only `>= 0` enforced) missed this simulation-side @:PROHIBIT. Since
  `amr_cpat_mar = ceil(buff_size/rr) + 1 <= buff_size` for buff_size >= 2, no VALID
  configuration produces a patch crossing the domain boundary. The implementation
  asserts the premise instead of handling the case: `@:ASSERT(amr_cpat_mar <=
  buff_size)` at shell-helper init.

## Expected payoff, stated as a bound

stepfill is bytes + skew (CMA-off control: bandwidth is real but not all of it). Bytes
cut ~71% on a 14.3%-of-wall family + pack/unpack kernels iterate ~71% fewer cells + ~71%
less PCIe. Ceiling ~10% of np=8 wall; realistic 4-8%. Judged by differenced same-node
arms against the 4% run-to-run variance (single runs are not evidence), and the
np-doubling ratio reported against the SOTA bar (AMReX 1.20x/1.15x, amrexs0-0822_1330).
