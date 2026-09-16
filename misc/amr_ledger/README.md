# AMR performance ledger: what was learned

This directory holds the records of the block-structured AMR performance campaign (July to September
2026). It is kept out of the built documentation on purpose: `docs/documentation/amr.md` says what the
AMR does and `amr_implementation.md` says how it is built; this directory says how we found out, what we
tried, and what turned out to be false. Nothing here is needed to use or modify the code.

## Where things are

| file | what it is |
|---|---|
| `status.md` | the last status page: the three scorecard statements, their evidence, what blocks the pull request |
| `action_plan.md` | the raw notebook: 165 dated entries, append-only, with retractions inline; the evidence behind every claim below |
| `endstate.md` | the architecture the work converged on: four pillars, the weak-scaling invariants W1 to W8 |
| `plan_based_exchange.md` | the exchange-wave design (families, keyed tags, plans, the audit oracle) |
| `block_batching.md` | why the fine advance is batched over a flat store, and the per-block swap cost that forced it |
| `multilevel.md`, `fine_distribution.md`, `per_level_distribution.md` | nesting rules, ownership and the space-filling-curve distribution |
| `regrid_gather_batching.md`, `stepfill_ring_clip.md` | two exchange-volume reductions and their measurements |
| `tax_review.md`, `slowness_analysis.md` | the early phase budgets and the causal model they led to |

The measurement harness (scripts, decks, pre-registrations, logs) is a separate repository, `amr-bench`,
next to the MFC checkouts on the cluster where the work was done.

## The result

On a healthy node, MFC's AMR overhead is at parity with AMReX on the same problem. Measured as the excess
of the AMR step over the code's own uniform step at the same fine resolution (the quantity that isolates
the AMR machinery from the base scheme): MFC 0.42 s/step, AMReX 0.38 to 0.39, ratio 1.08x, with a
measurement floor of about 0.03. Per cell-update, MFC's AMR inflates its own solver by 1.78x and AMReX's
inflates its own by 1.86x. MFC's AMR step is slower than AMReX's in absolute terms because its base scheme
(WENO5 with mapping and monotonicity preservation, HLLC, the 5-equation model) costs 1.7x per cell, which
is a scheme choice and not an AMR term. Weak scaling from 8 to 32 GPUs is 1.28x per doubling against
AMReX's 1.23x on matched rungs, with the growth in the base grid's cross-node halo rather than in the
regrid family.

## What the overhead is made of

About a fifth of the remaining excess is fine-solver inflation (ghost planes per block, batch padding, the
per-launch dispatch floor), about half is MPI wait, and the rest is real AMR work (regrid, restriction,
gather and reflux packing, seams, migration). The MPI wait is not a transport problem: roughly fourteen
rendezvous per step, and the wait at each is set by the same two ranks being slowest in every segment
between them. Merging or reordering rendezvous moves where the light ranks idle and does not shorten the
step. The wall is the heaviest rank's own work plus its own transfer floors.

## The levers that moved the wall, and why

- **Fixing the block store's growth ratchet** (2.3x). The store grew through the host on every regrid and
  never shrank; capped, in-place re-densified growth with early free of consumed slots fixed both the
  time and the device-memory blow-up.
- **Replacing per-box rendezvous with wave exchanges** (17 to 22 %). Every exchange family became one
  posted wave per stage with keyed tags, instead of one blocking send per box.
- **Batching the fine advance over the flat store** (about 12x on the AMR-local kernels). A per-block
  advance costs about one monolithic step regardless of block size, because launches and mapped-array
  descriptors are per launch, not per cell.
- **Device-resident wire pools** (73 ms/step). Under RDMA-capable MPI the exchange pools live on the
  device and MPI sends them by device address, which removes a synchronous host copy per box.
- **Block-local work arrays in the kernels** (69 ms/step). On amdflang a `private` fixed-size array
  costs a descriptor copy per launch; an array declared in a `block` inside the loop body costs none.
  Reverted afterwards: NVHPC rejects a block inside a parallel region and CCE OpenACC loses the
  block-local arrays from its present table, and the block form was not wanted in the code base. The
  gain is real on amdflang only and is available to anyone who accepts a compiler-specific kernel form.
- **Clipping reconstruction to the block interior** (12 ms/step of RHS). Transverse ghost planes were
  reconstructed and never read.

## The levers that did not, and why

- **Load balancing on cell counts, then on measured time.** Cells were already balanced to 1 %. The
  imbalance is per-rank launch count and block-shape diversity; moving blocks relocates launches and
  multiplies them. A time-feedback balancer made the light ranks heavier by more than the heavy ranks got
  lighter.
- **Fewer, larger batches.** Halving the batched calls doubled the per-call fixed cost: the cost is per
  member (fine-field copies, captures, ghost fills inside the batch), not per call.
- **Equalising tile shapes after clustering.** Shape diversity on the level-2 boxes cannot be trimmed
  after clustering without uncovering tagged cells.
- **Asynchronous kernels.** A `nowait` target region encountered outside a parallel region runs inline on
  this OpenMP stack; it overlaps nothing.
- **Kernel bodies as device routines.** The inlined body inherits the routine's register budget and the
  kernel runs 30 to 47 % slower per launch; and a `declare target` routine reads the never-updated device
  copy of a host-only module scalar, which is a GPU-only NaN with byte-identical CPU goldens.
- **Reordering reflux posts and drains, deferring sends, merging rendezvous.** Sign unresolved or zero:
  see the wall-is-the-heaviest-rank finding above.
- **Subcycling.** Parity with lock-step at matched fidelity; the earlier 2.8x came from a broken control.
- **Lowering precision.** Not a lever. The comparison with AMReX is at double precision by decision.

## Things that were believed and were false

- "AMR launches 173x more kernels, so launch count is the cost." Launch count was flat in the rank count
  and only 18 % of device time; the host wait behind it was the per-launch descriptor and dispatch floor.
- "Regrid is half the wall." An operating-point artefact of the store ratchet; after the fix it is 6 %.
- "The excess is MPI bandwidth." It is wait, and the wait is skew.
- "A flat, then rising, scaling ladder means regrid is O(P)." Both regrid paths were O(P) for a while (a
  window all-gather and a level-1 tree of global reductions), and both were fixed; after that the growth
  is the base grid's halo, which AMR does not own.
- "The node does not matter." The excess ratio for the same binary spans 1.1x to 1.9x across nodes of the
  same partition, because MFC's host-bound overhead is node-sensitive and AMReX's device-bound step is
  not. Half the campaign's flat week was sick nodes, dead InfiniBand ports pinned by the site's UCX
  configuration, and thermal throttling of one GPU package under back-to-back load.
- "The metric is stable." The excess is a difference of two noisy walls; the uniform arm was 96 % of the
  variance until both arms were lengthened to the same window, and the short window was biased as well
  as noisy.

## Rules that the mistakes bought

1. A performance claim needs the control and the treatment on the same node in the same allocation, the
   parent commit as the control, at least three reps per arm, and a stall detector on every arm.
2. Pre-register the prediction and the falsifier before submitting; declare any change to the statistic
   before reading the data.
3. Price a lever against the heaviest rank's own serial chain before building it. A lever that passes
   its mechanism gate but removes work off that chain does not move the wall.
4. A kernel or exchange change ships with byte-identical goldens (contraction pinned) on both the CPU and
   the GPU-versus-GPU pair; the CPU gate cannot see a device-only path.
5. Behaviour changes ride behind a default-off flag until measured; a flag whose one correct value is
   known is then removed, not kept.
6. Never wait on CI to decide; run the relevant goldens locally. Never parse aliased `ls`. Never kill by
   process-name pattern. One tree per concurrent job.

## What is left, as decisions rather than work

The batched advance was extended to MHD, IGR, Lagrangian bubbles, hypoelasticity, chemistry, the
6-equation model and prescribed-motion bodies, and the per-block fine advance was then retired from
this branch together with the physics only it served (subcycling, stretched or cylindrical grids,
QBMM, Euler bubbles, phase change, moving particle clouds). That code is correct and golden-tested
but about 1.3x slower per step; it lives on the branch `amr-per-block` (this branch plus one revert)
for a follow-on pull request. Level-0 tiling is a tested feature with no production user. The branch
bundles the AMR with kernel restructurings and unrelated physics and needs splitting before review.
