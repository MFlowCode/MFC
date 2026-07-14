# AMR on AMD flang (OpenMP target offload) — root cause and fix

**Status: RESOLVED** (2026-07-14, verified on AMD MI250X / gfx90a, ROCm flang, OpenMP offload).
The block-structured AMR GPU path now runs correctly on the AMD flang OpenMP-target-offload
backend, alongside the pre-existing NVHPC (OpenACC / OpenMP) and Cray CCE (OpenMP) support.

This note supersedes the original hypothesis (which pinned the failure on the empty AMD
`defaultmap` in `omp_macros.fpp`). Live reproduction on an MI250X **refuted** that: see below.

## Symptom (AMD `gpu-omp`, np ≥ 2, 2D+)
Four AMR tests failed only on the AMD offload lane, passing on every other backend:
- **Self-abort** — `244B1E42` (2D single-level subcycle np=2), `ADA042A2` (2D multi-level static
  np=2): NaN appears at a seam-region cell after a few steps and the NaN guard calls `MPI_Abort`.
- **Tolerance drift** — 2D dynamic-regrid cases (`660FFBFE` IGR, `B7704247` stretched grid): the
  field diverges beyond golden tolerance.

## Root cause (confirmed by isolation)
The empty AMD `defaultmap` and the Cray-only `ACC_SETUP_SFs`/`ACC_TEARDOWN_SFs` are **not** the
cause. Isolation on the MI250X shows the failure is **2D-and-higher, np ≥ 2 only**:
- `21C71558` (1D static np=1) — passes.
- `5EFB3277` (1D dynamic regrid np=2) — passes.
- `244B1E42` (2D subcycle np=2) — NaN → abort.

A backend-wide defaultmap/descriptor problem would break 1D and np=1 too; those pass. The one
code path that is 2D-and-higher-np≥2-specific is the **fine-fine seam halo**
`s_amr_fine_slice` (`src/simulation/m_amr.fpp`). It packed/unpacked on the host and moved the
near-seam slab with a `target update` map clause of a **strided array section of a doubly-nested
derived-type member**:

```
$:GPU_UPDATE(host='[amr_slots(slot)%q_cons(i)%sf(dlo:dhi, 0:fm(2), 0:fm(3))]')
```

flang miscomputes the offset/stride of that strided section (seam dim `d < num_dims`) in the map
clause. In 1D the section is contiguous, so flang gets it right — which is exactly why 1D np≥2
passes and 2D+ np≥2 does not. `LIBOMPTARGET_INFO=-1` on a failing run showed **no** mapping error
(no present-table miss): the update silently moves the wrong bytes, corrupting the seam ghosts → NaN.

The base-grid MPI halo (`s_mpi_sendrecv_variables_buffers`) never hits this: it **device-packs**
into a contiguous buffer and only ever transfers the contiguous buffer.

## Fix
Rewrite `s_amr_fine_slice` to mirror the base-grid halo — pack/unpack on the **device** straight
into the contiguous buffer, and move only that contiguous buffer host↔device (per-kernel
`copyin`/`copyout`, no strided map-clause section). Two consequences that the live run surfaced and
that the fix accounts for:
1. **`buf` must be mapped.** It is not `declare target`, and AMD's `default='present'` expands to an
   empty defaultmap, so the pack/unpack kernels map it explicitly (`copyout` on pack, `copyin` on
   unpack). Without this the kernel faults on a null `buf` device pointer.
2. **`amr_slots` is not `GPU_DECLARE`d.** Indexing the module array `amr_slots(slot)%q_cons(i)%sf`
   *inside a kernel* dereferences a null outer descriptor. The slot's `q_cons` is passed in as a
   `scalar_field` array argument and the kernel indexes the dummy — the same pattern every other AMR
   kernel already uses (e.g. `s_amr_fill_fine_ghosts`).

The change is numerically identical (same pack/unpack ordering and casts), so goldens are unchanged
on all backends; only the transfer mechanism differs.

## Verification (AMD MI250X, OpenMP offload, np=2)
All four previously-failing cases pass, no regression on the passing cases:
`244B1E42`, `ADA042A2`, `660FFBFE`, `B7704247` (fixed); `21C71558`, `5EFB3277` (still pass).

## Second, separate bug: AMR+IB fine-marker swap/restore hangs the offload runtime

Found while sweeping the suite. **Symptom:** AMR+IB tests hang (99% host CPU, GPU idle, no timestep
output) inside `s_ibm_swap_to_fine` / `s_ibm_restore_from_fine` (`m_ibm.fpp`), the AMR-fine IB
ghost-point park/pull kernels. `rocgdb` on the live hang shows the AMD offload runtime busy-looping in
a recursive `targetDataBegin → targetDataMapper → targetDataBegin` (cycling through `free`,
`SourceInfo::getSubstring`, `std::string::find`).

**Root cause:** those kernels access the device-resident `allocatable` derived-type arrays
`ghost_points` / `gp_park` with only `default='present'` — which on AMD emits **no** defaultmap, so
flang generates a map *entry* for them and lowers it to a per-element custom mapper (the same amdflang
failure the code already dodged for the whole-array `ghost_points` `GPU_UPDATE`). The runtime then
busy-loops processing that entry. Both an implicit `tofrom` map and an explicit `map(present,alloc:…)`
still create the entry and still hang; only **no entry at all** avoids it.

**Fix:** add an AMD-only `defaultmap(present:allocatable)` (via `extraOmpArgs`, gated on
`MFC_COMPILER == "LLVMFlang"`) to the four swap/restore kernels — asserts them present with no map
entry, so no mapper is generated. CCE already gets this via its `default='present'`; NVHPC/CPU
unchanged. `05A8C23C` is the only *multi-level* AMR+IB test (largest `gp_park`), which is why it
surfaced first; the fix is on the shared kernels so it covers all AMR+IB cases.

**Not minimally reproducible:** a standalone kernel implicitly mapping the same flat derived-type array
completes instantly even at MFC's element count and present-table size — the per-element-mapper path
only triggers in MFC's full offload context. Worth a flang/ROCm report with the `rocgdb` backtraces.

**Verification:** all AMR+IB tests pass on the MI250X: `05A8C23C` (multi-level, was hanging),
`2854A102`, `F980C769`, `7FC2F9F8`, `13945217`, `43AF9F25`.

## Reproduce / debug env
```
source ./mfc.sh load -c amdfund -m g   # must be chained (&&) with the command below in ONE shell
export OMP_TARGET_OFFLOAD=MANDATORY     # abort if offload cannot reach a GPU
export OFFLOAD_TRACK_NUM_KERNEL_LAUNCH_TRACES=8   # name the faulting kernel on a GPU memory fault
export LIBOMPTARGET_INFO=-1             # full host<->device mapping trace (very verbose)
./mfc.sh test --only 244B1E42 -- -b mpirun -n 2
```
For a *hang* (no error text), attach to the spinning rank to see where the runtime is stuck:
```
rocgdb -p $(pgrep -n simulation) -batch -ex 'bt 12'
```
