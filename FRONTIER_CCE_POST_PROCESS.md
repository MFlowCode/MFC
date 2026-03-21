# Frontier CCE Post-Process Segfault in 3D Cylindrical Tests

## Status

**Open issue.** The Frontier CCE CI jobs are `continue-on-error: true` so this does
not block merges, but the underlying bug should be investigated.

## Summary

Post-process segfaults on Frontier when compiled with the Cray CCE compiler for
3D Cylindrical coordinate tests. The same tests pass with gfortran (Ubuntu, macOS,
Phoenix) and with Cray on Frontier AMD. Only Cray CCE on Frontier triggers the crash.

## Affected Tests

All 5 are `3D -> Cylindrical -> *` tests:

| UUID     | Trace                                          |
|----------|-------------------------------------------------|
| 301B9153 | 3D -> Cylindrical -> Viscous -> weno_Re_flux   |
| 939D6718 | 3D -> Cylindrical -> Viscous -> weno_Re_flux (variant) |
| 09623DE3 | 3D -> Cylindrical -> Viscous -> weno_Re_flux (variant) |
| 07C33719 | 3D -> Cylindrical -> (variant)                  |
| 128954AD | 3D -> Cylindrical -> cfl_adap_dt (variant)      |

## Error

```
Post-processing a 29x29x29 case on 1 rank(s)
[  0%]  Saving        1 of 2 @ t_step =        0
srun: error: frontier10243: task 0: Segmentation fault (core dumped)
```

Post-process starts, allocates memory, begins saving timestep 0, and segfaults.
The simulation itself completes successfully for all these tests.

## Key Facts

1. **Compiler-specific**: Only fails with Cray CCE on Frontier. Passes with:
   - gfortran on Ubuntu (GitHub CI)
   - gfortran on macOS (GitHub CI)
   - gfortran on Phoenix (self-hosted CI)
   - Cray on Frontier AMD (`frontier_amd` template)

2. **Reproducible**: Rerunning the same Frontier CPU job produces the same 5 failures.

3. **Only in post_process**: The simulation (same binary, same grid, same compiler)
   completes successfully. Only the post_process binary crashes.

4. **Only 3D Cylindrical**: Other 3D tests (Cartesian geometry) pass post_process on
   Frontier CCE. The distinguishing parameters are `cyl_coord = T`, `grid_geometry = 3`,
   `bc_y%beg = -14` (axis BC), and `m = n = p = 29`.

5. **Pre-existing but silent on master**: On master, the `--test-all` path does not check
   `cmd.returncode` after running post_process, so post_process crashes are silently
   ignored and the test is reported as passing. PR #1331 added a returncode check
   (line 509-511 of `toolchain/mfc/test/test.py`) which surfaces the crash.

6. **Frontier is non-blocking**: The CI workflow has `continue-on-error: true` for
   Frontier CCE jobs (comment in `.github/workflows/test.yml` says "Frontier CCE
   compiler is periodically broken by toolchain updates").

## Relevant Code Paths

### Post-process entry
- `src/post_process/p_main.fpp` — main loop iterates over save files
- Line 37: `n_save = int(t_stop/t_save) + 1` — computes expected number of saves
- The crash happens during the save loop, likely in data I/O or variable conversion

### Cylindrical coordinate handling in post_process
- `src/post_process/m_global_parameters.fpp:979-984` — sets `grid_geometry = 3` for
  3D cylindrical (`cyl_coord .and. p > 0`)
- `src/post_process/m_global_parameters.fpp:770` — adds extra stress index for
  cylindrical: `if (cyl_coord) stress_idx%end = stress_idx%end + 1`

### Boundary conditions (related fix in PR #1332)
- `src/common/m_boundary_common.fpp` — the QBMM `present(pb_in)` fix in PR #1332
  addresses a different post_process crash (absent optional arguments), but the
  Cylindrical segfault is a separate issue. The `present()` guards are already in
  the threading branch and don't resolve the CCE Cylindrical crash.

### Files that handle cylindrical coordinates in post_process
- `src/post_process/m_global_parameters.fpp` (grid_geometry, stress_idx)
- `src/post_process/m_start_up.fpp` (cyl_coord in namelist)
- `src/post_process/m_mpi_proxy.fpp` (MPI broadcast of cyl_coord)
- `src/common/m_boundary_common.fpp` (axis BC handling, BC_AXIS = -14)

## Investigation Leads

1. **Array sizing with cylindrical stress indices**: The extra `stress_idx%end + 1`
   for cylindrical coordinates changes `sys_size`, which affects array allocations
   throughout post_process. A mismatch between how simulation and post_process compute
   `sys_size` for cylindrical cases could cause out-of-bounds access.

2. **Cray CCE optimization**: CCE may optimize differently from gfortran, exposing
   undefined behavior (uninitialized variables, out-of-bounds access) that gfortran
   happens to handle safely. Running with CCE debug flags (`-G0 -O0 -Rbcps`) might
   give a more informative error.

3. **Axis boundary condition**: The `bc_y%beg = -14` (BC_AXIS) triggers special
   handling in `m_boundary_common.fpp`. The axis BC in post_process for 3D cylindrical
   may have indexing issues specific to how CCE lays out memory.

4. **Grid size 29x29x29**: This specific grid is used by Cylindrical tests (not the
   base 24x24x24). The odd grid size combined with cylindrical geometry might trigger
   an edge case.

## How to Reproduce

On Frontier with CCE:
```bash
source ./mfc.sh load -c f -m c
./mfc.sh build -t pre_process simulation post_process -j 8
./mfc.sh test -a -j 1 -o 301B9153 -- -c frontier
```

Or run the full Cylindrical test set:
```bash
./mfc.sh test -a -j 1 -o "3D Cylindrical" -- -c frontier
```

## Context: PR #1331 (Test Speedup)

This issue was discovered during PR #1331 which speeds up the test suite 2.3x.
The PR added a returncode check for the `--test-all` post_process pass (previously
missing on master). This surfaced multiple pre-existing post_process crashes:

1. **QBMM absent optional args** — fixed in PR #1332 (`present(pb_in) .and. present(mv_in)`)
2. **Adaptive-dt missing restart files** — fixed in PR #1331 by clamping `t_stop = t_save`
3. **Frontier CCE Cylindrical segfault** — this issue (unfixed, non-blocking)

The returncode check is correct and should be kept. The Frontier CCE segfault is a
pre-existing bug that was always crashing silently.
