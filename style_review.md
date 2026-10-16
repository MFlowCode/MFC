# Style review: `allMachProjectionCore` against the existing code

**What I compared:** branch `allMachProjectionCore` at `0f880484`, against its merge-base with `upstream/master` (`47568a0f`).
**What I checked against:** existing modules in `src/simulation`, and the coding standards in
`docs/documentation/contributing.md`.

Nothing has been changed. This file is untracked; delete it when you're done.

The most significant findings come first.

## 1. Direction handling uses a different paradigm

**Status:** refactored, in uncommitted changes. Directions now use `m_weno`-style Fypp loops, and `cfx/cfy/cfz` follow the
`u_face` face indexing. Not yet built or tested.

- **Existing code:** sweeps permute bounds with `is1/is2/is3` and `irx/iry/irz`, or expand per direction with
  `#:for NORM_DIR`.
- **New code:** carries a runtime direction `d`, unit offsets `(oj, ok, ol)` from `s_unit`, and direction-generic helpers
  (`f_width`, `f_along`, `f_last`, `f_face_dist`, `f_bc_code`, `f_cf`). The helpers branch on `d` inside device code
  (`m_projection_helpers.fpp:22-101`, `m_projection.fpp:142-172`).
- **Face indexing is inconsistent within the new code:**
  - `u_face(j)` is the face between cells `j` and `j+1` (`m_projection.fpp:34`).
  - `cfx(j)` is the face between `j-1` and `j` (`:42`).
  - The multigrid conductances follow the `cfx` convention.
- **Redundant alias:** `d = id` (`:307`, `:397`, `:445`).

## 2. Terse names, some of which clash with MFC conventions

**Status:** renames applied (uncommitted), except `mg_smode`. `pflux` became `face_flux`. `rkc1..rkc4` became `rkc(4)`; in
`s_mg_kstep` the dot products are in `dots(5, 1)`, and `s1..s5` remain only as kernel reduction accumulators.

- **`vf` is a scalar face velocity** (`m_projection.fpp:274`). Elsewhere in MFC, the `_vf` suffix means a vector field.
- **`iv` is a momentum index** (`:193`, `:571`). In other modules `iv` is the velocity `int_bounds_info`.
- **`ar`/`al` mean alpha_rho/alpha** (`:439`, `:567`). Existing code uses `alpha_rho_K`/`alpha_K` or `alpha_rho`/`alpha`, and
  `al`/`ar` reads as left/right.
- **Other opaque names:**
  - `fbx/fby/fbz`, `cfx`, `gnd`, `kap`, `pk`, `fm`, `a_up`, `ar_up`, `dd`, `du`, `vn`, `tl`, `ga`, `gf`, `rc2`
  - `sl1 ... po3` (removed by the direction refactor)
  - `mg_d/e/f/r`, `mg_o`, `mg_om`, `mg_kt`
  - the `crs_*` family, `dk`
- **Mismatched docstring:** the `f_mg_ix` doc says "cell (i, j, k)", but the arguments are `xi, yi, zi`
  (`m_multigrid.fpp:277-278`). The multigrid also loops over `ii/jj/kk` where the rest of the code uses `j/k/l`.

### Proposed renames

**Where the proposals come from:**
- Existing MFC names where an equivalent exists:
  - `alpha_rho`, `alpha`, `gamma`, `pi_inf` (`m_sim_helpers`)
  - `Re_L/Re_R`, `Re_size_loc1/2`, `alpha_rho_L` (Riemann solvers)
  - `dyn_pres`
  - `pres_out` (`bc_x%pres_out`)
- Otherwise, a spelled-out name for what the variable holds.

**What's left alone:**
- Loop indices `i, j, k, l`.
- Textbook CG scalars (`alpha`, `beta`, `rz`, `bnorm`, `rnorm`).
- Names that already read clearly (`u_face`, `u0_face`, `p_stage`, `p_step0`, `rhostar`, `tau`, `face_bounds`, `seam_lo`).

#### `m_projection` module variables

| Current | Holds | Proposed |
|---|---|---|
| `cfx, cfy, cfz` | face conductances of the pressure operator | `cond_x, cond_y, cond_z` |
| `fbx, fby, fbz` | boundary code of each domain face | `bc_face_x, bc_face_y, bc_face_z` |
| `pk`, `pk_sf` | solved pressure with ghosts (the "search direction" in its comment is stale) | `p_solved`, `p_solved_sf` |
| `pguess` | initial guess for the solve, then its solution | `p_guess` |
| `dcoef` | Helmholtz diagonal D_c of the pressure system | `p_diag` |
| `bvec` | right-hand side of the pressure system | `p_sys_rhs` |
| `rhs_p` | pressure transport rate (reused as `s_mass_lag` scratch) | `p_rate` |
| `pflux` | upwind pressure flux, then face acceleration in the correction | `face_flux` (or split into two arrays) |
| `kap` | curvature (1) and \|grad alpha\| weight (2) | `curv` |
| `gnd` | grad(alpha_1) and its magnitude | `grad_alpha` |
| `wb_st` | well-balanced surface tension flag | `well_balanced_st` |
| `pout`, `has_pout` | outlet pressure per direction | `pres_out`, `has_pres_out` |

#### `m_projection` locals

| Current | Where | Holds | Proposed |
|---|---|---|---|
| `vf` | `s_projection_rhs` | face velocity (clashes with the `_vf` suffix) | `u_f` |
| `a_up`, `ar_up` | `s_projection_rhs` | upwind alpha, alpha_rho | `alpha_up`, `alpha_rho_up` |
| `fm` | `s_projection_rhs` | summed mass flux | `mass_flux` |
| `up_l` | `s_projection_rhs` | flow is left to right | `upwind_l` |
| `ar`, `al` | `heat`, `apply` | alpha_rho, alpha (`al`/`ar` read as left/right) | `alpha_rho`, `alpha` |
| `al`, `ar` | `face_props` | left and right alpha | `alpha_L`, `alpha_R` |
| `re_l`, `re_r` | `face_props` | interface Reynolds numbers | `Re_L`, `Re_R` |
| `rs1`, `rs2` | `face_props` | `Re_size` copies | `Re_size_loc1`, `Re_size_loc2` |
| `gam`, `pinf` | `heat`, `apply` | mixture Gamma, pi_inf | `gamma`, `pi_inf` |
| `rc2` | `apply` | rho c^2 | `rho_c2` |
| `dv` | `apply`, `f_div_uf` | div(u_face) in a cell | `div_u` |
| `ke` | `apply` | kinetic energy | `dyn_pres` |
| `acc` | `apply` | body-force acceleration | `accel` |
| `ga` | `apply` | tau times body-force acceleration | `dvel_bf` |
| `gf` | `apply` | pressure correction of a face velocity | `dvel_p` |
| `fix` | `apply` | face keeps its Dirichlet velocity | `dirichlet` |
| `code` | `apply` | boundary code of the face | `bc_code` |
| `iv` | `init_faces`, `dirichlet_faces`, `apply` | momentum equation index (`iv` is an `int_bounds_info` elsewhere) | `mom_idx` |
| `ia` | `s_compute_curvature` | alpha_1 equation index | `adv_idx` |
| `kv` | `s_compute_curvature` | curvature | `curv_c` |
| `e`, `g` | `face_codes`, init | domain side, ghost-layer count | `side`, `n_ghost` |
| `its`, `res` | `s_pressure_solve` | solve iterations and residual | `iters`, `res` |
| `c`, `dd` | `s_mass_lag` | lagged flux divergence of the variable, and of the velocity | `lag_div_q`, `lag_div_u` |
| `du`, `vn` | `s_mass_lag` | face velocity change, corrected face velocity | `du_face`, `u_f` |
| `tl`, `pre` | `s_mass_lag` | kernel copies of `tau`, `mom` | `tau_loc`, `mom_loc` |
| `iq`, `q`, `i1`, `i2` | `s_mass_lag` | variable counter, equation index, counter range | `lag_var`, `q_idx`, `var_beg`, `var_end` |
| `nalpha` | `s_mass_lag` | volume fractions to correct | `n_alpha` |

#### Projection functions and helper arguments

| Current | Holds | Proposed |
|---|---|---|
| `f_div_uf` | div(u_face) in a cell | `f_div_u_face` |
| `f_nrm`, `nd` | interface normal component from `gnd` | `f_interface_normal`, `n_d` |
| `f_normal(gd, g)` | grad_d(alpha)/\|grad alpha\| | `f_normal(grad_d, grad_mag)` |
| `f_face_gf` | pressure correction of a face velocity | `f_face_dvel_p` |
| `f_cond(ra, rb, ...)`, `kf` | face conductance | `f_cond(rho_a, rho_b, ...)`, `cond` |
| `f_capillary_accel(ka, kb, wa, wb, ca, cb, ra, rb, dist)` | curvatures, weights, colors, densities | `(kappa_a, kappa_b, w_a, w_b, c_a, c_b, rho_a, rho_b, dist)` |
| `f_face_capillary(..., jn, kn, ln, ...)` | the neighbor cell | `j2, k2, l2` (or keep) |
| `s_bounds` | sets loop bounds | `s_set_bounds` |

#### `m_multigrid*` storage (`m_multigrid_levels`)

| Current | Holds | Proposed |
|---|---|---|
| `mg_d` | level diagonals | `mg_diag` |
| `mg_kx, mg_ky, mg_kz` | low-face conductances | `mg_cond_x, mg_cond_y, mg_cond_z` |
| `mg_e` | cycle correction, plus the PCG search direction at `mg_poff` | `mg_corr` |
| `mg_f` | level right-hand side | `mg_rhs` |
| `mg_r` | level residual | `mg_res` |
| `mg_o` | flat index of cell (0, 0, 0) per level | `mg_origin` |
| `mg_sy, mg_sz` | y and z strides per level | `mg_stride_y, mg_stride_z` |
| `mg_off` | level start offset | `mg_lev_off` |
| `mg_poff` | offset of the PCG direction block | `mg_dir_off` |
| `mg_gx, mg_gy, mg_gz` | ghost width per direction | `mg_ghost_x, mg_ghost_y, mg_ghost_z` |
| `mg_nlev` | level count | `mg_num_lev` |
| `mg_nbr` | neighbor rank per side | `mg_nbr_rank` |

#### `m_multigrid`

| Current | Holds | Proposed |
|---|---|---|
| `rs, zs, qs` | PCG residual, preconditioned residual, A times direction | `pcg_res, pcg_prec, pcg_Ad` |
| `mg_kr, mg_kc, mg_kv` | K-cycle rhs, first direction, its image | `mg_k_rhs, mg_k_dir, mg_k_Adir` |
| `mg_kfull`, `mg_kspace` | levels between K-cycle steps (default, this solve) | `mg_k_spacing_dflt`, `mg_k_spacing` |
| `mg_blev` | level this solve stops at | `mg_bottom_lev` |
| `mg_om`, `om` | coarse-correction scale | `mg_cc_scale`, `cc_scale` |
| `dg`, `nb` | row diagonal, conductance-weighted neighbor sum | `diag`, `nb_sum` |
| `sh` | offset to the PCG direction | `dir_off` |
| `f_mg_ix(lv, xi, yi, zi)`, `ii, jj, kk` | cell indices | `j, k, l` |
| `cnx, cny, cnz` | coarse sizes | `nx_c, ny_c, nz_c` |
| `fi, fj, fk` | fine index of a coarse face | `j_f, k_f, l_f` |
| `a0, a1` (`b`, `c`) | child range per direction | `j_lo, j_hi` (`k_*`, `l_*`) |
| `sd, skx, sky, skz` | summed diagonal and conductances | `sum_diag, sum_cond_x, sum_cond_y, sum_cond_z` |
| `s1 ... s5` (`s_mg_kstep`) | fused dot products | named after each product, e.g. `c1_Ac1, c1_r, r_r, Ac1_r, Ac1_Ac1` |
| `cc`, `cd`, `a1` | K-step coefficients | `coef_c1, coef_d, alpha1` |
| `sl`, `sg` | local and global reduction vectors | `sums_loc`, `sums_glb` |
| `dk` | bottom cells' diagonal and conductance rows | `cell_rows` |
| `kl`, `ku` | K-cycle vector bounds | `k_beg`, `k_end` |
| `cl` | kernel copy of `color` | `color_loc` |

#### `m_multigrid_exchange`

| Current | Holds | Proposed |
|---|---|---|
| `mg_smode` | side action, magic values 0/1/3 (2 at unpack) | `mg_side_mode`, with named constants |
| `mg_boff` | side offset in the send (1) and receive (2) buffer | `mg_buf_off` |
| `mg_nmsg`, `mg_nlist` | message count, rank of each message | `mg_num_msg`, `mg_msg_rank` |
| `mg_sbeg, mg_slen, mg_rbeg, mg_rlen` | message offsets and lengths | `mg_send_beg, mg_send_len, mg_recv_beg, mg_recv_len` |
| `mg_sbuf, mg_rbuf` | send and receive buffers | `mg_send_buf, mg_recv_buf` |
| `q`, `qq` | side 1..6 | `side`, `side_nb` |
| `md` | side mode in this phase | `mode` |
| `nn`, `st` | cells and stride normal to the side | `n_norm`, `stride_norm` |
| `n1, n2, s1, s2` | tangential sizes and strides | `n_t1, n_t2, stride_t1, stride_t2` |
| `ex, ey` | x and y extents with ghosts | `ext_x, ext_y` |
| `lay`, `gho` | interior layer and ghost-layer index | `layer`, `ghost` |
| `lvl`, `ph` | kernel copies of `lv`, `phase` | `lv_loc`, `phase_loc` |

#### `m_multigrid_bottom`

`crs_` (short for coarse) becomes `bot_`, matching the module name.

| Current | Holds | Proposed |
|---|---|---|
| `crs_n`, `crs_tot` | bottom cells; cells over all bottom levels | `bot_ncells`, `bot_ncells_all` |
| `crs_nl`, `crs_ln`, `crs_loff`, `crs_ld` | level count, cells, offset and rank-grid dimensions per level | `bot_nlev`, `bot_lev_ncells`, `bot_lev_off`, `bot_lev_dims` |
| `crs_cnt`, `crs_disp`, `crs_sz` | cells, displacement and block dimensions per rank | `bot_rank_ncells`, `bot_rank_disp`, `bot_rank_dims` |
| `crs_ad`, `crs_ak`, `crs_aj` | stencil diagonal, conductances, neighbor cells | `bot_diag`, `bot_cond`, `bot_nbr` |
| `crs_agg`, `crs_co` | parent cell on the next level, rank-grid coordinates | `bot_parent`, `bot_coords` |
| `crs_l` | dense Cholesky factor | `bot_chol` |
| `crs_x`, `crs_b` | solution, right-hand side | `bot_x`, `bot_rhs` |
| `crs_exact`, `crs_cg` | solve flags | `bot_exact`, `bot_cg` |
| `f_crs_nbr`, `f_crs_ax`, `s_crs_cycle`, `s_crs_cg`, `s_crs_gs` | routines | `f_bot_nbr`, `f_bot_ax`, `s_bot_cycle`, `s_bot_cg`, `s_bot_gs` |
| `pc`, `ic`, `nd` | rank-grid or cell coordinates, side direction | `coords`, `cell`, `dir` |
| `row`, `rows` | packed local and gathered stencil rows | `rows_loc`, `rows_glb` |
| `jr`, `ci`, `cj`, `o`, `nn` | Cholesky column, parent cells, level offset and size | `col`, `parent_i`, `parent_j`, `lev_off`, `lev_n` |

## 3. Multigrid storage is flat 1-D with manual index arithmetic

- **Flat arrays:** all levels live in 1-D arrays (`mg_d, mg_kx, mg_e, ...`) addressed through `f_mg_ix` and hand-computed
  strides (`m_multigrid_levels.fpp:20-31`). MFC otherwise uses multi-dimensional arrays or `scalar_field`.
- **Hand-coded exchange:** `s_mg_sides` builds strides and offsets by hand (`m_multigrid_exchange.fpp:253-284`).
- **Integers sent as reals:** `s_mg_bottom_build` packs integer neighbor indices into a real array so they can travel through
  one `MPI_ALLGATHERV`, then recovers them with `nint` (`m_multigrid_bottom.fpp:134-150`).
  **Kept:** separating them would cost a second `MPI_ALLGATHERV`. The level-count reductions in `m_multigrid_levels` now use
  the new `s_mpi_allreduce_integer_max` (same call count).

## 4. Shared public mutable state

- **Public by default:** `m_multigrid_levels` has no `private; public ::` line. Everything is public except `f_seam_rank`
  (`:18`), which inverts the module convention. Three modules read and write its `mg_*` arrays. This conflicts with "No new
  global state; private helpers stay inside their defining module" in `contributing.md`.
- **Solver diagnostics in global parameters:** `proj_pcg_iters` and `proj_pcg_res` are written by `m_projection` and read by
  `m_data_output` (`m_global_parameters.fpp:24-25`).

## 5. Semicolon-chained statements are about four times as dense

- **Density:** about 46 per 1000 lines in the new modules, against about 12 per 1000 in the existing simulation code.
- **Worst cases:**
  - A continuation line that starts with a semicolon: `& ; om = mg_om(lv)` (`m_multigrid.fpp:776-777`).
  - `r = mg_nbr(q); if (r < 0) return` (`m_multigrid_bottom.fpp:104`, `:110`).
  - One-line `case (1); f_side_size = ...` (`m_multigrid_exchange.fpp:236-240`).
  - Long assignment chains (`m_multigrid.fpp:365`, `:606`, `:619-620`; `m_multigrid_bottom.fpp:62`, `:66-67`).

## 6. Function declarations are inconsistent

- **Two forms are mixed:** 11 functions use the typed-prefix form (`real(wp) function f_cf`, `integer function f_bc_code`,
  `f_width`, `f_face_dist`, `f_last`, `f_along`, `f_wall_face`, `f_opposite`, `f_side_size`, `f_crs_nbr`, `f_crs_ax`); the
  rest use `result()`. Existing simulation code favours `result()` by about 25 to 6. `f_crs_nbr` and `f_crs_ax` use both
  forms at once (`integer function ... result(g)`).
- **`pure` is inconsistent:**
  - Not `pure`: `f_width`, `f_face_dist`, `f_cf`, `f_bc_code`, `f_div_uf`.
  - `pure`: their siblings (`f_last`, `f_along`, `f_cond`).
- **Blank lines are inconsistent:** usually a blank line separates the signature from `GPU_ROUTINE`, but not in `f_cf`
  (`m_projection.fpp:159-160`), `f_width` (`m_projection_helpers.fpp:72-73`), `f_face_dist` (`:88-89`), or `f_crs_ax`
  (`m_multigrid_bottom.fpp:284-285`).
- **Module `@brief` placement:** in `m_projection.fpp:5-8` it comes before the `#:include` lines. The other five new files,
  and most existing ones, put it directly above `module`.

## 7. Comments are long prose and carry history

- **Paragraph-length blocks that derive or justify the method:**
  - `s_mass_lag` (8 lines, `m_projection.fpp:749-756`)
  - `s_mg_set_operator` (5 lines, `m_multigrid.fpp:329-333`)
  - `s_mg_omega` (4-line header plus a 5-line block inside, `:429-432`, `:461-465`)
  - `m_multigrid_bottom` module variables (7 lines, `:24-30`)
  - `m_multigrid` header (6 lines)

  `CLAUDE.md` asks for comments "as short as possible".
- **History and debugging notes:**
  - "(seen at 6e10)" (`m_projection.fpp:319-321`)
  - "fitted to TGV iteration counts from Mach 0.1 to 0.001" (`m_multigrid.fpp:432`)
  - "(NVHPC 26.3)" (`cmake/MFCTargets.cmake`)
- **Machine-specific notes on defaults:** "Frontier, 300^3 cells per rank" (`m_global_parameters.fpp:414`, `:416`).
- **Inline citations:** "Fuster & Popinet, JCP 374, 2018" and "Notay & Vassilevski, NLAA 15, 2008". The docs use
  `references.bib` keys.

## 8. Magic numbers where existing code uses named constants

- MPI tag `7101` (`m_multigrid_exchange.fpp:198`, `:209`).
- Side-action modes `0/1/2/3` (`:130`, `:153`, `:264`, `:277-282`).
- Side encoding `q = 1..6`, decoded with `(q + 1)/2` and `mod(q, 2)`, in three files.
- Packed row layouts `dk(7, :)` and `row(13, :)` with stride-2 slices (`m_multigrid_bottom.fpp:134`).

## 9. The workaround for dummies inside kernels is applied unevenly

- **The workaround:** several routines copy dummy arguments to locals before a kernel, with the comment "a dummy may alias a
  host variable":
  - `tl = tau; pre = mom` (`m_projection.fpp:766`)
  - `lvl = lv; ph = phase` (`m_multigrid_exchange.fpp:253`)
  - `cl = color` (`m_multigrid.fpp:689`)
- **Where it's skipped:** `lv` is still used directly inside the kernels of `s_mg_smooth`, `s_mg_restrict` and
  `s_mg_prolong`, so either the copies are unnecessary or those kernels are exposed.
- **Similar copies elsewhere:** `sl1 ... po3` (`m_projection.fpp:839-840`), `rs1/rs2` (`:395`) and `po` (`:919`, `:945`).

## 10. Constants are placed inconsistently

- **In `src/common/m_constants.fpp`:** the multigrid constants, so they compile into all three targets although only
  `simulation` uses them.
- **As module parameters:**
  - `proj_ramp_default` (`m_time_steppers.fpp:54`)
  - `capillary_alpha_min` (`m_sim_helpers.fpp:25`)

## 11. `HOST_PAGE_LOCK` macro (`src/common/include/macros.fpp:53-88`)

- **Raw CUDA call:** it binds the CUDA runtime directly through `bind(C)` behind
  `#if defined(MFC_GPU) && defined(__PGI) && ...`.
- **Failure handling:** it reports failures with `print` and doesn't abort. The hard rules ask for `s_mpi_abort`. A warning
  may be what you intend here, but it's the only `print` in the macro file.

## 12. Changes to existing behaviour (scope rather than style, but reviewers will flag them)

- **`run_time.inf` format changes for every run** (`m_data_output.fpp:115` onward, and the summary near `:1946`):
  - Column spacing goes from `13X` to `2X`, and values from `F10.6` to `ES10.3`.
  - The ICFL header becomes `AdvCFL Max` under the projection.

  `contributing.md` asks that output-format changes be flagged.
- **The explicit solver's adaptive `dt` changes:**
  - **VCFL:** now `(4/3 mu + mu_b)/rho * sum(1/dx^2)`, previously `max(1/(rho Re))/min(dx)^2`.
  - **CCFL:** now uses Brackbill's capillary density, previously the mixture density, and is zero in single-phase cells.
  - **TCFL:** now `sum(1/dx^2)`, previously `1/min(dx)^2`.

  Non-projection `cfl_adap_dt`/`cfl_const_dt` cases with viscosity, surface tension or conduction get a different `dt`.
- **Lost docstring:** a blank line now sits between the doc comment and `s_compute_stability_from_dt`
  (`m_sim_helpers.fpp:195-197`), so Doxygen no longer attaches the comment.
- **Refactor of existing code:** `s_compute_directional_rhs`'s bounds setup moved into the new `s_set_face_bounds`
  (`m_rhs.fpp:1061`).

## Matches existing conventions (no action needed)

- File headers, `#:include` lines, `implicit none` and `private; public ::`, except in `m_multigrid_levels`.
- `s_`/`f_` prefixes and `m_` module names.
- GPU macros only, with no raw pragmas; `@:ALLOCATE`/`@:DEALLOCATE` pairs on device arrays; the AMD `dimension(3)` guard.
- `BC_*` named constants; `eqn_idx` for every equation index; `wp`/`stp` kinds with explicit `real(..., wp)` conversions.
- Overall comment density (about 10% of non-blank lines, against about 8% in existing code). Only the length of individual
  blocks stands out.
- Raw `allocate` for host-only arrays and direct `MPI_*` calls, which `m_ibm`, `m_bubbles_EL` and `m_start_up` also use.
