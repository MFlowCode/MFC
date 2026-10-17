# Style review: `allMachProjectionCore` against the existing code

Open items only; completed and outdated ones have been removed. Compared against existing modules in `src/simulation` and
`docs/documentation/contributing.md`. This file is untracked; delete it when you're done.

## 1. Shared public mutable state

- **Public by default:** `m_multigrid_levels` has no `private; public ::` line. Everything is public except `f_seam_rank`
  (`:18`), which inverts the module convention. Three modules read and write its `mg_*` arrays. This conflicts with "No new
  global state; private helpers stay inside their defining module" in `contributing.md`.
- **Solver diagnostics in global parameters:** `proj_pcg_iters` and `proj_pcg_res` are written by `m_projection` and read by
  `m_data_output` (`m_global_parameters.fpp:24-25`).

## 2. Function declarations are inconsistent

- **Two forms are mixed:** five functions use the typed-prefix form (`f_wall_face`, `f_opposite`, `f_side_size`, `f_bot_nbr`,
  `f_bot_ax`); the rest use `result()`. Existing simulation code favours `result()` by about 25 to 6. `f_bot_nbr` and
  `f_bot_ax` use both forms at once (`integer function ... result(g)`).
- **`pure` is inconsistent:** `f_div_u_face`, `f_face_capillary` and `f_interface_normal` aren't `pure`; the helpers in
  `m_projection_helpers` are.
- **Blank line:** `f_bot_ax` has none between its signature and its declarations, unlike the other functions.
- **Module `@brief` placement:** in `m_projection.fpp:5` it comes before the `#:include` lines. The other new files, and most
  existing ones, put it directly above `module`.

## 3. Comments are long prose and carry history

- **Paragraph-length blocks that derive or justify the method:**
  - `s_mass_lag` header (8 lines)
  - `s_mg_set_operator` header (5 lines)
  - `s_mg_omega` (4-line header plus a 5-line block inside)
  - `m_multigrid_bottom` module variables (7 lines)
  - `m_multigrid` header (6 lines)

  `CLAUDE.md` asks for comments "as short as possible".
- **History and debugging notes:**
  - "(seen at 6e10)" (`m_projection.fpp:310`)
  - "fitted to TGV iteration counts from Mach 0.1 to 0.001" (`m_multigrid.fpp:417`)
  - "mg_e/mg_f (NVHPC 26.3)" (`cmake/MFCTargets.cmake:75`); `mg_e/mg_f` no longer exist either
- **Machine-specific notes on defaults:** "Frontier, 300^3 cells per rank" (`m_global_parameters.fpp:414`, `:416`).
- **Inline citations:** "Fuster & Popinet, JCP 374, 2018" (`m_projection.fpp:5`) and "Notay & Vassilevski, NLAA 15, 2008"
  (`m_multigrid.fpp:582`). The docs use `references.bib` keys.

## 4. Magic numbers where existing code uses named constants

- Side-action modes `0/1/3` (2 at unpack) in `mg_smode`; proposed `mg_side_mode` with named constants.
- Side encoding `q = 1..6`, decoded with `(q + 1)/2` and `mod(q, 2)`, in three files.
- Packed row layouts `cell_rows(7, :)` and `rows_loc(13, :)` with stride-2 slices (`s_mg_bottom_build`).

## 5. The workaround for dummies inside kernels is applied unevenly

- **The workaround:** two routines copy dummy arguments to locals before a kernel ("a dummy may alias a host variable"):
  `tau_loc = tau; mom_loc = mom` (`m_projection.fpp:790`) and `lv_loc = lv; phase_loc = phase` (`m_multigrid_exchange.fpp:159`).
- **Where it's skipped:** kernels index `mg_lev(lv)` with the dummy `lv` directly (smooth, restrict, prolong, omega,
  set_operator). This passes on NVHPC OpenACC but is untested on Cray/AMD OpenMP, so either the copies are unnecessary or
  those kernels are exposed.

## 6. Changes to existing behaviour (scope rather than style, but reviewers will flag them)

- **`run_time.inf` format changes for every run** (`m_data_output.fpp`, header and summary):
  - Column spacing goes from `13X` to `2X`, and values from `F10.6` to `ES10.3`.
  - The ICFL header becomes `AdvCFL Max` under the projection.

  `contributing.md` asks that output-format changes be flagged.
- **The explicit solver's adaptive `dt` changes:**
  - **VCFL:** now `(4/3 mu + mu_b)/rho * sum(1/dx^2)`, previously `max(1/(rho Re))/min(dx)^2`.
  - **CCFL:** now uses Brackbill's capillary density, previously the mixture density, and is zero in single-phase cells.
  - **TCFL:** now `sum(1/dx^2)`, previously `1/min(dx)^2`.

  Non-projection `cfl_adap_dt`/`cfl_const_dt` cases with viscosity, surface tension or conduction get a different `dt`.
- **Lost docstring:** a blank line sits between the doc comment and `s_compute_stability_from_dt`
  (`m_sim_helpers.fpp:194-196`), so Doxygen doesn't attach the comment.
- **Refactor of existing code:** `s_compute_directional_rhs`'s bounds setup moved into the new `s_set_face_bounds`
  (`m_rhs.fpp:1061`).
