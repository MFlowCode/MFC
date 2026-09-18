# Fourier Heat Conduction in the Energy Equation

Date: 2026-09-17
Status: approved design, not yet implemented

## Goal

Add Fourier heat conduction, $\nabla\cdot(k\nabla T)$, to MFC's mixture energy
equation for ordinary single- and multi-fluid Navier-Stokes runs. The governing
system becomes

$$\partial_t(\rho E) + \nabla\cdot\big((\rho E + p)\mathbf{u}\big)
  = \nabla\cdot(\boldsymbol\tau\cdot\mathbf{u}) + \nabla\cdot(k\nabla T).$$

Target problems are thermal boundary layers, heated or cooled walls, and
shock-heated gas. Interfacial heat transfer, bubble thermal coupling, and
conjugate heat transfer with immersed boundaries are explicitly out of scope.

## Current state

Conduction already exists in the tree, but only for reacting flow.
`s_compute_chemistry_diffusion_flux` (`src/common/m_chemistry.fpp:321`) computes
a face-averaged conductivity, a face temperature difference, and adds
`lambda_Cell*dT_dxi` to the energy source flux
(`src/common/m_chemistry.fpp:441-467`). It is gated behind `chemistry` and
`chem_params%diffusion`, and draws $\lambda$ from Cantera's mixture-averaged
transport.

Three things are missing for the general stiffened-gas path:

1. Temperature. The `q_T_sf` field exists (`src/common/m_derived_types.fpp:166`)
   but is filled only when `chemistry` is on
   (`src/simulation/m_start_up.fpp:948`).
2. A conductivity input. `physical_parameters` has no thermal conductivity
   member. Note `%K` is taken by the Herschel-Bulkley consistency index.
3. A thermal time step constraint, and the cylindrical geometric source.

Characteristic boundary conditions need no work: `m_cbc.fpp` only carries
`flux_src` over the advection indices (`src/simulation/m_cbc.fpp:131`), so the
energy source flux bypasses CBC exactly as the viscous work term already does.

## Approach

Add the conduction flux as a face-centered central difference into
`flux_src_vf(eqn_idx%E)`, generalizing the chemistry discretization to the
stiffened-gas EOS.

$\nabla\cdot(k\nabla T)$ contains no cross-derivatives, so a direction-split
face difference is exact for it, conservative, and second-order accurate. The
alternative -- threading `T` through the reconstructed-gradient machinery in
`s_get_viscous` so it rides along with the velocity gradients -- would require
growing the `iv` ranges, the reconstruction buffers, and the `rs{x,y,z}` arrays
to hold a variable that is not a member of `q_prim_vf`. That is large churn and
extra memory for no accuracy gain on this term. A cell-centered RHS source built
from `s_compute_fd_gradient` was also rejected: it is not in flux form and so
not discretely conservative.

No new halo exchange is needed. `q_prim_vf` buffers are already exchanged, so
temperature evaluated over `idwbuff` is correct in ghost cells for free.

## Design

### 1. Inputs and gating

Add `real(wp) :: k_therm` to `physical_parameters`
(`src/common/m_derived_types.fpp`), documented as thermal conductivity. The name
avoids the taken `%K`.

Register it once in `toolchain/mfc/params/definitions.py` beside the `Re(1)` and
`Re(2)` entries (around line 962), under a new `heat_conduction` group rather
than the existing `viscosity` group, with math symbol `\f$k_k\f$`. Add the
one-line description in `toolchain/mfc/params/descriptions.py`.

Add a `heat_conduction` logical to the simulation `m_global_parameters.fpp`,
defaulted `.false.` and set the same way `shear_stress` is (`src/simulation/m_global_parameters.fpp:831-836`): count
fluids with `k_therm > 0`, and if any, set the flag. It is independent of
`viscous`, so conduction runs in an otherwise inviscid case. Add it to the
`GPU_DECLARE`/`GPU_UPDATE` lists alongside `shear_stress` and `bulk_stress`, and
to the MPI broadcast list in `m_mpi_proxy.fpp`. It is a derived flag, not a
namelist entry, so it needs no per-target registration; only `k_therm` is read
from the case file, and it goes in whichever per-target parameter lists in
`definitions.py` already carry the other `fluid_pp` members.

Conductivities are stored per fluid in a device-resident
`fluid_k_therm(1:num_fluids)` array built at initialization, mirroring
`fluid_inv_re` (`src/simulation/m_global_parameters.fpp:878-879`).

### 2. Temperature

`m_phase_change.fpp:278` already defines the thermal-equilibrium mixture
temperature for stiffened gas:

$$T = \frac{\rho e + p - \sum_i \alpha_i\rho_i q_{v,i}}
             {\sum_i \alpha_i\rho_i c_{v,i} n_i},$$

with $n_i$ the isentrope exponent (`isentrope_n`) so that
$c_{p,i} = n_i c_{v,i}$. Extract this into a single shared function in
`m_variables_conversion.fpp`,

```fortran
function f_mixture_temperature(rho_e, pres, alpha_rho) result(T)
```

marked `$:GPU_ROUTINE(parallelism='[seq]')`, and call it from both
`s_infinite_pt_relaxation_k` and the new conduction path. One definition, no
second copy of the formula.

Per-phase temperatures are deliberately not used. The 5-equation system carries
a single mixture energy equation, so there is no per-phase energy to distribute
a per-phase conduction flux into. Thermal equilibrium within a cell is the
closure consistent with the model already in place, and it is the same closure
the phase-change module assumes.

Populate `q_T_sf` over `idwbuff` when conduction is active by widening the gate
at `src/simulation/m_start_up.fpp:948` from `if (chemistry)` to
`if (chemistry .or. heat_conduction)`, and by filling it from
`f_mixture_temperature` in the non-chemistry branch. The field must be refreshed
each RHS evaluation, not only at startup, so the fill belongs next to the
existing primitive-variable conversion in `m_rhs.fpp` rather than only in
start-up.

This makes `cv` a required input for any fluid with `k_therm > 0`.

### 3. Conduction flux

New module `src/simulation/m_conduction.fpp` exporting

```fortran
subroutine s_compute_conduction_source_flux(q_prim_vf, q_T_sf, flux_src_vf, norm_dir)
```

`m_riemann_state.fpp` is already 1167 lines; a separate module keeps the
conduction physics in one readable place.

For each face in the `norm_dir` direction:

- $\bar\alpha_i$ is the arithmetic average of the two adjacent cell volume
  fractions, clamped to $[0,1]$ -- the same treatment the non-Newtonian branch
  gives `alpha_avg` (`src/simulation/m_riemann_state.fpp:1085-1090`), since raw
  cell-centered alphas over- and undershoot near interfaces.
- $k_f = \sum_i \bar\alpha_i k_i$.
- $\partial T/\partial x_n = (T_R - T_L)/\Delta$, with $\Delta$ the
  cell-center-to-cell-center spacing in `norm_dir`.
- `flux_src_vf(eqn_idx%E)%sf(...) -= k_f * dT_dn`.

The loop uses the same `$:GPU_PARALLEL_LOOP(collapse=3, ...)` structure and
`idx_right_phys` indexing as `s_compute_cartesian_viscous_source_flux`.

Call it from the two sites that already call `s_compute_viscous_source_flux`:
`src/simulation/m_riemann_solver_hll.fpp:717,727` and
`src/simulation/m_riemann_solver_hllc.fpp:1612,1622`, gated on
`heat_conduction`. Because it accumulates into `flux_src_vf`, the existing
`flux_src` differencing in `m_rhs.fpp` carries it into the RHS with no further
change.

### 4. Cylindrical geometry

The axisymmetric form contributes a geometric source
$\tfrac{1}{r}k\,\partial_r T$ beyond the Cartesian divergence. `m_rhs.fpp`
already handles the analogous viscous geometric source through `tau_Re_vf`, with
special axis treatment at `src/simulation/m_rhs.fpp:1820-1909`. Add the
conduction contribution into `tau_Re_vf(eqn_idx%E)` in that same cylindrical
branch rather than introducing a parallel array, so the existing axis reflection
(`rhs_vf(i) += (tau_Re_vf(i)(j,-1,l) - tau_Re_vf(i)(j,1,l))/(y_cc(1) - y_cc(-1))`)
applies unchanged. This requires `tau_Re_vf` to be allocated when
`heat_conduction` is on even if `viscous` is off -- widen the allocation
condition at `src/simulation/m_rhs.fpp:350`.

### 5. Time step constraint

Conduction imposes $\Delta t \le \mathrm{CFL}\,\Delta x^2 \rho c_v / k$.

`s_compute_dt_from_cfl` (`src/simulation/m_sim_helpers.fpp:183`) fills
`max_dt(1:3)` for the inviscid, viscous, and capillary limits. Widen it to
`max_dt(1:4)` with a thermal slot, following the viscous branch's structure
including the `grid_geometry == 3` filtered-dtheta case. Widen
`dt_candidates_loc`/`dt_candidates_glb` in
`src/simulation/m_time_steppers.fpp:737-740` from 4 to 5 entries, and
`dt_limiter_names` (`src/simulation/m_sim_helpers.fpp:21`) from
`('ICFL','VCFL','CCFL','COLL')` to include `'TCFL'`, so run-time info reports
when conduction is the binding constraint.

Add the matching `tcfl` diagnostic to `s_compute_stability_from_dt`
(`src/simulation/m_sim_helpers.fpp:108`) and its reporting in
`src/simulation/m_data_output.fpp:165-246`. That reporting path runs through
`s_mpi_reduce_stability_criteria_extrema`, whose positional argument list is
already ten items long; since this change has to touch it, convert its max-
reduced and min-reduced scalars into two small arrays instead of appending an
eleventh argument.

### 6. Input validation

In `src/simulation/m_checker.fpp` and the shared `m_checker_common.fpp`, for any
fluid with `k_therm > 0`:

- require `cv > 0`, since the mixture temperature is undefined without it;
- require `eos` to be stiffened gas or ideal gas. Mie-Gruneisen, JWL, and Vinet
  carry their own reference-temperature paths (`mg_t0`, `jwl_t0`, `vinet_t0`)
  and are out of scope;
- reject `k_therm < 0`.

Hard error when `heat_conduction` is combined with `igr`, which has its own RHS
assembly and viscous branches in roughly ten places in `m_igr.fpp` and is
deferred to a follow-up.

Mirror the same checks in `toolchain/mfc/case_validator.py` so bad cases fail
before a job is submitted.

### 7. Verification

Three levels, in increasing order of what they actually catch:

1. **Analytic convergence.** A 1D case at uniform pressure and velocity with a
   Gaussian initial temperature, compared against the exact heat-kernel
   solution. Refine the mesh and confirm second-order convergence in the
   $L_2$ error. This is the test that catches a wrong coefficient, a factor of
   two, or a sign error -- the golden-file tests cannot.
2. **Golden tests.** Two cases in `toolchain/mfc/test/cases.py`: a 1D Cartesian
   conduction case and a 2D axisymmetric one exercising the cylindrical source.
3. **Regression.** Confirm that a case with no `k_therm` set produces
   bit-identical output to master, i.e. the feature is inert when off.

`T_wrt` already exists in post-process (`src/post_process/m_start_up.fpp:549`)
and is reused to dump the temperature field for these checks.

Add one `examples/` case, and document the new parameter in
`docs/documentation/case.md` and the term itself in
`docs/documentation/equations.md`.

## Sequencing

1. `k_therm` input plumbing and `heat_conduction` gate, inert.
2. `f_mixture_temperature` extraction and `q_T_sf` population, with the
   phase-change call site switched over and existing phase-change tests green.
3. Cartesian conduction flux plus the analytic convergence test.
4. Thermal time step constraint and reporting.
5. Cylindrical geometric source plus its golden test.
6. Validation, docs, example.

Each step leaves the tree building and the test suite green.
