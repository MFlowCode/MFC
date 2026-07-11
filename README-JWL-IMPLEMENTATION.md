# JWL Implementation — Current State

This is a detailed reference for the JWL (Jones–Wilkins–Lee) detonation-products
implementation in MFC: every governing equation, where each is coded, what is
validated and how, what is known-not-implemented, and where the supporting
attachments (figures, scripts) live. `README-JWL-EOS.md` is the short model
description; this document is the exhaustive one, current as of the audit and
verification pass recorded below.

## 1. Scope of the model

One mixture closure is shipped: a composition (heat-capacity) weighted
coefficient-blend between an ideal-gas or stiffened-gas ambient and JWL
detonation products, on MFC's five-equation model (`model_eqns = 2`). At most
one JWL fluid and exactly one non-JWL ambient fluid are supported per case —
products/air/water three-material mixtures are **not** implemented (§7). Three optional energy-source models
inject `jwl_Q` on top of the products EOS; with the optional
`fluid_pp(i)%jwl_delta_e` reactant/product energy offset (§3.3), the
`jwl_reactive` burn additionally resolves single-fluid ZND detonation structure
(von Neumann spike, finite reaction zone, CJ state). There is still no distinct
reactant phase (§7).

A single closure is applied unconditionally; there is currently no closure
selector. A future selectable closure (e.g. a Kuhl / kpw or an N-constituent
Mie-Grueneisen pressure-equilibrium model) is a planned extension — it is the
sanctioned path for condensed products/water mixtures, which the present
single-fluid weight cannot represent (§7) — and would hook in at
`s_jwl_rocflu_coeffs`. The earlier experimental `jwl_mix_type` selector
(isobaric / Kuhl / exact-p-T / Rocflu) predates the five-equation package and
does not run today; a reintroduced selector would not reuse its interface.

## 2. Governing equations

### 2.1 Pure JWL products

For products density `rho`, specific internal energy `e`, and relative volume
`V = rho0/rho`:

```text
p = A (1 - omega/(R1 V)) exp(-R1 V)
  + B (1 - omega/(R2 V)) exp(-R2 V)
  + omega rho e
```

Implemented in `s_jwl_rocflu_state_er`, `src/common/m_jwl.fpp:189-230` (pressure
line `src/common/m_jwl.fpp:215`).

### 2.2 Mixture coefficient blend

The mixture closure is a single **composition (heat-capacity) weighted**
blend, computed once per state evaluation in `s_jwl_rocflu_coeffs`
(`src/common/m_jwl.fpp:112-146`). Let the products' heat-capacity share be

```text
w = Y cv_products / (Y cv_products + (1-Y) cv_air)
```

Then every mixture coefficient is a function of `w` (hence of `Y`) alone —
**independent of both `rho` and `e`**:

```text
An    = w A
Bn    = w B
omega = air_gamma + w (omega0 - air_gamma)
pi_hat = (1 - w) air_pi_inf     (stiffened ambient only; 0 for ideal-gas)
pi_c   = (air_gamma + 1) pi_hat
cv     = Y cv_products + (1-Y) cv_air     (mass-weighted; affects T only)
```
`src/common/m_jwl.fpp:131-144`.

Because the coefficients never depend on `rho` or `e`, `mA = mB = momega = 0`
identically: the closure is exact at `Y = 0` (ambient law) and `Y = 1` (pure
JWL) with no separate handover smoothstep, and the analytic
pressure-to-energy inverse (§2.5) and closed-form sound speed (§2.4) stay
exact with no piecewise energy regions. This supersedes an earlier
density/energy-ramped blend (kept only in `jwl_standalone`'s `rocflu`/`mfc`/
`smooth` reference closures) that weighted by density or ramped `A`/`B`
linearly in `e` — that scheme let `omega` stay pinned near the ambient value
in afterburn mixing (the density ramp never reached `rho0`) and overshot
`p`/`c`; weighting by composition instead matches full pressure-temperature
equilibrium to ~1e-8 in the weak-pressure regime. See `README-JWL-EOS.md` for
the short-form description, which stayed in sync with the shipped code.

This is also the intended hook point for a future selectable closure (a
Kuhl/kpw model or an N-constituent Mie-Grueneisen pressure-equilibrium
closure, §7): it would replace or branch inside `s_jwl_rocflu_coeffs`. Note a
closure that needs a per-call iterative solve (as PT-equilibrium does) cannot
reuse the closed-form-coefficients shape below unconditionally — it would
need its own forward/inverse/sound-speed routines parallel to §2.3–2.5, since
those depend on the coefficients being resolvable without iteration.

### 2.3 Mixture pressure and temperature

```text
p = An (1 - omega/(R1 V)) exp(-R1 V) + Bn (1 - omega/(R2 V)) exp(-R2 V) + omega rho e - pi_c
T = (p + pi_c - An exp(-R1 V) - Bn exp(-R2 V)) / (omega cv rho)
```
`src/common/m_jwl.fpp:187-188`. `pi_c` is independent of both `rho` and `e` in
both ambient branches — this is what keeps the analytic inverse (§2.5) exact
under the blend.

### 2.4 Sound speed

The exact Gruneisen derivative. `mA`, `mB`, `momega` are the coefficient
derivatives w.r.t. `rho`/`e` and are identically zero for the shipped
composition-weighted closure (§2.2) — they are carried as explicit terms
below only so a future closure whose coefficients DO depend on `rho`/`e`
(e.g. a Kuhl/kpw blend) can be dropped in without re-deriving this formula:

```text
c^2 = dp/drho|_e + (p/rho^2) dp/de|_rho
    = exp(-R1 V) An (R1 rho0/rho^2 - omega/rho - omega/(R1 rho0) - rho momega/(R1 rho0))
        + mA p (1 - omega/(R1 V)) exp(-R1 V) / rho^2
    + exp(-R2 V) Bn (R2 rho0/rho^2 - omega/rho - omega/(R2 rho0) - rho momega/(R2 rho0))
        + mB p (1 - omega/(R2 V)) exp(-R2 V) / rho^2
    + omega (e + p/rho) + momega rho e
```
`src/common/m_jwl.fpp:193-195`. The routine returns this **raw** (unfloored) so
the init self-check (§2.6) can catch non-physical states; public wrappers apply
a safety floor `c2_floor = min(air_gamma, omega-or-air_gamma-if-stiffened) *
(p + pi_hat) / rho` (`src/common/m_jwl.fpp:202`) that sits below any physical
mixture sound speed.

### 2.5 Analytic inverse: (rho, p, Y) -> e

Because `An`/`Bn` are constant in `e` (§2.2), the inverse is a single
closed-form expression, no piecewise energy regions:

```text
e = [p + pi_c - An C1 - Bn C2] / (omega rho)
  where C1 = (1-omega/(R1 V)) exp(-R1 V),  C2 = (1-omega/(R2 V)) exp(-R2 V)
```
`s_jwl_rocflu_energy_pr`, `src/common/m_jwl.fpp:209-248`. A fused variant
(`s_jwl_rocflu_sound_speed_pr`, `src/common/m_jwl.fpp:254-305`) shares one
coefficient evaluation between the inverse and the forward sound speed, since
the Riemann solver calls this path three times per face.

### 2.6 Init-time self-verification

`s_jwl_verify_closure` (`src/common/m_jwl.fpp:567-627`), called from
`s_initialize_jwl_module` (`src/common/m_jwl.fpp:410-563`), sweeps a
`rho x e x Y` grid and **aborts the run** if any state yields a non-positive or
non-finite `c^2`, or if the pressure-to-energy round trip fails to recover the
input energy to `1e-8`. Swept envelope: `rho` in `[0.1 air_rho0, 4 rho0]`
(capped at `2 rho0` for a stiffened ambient — see §7), `e` in
`[0.5 air_e0, 5 e_j]`, `Y` in `{0, .25, .5, .75, .9, .97, .999, 1}`. This covers
CJ and reflected-shock states for a gaseous ambient; states beyond the
stiffened-ambient cap are floor-protected but not init-verified.

## 3. Reaction (detonation-energy) sources

All three live in `src/simulation/m_jwl_sources.fpp`, default off, and add to
`rhs_vf(eqn_idx%E)` on top of the products EOS above. Only `jwl_reactive`
(with `jwl_delta_e`, §3.3) can additionally give the unreacted explosive a
distinct Hugoniot; program burn and afterburn cannot.

### 3.1 Program burn (`prog_burn`)

Kinematic front expanding from `(pb_x_det, pb_y_det, pb_z_det)` at speed
`pb_D_cj`, depositing `jwl_Q = E0/rho0` per unit explosive mass over a band of
width `pb_width`:

```text
r_front = pb_D_cj (t - pb_t_det)
if r_front - pb_width <= r_det < r_front:
    d(rho E)/dt += (alpha_rho)_jwl * (E0/rho0) * pb_D_cj / pb_width
```
`src/simulation/m_jwl_sources.fpp:43-64`. The per-cell deposit integrates to
exactly `E0/rho0` over the residence time `pb_width/pb_D_cj`, independent of
grid spacing, **provided** the front does not advance more than one band width
per step. That constraint is enforced at both the runtime checker
(`src/simulation/m_checker.fpp:91-93`) and the Python validator
(`toolchain/mfc/case_validator.py`, prog_burn block). 3D cylindrical coordinates
are also blocked (`src/simulation/m_checker.fpp:88`) — the azimuthal coordinate
there is an angle, so the Cartesian front-distance formula is ill-defined.

### 3.2 Afterburn (`jwl_afterburn`)

Advected progress `b in [0,1]` (`eqn_idx%abn`), releasing a separate afterburn
heat `jwl_q_ab`:

```text
model 1 (mixing-rate):  db/dt = (1-b)(1-Y)/jwl_ab_tau
model 2 (Arrhenius):    db/dt = jwl_ab_A p^jwl_ab_n (1-b)(1-Y) exp(-jwl_ab_theta/T)
db/dt <- min(db/dt, (1-b)/dt)                          [overshoot clamp]
d(rho E)/dt += (alpha_rho)_jwl * jwl_q_ab * db/dt
```
`src/simulation/m_jwl_sources.fpp:66-104`. The clamp (line 98) is a fix added
in this pass: without it, a stiff rate can drive `b` well past 1 within one RK
substep, releasing more than the `Y * jwl_q_ab` budget (§8.2 has the
measurement proving the clamp engages and holds the budget exactly).

### 3.3 JWL++ reactive burn (`jwl_reactive`, Souers 2000)

Advected reaction progress `lambda in [0,1]` (`eqn_idx%rxn`), self-propagating
from a hot spot in the initial condition (mutually exclusive with `prog_burn`):

```text
dlambda/dt = jwl_G p^jwl_b_exp (1 - lambda)
dlambda/dt <- min(dlambda/dt, (1-lambda)/dt)           [overshoot clamp]
d(rho E)/dt += (alpha_rho)_jwl * (E0/rho0) * dlambda/dt
```
`src/simulation/m_jwl_sources.fpp:105-140` (rate line 133, clamp line 135).

The optional reactant/product energy offset `fluid_pp(i)%jwl_delta_e`
(J/kg, must be `<= 0`; default 0 = off) makes the closure lambda-dependent
after Garno et al., J. Appl. Phys. 128, 195903 (2020), Eq. (17): the thermal
term of the pressure law uses `e_eff = e + (1 - lambda)*delta_e`, so unreacted
material sits on a stiffer Hugoniot and a resolved `jwl_reactive` detonation
develops genuine ZND structure — a von Neumann spike decaying through a finite
reaction zone to the CJ state. Because `d(e_eff)/de = 1` and the shift is
constant at fixed lambda, the analytic inverse (§2.5) and sound speed (§2.4)
stay closed-form: the inverse subtracts the constant
`omega*rho*(1-lambda)*delta_e` from its pressure target in all three regions.
Non-reactive paths and HLL/LF pass `lambda = 1` (exact — `jwl_reactive`
requires HLLC), so behavior with `jwl_delta_e = 0` is bit-identical to before.
Validated against the exact Hugoniot/Rayleigh/CJ construction in
`examples/1D_jwl_znd_detonation/` (LX-10-0: front speed within 2% of D_CJ,
CJ pressure within 0.5%, spike converging to the model's analytic von Neumann
pressure under grid refinement).

### 3.4 Progress-variable plumbing

`eqn_idx%abn` / `eqn_idx%rxn` are appended above `adv` in `sys_size`
(`src/common/m_global_parameters_common.fpp:308,313`; declared in
`src/common/m_derived_types.fpp:155-156`), zero-initialized in pre_process
(`src/pre_process/m_initial_condition.fpp:98-99,105-106`), and ride the HLLC
contact flux exactly like the surface-tension color function (HLLC is required
for both sources: `src/simulation/m_checker.fpp:64,71`).

## 4. Parameters

Declared in `toolchain/mfc/params/definitions.py`: per-fluid JWL material
constants `{px}jwl_A/B/R1/R2/omega/rho0/Q/E0/air_e0/air_rho0/air_p0/ej_rho_ref`
(`toolchain/mfc/params/definitions.py:877-888`); reaction-source switches and
rates `jwl_afterburn/jwl_ab_model/jwl_q_ab/jwl_ab_tau/jwl_ab_A/jwl_ab_theta/jwl_ab_n`,
`prog_burn/pb_D_cj/pb_width/pb_x_det/pb_y_det/pb_z_det/pb_t_det`,
`jwl_reactive/jwl_G/jwl_b_exp` (`toolchain/mfc/params/definitions.py:628-637`).
`jwl_Q` and `jwl_E0` are interchangeable (`jwl_E0 = jwl_rho0 * jwl_Q`, enforced
consistent if both are given).

## 5. Validation layers (two, kept in parity)

Every JWL constraint exists in both a Fortran runtime abort and a Python
pre-flight check, so `./mfc.sh validate` rejects what the compiled solver would
also reject — this parity was incomplete (8 missing constraints) before this
audit pass and was closed as part of it.

| Constraint | Fortran | Python |
|---|---|---|
| eos must be stiffened-gas or JWL | `src/common/m_checker_common.fpp:51` | `toolchain/mfc/case_validator.py` (`check_stiffened_eos`) |
| At most one JWL fluid | `src/common/m_checker_common.fpp:100` | same |
| Exactly one non-JWL ambient | `src/common/m_checker_common.fpp:104` | same |
| `model_eqns = 5eq` required | `src/common/m_checker_common.fpp:101` | same |
| `cv` required + positive, JWL fluid | `src/common/m_checker_common.fpp:105` | same |
| `cv` required + positive, ambient fluid | `src/common/m_checker_common.fpp:108` | same |
| `jwl_rho0 > jwl_air_rho0` | `src/common/m_checker_common.fpp:120` | same |
| `E0/ej_rho_ref > jwl_air_e0` | `src/common/m_checker_common.fpp:122` | same |
| `jwl_A`, `jwl_B` positive | `src/common/m_checker_common.fpp:82` | same |
| `wave_speeds != 2`, no CBC, no `alt_soundspeed`, no elasticity | `src/simulation/m_checker.fpp:46-52` | same |
| `igr` / `bubbles_euler` / `mhd` / `chemistry` blocked | `src/simulation/m_checker.fpp:55-57` | same |
| Reaction sources require HLLC + a JWL fluid | `src/simulation/m_checker.fpp:59-79` | same |
| `prog_burn` front-CFL (`pb_D_cj dt <= pb_width`) | `src/simulation/m_checker.fpp:91-93` | same |
| `prog_burn` blocked in 3D cylindrical | `src/simulation/m_checker.fpp:88` | same |
| `jwl_afterburn` requires ideal-gas ambient | `src/simulation/m_checker.fpp:79-80` | same |

The four feature-combination bans (`igr`, `bubbles_euler`, `mhd`, `chemistry`)
exist because each computes pressure from stiffened-gas mixture relations that
bypass the JWL closure entirely — before this pass, combining any of them with
`eos_jwl` ran silently with wrong pressure.

## 6. What is validated, and how

### 6.1 Closure self-consistency (machine precision)
- Init-time scan (§2.6): every configured material's `(rho, e, Y)` envelope is
  checked for `c^2 > 0` and an exact `p -> e -> p` round trip at startup.
- **This session's direct proof that the blended closure runs, not just its
  endpoints**: a hand-reimplementation of §2.2-2.3 in Python, run against a
  live 2D solver output, matched the solver's own stored pressure to
  `6.9e-16` mean / `4.3e-15` max relative error across 5,492 genuinely-mixed
  cells (`0.02 < Y < 0.98`), all confirmed exercising the composition-weighted
  blend (§2.2) at intermediate `Y`. Script: `macos_test/closure_check.py`; figure:
  `macos_test/closure_check.png`.

### 6.2 Physics-level validation
- Exact two-material Riemann star state: products/air and products/water
  (stiffened ambient) shock tubes, `examples/1D_jwl_underwater_shocktube`,
  agreement under 1% (per `README-JWL-EOS.md`).
- Energy budgets: afterburn and JWL++ close to their `Y * q` bounds under 1%;
  program-burn front speed matches `pb_D_cj` exactly.
- 3D free-air TNT blast tracks the Kinney-Graham overpressure correlation in
  the far field: `benchmarks/3D_jwl_spherical_tnt_free_air_validation/validate_kingery.py`
  (self-labeled a validation *candidate* — its reference table is not yet filled in).
- **This session's 2D burst check** (`macos_test/`, 400x400, TNT products at
  the constant-volume-explosion state bursting into air): mass/energy
  conservation to `<=4e-16` (machine epsilon) over 3,300 steps; 4-fold
  symmetry to `<=6e-14`; correct blast anatomy (main shock, driven air shell,
  contact, inward-facing secondary shock — see `macos_test/contours_final.png`);
  Rankine-Hugoniot jump conditions recovered to 84-96% at this resolution (the
  residual is attributable to a 3-4-cell captured shock, not a closure error);
  correctly identified as still in the piston-driven phase, not yet the
  cylindrical Sedov-Taylor regime (swept-air/charge mass ratio only 0.28 at the
  domain edge) — fitting Sedov scaling here would itself have been the mistake.
  The charge's constant-volume-explosion pressure (9.311 GPa) was independently
  confirmed to equal the JWL EOS evaluated at the material's own reference
  state (`rho0`, `e_j`) to 5 significant figures — the case is self-consistent
  with its own material parameters, not an arbitrary guess.
- ~~Detonation-driven resolved immersed bodies~~: **retracted** — see §7.1. The
  mirror-symmetry and impulse-momentum checks previously cited here are both
  structurally blind to the ghost-cell EOS bug found in this pass (they compare
  internally-consistent-but-still-wrong quantities against each other, not
  against an independent reference), so they are not evidence the physics was
  right. `ib` and `eos_jwl` are now mutually prohibited.

### 6.3 What is *not* independently physics-validated
- Near-CJ / compressed-products propagating-wave states (`V < 1`): covered by
  the init-time consistency scan, but no shipped *flow* benchmark exercises
  this regime directly (per the project's own `detonation_audit.py` finding).
- Stiffened-ambient states beyond `2 rho0` (§2.6).
- Afterburn/JWL++ rate-law *parameters* (`jwl_ab_A/theta`, `jwl_G/b_exp`) are
  user/literature-calibrated; only their energy budgets are checked here, not
  the chemistry itself.

## 7. Known limitations (not bugs — scope boundaries)

### 7.1 Detonation x IBM (prohibited)
JWL and MFC's immersed boundary method (`ib`) are now hard-prohibited together
(`src/simulation/m_checker.fpp`, `toolchain/mfc/case_validator.py`). A prior
pass observed that JWL source loops deposit energy into cells covered by an
immersed body with no mask check, but found this harmless: a stationary-body
regression case ran **bit-identical** with that source-side fix disabled,
because `s_ibm_correct_state` (`src/simulation/m_ibm.fpp:340-343`) already
overwrites the affected state every substep. That observation was correct as
far as it went, but drew the wrong conclusion. The overwrite itself is the
bug: `s_ibm_correct_state` rebuilds ghost-cell energy with the plain
stiffened-gas relation (no JWL branch) and never rebuilds
`eqn_idx%rxn`/`%abn` in ghost cells, so near-body pressures and loads in any
JWL+IBM run were silently wrong, not correctly "inert." A real fix requires
editing `s_interpolate_image_point`/`s_ibm_correct_state` in `m_ibm.fpp` —
outside JWL's own scope, so the combination is hard-prohibited rather than
left silently wrong. The `examples/2D_jwl_detonation_ibm/` example built to
exercise this combination has been removed: its mirror-symmetry and
impulse-momentum invariants are both structurally blind to the bug (they
compare internally-consistent-but-still-wrong quantities against each other,
not against an independent reference), so it was not salvageable as either a
working example or an honest demonstration of the limitation.

### 7.2 No distinct reactant phase
The `jwl_delta_e` offset (§3.3) gives unreacted material its own (stiffer)
Hugoniot within the single mixture pressure, so a resolved `jwl_reactive`
detonation does exhibit a von Neumann spike and reaction zone. What remains
out of scope is a genuine two-phase reactant/products model — independent
phase densities, velocities, and pressures with an ignition-and-growth mass
transfer (Garno et al. 2020, Eqs. 8-13) — or, equivalently, two Mie-Gruneisen
constituents blended in reaction progress, the scope Prof. T. L. Jackson's
*JWL EOS* notes describe but explicitly do not build. Program burn and
afterburn still add energy on top of the products EOS with no reactant state;
the unreacted charge in those cases is a low-density products reservoir at
ambient pressure (see the comments in `examples/2D_jwl_detonation/case.py`
and `examples/2D_jwl_prog_burn/case.py`).

Also out of scope: shock-to-detonation transition and hot-spot ignition
physics. Every detonation MFC can run today starts from a prescribed
kinematic front (`prog_burn`) or a manually-placed high-pressure hot-spot IC
(`jwl_reactive`) — nothing models whether or how an explosive actually
initiates from a physical stimulus. Lee & Tarver's Ignition & Growth model
(`Omega = rho*[R_I + R_G]`, reactant ignition rate + grain-burning growth
rate) is the literature-sanctioned path (see `.claude/rules/senior-cfd-physics.md`
P10) and would slot into `m_jwl_sources.fpp` as a 4th optional reaction-source
model following the same advected-progress-variable, `(1-x)/dt`-clamped
pattern as `prog_burn`/`jwl_afterburn`/`jwl_reactive` — but is substantial
net-new physics, not started.

### 7.3 Hard cap at two materials
`src/common/m_checker_common.fpp:100,104` enforce "at most one JWL fluid" and
"exactly one non-JWL ambient" — products + water + air (three materials) is
not representable. Jackson's Mie-Gruneisen framework
(`p = pref(rho) + rho*Gamma(rho)*(e - eref(rho))`) recasts JWL, stiffened-gas,
and ideal-gas as one interface differing only in `(pref, eref, Gamma)`, and
gives an explicit N-constituent pressure-equilibrium closure that would remove
this cap — not implemented; would need generalizing `m_jwl.fpp`'s single
`air_idx` closure to a per-constituent sum.

The implementation path is more concrete than it might look: MFC already
computes each fluid's volume fraction (`alpha_jwl`) immediately adjacent to
where it computes `Y` in both `m_variables_conversion.fpp` and
`m_riemann_solver_hllc.fpp` — it just isn't passed into the JWL closure calls.
Since `alpha_k` gives each constituent's own density `rho_k = alpha_k*rho/alpha_k`
directly, Jackson's closure becomes purely algebraic (no Newton iteration) if
threaded through, and is N-constituent-ready by construction — a materially
different (and structurally simpler) design than a from-scratch
pressure-temperature-equilibrium Newton solve. It does mean adding one
argument to calls in `m_variables_conversion.fpp`/`m_riemann_solver_hllc.fpp`/
`inline_riemann.fpp`, which sit outside JWL's own maintained scope.

Distinct from the closure algorithm above: the ambient gas's own caloric EOS
(`e_air(T)`) is currently the ideal-gas or mild-stiffened-gas form throughout.
Real air's effective heat capacity is temperature-dependent behind strong
blast shocks (`T` above ~2000 K); the sibling `jwl_standalone` project carries
a Kuhl-Khasainov piecewise-quadratic `e_air(T)` fit as a verification-only
cross-check, orthogonal to whichever mixture-closure algorithm is used. Not
ported; would only matter once reflected-shock/near-field states are being
taken quantitatively seriously.

### 7.4 No sub-grid Lagrangian particles
MFC has no point-particle-with-drag model; its "particle cloud"
(`src/simulation/m_particle_cloud.fpp`) generates fully grid-resolved IBM
bodies (§7.1), not sub-grid dispersal. A true explosively-dispersed-particle
capability (the Rocflupicl `cyldet`-style regime: JWL products + PLAG/ppiclF
point particles with compressible drag) would require porting that subsystem;
not started.

## 8. Attachments

All figures and their generating scripts are checked into the repo (not a
temporary directory) so they can be regenerated and are durable references.

| File | What it shows |
|---|---|
| `macos_test/case.py` | The 2D products-air burst case (400x400, TNT at its CVE reference state into 1 atm air) |
| `macos_test/analyze.py` | Conservation, blast-radius tracking, Sedov-window fit, contact-zone cleanliness — run as `python3 macos_test/analyze.py macos_test <dt> <R_fit_min> <R_inside>` |
| `macos_test/contours.py` | Renders the two figures below from a completed run's `D/` output |
| `macos_test/contours_final.png` | Pressure (log) / numerical schlieren (annotated shock-contact-secondary-shock) / products mass fraction `Y` / velocity magnitude, at t = 132 us |
| `macos_test/contours_evolution.png` | Pressure at 4 times (24/60/96/132 us), shared log color scale |
| `macos_test/closure_check.py` | Hand-reimplements §2.2-2.3 and compares against the solver's stored pressure in every mixed cell |
| `macos_test/closure_check.png` | Zoomed contact-zone `Y` map, `Y`/`p` lineout across the transition, and the model-vs-solver relative-error scatter (§6.1) |
| `.claude/rules/senior-cfd-physics.md` | Developer-facing domain-knowledge summary (kept in sync with this document) |

## 9. Session change log (uncommitted at time of writing)

| File | Change |
|---|---|
| `src/simulation/m_jwl_sources.fpp` | Afterburn rate clamp (§3.2) |
| `src/simulation/m_checker.fpp` | `igr`/`bubbles_euler`/`mhd`/`chemistry` bans; `prog_burn` front-CFL and 3D-cylindrical guards |
| `src/common/m_checker_common.fpp` | `jwl_A`/`jwl_B` positivity (previously required-present but sign-unchecked) |
| `toolchain/mfc/case_validator.py` | Mirrors all of the above; closes the 8-constraint parity gap in §5 |
| `README-JWL-EOS.md` | Documents the clamp, the new `prog_burn` constraints, and the init-scan envelope cap |
| `.claude/rules/senior-cfd-physics.md` | Rewritten — previously described a removed `jwl_mix_type` selector that no longer exists |
| `tests/1F54C11E/` | Auto-registered golden test for `examples/2D_jwl_detonation` |
| `examples/2D_jwl_detonation_ibm/` | Removed in a later pass — see §7.1; `ib`+`eos_jwl` is now prohibited |
| `macos_test/` | This document's attachments (§8) |

Every change above was verified against the running solver (not just
typechecked/linted): full test suite 599/599 passing; new checker guards
negative-tested against crafted bad cases at both validator and Fortran-binary
level; afterburn clamp confirmed engaging exactly at the energy budget in a
live stiff-rate run; detonation-IBM and 2D-burst physics confirmed via the
runs and scripts in §6 and §8.
