# EOS Family Registry Implementation Plan

**Goal:** One `EOS_FAMILIES` registry drives every mechanical EOS site; adding a family touches the registry plus the three genuinely per-family artefacts.

**Spec:** registry-design-STAGED.md (to become docs/superpowers/specs/2026-09-12-eos-family-registry-design.md)

**Branch:** module-eos, as separate commits after the extraction commits.

## Global Constraints

- Bit-for-bit goldens: `./mfc.sh test` = 713 passed / 0 failed / 39 skipped / 752 total. Never `--generate`.
- `m_eos` stays a leaf (m_derived_types, m_global_parameters_common, m_constants). A generated `#:include` is textual and adds no `use`.
- Build/test only via `./mfc.sh`; sticky lock in `build/`.
- Commits use `--no-verify`; precheck must show ONLY the 3 known doc-reference failures.
- Never `git add -A`. Scratch files go to /tmp or the gitignored workspace.

## Verification strategy

Python-only tasks (R1, R2, R5) run the toolchain unit tests, which are fast — no 2-hour suite.
R3 and R4 touch Fortran and need a build. R4's generated predicates are held to a stronger
bar than a test run: the generator must emit text byte-identical to today's hand-written
bodies, which makes the change provably behaviour-preserving without a suite. The full suite
runs once, at R6.

---

### Task R1: The registry, and definitions.py

**Files:** create `toolchain/mfc/params/eos_families.py`; modify `toolchain/mfc/params/definitions.py:914-956`.

One frozen dataclass per family carrying: `value`, `suffix` (Fortran constant name after `eos_`),
`label`, `prefix`, `required` and `optional` as ordered (name, math-symbol) pairs,
`state_dependent: bool`, `isentropic_reference: bool`, `coefficients_fn` (the `eos.py` name).

Then derive, replacing the literals in place:
- `_EOS_NAMES` -> `{f.suffix: f.value for f in EOS_FAMILIES}`
- `_EOS_VALUE_LABELS` -> `{f.value: f.label for f in EOS_FAMILIES}`
- `"choices"` -> `sorted(f.value for f in EOS_FAMILIES)`
- the three `mg_*`/`jwl_*`/`vinet_*` registration blocks -> one loop over the registry

Verify: `_EOS_NAMES`, `_EOS_VALUE_LABELS` and `choices` must equal their previous literal values
exactly, and the set of registered `fluid_pp(1)%*` parameter names must be unchanged. Assert this
with a test that hard-codes the old values — the point is to prove the refactor changed nothing.

Run the toolchain unit tests. Commit.

### Task R2: case_validator.py

**Files:** modify `toolchain/mfc/case_validator.py` at 969-974, 1004-1007, 1067-1073, 2050.

- `families` -> `{f.value: (f.prefix, tuple(n for n, _ in f.required)) for f in EOS_FAMILIES if f.state_dependent}`
- `optional` -> the same shape from `f.optional`
- `_check_initial_states_inside_eos`'s dispatch -> resolve `f.coefficients_fn` from the registry
- both hard-coded state-dependent trios -> `{f.value for f in EOS_FAMILIES if f.state_dependent}`

Verify: the existing validator tests pass unchanged, and the derived `families`/`optional` dicts
equal their previous literals. Run the toolchain unit tests. Commit.

### Task R3: Generate the Fortran eos_* constants

**Files:** modify `toolchain/mfc/params/generators/fortran_gen.py:257-272`; modify `src/common/m_constants.fpp:117-124`.

Add an optional `"fortran_prefix"` to the `fluid_pp(1)%eos` constraint entry. In
`generate_constants_fpp`, emit `<fortran_prefix>_<name>` for a compound key that carries one
instead of skipping it. Delete the hand-written `eos_*` block from `m_constants.fpp`.

Verify: the generated `generated_constants.fpp` must contain exactly the five constants with the
same values the hand-written block had. `test_eos_selector.py::test_fortran_and_python_enums_agree`
must still pass (it now proves the generator, not a hand-copy). Build. Commit.

### Task R4: Generate the two predicates

**Files:** `toolchain/mfc/params/generators/fortran_gen.py` (new emitter); `src/common/m_eos.fpp` (`f_is_state_dependent`, `f_has_isentropic_reference`).

Emit a `generated_eos.fpp` include holding the two predicate bodies as `select case` blocks over
the registry. `m_eos` consumes it with `#:include`.

**The bar:** before switching `m_eos` over, assert the generated text is byte-identical in meaning
to the current bodies — same families in `f_is_state_dependent` (mie_gruneisen, jwl, vinet), same
condition in `f_has_isentropic_reference` (jwl or vinet, AND `gruneisen_a == 0`). Note the second
predicate has a runtime term that is NOT a family property; the registry supplies only the family
set, and the `gruneisen_a` test stays in the Fortran. Getting this wrong silently changes which
fluids take the closed-form isentrope.

Build. Commit.

### Task R5: The agreement tests

**Files:** create `toolchain/mfc/params_tests/test_eos_families.py`.

Four tests, each closing one leg of the web the spec maps:
1. Every registry parameter has a matching field on `physical_parameters` in `m_derived_types.fpp`.
2. Every required parameter of a family appears in that family's `case` body in `s_initialize_eos_module`. **This is the test that catches "validates cleanly but computes with dflt_real".**
3. Every `coefficients_fn` named by the registry resolves in `toolchain/mfc/eos.py`.
4. Every family's `suffix` appears as a `case (eos_<suffix>)` in `s_reference_curve`.

Verify each test FAILS when its target is broken — introduce the break, observe the failure, revert.
A test that cannot fail is worse than no test. Run the toolchain unit tests. Commit.

### Task R6: Docs and final verification

Update `docs/documentation/contributing.md`'s "How to Add an Equation of State" to describe the
registry as the entry point and name the three artefacts that stay hand-written.

Then: clean build; full suite = 713/0/39/752; precheck = the 3 known failures only;
`git diff --stat upstream/master -- tests/` empty; `git status --porcelain` empty.

---

## Pre-flight findings (controller, verified before execution)

**R3 is mechanically ready.** `src/common/m_constants.fpp:131` already carries
`#:include 'generated_constants.fpp'`, and `cmake/ParamsCodegen.cmake` emits that file per
target. So deleting the hand-written `eos_*` block at `m_constants.fpp:117-124` and teaching
`generate_constants_fpp` to emit compound-key constants is the whole change. No new plumbing.

**Direction of truth, confirmed.** `definitions.py` reads Fortran constants via `_fc()` only for
ARRAY BOUNDS (`num_fluids_max`, `num_patches_max`, ...) where Fortran is authoritative. For
`eos_*` the direction is the opposite: Python's `_EOS_NAMES` is the source and `m_constants.fpp`
holds a hand copy. Generating it therefore removes a copy rather than inverting a dependency.

**R3 breaks an existing test, and must.** `params_tests/test_eos_selector.py::test_fortran_and_python_enums_agree`
calls `get_fortran_constants()`, which parses `src/common/m_constants.fpp` TEXTUALLY
(`namelist_parser.py:95`) and does not follow `#:include`. Once the constants are generated that
test finds nothing and fails. It cannot simply be pointed at the generated file either: asserting
the generator agrees with the registry it was generated from is circular and proves nothing.

Replace it with two tests that are not circular:
1. `generate_constants_fpp()` emits exactly one `eos_<suffix> = <value>` per registry family, with
   the values the hand-written block had. This tests the emitter's compound-key naming, which is
   the part that can actually break.
2. **Every `eos_*` symbol referenced anywhere under `src/` is emitted by the generator.** Grep the
   Fortran for `eos_[a-z_]+` identifiers, subtract the generated set, and require the remainder to
   be empty. This catches a family referenced in Fortran but absent from the registry — a build
   break today, and after R4 a silent one, since a missing `case` arm falls through rather than
   failing to compile. Add it to R5's list as test 5.

**Also check** `params_tests/test_fortran_gen.py:300`, which enumerates the generated files and
will see the new constants.

**R5 test 2 as first written is WRONG — corrected.** `s_initialize_eos_module` has two layers:
an UNCONDITIONAL block reading `mg_c0, mg_s, mg_s2, mg_s3, jwl_a, jwl_b, jwl_r1, jwl_r2,
vinet_k0, vinet_k0p` for every fluid regardless of family, and a `select case` that dispatches
only `rho0, t0, gruneisen0, gruneisen_a`. So "every required parameter appears in that family's
case body" would fail on `mg_c0`, `jwl_a` and the rest, which are required but read outside the
case. Correct formulation: **every required parameter of every family is read somewhere in
`s_initialize_eos_module`**. That still catches the failure the goal describes — a parameter the
validator demands but the init never reads — without false-flagging the unconditional block.

**Design refinement, and it is what actually reaches the goal's "done when".** The registry should
carry, per family, a `param -> eos_coeffs field` mapping (`mg_c0 -> c0`, `jwl_a -> a`,
`mg_rho0 -> rho0`, ...). With that, BOTH init layers are generatable: the unconditional block and
the per-family case arms. Without it, adding a tenth family still means hand-editing
`s_initialize_eos_module` twice, and the goal's success criterion is not met. This raises R1's
scope slightly (richer registry entries) and adds a generated-init step to R4, but it is the
difference between "fewer edits" and the stated target.

Note for whoever implements R4: `case default` currently sits BETWEEN `case (eos_jwl)` and
`case (eos_vinet)` in the init's second select. Legal Fortran, unusual order. A generator will
naturally emit it last; that is a behaviour-preserving change but it WILL show up in the diff, so
call it out rather than letting a reviewer wonder.

## Verified registry content (extracted from the code, not from memory)

The `param -> eos_coeffs field` mapping R1 must encode, read out of
`s_initialize_eos_module` directly:

| Family | Parameter | eos_coeffs field | Layer |
| --- | --- | --- | --- |
| mie_gruneisen | mg_c0 | c0 | unconditional |
| mie_gruneisen | mg_s | s | unconditional |
| mie_gruneisen | mg_s2 | s2 | unconditional |
| mie_gruneisen | mg_s3 | s3 | unconditional |
| mie_gruneisen | mg_rho0 | rho0 | case |
| mie_gruneisen | mg_t0 | t0 | case |
| mie_gruneisen | mg_gruneisen | gruneisen0 | case |
| mie_gruneisen | mg_gruneisen_a | gruneisen_a | case |
| jwl | jwl_a | a | unconditional |
| jwl | jwl_b | b | unconditional |
| jwl | jwl_r1 | r1 | unconditional |
| jwl | jwl_r2 | r2 | unconditional |
| jwl | jwl_rho0 | rho0 | case |
| jwl | jwl_t0 | t0 | case |
| jwl | jwl_omega | gruneisen0 | case |
| vinet | vinet_k0 | k0 | unconditional |
| vinet | vinet_k0p | k0p | unconditional |
| vinet | vinet_rho0 | rho0 | case |
| vinet | vinet_t0 | t0 | case |
| vinet | vinet_gruneisen | gruneisen0 | case |
| vinet | vinet_gruneisen_a | gruneisen_a | case |

Three things in that table are NOT a plain parameter copy, and a generator that assumes
uniformity will get them wrong:

1. **JWL has no `gruneisen_a` parameter.** Its case arm assigns the literal `0._wp`. The registry
   entry needs a constant, not a parameter name. (It is also why `f_has_isentropic_reference`
   still tests `gruneisen_a == 0` at runtime rather than trusting the family.)
2. **`mu_max` is computed, not copied** — `f_hugoniot_compression_limit(mg_c0, mg_s, mg_s2, mg_s3)`,
   evaluated unconditionally for every fluid. It is a derived field with no parameter behind it.
3. **`case default`** assigns `rho0/t0/gruneisen0 = dflt_real` and `gruneisen_a = 0._wp` — the
   non-state-dependent families' arm. The generator must emit it, and it is what makes a missing
   registry entry fail SILENTLY (dflt_real coefficients) rather than loudly. This is precisely the
   failure mode the goal describes, and it is why R5's tests matter more than the generation does.

So an `EOS_FAMILIES` entry needs, per field: either a parameter name, a literal, or a marker that
it is computed. Not a flat `param -> field` dict.
