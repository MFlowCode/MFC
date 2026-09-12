# Collapsing the EOS family registry

Date: 2026-09-12
Depends on: the `m_eos` extraction (branch `module-eos`)
Status: design

## The actual duplication

The goal named eight edit sites. Reading the code turns up sixteen, and the enum
itself is stored three times.

**Python**

| # | Site | Holds |
| --- | --- | --- |
| 1 | `definitions.py:914` `_EOS_NAMES` | name -> value |
| 2 | `definitions.py:915` `_EOS_VALUE_LABELS` | value -> display label |
| 3 | `definitions.py:923` `"choices": [1, 2, 3, 4, 5]` | the value set, a third copy |
| 4 | `definitions.py:927-956` | three per-family parameter blocks with Doxygen math symbols |
| 5 | `case_validator.py:1004-1006` `families` | value -> (prefix, required parameters) |
| 6 | `case_validator.py:1007` `optional` | prefix -> optional parameters |
| 7 | `case_validator.py:969-974` | value -> the `eos.py` coefficient function |
| 8 | `case_validator.py:1067-1073` | the state-dependent trio, hard-coded |
| 9 | `case_validator.py:2050` | the state-dependent trio again |
| 10 | `eos.py` | the Python mirror of each reference curve |

**Fortran**

| # | Site | Holds |
| --- | --- | --- |
| 11 | `m_constants.fpp:119-124` | the `eos_*` parameters |
| 12 | `m_derived_types.fpp` | the family's fields on `physical_parameters` |
| 13 | `m_eos.fpp` `s_initialize_eos_module` | `select case` mapping user parameters into `eos_coeffs` |
| 14 | `m_eos.fpp` `s_reference_curve` | `select case` computing the curve |
| 15 | `m_eos.fpp` `f_is_state_dependent` | a `.or.` chain |
| 16 | `m_eos.fpp` `f_has_isentropic_reference` | another `.or.` chain |

Nothing enforces agreement between any pair of these except one test,
`params_tests/test_eos_selector.py::test_fortran_and_python_enums_agree`, which checks
site 11 against site 1 and nothing else.

## Why the Fortran constants are hand-written

Not an oversight. `generate_constants_fpp` (`params/generators/fortran_gen.py:257`)
skips any registry key containing `%` or `(`:

```python
# Compound keys (fluid_pp(1)%eos) do not form valid Fortran identifiers, so their
# constants are hand-written in m_constants.fpp and guarded by test_eos_selector.py.
if "%" in param or "(" in param:
    continue
```

The obstacle is only the emitted *name*: the constants must be `eos_mie_gruneisen`,
not `fluid_pp(1)%eos_mie_gruneisen`. That is a one-line fix, not a redesign.

## Target

One `EOS_FAMILIES` registry in `toolchain/mfc/params/eos_families.py`, one entry per
family, from which everything mechanical is derived.

Each entry carries: the selector value; the Fortran constant suffix; the display label;
the parameter prefix (`mg`, `jwl`, `vinet`); required and optional parameters with their
Doxygen math symbols; whether the family is state-dependent; whether its reference curve
is itself an isentrope; and the name of its `eos.py` coefficient function.

Derived from it:

- **Site 3** — `choices` becomes `sorted(f.value for f in EOS_FAMILIES)`.
- **Sites 1, 2** — `_EOS_NAMES` and `_EOS_VALUE_LABELS` become comprehensions.
- **Site 4** — the three parameter blocks become one loop over the registry.
- **Sites 5, 6** — `families` and `optional` become comprehensions.
- **Sites 7, 8, 9** — the dispatch and both hard-coded trios read the registry.
- **Site 11** — `generate_constants_fpp` gains an optional `fortran_prefix` on the
  constraint entry and emits `eos_<suffix>` for compound keys instead of skipping them.
- **Sites 15, 16** — both predicates become generated `select case` bodies emitted into
  a new `generated_eos.fpp` include, consumed by `m_eos`.

## What stays hand-written

- **Site 14, `s_reference_curve`** — genuinely per-family mathematics. Its `case`
  labels come from the registry; its bodies do not.
- **Site 10, `eos.py`** — the Python mirror the tests check the solver against. A
  generated mirror would test the generator, not the solver. The registry names the
  function; a test asserts every family's name resolves.
- **Site 12, `m_derived_types.fpp`** — the `physical_parameters` fields. Generating a
  derived type's members is a larger change to the codegen than this work justifies,
  and the fields are read by name in `s_initialize_eos_module` anyway. A test asserts
  every registry parameter has a matching field.
- **Site 13, `s_initialize_eos_module`** — the mapping is per-family and short. Its
  `case` labels come from the registry.

## Closing the gap the goal names

"Validates cleanly but computes with `dflt_real` coefficients" happens when a family's
required-parameter set (site 5) and the fields the init actually reads (site 13)
disagree. With site 5 generated from the registry, a test can assert the third leg:
every required parameter of every family appears in `s_initialize_eos_module`'s body
for that family's `case`. That is a grep-level check, and it is the one that would have
caught the failure the goal describes.

## Constraints

- Bit-for-bit goldens: `./mfc.sh test` reproduces 713 passed / 0 failed / 39 skipped /
  752 total. Never `--generate`.
- `m_eos` stays a leaf: `m_derived_types`, `m_global_parameters_common`, `m_constants`.
  A generated include is textual and does not add a `use`.
- Build and test only through `./mfc.sh`; it leaves a sticky lock in `build/`.
- Short routines, terse comments.
- Separate branch off `module-eos`, separate PR, from the fork, following the template,
  stating it was made with Claude Code.

## Done when

Adding a tenth family means one `EOS_FAMILIES` entry, one `s_reference_curve` case, one
`eos.py` function, and its `physical_parameters` fields — with a test failing if any of
those four disagree. Every other site follows.
