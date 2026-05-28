# Fortran Param Codegen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Auto-generate Fortran namelist declarations and simple scalar variable declarations from `definitions.py`, reducing new-parameter additions from 7 manual file edits to ≤2.

**Architecture:** A new `fortran_gen.py` generator (fitting the existing generator pattern in `params/generators/`) reads `definitions.py` + a new `namelist_targets.py` to produce six `.fpp` include files — one namelist fragment and one declarations fragment per target (pre/sim/post). These are checked into `src/common/include/`, regenerated via `./mfc.sh generate`, and verified up-to-date by a new `precheck.sh` step.

**Tech Stack:** Python 3 (existing toolchain), Fypp preprocessor (for `#:include`), Fortran namelists, MFC's existing `generate.py` check-or-write infrastructure.

---

## Background: Current Pain

Adding one simple scalar used by all three targets currently requires edits in 7 places:
1. `toolchain/mfc/params/definitions.py`
2-4. `src/{pre_process,simulation,post_process}/m_global_parameters.fpp` (variable declaration)
5-7. `src/{pre_process,simulation,post_process}/m_start_up.fpp` (namelist entry)

After this plan: **1–2 places** — just `definitions.py` + `namelist_targets.py`.

---

## File Map

### New files
| File | Responsibility |
|------|----------------|
| `toolchain/mfc/params/namelist_targets.py` | `NAMELIST_VARS` dict (namelist_var → targets) + `CASE_OPT_EXCLUDE` set |
| `toolchain/mfc/params/generators/fortran_gen.py` | Generator: produces namelist + decls `.fpp` content per target |
| `src/common/include/generated_namelist_pre.fpp` | Generated: pre_process namelist fragment |
| `src/common/include/generated_namelist_sim.fpp` | Generated: simulation namelist fragment (with Fypp CASE_OPT guard) |
| `src/common/include/generated_namelist_post.fpp` | Generated: post_process namelist fragment |
| `src/common/include/generated_decls_pre.fpp` | Generated: pre_process simple scalar declarations |
| `src/common/include/generated_decls_sim.fpp` | Generated: simulation simple scalar declarations |
| `src/common/include/generated_decls_post.fpp` | Generated: post_process simple scalar declarations |

### Modified files
| File | Change |
|------|--------|
| `toolchain/mfc/params/schema.py` | Add `str_len: str = "name_len"` to `ParamDef` |
| `toolchain/mfc/params/definitions.py` | Set `str_len="path_len"` on `case_dir` |
| `toolchain/mfc/generate.py` | Register six new generated files |
| `toolchain/bootstrap/precheck.sh` | Add check 7/7: `./mfc.sh generate --check` |
| `src/pre_process/m_start_up.fpp:77-87` | Replace namelist block with `#:include` |
| `src/simulation/m_start_up.fpp:85-115` | Replace namelist block with `#:include` |
| `src/post_process/m_start_up.fpp:62-73` | Replace namelist block with `#:include` |
| `src/pre_process/m_global_parameters.fpp` | Remove generated scalars, add `#:include` after implicit none |
| `src/simulation/m_global_parameters.fpp` | Remove generated scalars, add `#:include` after implicit none |
| `src/post_process/m_global_parameters.fpp` | Remove generated scalars, add `#:include` after implicit none |

---

## Key Design Decisions

**`namelist_var` derivation (no new ParamDef field needed):**
```python
import re

def get_namelist_var(param_name: str) -> str:
    # fluid_pp(1)%gamma -> fluid_pp
    m = re.match(r'^([a-zA-Z_]\w*)\(\d', param_name)
    if m:
        return m.group(1)
    # bc_x%beg, lag_params%foo -> bc_x, lag_params
    if '%' in param_name:
        return param_name.split('%')[0]
    # Simple scalar: m, n, dt, model_eqns
    return param_name
```

**Simple scalars** (candidates for declaration generation) = params with no `%` or `(` in their name, with `ParamType` of INT, REAL, LOG, or STR.

**`NAMELIST_TARGETS` structure** (one source of truth for what goes where):
```python
NAMELIST_VARS: Dict[str, Set[str]] = {
    "m": {"pre", "sim", "post"},
    "bc_x": {"pre", "sim", "post"},
    "patch_icpp": {"pre"},
    "run_time_info": {"sim"},
    ...
}
CASE_OPT_EXCLUDE: Set[str] = {
    "nb", "mapped_weno", "wenoz", "teno", "wenoz_q", "weno_order",
    "num_fluids", "mhd", "relativity", "igr_order", "viscous",
    "igr_iter_solver", "igr", "igr_pres_lim", "recon_type", "muscl_order", "muscl_lim",
}
```

**Generated namelist format** (for simulation):
```fortran
! AUTO-GENERATED — do not edit. Regenerate with: ./mfc.sh generate
namelist /user_inputs/ m, n, p, dt, model_eqns, bc_x, bc_y, bc_z, &
    & fluid_pp, bub_pp, ...
#:if not MFC_CASE_OPTIMIZATION
    & nb, mapped_weno, wenoz, weno_order, ...
#:endif
    & last_var
```

**Generated declarations format**:
```fortran
! AUTO-GENERATED — do not edit. Regenerate with: ./mfc.sh generate
integer  :: model_eqns
real(wp) :: dt
logical  :: cyl_coord
character(LEN=path_len) :: case_dir
```

---

## Task 1: Extend ParamDef for STR length

**Files:**
- Modify: `toolchain/mfc/params/schema.py`
- Modify: `toolchain/mfc/params/definitions.py`

- [ ] **Step 1: Write a failing test in `toolchain/tests/params/test_schema.py`**

```python
def test_paramdef_str_len_default():
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="foo", param_type=ParamType.STR)
    assert p.str_len == "name_len"

def test_paramdef_str_len_override():
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="case_dir", param_type=ParamType.STR, str_len="path_len")
    assert p.str_len == "path_len"
```

- [ ] **Step 2: Run to verify failure**

```bash
cd /path/to/mfc
python3 -m pytest toolchain/tests/params/test_schema.py::test_paramdef_str_len_default -v
```
Expected: `AttributeError: str_len`

- [ ] **Step 3: Add `str_len` to `ParamDef` in `schema.py`**

In `toolchain/mfc/params/schema.py`, add to the `@dataclass class ParamDef:` block after `math_symbol`:
```python
str_len: str = "name_len"
# For STR type: Fortran character length constant. Default "name_len"; set "path_len" for case_dir.
```

- [ ] **Step 4: Set `str_len` on `case_dir` in `definitions.py`**

In `toolchain/mfc/params/definitions.py`, find the line:
```python
    _r("case_dir", STR)
```
Replace with:
```python
    _r("case_dir", STR, str_len="path_len")
```

This requires updating `_r()` to accept and forward `str_len`. Find `_r()` at line ~818:
```python
def _r(name, ptype, tags=None, desc=None, hint=None, math=None):
```
Replace with:
```python
def _r(name, ptype, tags=None, desc=None, hint=None, math=None, str_len=None):
    ...
    REGISTRY.register(
        ParamDef(
            ...
            str_len=str_len if str_len is not None else "name_len",
        )
    )
```

- [ ] **Step 5: Run the tests**

```bash
python3 -m pytest toolchain/tests/params/test_schema.py -v
```
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add toolchain/mfc/params/schema.py toolchain/mfc/params/definitions.py
git commit -m "feat(params): add str_len to ParamDef for STR character length"
```

---

## Task 2: Create `namelist_targets.py`

**Files:**
- Create: `toolchain/mfc/params/namelist_targets.py`

This file is the authoritative source for which Fortran namelist variable belongs to which target(s). It is derived by reading the existing three namelist blocks; after this plan, those blocks will be generated from it.

- [ ] **Step 1: Write a failing test**

Create `toolchain/tests/params/test_namelist_targets.py`:
```python
def test_common_vars_in_all_targets():
    from mfc.params.namelist_targets import NAMELIST_VARS
    for var in ["m", "n", "p", "bc_x", "bc_y", "bc_z", "model_eqns", "cyl_coord", "fluid_pp"]:
        assert {"pre", "sim", "post"}.issubset(NAMELIST_VARS.get(var, set())), \
            f"{var!r} not marked for all targets"

def test_sim_only_vars():
    from mfc.params.namelist_targets import NAMELIST_VARS
    for var in ["run_time_info", "dt", "riemann_solver", "acoustic", "probe"]:
        targets = NAMELIST_VARS.get(var, set())
        assert "sim" in targets, f"{var!r} not marked for sim"
        assert "pre" not in targets, f"{var!r} incorrectly marked for pre"
        assert "post" not in targets, f"{var!r} incorrectly marked for post"

def test_pre_only_vars():
    from mfc.params.namelist_targets import NAMELIST_VARS
    for var in ["old_grid", "old_ic", "patch_icpp", "simplex_params"]:
        targets = NAMELIST_VARS.get(var, set())
        assert "pre" in targets, f"{var!r} not marked for pre"
        assert "sim" not in targets, f"{var!r} incorrectly in sim"

def test_post_only_vars():
    from mfc.params.namelist_targets import NAMELIST_VARS
    for var in ["format", "sim_data", "lag_header", "output_partial_domain"]:
        targets = NAMELIST_VARS.get(var, set())
        assert "post" in targets, f"{var!r} not marked for post"
        assert "sim" not in targets, f"{var!r} incorrectly in sim"

def test_case_opt_exclude_vars():
    from mfc.params.namelist_targets import CASE_OPT_EXCLUDE
    for var in ["nb", "mapped_weno", "wenoz", "weno_order", "num_fluids"]:
        assert var in CASE_OPT_EXCLUDE
```

- [ ] **Step 2: Run to verify failure**

```bash
python3 -m pytest toolchain/tests/params/test_namelist_targets.py -v
```
Expected: `ModuleNotFoundError`

- [ ] **Step 3: Create `namelist_targets.py`**

Create `toolchain/mfc/params/namelist_targets.py` with the complete content:

```python
"""
Namelist target mapping for Fortran codegen.

NAMELIST_VARS maps each Fortran namelist variable (struct root or simple scalar)
to the set of MFC executables whose namelist it appears in.

CASE_OPT_EXCLUDE is the set of simulation namelist variables excluded under
MFC_CASE_OPTIMIZATION (they become compile-time constants instead).

When adding a new parameter:
  1. Add to definitions.py (type, constraints, etc.)
  2. Add the namelist root variable to NAMELIST_VARS with its target set
  3. Run ./mfc.sh generate to regenerate the .fpp files
"""

from typing import Dict, Set

# All three targets
_ALL = {"pre", "sim", "post"}
_PRE_SIM = {"pre", "sim"}
_SIM_POST = {"sim", "post"}

NAMELIST_VARS: Dict[str, Set[str]] = {
    # --- Grid (all targets) ---
    "m": _ALL,
    "n": _ALL,
    "p": _ALL,
    "cyl_coord": _ALL,
    "x_domain": {"pre", "sim"},
    "y_domain": {"pre", "sim"},
    "z_domain": {"pre", "sim"},
    "x_output": {"post"},
    "y_output": {"post"},
    "z_output": {"post"},
    # --- Grid stretching (pre + sim) ---
    "stretch_x": {"pre"},
    "stretch_y": {"pre"},
    "stretch_z": {"pre"},
    "a_x": {"pre"},
    "a_y": {"pre"},
    "a_z": {"pre"},
    "x_a": {"pre", "sim"},
    "y_a": {"pre", "sim"},
    "z_a": {"pre", "sim"},
    "x_b": {"pre", "sim"},
    "y_b": {"pre", "sim"},
    "z_b": {"pre", "sim"},
    "loops_x": {"pre"},
    "loops_y": {"pre"},
    "loops_z": {"pre"},
    # --- Time (sim) ---
    "dt": {"sim"},
    "t_step_start": _ALL,
    "t_step_stop": {"sim", "post"},
    "t_step_save": {"sim", "post"},
    "t_step_print": {"sim"},
    "t_step_old": {"pre", "sim"},
    "time_stepper": {"sim"},
    "t_stop": {"sim", "post"},
    "t_save": {"sim", "post"},
    "cfl_target": {"sim", "post"},
    "cfl_adap_dt": _ALL,
    "cfl_const_dt": _ALL,
    "n_start": _ALL,
    "n_start_old": {"pre"},
    "adap_dt": {"sim"},
    "adap_dt_tol": {"sim"},
    "adap_dt_max_iters": {"sim"},
    # --- Physics model (all targets) ---
    "model_eqns": _ALL,
    "num_fluids": {"pre", "post"},
    "mpp_lim": _ALL,
    "relax": _ALL,
    "relax_model": _ALL,
    "palpha_eps": _ALL,
    "ptgalpha_eps": _ALL,
    # --- WENO / reconstruction (pre has weno_order; sim has it under CASE_OPT) ---
    "weno_order": {"pre", "post"},
    "weno_eps": {"sim"},
    "teno_CT": {"sim"},
    "wenoz_q": {"sim"},
    "mp_weno": {"sim"},
    "weno_avg": {"sim"},
    "weno_Re_flux": {"sim"},
    "null_weights": {"sim"},
    "muscl_eps": {"sim"},
    "recon_type": {"pre", "post"},
    "muscl_order": {"pre", "post"},
    "muscl_lim": {"post"},
    "int_comp": {"sim"},
    "ic_eps": {"sim"},
    "ic_beta": {"sim"},
    # --- Riemann solver (sim only) ---
    "riemann_solver": {"sim"},
    "wave_speeds": {"sim"},
    "avg_state": {"sim", "post"},
    "low_Mach": {"sim"},
    # --- MHD (all targets) ---
    "mhd": {"pre", "post"},
    "hyper_cleaning": _ALL,
    "hyper_cleaning_speed": {"sim"},
    "hyper_cleaning_tau": {"sim"},
    "Bx0": _ALL,
    # --- BCs (all targets) ---
    "bc_x": _ALL,
    "bc_y": _ALL,
    "bc_z": _ALL,
    "num_bc_patches": _ALL,
    "patch_bc": {"pre"},
    # --- ICs (pre only) ---
    "num_patches": {"pre"},
    "patch_icpp": {"pre"},
    # --- Fluid properties (all) ---
    "fluid_pp": _ALL,
    "bub_pp": _ALL,
    "rhoref": _ALL,
    "pref": _ALL,
    # --- Bubbles (all) ---
    "bubbles_euler": _ALL,
    "bubbles_lagrange": _ALL,
    "R0ref": _ALL,
    "nb": {"pre", "post"},
    "polytropic": _ALL,
    "thermal": _ALL,
    "Ca": _ALL,
    "Web": _ALL,
    "Re_inv": _ALL,
    "polydisperse": _ALL,
    "poly_sigma": _ALL,
    "qbmm": _ALL,
    "sigma": _ALL,
    "adv_n": _ALL,
    "bubble_model": {"sim"},
    "sigR": {"pre", "post"},
    "sigV": {"pre"},
    "dist_type": {"pre"},
    "rhoRV": {"pre"},
    "lag_params": {"sim"},
    # --- Lagrangian output (post) ---
    "lag_header": {"post"},
    "lag_txt_wrt": {"post"},
    "lag_db_wrt": {"post"},
    "lag_id_wrt": {"post"},
    "lag_pos_wrt": {"post"},
    "lag_pos_prev_wrt": {"post"},
    "lag_vel_wrt": {"post"},
    "lag_rad_wrt": {"post"},
    "lag_rvel_wrt": {"post"},
    "lag_r0_wrt": {"post"},
    "lag_rmax_wrt": {"post"},
    "lag_rmin_wrt": {"post"},
    "lag_dphidt_wrt": {"post"},
    "lag_pres_wrt": {"post"},
    "lag_mv_wrt": {"post"},
    "lag_mg_wrt": {"post"},
    "lag_betaT_wrt": {"post"},
    "lag_betaC_wrt": {"post"},
    # --- Elasticity (all) ---
    "hypoelasticity": _ALL,
    "hyperelasticity": _ALL,
    # --- Surface tension (all) ---
    "surface_tension": _ALL,
    # --- Relativity (all) ---
    "relativity": _ALL,
    # --- Immersed boundaries (all) ---
    "ib": _ALL,
    "num_ibs": _ALL,
    "patch_ib": {"pre", "sim"},
    "collision_model": {"sim"},
    "coefficient_of_restitution": {"sim"},
    "collision_time": {"sim"},
    "ib_coefficient_of_friction": {"sim"},
    "ib_state_wrt": {"sim", "post"},
    # --- Continuum damage (all) ---
    "cont_damage": _ALL,
    "tau_star": {"sim"},
    "cont_damage_s": {"sim"},
    "alpha_bar": {"sim"},
    # --- IGR (all) ---
    "igr": {"pre", "post"},
    "igr_order": {"pre", "post"},
    "down_sample": _ALL,
    # --- Probes (sim) ---
    "probe_wrt": {"sim"},
    "num_probes": {"sim"},
    "probe": {"sim"},
    "integral_wrt": {"sim"},
    "num_integrals": {"sim"},
    "integral": {"sim"},
    "fd_order": {"sim", "post"},
    # --- Acoustic sources (sim) ---
    "acoustic_source": {"sim"},
    "num_source": {"sim"},
    "acoustic": {"sim"},
    # --- Chemistry (sim) ---
    "chem_params": {"sim"},
    # --- Body forces (sim) ---
    "bf_x": {"sim"},
    "bf_y": {"sim"},
    "bf_z": {"sim"},
    "k_x": {"sim"},
    "k_y": {"sim"},
    "k_z": {"sim"},
    "w_x": {"sim"},
    "w_y": {"sim"},
    "w_z": {"sim"},
    "p_x": {"sim"},
    "p_y": {"sim"},
    "p_z": {"sim"},
    "g_x": {"sim"},
    "g_y": {"sim"},
    "g_z": {"sim"},
    # --- Viscous (pre) ---
    "viscous": {"pre"},
    # --- Output (all) ---
    "precision": _ALL,
    "parallel_io": _ALL,
    "file_per_process": _ALL,
    "prim_vars_wrt": {"sim", "post"},
    "cons_vars_wrt": {"post"},
    "run_time_info": {"sim"},
    "fft_wrt": _ALL,
    "pi_fac": {"pre", "sim"},
    # --- Post-process output ---
    "format": {"post"},
    "output_partial_domain": {"post"},
    "rho_wrt": {"post"},
    "E_wrt": {"post"},
    "pres_wrt": {"post"},
    "c_wrt": {"post"},
    "omega_wrt": {"post"},
    "qm_wrt": {"post"},
    "liutex_wrt": {"post"},
    "schlieren_wrt": {"post"},
    "schlieren_alpha": {"post"},
    "gamma_wrt": {"post"},
    "heat_ratio_wrt": {"post"},
    "pi_inf_wrt": {"post"},
    "pres_inf_wrt": {"post"},
    "alpha_rho_wrt": {"post"},
    "mom_wrt": {"post"},
    "vel_wrt": {"post"},
    "flux_wrt": {"post"},
    "alpha_wrt": {"post"},
    "cf_wrt": {"post"},
    "chem_wrt_T": {"post"},
    "chem_wrt_Y": {"post"},
    "alt_soundspeed": {"sim", "post"},
    "mixture_err": {"sim", "post"},
    "flux_lim": {"post"},
    "sim_data": {"post"},
    "alpha_rho_e_wrt": {"post"},
    "G": {"post"},
    # --- Pre-process IC perturbations ---
    "perturb_flow": {"pre"},
    "perturb_flow_fluid": {"pre"},
    "perturb_flow_mag": {"pre"},
    "perturb_sph": {"pre"},
    "perturb_sph_fluid": {"pre"},
    "fluid_rho": {"pre"},
    "mixlayer_vel_profile": {"pre"},
    "mixlayer_vel_coef": {"pre"},
    "mixlayer_perturb": {"pre"},
    "mixlayer_perturb_nk": {"pre"},
    "mixlayer_perturb_k0": {"pre"},
    "pre_stress": {"pre"},
    "elliptic_smoothing": {"pre"},
    "elliptic_smoothing_iters": {"pre"},
    "simplex_perturb": {"pre"},
    "simplex_params": {"pre"},
    # --- Pre-process restart ---
    "old_grid": {"pre"},
    "old_ic": {"pre"},
    # --- Sim-specific physics ---
    "rdma_mpi": {"sim"},
    "alf_factor": {"sim"},
    "num_igr_iters": {"sim"},
    "num_igr_warm_start_iters": {"sim"},
    "igr_iter_solver": {"sim"},
    "igr_pres_lim": {"sim"},
    "nv_uvm_out_of_core": {"sim"},
    "nv_uvm_igr_temps_on_gpu": {"sim"},
    "nv_uvm_pref_gpu": {"sim"},
    # --- Logistics ---
    "case_dir": _ALL,
}

# Variables excluded from the sim namelist when MFC_CASE_OPTIMIZATION is active
# (they become compile-time integer/logical parameters instead).
# Must all be present in NAMELIST_VARS with "sim" in their target set.
CASE_OPT_EXCLUDE: Set[str] = {
    "nb",
    "mapped_weno",
    "wenoz",
    "teno",
    "wenoz_q",
    "weno_order",
    "num_fluids",
    "mhd",
    "relativity",
    "igr_order",
    "viscous",
    "igr_iter_solver",
    "igr",
    "igr_pres_lim",
    "recon_type",
    "muscl_order",
    "muscl_lim",
}

# Add CASE_OPT_EXCLUDE vars to NAMELIST_VARS for sim target
# (they appear in the namelist when NOT using case optimization)
for _v in CASE_OPT_EXCLUDE:
    if _v not in NAMELIST_VARS:
        NAMELIST_VARS[_v] = {"sim"}
    else:
        NAMELIST_VARS[_v].add("sim")
```

- [ ] **Step 4: Run tests**

```bash
python3 -m pytest toolchain/tests/params/test_namelist_targets.py -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/params/namelist_targets.py toolchain/tests/params/test_namelist_targets.py
git commit -m "feat(params): add namelist_targets.py with target mapping for all namelist vars"
```

---

## Task 3: Write `fortran_gen.py`

**Files:**
- Create: `toolchain/mfc/params/generators/fortran_gen.py`
- Test: `toolchain/tests/params/test_fortran_gen.py`

The generator reads `definitions.py` (via `REGISTRY`) + `namelist_targets.py` and produces:
- **Namelist fragment**: `namelist /user_inputs/ var1, var2, ...` with Fypp case-opt guard for sim
- **Decls fragment**: `integer :: foo` / `real(wp) :: bar` for each simple scalar in that target

- [ ] **Step 1: Write failing tests**

Create `toolchain/tests/params/test_fortran_gen.py`:

```python
import pytest


def test_get_namelist_var_simple():
    from mfc.params.generators.fortran_gen import get_namelist_var
    assert get_namelist_var("m") == "m"
    assert get_namelist_var("dt") == "dt"
    assert get_namelist_var("cyl_coord") == "cyl_coord"


def test_get_namelist_var_indexed_family():
    from mfc.params.generators.fortran_gen import get_namelist_var
    assert get_namelist_var("fluid_pp(1)%gamma") == "fluid_pp"
    assert get_namelist_var("patch_icpp(3)%geometry") == "patch_icpp"
    assert get_namelist_var("patch_ib(100)%radius") == "patch_ib"


def test_get_namelist_var_struct_member():
    from mfc.params.generators.fortran_gen import get_namelist_var
    assert get_namelist_var("bc_x%beg") == "bc_x"
    assert get_namelist_var("bc_y%grcbc_in") == "bc_y"
    assert get_namelist_var("lag_params%solver_approach") == "lag_params"
    assert get_namelist_var("chem_params%diffusion") == "chem_params"
    assert get_namelist_var("bub_pp%R0ref") == "bub_pp"


def test_fortran_type_for_int():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="model_eqns", param_type=ParamType.INT)
    assert fortran_type_decl(p) == "integer"


def test_fortran_type_for_real():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="dt", param_type=ParamType.REAL)
    assert fortran_type_decl(p) == "real(wp)"


def test_fortran_type_for_log():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="cyl_coord", param_type=ParamType.LOG)
    assert fortran_type_decl(p) == "logical"


def test_fortran_type_for_str():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType
    p = ParamDef(name="case_dir", param_type=ParamType.STR, str_len="path_len")
    assert fortran_type_decl(p) == "character(LEN=path_len)"


def test_generate_namelist_contains_common_vars():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp
    for target in ("pre", "sim", "post"):
        content = generate_namelist_fpp(target)
        for var in ("m", "n", "p", "bc_x", "bc_y", "bc_z", "fluid_pp", "case_dir"):
            assert var in content, f"{var!r} missing from {target} namelist"


def test_sim_namelist_has_case_opt_guard():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp
    content = generate_namelist_fpp("sim")
    assert "#:if not MFC_CASE_OPTIMIZATION" in content
    assert "weno_order" in content
    assert "num_fluids" in content


def test_pre_namelist_has_patch_icpp():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp
    content = generate_namelist_fpp("pre")
    assert "patch_icpp" in content
    assert "run_time_info" not in content


def test_post_namelist_has_output_vars():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp
    content = generate_namelist_fpp("post")
    assert "sim_data" in content
    assert "lag_header" in content
    assert "patch_icpp" not in content


def test_generate_decls_contains_simple_scalars():
    from mfc.params.generators.fortran_gen import generate_decls_fpp
    for target in ("pre", "sim", "post"):
        content = generate_decls_fpp(target)
        assert "integer" in content
        assert "real(wp)" in content
        assert "logical" in content


def test_generate_decls_has_dt_for_sim():
    from mfc.params.generators.fortran_gen import generate_decls_fpp
    content = generate_decls_fpp("sim")
    assert "real(wp) :: dt" in content


def test_generate_decls_no_percent_vars():
    from mfc.params.generators.fortran_gen import generate_decls_fpp
    for target in ("pre", "sim", "post"):
        content = generate_decls_fpp(target)
        # Derived type members are NOT declared — only their struct root
        assert "bc_x%beg" not in content
        assert "fluid_pp(1)" not in content


def test_generate_decls_has_case_dir():
    from mfc.params.generators.fortran_gen import generate_decls_fpp
    for target in ("pre", "sim", "post"):
        content = generate_decls_fpp(target)
        assert "character(LEN=path_len) :: case_dir" in content
```

- [ ] **Step 2: Run to verify failure**

```bash
python3 -m pytest toolchain/tests/params/test_fortran_gen.py -v
```
Expected: `ModuleNotFoundError`

- [ ] **Step 3: Implement `fortran_gen.py`**

Create `toolchain/mfc/params/generators/fortran_gen.py`:

```python
"""
Fortran parameter code generator.

Generates namelist fragments and simple scalar declaration fragments from
definitions.py + namelist_targets.py. Output goes to src/common/include/.

Usage: called from generate.py. Not invoked directly.
"""

import re
from typing import Optional

from ..namelist_targets import CASE_OPT_EXCLUDE, NAMELIST_VARS
from ..registry import REGISTRY
from ..schema import ParamDef, ParamType

import mfc.params.definitions  # noqa: F401 — triggers registry population


_HEADER = """\
! AUTO-GENERATED — do not edit directly.
! Source: toolchain/mfc/params/definitions.py + namelist_targets.py
! Regenerate: ./mfc.sh generate
!
"""

_FORTRAN_INDENT = "    "


def get_namelist_var(param_name: str) -> str:
    """Return the Fortran namelist variable root for a parameter name.

    fluid_pp(1)%gamma -> fluid_pp
    bc_x%beg          -> bc_x
    m                 -> m
    """
    # Indexed family: fluid_pp(1)%gamma -> fluid_pp
    m = re.match(r"^([a-zA-Z_]\w*)\(", param_name)
    if m and "%" in param_name:
        return m.group(1)
    # Struct member: bc_x%beg -> bc_x
    if "%" in param_name:
        return param_name.split("%")[0]
    # Simple scalar
    return param_name


def fortran_type_decl(param: ParamDef) -> str:
    """Return the Fortran type keyword for a parameter."""
    mapping = {
        ParamType.INT: "integer",
        ParamType.REAL: "real(wp)",
        ParamType.LOG: "logical",
        ParamType.ANALYTIC_INT: "integer",
        ParamType.ANALYTIC_REAL: "real(wp)",
    }
    if param.param_type == ParamType.STR:
        return f"character(LEN={param.str_len})"
    return mapping[param.param_type]


def _is_simple_scalar(param_name: str) -> bool:
    """Return True if this param is a simple scalar (no % or () in name)."""
    return "%" not in param_name and "(" not in param_name


def _namelist_vars_for_target(target: str) -> list:
    """Return ordered list of namelist variables for a target."""
    result = []
    seen = set()
    for var, targets in NAMELIST_VARS.items():
        if target in targets and var not in seen:
            result.append(var)
            seen.add(var)
    return result


def generate_namelist_fpp(target: str) -> str:
    """Generate the namelist /user_inputs/ fragment for a target.

    For sim, wraps CASE_OPT_EXCLUDE vars in a Fypp guard.
    """
    assert target in ("pre", "sim", "post"), f"Unknown target: {target!r}"

    all_vars = _namelist_vars_for_target(target)

    if target == "sim":
        normal_vars = [v for v in all_vars if v not in CASE_OPT_EXCLUDE]
        opt_vars = [v for v in all_vars if v in CASE_OPT_EXCLUDE]
    else:
        normal_vars = all_vars
        opt_vars = []

    lines = [_HEADER]

    def _fmt_varlist(varlist: list, continuation: bool = False) -> list:
        """Format a list of variable names as Fortran continuation lines."""
        out = []
        chunk = []
        for v in varlist:
            chunk.append(v)
            if len(", ".join(chunk)) > 70:
                prefix = f"{_FORTRAN_INDENT}& " if continuation else f"{_FORTRAN_INDENT}"
                out.append(prefix + ", ".join(chunk[:-1]) + ", &")
                chunk = [chunk[-1]]
                continuation = True
        if chunk:
            prefix = f"{_FORTRAN_INDENT}& " if continuation else f"{_FORTRAN_INDENT}"
            out.append(prefix + ", ".join(chunk))
        return out

    if not opt_vars:
        var_lines = _fmt_varlist(normal_vars)
        lines.append(f"namelist /user_inputs/ {var_lines[0].lstrip()}")
        if len(var_lines) > 1:
            lines[-1] += ", &"
            lines.extend(v + (", &" if i < len(var_lines) - 2 else "") for i, v in enumerate(var_lines[1:]))
    else:
        # Sim: normal vars first, then CASE_OPT guard, then trailing vars
        # Emit normal vars up to end, then the guard block
        var_lines = _fmt_varlist(normal_vars)
        lines.append(f"namelist /user_inputs/ {var_lines[0].lstrip()}, &")
        for v in var_lines[1:]:
            lines.append(v + ", &")

        lines.append(f"#:if not MFC_CASE_OPTIMIZATION")
        opt_lines = _fmt_varlist(opt_vars, continuation=True)
        for i, v in enumerate(opt_lines):
            lines.append(v + (", &" if i < len(opt_lines) - 1 else ", &"))
        lines.append(f"#:endif")
        # Need a placeholder last entry — use a trailing comment-safe idiom
        # Actually Fortran namelist doesn't need a sentinel, just remove trailing comma
        # Fix: the last normal line had a trailing ", &" — correct the last line
        # Remove trailing ", &" from the last real entry
        lines[-1] = lines[-1].rstrip(", &")

    return "\n".join(lines) + "\n"


def generate_decls_fpp(target: str) -> str:
    """Generate simple scalar Fortran declarations for a target.

    Only generates declarations for params that are:
    - Simple scalars (no % or () in name)
    - In NAMELIST_VARS for this target
    - Registered in definitions.py
    """
    assert target in ("pre", "sim", "post"), f"Unknown target: {target!r}"

    vars_for_target = set(_namelist_vars_for_target(target))
    lines = [_HEADER]

    all_params = REGISTRY.all_params
    for name in sorted(vars_for_target):
        if not _is_simple_scalar(name):
            continue
        param = all_params.get(name)
        if param is None:
            continue
        ftype = fortran_type_decl(param)
        lines.append(f"{ftype} :: {name}")

    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run tests**

```bash
python3 -m pytest toolchain/tests/params/test_fortran_gen.py -v
```
Expected: all PASS (fix any off-by-one in namelist formatting if needed)

- [ ] **Step 5: Commit**

```bash
git add toolchain/mfc/params/generators/fortran_gen.py toolchain/tests/params/test_fortran_gen.py
git commit -m "feat(params): add fortran_gen.py to generate namelist and declaration .fpp files"
```

---

## Task 4: Wire into `generate.py`

**Files:**
- Modify: `toolchain/mfc/generate.py`

The generator produces six files in `src/common/include/`. They are checked in and verified by `generate --check`.

- [ ] **Step 1: Write a failing test**

In `toolchain/tests/params/test_fortran_gen.py`, add:

```python
def test_generate_function_returns_six_paths():
    from mfc.params.generators.fortran_gen import get_generated_files
    from pathlib import Path
    files = get_generated_files()
    assert len(files) == 6
    names = {p.name for p, _ in files}
    assert "generated_namelist_pre.fpp" in names
    assert "generated_decls_sim.fpp" in names
```

- [ ] **Step 2: Add `get_generated_files()` to `fortran_gen.py`**

At the bottom of `fortran_gen.py`, add:

```python
from pathlib import Path
from ..common import MFC_ROOT_DIR  # already used elsewhere in the toolchain


def get_generated_files() -> list:
    """Return list of (output_path, content) tuples for all six generated files."""
    include_dir = Path(MFC_ROOT_DIR) / "src" / "common" / "include"
    return [
        (include_dir / "generated_namelist_pre.fpp",  generate_namelist_fpp("pre")),
        (include_dir / "generated_namelist_sim.fpp",  generate_namelist_fpp("sim")),
        (include_dir / "generated_namelist_post.fpp", generate_namelist_fpp("post")),
        (include_dir / "generated_decls_pre.fpp",     generate_decls_fpp("pre")),
        (include_dir / "generated_decls_sim.fpp",     generate_decls_fpp("sim")),
        (include_dir / "generated_decls_post.fpp",    generate_decls_fpp("post")),
    ]
```

Note: `MFC_ROOT_DIR` is defined in `toolchain/mfc/common.py`. Check the import path before using.

- [ ] **Step 3: Register in `generate.py`**

In `toolchain/mfc/generate.py`, in the `generate()` function, add after the existing `files = [...]` list:

```python
    from .params.generators.fortran_gen import get_generated_files
    files += get_generated_files()
```

- [ ] **Step 4: Run `./mfc.sh generate` to produce the six files**

```bash
./mfc.sh generate
```
Expected: six new/updated `.fpp` files in `src/common/include/`

Inspect the output:
```bash
cat src/common/include/generated_namelist_sim.fpp | head -30
cat src/common/include/generated_decls_pre.fpp | head -20
```

- [ ] **Step 5: Verify content matches the existing namelists**

Manually cross-check that every variable in the current `m_start_up.fpp` namelist blocks appears in the corresponding generated file. Key checks:
- `generated_namelist_sim.fpp` should contain `weno_order` inside the `#:if not MFC_CASE_OPTIMIZATION` guard
- `generated_namelist_pre.fpp` should contain `patch_icpp`, `simplex_params`, `old_grid`
- `generated_namelist_post.fpp` should contain `sim_data`, `lag_header`, `format`

- [ ] **Step 6: Run tests**

```bash
python3 -m pytest toolchain/tests/params/test_fortran_gen.py -v
```
Expected: all PASS

- [ ] **Step 7: Commit the generated files and generate.py change**

```bash
git add src/common/include/generated_namelist_*.fpp src/common/include/generated_decls_*.fpp
git add toolchain/mfc/generate.py toolchain/mfc/params/generators/fortran_gen.py
git commit -m "feat(params): generate Fortran namelist and declaration .fpp files from definitions.py"
```

---

## Task 5: Replace namelists in `m_start_up.fpp` (all three targets)

**Files:**
- Modify: `src/pre_process/m_start_up.fpp:77-87`
- Modify: `src/simulation/m_start_up.fpp:85-115`
- Modify: `src/post_process/m_start_up.fpp:62-73`

This is the first Fortran change. Do one target at a time and build between each to catch errors early.

### Pre-process

- [ ] **Step 1: In `src/pre_process/m_start_up.fpp`, replace lines 77–87**

Current (lines 77–87):
```fortran
        namelist /user_inputs/ case_dir, old_grid, old_ic, t_step_old, t_step_start, m, n, p, x_domain, y_domain, z_domain, &
            & stretch_x, stretch_y, stretch_z, a_x, a_y, a_z, x_a, y_a, z_a, x_b, y_b, z_b, model_eqns, num_fluids, mpp_lim, &
            ...
            & simplex_perturb, simplex_params, fft_wrt
```

Replace with:
```fortran
        #:include 'generated_namelist_pre.fpp'
```

- [ ] **Step 2: Build pre_process only**

```bash
./mfc.sh build -t pre_process -j 8
```
Expected: successful compilation. If there are "Undefined variable" errors, a namelist variable is missing from `generated_namelist_pre.fpp` — add it to `NAMELIST_VARS` in `namelist_targets.py` with target `"pre"`, then re-run `./mfc.sh generate`.

### Simulation

- [ ] **Step 3: In `src/simulation/m_start_up.fpp`, replace lines 85–115**

Current (lines 85–115):
```fortran
        namelist /user_inputs/ case_dir, run_time_info, m, n, p, dt, &
            t_step_start, t_step_stop, t_step_save, t_step_print, &
            ...
            & int_comp, ic_eps, ic_beta, nv_uvm_out_of_core, nv_uvm_igr_temps_on_gpu, nv_uvm_pref_gpu, down_sample, fft_wrt
```

Replace with:
```fortran
        #:include 'generated_namelist_sim.fpp'
```

- [ ] **Step 4: Build simulation**

```bash
./mfc.sh build -t simulation -j 8
```
Expected: successful. Errors → check NAMELIST_VARS for the missing variable.

### Post-process

- [ ] **Step 5: In `src/post_process/m_start_up.fpp`, replace lines 62–73**

Current (lines 62–73):
```fortran
        namelist /user_inputs/ case_dir, m, n, p, t_step_start, t_step_stop, t_step_save, model_eqns, num_fluids, mpp_lim, &
            ...
            & alpha_rho_e_wrt, ib_state_wrt
```

Replace with:
```fortran
        #:include 'generated_namelist_post.fpp'
```

- [ ] **Step 6: Build post_process**

```bash
./mfc.sh build -t post_process -j 8
```
Expected: successful.

- [ ] **Step 7: Build all three targets together**

```bash
./mfc.sh build -j 8
```
Expected: all three compile clean.

- [ ] **Step 8: Run a subset of tests to verify runtime behavior**

```bash
./mfc.sh test --only 1D -j 8
```
Expected: all pass.

- [ ] **Step 9: Commit**

```bash
git add src/pre_process/m_start_up.fpp src/simulation/m_start_up.fpp src/post_process/m_start_up.fpp
git commit -m "refactor(fortran): replace hand-written namelists with generated includes"
```

---

## Task 6: Replace declarations in `m_global_parameters.fpp`

**Files:**
- Modify: `src/simulation/m_global_parameters.fpp`
- Modify: `src/pre_process/m_global_parameters.fpp`
- Modify: `src/post_process/m_global_parameters.fpp`

### Strategy

For each target:
1. Add `#:include 'generated_decls_{target}.fpp'` right after `implicit none`
2. Remove every declaration line for a variable that will be generated
3. Build and fix errors (duplicate declaration → remove from file; missing declaration → add to generator)

**Which declarations get removed?** Any `integer :: foo`, `real(wp) :: bar`, `logical :: baz`, `character(...) :: qux` line where `foo`/`bar`/`baz`/`qux` is a simple scalar in `NAMELIST_VARS` for this target.

**Which stay?** Everything else: derived types, allocatables, GPU_DECLARE macros, internal variables not in definitions.py, case optimization conditionals.

### Simulation (`m_global_parameters.fpp`)

- [ ] **Step 1: Add include after `implicit none` (line ~19)**

After line `implicit none` in `src/simulation/m_global_parameters.fpp`, add:
```fortran
    #:include 'generated_decls_sim.fpp'
```

- [ ] **Step 2: Build — expect duplicate declaration errors**

```bash
./mfc.sh build -t simulation -j 8 2>&1 | grep -i "duplicate\|redeclared\|already"
```

- [ ] **Step 3: For each duplicate, remove the declaration from the file**

Systematically remove every line of the form `integer :: foo`, `real(wp) :: bar`, `logical :: baz` where `foo`/`bar`/`baz` is in `generated_decls_sim.fpp`. Leave compound lines (e.g., `integer :: m, n, p`) if only some of those vars are generated — split them if needed.

Also remove the Fypp loop at lines ~183-187 that generates `k_x, w_x, ...` etc., since those are now in `generated_decls_sim.fpp`.

- [ ] **Step 4: Build simulation cleanly**

```bash
./mfc.sh build -t simulation -j 8
```
Expected: no errors.

- [ ] **Step 5: Run simulation tests**

```bash
./mfc.sh test --only 1D -j 8
```
Expected: pass.

### Pre-process (`m_global_parameters.fpp`)

- [ ] **Step 6: Add include after `implicit none` in pre_process file**

```fortran
    #:include 'generated_decls_pre.fpp'
```

- [ ] **Step 7: Build and remove duplicates**

```bash
./mfc.sh build -t pre_process -j 8 2>&1 | grep -i "duplicate\|redeclared"
```
Remove listed duplicates from the file.

- [ ] **Step 8: Build cleanly**

```bash
./mfc.sh build -t pre_process -j 8
```

### Post-process (`m_global_parameters.fpp`)

- [ ] **Step 9: Add include after `implicit none` in post_process file**

```fortran
    #:include 'generated_decls_post.fpp'
```

- [ ] **Step 10: Build and remove duplicates**

```bash
./mfc.sh build -t post_process -j 8 2>&1 | grep -i "duplicate\|redeclared"
```
Remove listed duplicates from the file.

- [ ] **Step 11: Build all three**

```bash
./mfc.sh build -j 8
```
Expected: clean.

- [ ] **Step 12: Run full 1D + 2D tests**

```bash
./mfc.sh test --only 1D 2D -j 8
```
Expected: pass.

- [ ] **Step 13: Commit**

```bash
git add src/simulation/m_global_parameters.fpp src/pre_process/m_global_parameters.fpp src/post_process/m_global_parameters.fpp
git commit -m "refactor(fortran): replace hand-written scalar decls with generated includes"
```

---

## Task 7: Add precheck step

**Files:**
- Modify: `toolchain/bootstrap/precheck.sh`

- [ ] **Step 1: Add check 7/7 to `precheck.sh`**

In `toolchain/bootstrap/precheck.sh`, find the block for check 6/6 (parameter docs). After the block that waits for `PID_PARAM_DOCS`, add a new parallel job for the generate check:

```bash
# Generated Fortran files up-to-date check
(
    if ./mfc.sh generate --check > /dev/null 2>&1; then
        echo "0" > "$TMPDIR_PC/generate_exit"
    else
        echo "1" > "$TMPDIR_PC/generate_exit"
    fi
) &
PID_GENERATE=$!
```

Then update the results-collection section. Change `6/6` references to `7/7` for the last check, and change the existing `6/6` block:

Before the existing `6/6` log line:
```bash
log "[$CYAN 6/6$COLOR_RESET] Checking$MAGENTA parameter docs$COLOR_RESET..."
```
Change to:
```bash
log "[$CYAN 6/7$COLOR_RESET] Checking$MAGENTA parameter docs$COLOR_RESET..."
```

Then add after the parameter docs block:

```bash
wait $PID_GENERATE
log "[$CYAN 7/7$COLOR_RESET] Checking$MAGENTA generated Fortran files$COLOR_RESET..."
GENERATE_RC=$(cat "$TMPDIR_PC/generate_exit" 2>/dev/null || echo "1")
if [ "$GENERATE_RC" = "0" ]; then
    ok "Generated Fortran files are up to date."
else
    error "Generated Fortran files are out of date. Run$MAGENTA ./mfc.sh generate$COLOR_RESET to update."
    FAILED=1
fi
```

Also update the first log line from `(same checks as CI lint-gate)` wording and `1/6` references to `1/7`.

- [ ] **Step 2: Test the new check**

```bash
./mfc.sh precheck -j 8
```
Expected: `[7/7] Checking generated Fortran files... OK`

- [ ] **Step 3: Verify it catches a stale file**

Edit one line in `src/common/include/generated_namelist_sim.fpp` (add a comment), then run:
```bash
./mfc.sh precheck -j 8
```
Expected: `[7/7] Checking generated Fortran files... FAIL`. Revert the test change.

- [ ] **Step 4: Commit**

```bash
git add toolchain/bootstrap/precheck.sh
git commit -m "feat(precheck): add generated Fortran files up-to-date check (7/7)"
```

---

## Task 8: Final verification — full test suite

- [ ] **Step 1: Run precheck**

```bash
./mfc.sh precheck -j 8
```
Expected: all 7/7 checks pass.

- [ ] **Step 2: Build with case optimization to verify CASE_OPT guard works**

```bash
./mfc.sh build -t simulation --case-optimization -i examples/3d_taylor_green_vortex/case.py -j 8
```
Expected: successful compilation. The `#:if not MFC_CASE_OPTIMIZATION` guard in `generated_namelist_sim.fpp` must exclude the case-opt vars from the namelist at compile time.

- [ ] **Step 3: Run full test suite**

```bash
./mfc.sh test -j 8
```
Expected: same pass/fail count as before this change (zero new failures).

- [ ] **Step 4: Verify new-parameter workflow end-to-end**

Add a dummy parameter to prove the 2-location workflow works:

In `definitions.py`, add:
```python
_r("test_codegen_param", INT)
```

In `namelist_targets.py`, add to `NAMELIST_VARS`:
```python
"test_codegen_param": {"sim"},
```

Run `./mfc.sh generate` and verify `generated_namelist_sim.fpp` and `generated_decls_sim.fpp` contain `test_codegen_param`.

Then revert both additions (this was just a smoke test):
```bash
git diff toolchain/mfc/params/definitions.py toolchain/mfc/params/namelist_targets.py
git checkout -- toolchain/mfc/params/definitions.py toolchain/mfc/params/namelist_targets.py
./mfc.sh generate
```

---

## Test Plan

### Unit tests (run via `./mfc.sh lint`)

| Test file | What it covers |
|-----------|----------------|
| `toolchain/tests/params/test_schema.py` | `str_len` field on `ParamDef` |
| `toolchain/tests/params/test_namelist_targets.py` | `NAMELIST_VARS` coverage per target; `CASE_OPT_EXCLUDE` |
| `toolchain/tests/params/test_fortran_gen.py` | `get_namelist_var`, `fortran_type_decl`, namelist/decl content per target |

### Integration tests

| Test | Command | Pass criterion |
|------|---------|----------------|
| Build all targets | `./mfc.sh build -j 8` | Zero compile errors |
| Case-optimized build | `./mfc.sh build -t simulation --case-optimization -i examples/3d_taylor_green_vortex/case.py -j 8` | Compiles; CASE_OPT vars absent from runtime namelist |
| 1D regression tests | `./mfc.sh test --only 1D -j 8` | Same pass count as before |
| 2D regression tests | `./mfc.sh test --only 2D -j 8` | Same pass count as before |
| Full test suite | `./mfc.sh test -j 8` | Zero new failures |

### CI / precheck tests

| Test | Command | Pass criterion |
|------|---------|----------------|
| Generated files up to date | `./mfc.sh generate --check` | Exit 0 |
| Full precheck | `./mfc.sh precheck -j 8` | All 7/7 checks pass |
| Stale file detection | Modify a generated `.fpp`, run `./mfc.sh generate --check` | Exit 1 |

### Regression proof: new-parameter workflow

**Before**: adding `test_codegen_param` to simulation only required editing 3 files (definitions.py, simulation/m_global_parameters.fpp, simulation/m_start_up.fpp).

**After**: editing only `definitions.py` + `namelist_targets.py` + running `./mfc.sh generate` is sufficient. This is verified in Task 8 Step 4.
