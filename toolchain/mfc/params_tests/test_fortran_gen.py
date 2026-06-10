def test_get_namelist_var_simple():
    from mfc.params.generators.fortran_gen import get_namelist_var

    assert get_namelist_var("m") == "m"
    assert get_namelist_var("dt") == "dt"


def test_get_namelist_var_indexed_family():
    from mfc.params.generators.fortran_gen import get_namelist_var

    assert get_namelist_var("fluid_pp(1)%gamma") == "fluid_pp"
    assert get_namelist_var("patch_icpp(3)%geometry") == "patch_icpp"


def test_get_namelist_var_struct_member():
    from mfc.params.generators.fortran_gen import get_namelist_var

    assert get_namelist_var("bc_x%beg") == "bc_x"
    assert get_namelist_var("lag_params%solver_approach") == "lag_params"


def test_fortran_type_int():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType

    assert fortran_type_decl(ParamDef(name="x", param_type=ParamType.INT)) == "integer"


def test_fortran_type_real():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType

    assert fortran_type_decl(ParamDef(name="x", param_type=ParamType.REAL)) == "real(wp)"


def test_fortran_type_real_storage_precision():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType

    p = ParamDef(name="x", param_type=ParamType.REAL, storage_precision=True)
    assert fortran_type_decl(p) == "real(stp)"
    a = ParamDef(name="x", param_type=ParamType.ANALYTIC_REAL, storage_precision=True)
    assert fortran_type_decl(a) == "real(stp)"


def test_fortran_type_log():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType

    assert fortran_type_decl(ParamDef(name="x", param_type=ParamType.LOG)) == "logical"


def test_fortran_type_str():
    from mfc.params.generators.fortran_gen import fortran_type_decl
    from mfc.params.schema import ParamDef, ParamType

    p = ParamDef(name="case_dir", param_type=ParamType.STR, str_len="path_len")
    assert fortran_type_decl(p) == "character(LEN=path_len)"


def test_namelist_contains_common_vars():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp

    for target in ("pre", "sim", "post"):
        c = generate_namelist_fpp(target)
        for v in ("m", "n", "p", "bc_x", "case_dir", "fluid_pp"):
            assert v in c, f"{v!r} missing from {target} namelist"


def test_sim_namelist_case_opt_guard():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp

    c = generate_namelist_fpp("sim")
    # Case-opt guard: two complete namelist statements wrapped in #:if/#:else/#:endif
    assert "#:if MFC_CASE_OPTIMIZATION" in c
    assert "#:else" in c
    assert "#:endif" in c
    assert "weno_order" in c
    assert "num_fluids" in c
    # No dangling continuation before the #:if block or after #:else
    lines = c.splitlines()
    for i, line in enumerate(lines):
        if line.strip() in ("#:if MFC_CASE_OPTIMIZATION", "#:endif"):
            assert not lines[i - 1].rstrip().endswith("&"), f"Dangling & before {line!r}"


def test_pre_namelist_has_patch_icpp():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp

    c = generate_namelist_fpp("pre")
    assert "patch_icpp" in c
    assert "run_time_info" not in c


def test_post_namelist_has_sim_data():
    from mfc.params.generators.fortran_gen import generate_namelist_fpp

    c = generate_namelist_fpp("post")
    assert "sim_data" in c
    assert "patch_icpp" not in c


def test_decls_contains_simple_scalars():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    for target in ("pre", "sim", "post"):
        c = generate_decls_fpp(target)
        assert "integer" in c
        assert "real(wp)" in c
        assert "logical" in c


def test_decls_dt_for_sim():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    # Column-aligned output: type padded to 24 chars before '::'
    assert "real(wp)                :: dt" in generate_decls_fpp("sim")


def test_decls_no_percent_vars():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    for target in ("pre", "sim", "post"):
        c = generate_decls_fpp(target)
        assert "bc_x%beg" not in c
        assert "fluid_pp(1)" not in c


def test_decls_case_dir():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    for target in ("pre", "sim", "post"):
        assert "character(LEN=path_len) :: case_dir" in generate_decls_fpp(target)


def test_decls_array_dims():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    post = generate_decls_fpp("post")
    assert "dimension(num_fluids_max)" in post and ":: alpha_wrt" in post
    assert "dimension(num_fluids_max)" in post and ":: alpha_rho_wrt" in post
    assert "dimension(3)" in post and ":: mom_wrt" in post
    # fluid_pp is now emitted via TYPED_DECLS; bc_x still must not appear
    assert ":: fluid_pp" in post
    assert ":: bc_x" not in post

    pre = generate_decls_fpp("pre")
    assert "dimension(num_fluids_max)" in pre and ":: fluid_rho" in pre


def test_generate_decls_emits_typed_declarations():
    from mfc.params.generators.fortran_gen import generate_decls_fpp

    pre = generate_decls_fpp("pre")
    sim = generate_decls_fpp("sim")
    post = generate_decls_fpp("post")

    # patch_icpp is pre-only
    assert "patch_icpp" in pre and "patch_icpp" not in post and "patch_icpp" not in sim

    # fluid_pp appears in all three targets
    for out in (pre, sim, post):
        assert "fluid_pp" in out

    # acoustic is sim-only
    assert "acoustic" in sim and "acoustic" not in pre

    # Exact line shape (harvested from the manual declaration it replaces).
    # Note: _ARRAY_DECL_COL=36; the type+dim string is longer so no padding space before ::
    assert "type(physical_parameters), dimension(num_fluids_max) :: fluid_pp" in sim

    # chem_params is sim-only with GPU declare
    assert "chem_params" in sim and "chem_params" not in pre and "chem_params" not in post
    assert "$:GPU_DECLARE(create='[chem_params]')" in sim

    # bub_pp is in all three targets (derived-type scalar, dim=None)
    for out in (pre, sim, post):
        assert "bub_pp" in out

    # lag_params is sim-only; the generator owns its GPU_DECLARE (removed from the
    # grouped list in m_global_parameters.fpp by the Task-2 Fortran edit)
    assert "lag_params" in sim and "lag_params" not in pre
    for gpu_var in ("patch_ib", "ib_airfoil", "acoustic", "lag_params", "chem_params"):
        assert f"$:GPU_DECLARE(create='[{gpu_var}]')" in sim

    # Doxygen descriptions from TYPED_DECLS are emitted as inline comments
    assert ":: fluid_pp !< Per-fluid stiffened-gas EOS parameters, Reynolds numbers, and shear modulus" in sim

    # simplex_params is pre-only
    assert "simplex_params" in pre and "simplex_params" not in sim

    # No TYPED_DECLS name emitted twice
    for name in ("fluid_pp", "bub_pp", "chem_params", "lag_params"):
        assert sim.count(f":: {name}") == 1, f"{name!r} emitted more than once in sim"

    # TYPED_DECLS keys must not overlap with FORTRAN_ARRAY_DIMS
    from mfc.params.definitions import FORTRAN_ARRAY_DIMS, TYPED_DECLS

    assert not (set(TYPED_DECLS) & set(FORTRAN_ARRAY_DIMS)), "TYPED_DECLS and FORTRAN_ARRAY_DIMS share keys"


def test_check_target_raises_on_bad_target():
    import pytest

    from mfc.params.generators.fortran_gen import generate_decls_fpp, generate_namelist_fpp

    with pytest.raises(ValueError, match="Unknown target"):
        generate_namelist_fpp("bad")
    with pytest.raises(ValueError, match="Unknown target"):
        generate_decls_fpp("bad")


def test_get_generated_files_returns_nine():
    from pathlib import Path

    from mfc.params.generators.fortran_gen import get_generated_files

    files = get_generated_files(Path("/build"))
    assert len(files) == 9
    paths = [str(p) for p, _ in files]
    assert any("pre_process/generated_namelist.fpp" in p for p in paths)
    assert any("simulation/generated_decls.fpp" in p for p in paths)
    assert any("post_process/generated_namelist.fpp" in p for p in paths)


def test_generate_constants_fpp_content():
    from mfc.params.generators.fortran_gen import generate_constants_fpp

    out = generate_constants_fpp()
    assert "AUTO-GENERATED" in out
    assert "integer, parameter :: riemann_solver_hllc = 2" in out
    assert "integer, parameter :: model_eqns_5eq = 2" in out
    assert "integer, parameter :: time_stepper_rk3 = 3" in out
    # Deterministic: two calls produce identical output
    assert out == generate_constants_fpp()


def test_get_generated_files_includes_constants():
    from pathlib import Path

    from mfc.params.generators.fortran_gen import get_generated_files

    files = get_generated_files(Path("/tmp/x"))
    names = {p.name for p, _ in files}
    assert names == {"generated_namelist.fpp", "generated_decls.fpp", "generated_constants.fpp"}
    assert len(files) == 9
