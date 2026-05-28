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
    assert "#:if not MFC_CASE_OPTIMIZATION" in c
    assert "weno_order" in c
    assert "num_fluids" in c

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
    assert "real(wp) :: dt" in generate_decls_fpp("sim")

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

def test_get_generated_files_returns_six():
    from pathlib import Path

    from mfc.params.generators.fortran_gen import get_generated_files
    files = get_generated_files(Path("/tmp"))
    assert len(files) == 6
    names = {p.name for p, _ in files}
    assert "generated_namelist_pre.fpp" in names
    assert "generated_decls_sim.fpp" in names
