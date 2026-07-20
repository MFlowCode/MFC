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


def test_get_generated_files_returns_fifteen():
    from pathlib import Path

    from mfc.params.generators.fortran_gen import get_generated_files

    files = get_generated_files(Path("/build"))
    assert len(files) == 15
    paths = [str(p) for p, _ in files]
    assert any("pre_process/generated_namelist.fpp" in p for p in paths)
    assert any("simulation/generated_decls.fpp" in p for p in paths)
    assert any("post_process/generated_namelist.fpp" in p for p in paths)
    assert any("simulation/generated_case_opt_decls.fpp" in p for p in paths)
    assert any("pre_process/generated_bcast.fpp" in p for p in paths)
    assert any("simulation/generated_bcast.fpp" in p for p in paths)
    assert any("post_process/generated_bcast.fpp" in p for p in paths)


def test_generate_constants_fpp_content():
    from mfc.params.generators.fortran_gen import generate_constants_fpp

    out = generate_constants_fpp()
    assert "AUTO-GENERATED" in out
    assert "integer, parameter :: riemann_solver_hllc = 2" in out
    assert "integer, parameter :: model_eqns_5eq = 2" in out
    assert "integer, parameter :: time_stepper_rk3 = 3" in out
    # Deterministic: two calls produce identical output
    assert out == generate_constants_fpp()


def test_get_generated_files_includes_bcast():
    from pathlib import Path

    from mfc.params.generators.fortran_gen import get_generated_files

    files = get_generated_files(Path("/tmp/x"))
    names = {p.name for p, _ in files}
    assert names == {
        "generated_namelist.fpp",
        "generated_decls.fpp",
        "generated_constants.fpp",
        "generated_case_opt_decls.fpp",
        "generated_bcast.fpp",
    }
    assert len(files) == 15


def test_generate_case_opt_decls_fpp():
    from mfc.params.generators.fortran_gen import generate_case_opt_decls_fpp

    out = generate_case_opt_decls_fpp()

    # Must have the case-opt guard structure
    assert "#:if MFC_CASE_OPTIMIZATION" in out
    assert "#:else" in out
    assert "#:endif" in out

    # Representative registry-driven parameter lines (#:if branch)
    assert "integer, parameter :: weno_order = ${weno_order}$" in out
    assert "integer, parameter :: num_fluids = ${num_fluids}$" in out
    assert "logical, parameter :: mapped_weno = (${mapped_weno}$ /= 0)" in out
    assert "real(wp), parameter :: wenoz_q = ${wenoz_q}$" in out

    # #:else branch variable forms
    assert "integer                 :: weno_order" in out
    assert "logical                 :: mapped_weno" in out
    assert "real(wp)                :: wenoz_q" in out

    # case.py-computed extras
    assert "integer, parameter :: num_dims = ${num_dims}$" in out
    assert "integer, parameter :: weno_polyn = ${weno_polyn}$" in out
    assert "logical, parameter :: wenojs = (${wenojs}$ /= 0)" in out
    assert "integer                 :: num_dims" in out
    assert "integer                 :: weno_polyn" in out
    assert "logical                 :: wenojs" in out

    # nb must NOT be in this block (it lives in a separate block)
    assert ":: nb" not in out

    # All CASE_OPT_PARAMS (minus nb) must appear in both branches
    from mfc.params.definitions import CASE_OPT_PARAMS

    for name in CASE_OPT_PARAMS - {"nb"}:
        assert f":: {name}" in out, f"{name!r} missing from generated_case_opt_decls"

    # AUTO-GENERATED header present
    assert "AUTO-GENERATED" in out


# ── generate_bcast_fpp tests ──────────────────────────────────────────────────


def test_generate_bcast_fpp_header_and_case_dir():
    """All targets emit the AUTO-GENERATED header and case_dir broadcast."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    for target in ("pre", "sim", "post"):
        out = generate_bcast_fpp(target)
        assert "AUTO-GENERATED" in out, f"{target}: missing header"
        assert "call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)" in out


def test_generate_bcast_fpp_class_a_int_scalars():
    """Class-(a) integer scalars are broadcast for all targets."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    pre = generate_bcast_fpp("pre")
    sim = generate_bcast_fpp("sim")
    post = generate_bcast_fpp("post")

    # m, n, p are in all three targets
    for out, target in [(pre, "pre"), (sim, "sim"), (post, "post")]:
        assert "call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: m missing"
        assert "call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: n missing"
        assert "call MPI_BCAST(p, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: p missing"

    # precision is in all three
    for out, target in [(pre, "pre"), (sim, "sim"), (post, "post")]:
        assert "call MPI_BCAST(precision, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: precision missing"


def test_generate_bcast_fpp_class_a_log_scalars():
    """Class-(a) logical scalars are broadcast per target."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    pre = generate_bcast_fpp("pre")
    sim = generate_bcast_fpp("sim")
    post = generate_bcast_fpp("post")

    # bubbles_euler is in all three
    for out, target in [(pre, "pre"), (sim, "sim"), (post, "post")]:
        assert "call MPI_BCAST(bubbles_euler, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: bubbles_euler missing"

    # run_time_info is sim-only
    assert "run_time_info" in sim
    assert "run_time_info" not in pre
    assert "run_time_info" not in post

    # old_grid is pre-only
    assert "old_grid" in pre
    assert "old_grid" not in sim
    assert "old_grid" not in post


def test_generate_bcast_fpp_class_a_real_scalars():
    """Class-(a) real scalars (mpi_p) are broadcast per target."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")
    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")

    # dt is sim-only REAL
    assert "call MPI_BCAST(dt, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in sim
    assert "call MPI_BCAST(dt, " not in pre and "call MPI_BCAST(dt, " not in post

    # pref is in all three
    for out, target in [(pre, "pre"), (sim, "sim"), (post, "post")]:
        assert "call MPI_BCAST(pref, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out, f"{target}: pref missing"


def test_generate_bcast_fpp_class_a_str_scalars():
    """Class-(a) character scalars use the len() MPI_CHARACTER form, never mpi_p."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    pre = generate_bcast_fpp("pre")
    sim = generate_bcast_fpp("sim")
    post = generate_bcast_fpp("post")

    # files_dir/file_extension are pre-only STR scalars
    for name in ("files_dir", "file_extension"):
        assert f"call MPI_BCAST({name}, len({name}), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)" in pre, f"pre: {name} missing"
        assert f"call MPI_BCAST({name}, " not in sim and f"call MPI_BCAST({name}, " not in post
        assert f"call MPI_BCAST({name}, 1, mpi_p" not in pre


def test_generate_bcast_fpp_case_opt_guard_sim():
    """Sim case-opt scalars are wrapped in #:if not MFC_CASE_OPTIMIZATION."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    out = generate_bcast_fpp("sim")
    assert "#:if not MFC_CASE_OPTIMIZATION" in out
    assert "#:endif" in out

    # weno_order and num_fluids are CASE_OPT_PARAMS — must be inside the guard
    lines = out.splitlines()
    guard_start = next(i for i, ln in enumerate(lines) if "#:if not MFC_CASE_OPTIMIZATION" in ln)
    guard_end = next(i for i, ln in enumerate(lines) if i > guard_start and "#:endif" in ln)
    guard_body = "\n".join(lines[guard_start:guard_end])
    assert "weno_order" in guard_body
    assert "num_fluids" in guard_body
    assert "mapped_weno" in guard_body

    # muscl_eps IS inside the case-opt guard (it is sim-only, non-CASE_OPT_PARAM,
    # so it appears in the real-scalars section, which is outside the case-opt block)
    assert "call MPI_BCAST(muscl_eps, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out


def test_generate_bcast_fpp_case_opt_not_in_pre_post():
    """CASE_OPT_PARAMS that apply to pre/post are emitted unconditionally there."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")

    # weno_order and nb are in pre and post (not only sim) — emitted unconditionally
    assert "call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in pre
    assert "call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in post
    assert "#:if not MFC_CASE_OPTIMIZATION" not in pre
    assert "#:if not MFC_CASE_OPTIMIZATION" not in post


def test_generate_bcast_fpp_fluid_pp_loop():
    """fluid_pp member-loop is emitted for all targets; Re(1) count=2 only for sim."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    for target in ("pre", "sim", "post"):
        out = generate_bcast_fpp(target)
        assert "do i = 1, num_fluids_max" in out, f"{target}: fluid_pp loop missing"
        assert "call MPI_BCAST(fluid_pp(i)%gamma, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out
        assert "call MPI_BCAST(fluid_pp(i)%cv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out
        # Herschel-Bulkley members (#1545) — dropping any is a multi-rank regression
        for hb in ("K", "nn", "tau0", "hb_m", "mu_min", "mu_max", "mu_bulk"):
            assert f"call MPI_BCAST(fluid_pp(i)%{hb}, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out, f"missing HB member {hb}"
        assert "call MPI_BCAST(fluid_pp(i)%non_newtonian, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in out

    sim = generate_bcast_fpp("sim")
    assert "call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)" in sim

    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")
    assert "fluid_pp(i)%Re(1)" not in pre
    assert "fluid_pp(i)%Re(1)" not in post


def test_generate_bcast_fpp_bub_pp_under_guard():
    """bub_pp members are emitted under a bubbles guard for all targets."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    for target in ("pre", "sim", "post"):
        out = generate_bcast_fpp(target)
        assert "if (bubbles_euler .or. bubbles_lagrange) then" in out, f"{target}: bub_pp guard missing"
        assert "call MPI_BCAST(bub_pp%R0ref, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out
        assert "call MPI_BCAST(bub_pp%gam_g, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in out


def test_generate_bcast_fpp_lag_chem_sim_only():
    """lag_params and chem_params broadcasts are emitted for sim only."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")
    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")

    # lag_params under bubbles_lagrange guard
    assert "if (bubbles_lagrange) then" in sim
    assert "call MPI_BCAST(lag_params%solver_approach, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)" in sim
    assert "lag_params" not in pre
    assert "lag_params" not in post

    # chem_params under chemistry guard
    assert "if (chemistry) then" in sim
    assert "call MPI_BCAST(chem_params%diffusion, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in sim
    assert "chem_params" not in pre
    assert "chem_params" not in post


def test_generate_bcast_fpp_fortran_array_dims():
    """FORTRAN_ARRAY_DIMS arrays are broadcast with the correct dim and type."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")

    # fluid_rho is pre-only, REAL → mpi_p
    assert "call MPI_BCAST(fluid_rho(1), num_fluids_max, mpi_p, 0, MPI_COMM_WORLD, ierr)" in pre
    assert "fluid_rho" not in post

    # alpha_wrt, mom_wrt are post-only LOGICAL → MPI_LOGICAL
    assert "call MPI_BCAST(alpha_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in post
    assert "call MPI_BCAST(mom_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in post

    # schlieren_alpha is post-only REAL → mpi_p
    assert "call MPI_BCAST(schlieren_alpha(1), num_fluids_max, mpi_p, 0, MPI_COMM_WORLD, ierr)" in post


def test_generate_bcast_fpp_chem_wrt_Y_latent_bug_fix():
    """chem_wrt_Y is broadcast for post (latent bug fix: previously missing)."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    post = generate_bcast_fpp("post")
    # chem_wrt_Y is FORTRAN_ARRAY_DIMS with dim=num_species, type=LOG
    assert "call MPI_BCAST(chem_wrt_Y(1), num_species, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)" in post

    # Must not appear in pre or sim (post-only)
    pre = generate_bcast_fpp("pre")
    sim = generate_bcast_fpp("sim")
    assert "chem_wrt_Y" not in pre
    assert "chem_wrt_Y" not in sim


def test_generate_bcast_fpp_excludes_manual_residue():
    """Class-(c) vars and excluded vars must not appear in generated output."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    for target in ("pre", "sim", "post"):
        out = generate_bcast_fpp(target)
        # These are non-namelist and must stay manual
        assert "m_glb" not in out, f"{target}: m_glb should be manual"
        assert "n_glb" not in out, f"{target}: n_glb should be manual"
        assert "p_glb" not in out, f"{target}: p_glb should be manual"

    sim = generate_bcast_fpp("sim")
    # shear_stress, bulk_stress, bodyForces are derived (non-namelist)
    assert "shear_stress" not in sim
    assert "bulk_stress" not in sim
    assert "bodyForces" not in sim


def test_generate_bcast_fpp_deterministic():
    """Calling generate_bcast_fpp twice produces identical output (sorted, stable)."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    for target in ("pre", "sim", "post"):
        assert generate_bcast_fpp(target) == generate_bcast_fpp(target), f"{target}: not deterministic"


def test_generate_bcast_fpp_bad_target():
    """Unknown target raises ValueError."""
    import pytest

    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    with pytest.raises(ValueError, match="Unknown target"):
        generate_bcast_fpp("bad")


def test_generate_bcast_fpp_muscl_eps_now_broadcast():
    """muscl_eps is broadcast for sim (latent-bug fix: derivation is rank-0-only).

    Previously excluded via _BCAST_EXCLUDE; every multi-rank MUSCL run had
    rank-divergent muscl_eps because f_is_default() only fires on rank 0.
    """
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")
    assert "call MPI_BCAST(muscl_eps, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)" in sim

    # muscl_eps is sim-only; must not appear in pre/post
    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")
    assert "muscl_eps" not in pre
    assert "muscl_eps" not in post


def test_generate_bcast_fpp_fluid_pp_registry_walk():
    """fluid_pp emitter walks registry; no dead members (mul0/ss/pv/gamma_v/M_v/mu_v/k_v/cp_v/D_v removed upstream).

    Re(1) count=2 is sim-only; all other registered REAL members appear in all three
    targets.
    """
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")
    pre = generate_bcast_fpp("pre")
    post = generate_bcast_fpp("post")

    # All registered members appear in every target
    for mem in ("gamma", "pi_inf", "cv", "qv", "qvp", "G"):
        for out, t in [(sim, "sim"), (pre, "pre"), (post, "post")]:
            assert f"fluid_pp(i)%{mem}" in out, f"{t}: fluid_pp(i)%{mem} missing"

    # Re(1) count=2 is sim-only
    assert "fluid_pp(i)%Re(1)" in sim
    assert "fluid_pp(i)%Re(1)" not in pre
    assert "fluid_pp(i)%Re(1)" not in post

    # Dead members must not appear
    for dead in ("mul0", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "cp_v", "D_v"):
        for out, t in [(sim, "sim"), (pre, "pre"), (post, "post")]:
            assert f"fluid_pp(i)%{dead}" not in out, f"{t}: dead member fluid_pp(i)%{dead} present"


def test_generate_bcast_fpp_fluid_pp_member_datatypes():
    # The HB merge added a LOGICAL fluid_pp member; datatypes must come from the
    # registry, not be assumed REAL.
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")
    assert "call MPI_BCAST(fluid_pp(i)%non_newtonian, 1, MPI_LOGICAL," in sim
    assert "call MPI_BCAST(fluid_pp(i)%tau0, 1, mpi_p," in sim
    assert "non_newtonian, 1, mpi_p" not in sim


def test_generate_bcast_fpp_lag_params_registry_walk():
    """lag_params emitter walks registry; no dead members (T0/Thost/c0/rho0/x0 removed upstream)."""
    from mfc.params.generators.fortran_gen import generate_bcast_fpp

    sim = generate_bcast_fpp("sim")

    # All registered members must appear
    for mem in ("solver_approach", "cluster_type", "smooth_type", "nBubs_glb"):
        assert f"lag_params%{mem}" in sim, f"lag_params%{mem} missing from sim"
    for mem in ("heatTransfer_model", "massTransfer_model", "pressure_corrector", "write_bubbles", "write_bubbles_stats"):
        assert f"lag_params%{mem}" in sim, f"lag_params%{mem} missing from sim"
    for mem in ("epsilonb", "charwidth", "valmaxvoid"):
        assert f"lag_params%{mem}" in sim, f"lag_params%{mem} missing from sim"

    # Dead members must not appear
    for dead in ("T0", "Thost", "c0", "rho0", "x0"):
        assert f"lag_params%{dead}" not in sim, f"dead member lag_params%{dead} present in sim"


def test_mpi_proxy_residue_pins_wall_velocity_and_bc_datatypes():
    """The vb/ve wall-velocity broadcasts and integer BC datatypes live in
    hand-written residue (not codegen); pin them so an edit or merge conflict
    that drops them fails loudly."""
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[3]
    for target in ("pre_process", "post_process"):
        src = (root / "src" / target / "m_mpi_proxy.fpp").read_text()
        assert "bc_${DIM}$%vb${DIR}$" in src, f"{target}: vb broadcasts missing"
        assert "bc_${DIM}$%ve${DIR}$" in src, f"{target}: ve broadcasts missing"
        assert "'bc_x%beg', 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end']" in src
        seg = src.split("'bc_z%beg', 'bc_z%end']", 1)[1]
        assert "MPI_INTEGER" in seg.split("#:endfor")[0], f"{target}: BC codes not MPI_INTEGER"
