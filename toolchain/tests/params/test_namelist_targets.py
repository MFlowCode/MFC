def test_common_vars_in_all_targets():
    from mfc.params.namelist_targets import NAMELIST_VARS

    for var in ["m", "n", "p", "bc_x", "bc_y", "bc_z", "model_eqns", "cyl_coord", "fluid_pp", "case_dir"]:
        assert {"pre", "sim", "post"}.issubset(NAMELIST_VARS.get(var, set())), f"{var!r} not marked for all targets"


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


def test_case_opt_exclude_vars_also_in_sim_namelist():
    from mfc.params.namelist_targets import CASE_OPT_EXCLUDE, NAMELIST_VARS

    for var in CASE_OPT_EXCLUDE:
        targets = NAMELIST_VARS.get(var, set())
        assert "sim" in targets, f"CASE_OPT_EXCLUDE var {var!r} must also be in sim namelist"
