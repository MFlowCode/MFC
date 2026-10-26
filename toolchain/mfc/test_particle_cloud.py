import pytest

from mfc.case_validator import CaseConstraintError, CaseValidator


def _valid_cloud_params():
    return {
        "m": 31,
        "n": 31,
        "p": 31,
        "model_eqns": 2,
        "num_fluids": 1,
        "num_patches": 1,
        "t_step_start": 0,
        "t_step_stop": 1,
        "t_step_save": 1,
        "dt": 1.0e-6,
        "x_domain%beg": -1.0,
        "x_domain%end": 1.0,
        "y_domain%beg": -1.0,
        "y_domain%end": 1.0,
        "z_domain%beg": -1.0,
        "z_domain%end": 1.0,
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "bc_z%beg": -1,
        "bc_z%end": -1,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fd_order": 2,
        "ib": "T",
        "num_ibs": 0,
        "num_particle_clouds": 1,
        "particle_cloud(1)%cloud_geometry": 20,
        "particle_cloud(1)%packing_method": 1,
        "particle_cloud(1)%x_centroid": 0.0,
        "particle_cloud(1)%y_centroid": 0.0,
        "particle_cloud(1)%z_centroid": 0.0,
        "particle_cloud(1)%num_particles": 8,
        "particle_cloud(1)%radius": 0.04,
        "particle_cloud(1)%mass": 1.0,
        "particle_cloud(1)%min_spacing": 0.01,
        "particle_cloud(1)%shell_inner_radius": 0.1,
        "particle_cloud(1)%shell_outer_radius": 0.3,
    }


def test_hemi_shell_validator_accepts_feasible_shell():
    CaseValidator(_valid_cloud_params()).validate("simulation")


def test_hemi_shell_validator_rejects_negative_inner_radius():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_inner_radius": -0.01}
    with pytest.raises(CaseConstraintError, match="shell_inner_radius >= 0"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_thinner_than_particle_diameter():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_outer_radius": 0.17}
    with pytest.raises(CaseConstraintError, match="shell_outer_radius > shell_inner_radius"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_outside_domain():
    params = {**_valid_cloud_params(), "particle_cloud(1)%x_centroid": 0.9, "particle_cloud(1)%shell_outer_radius": 0.3}
    with pytest.raises(CaseConstraintError, match="x-extent must lie within x_domain"):
        CaseValidator(params).validate("simulation")


def _el_cloud_params():
    return {**_valid_cloud_params(), "ib": "F", "bubbles_lagrange": "T", "lag_params%nBubs_glb": 8, "particle_cloud(1)%particle_type": 2}


def test_el_cloud_validator_accepts_shell():
    params = _el_cloud_params()
    del params["particle_cloud(1)%mass"]
    CaseValidator(params).validate("simulation")


def test_el_cloud_validator_requires_bubbles_lagrange():
    params = {**_el_cloud_params(), "bubbles_lagrange": "F"}
    with pytest.raises(CaseConstraintError, match="requires bubbles_lagrange"):
        CaseValidator(params).validate("simulation")


def test_el_cloud_validator_requires_enough_nbubs_glb():
    params = {**_el_cloud_params(), "lag_params%nBubs_glb": 7}
    with pytest.raises(CaseConstraintError, match="nBubs_glb must be >="):
        CaseValidator(params).validate("simulation")


def test_el_cloud_validator_rejects_custom_input_path():
    params = {**_el_cloud_params(), "lag_params%input_path": "bubbles.dat"}
    with pytest.raises(CaseConstraintError, match="input_path must be left at its default"):
        CaseValidator(params).validate("simulation")


def test_ib_cloud_validator_requires_ib():
    params = {**_valid_cloud_params(), "ib": "F"}
    with pytest.raises(CaseConstraintError, match="requires ib to be enabled"):
        CaseValidator(params).validate("simulation")


def test_file_cloud_validator_requires_filepath():
    params = {**_valid_cloud_params(), "particle_cloud(1)%cloud_geometry": 21}
    del params["particle_cloud(1)%packing_method"]
    with pytest.raises(CaseConstraintError, match="requires filepath"):
        CaseValidator(params).validate("simulation")
    params["particle_cloud(1)%filepath"] = "centres.dat"
    CaseValidator(params).validate("simulation")


@pytest.mark.parametrize("name, code", [("IB", 1), ("ib", 1), ("EL", 2), ("el", 2)])
def test_particle_type_accepts_names(name, code):
    from mfc.case import Case

    assert Case({"particle_cloud(3)%particle_type": name}).params["particle_cloud(3)%particle_type"] == code


def test_pdf_validator_accepts_rejection_shell():
    params = {**_valid_cloud_params(), "particle_cloud(1)%pdf": "1/((x-xc)**2 + (y-yc)**2 + (z-zc)**2)"}
    CaseValidator(params).validate("simulation")


def test_pdf_validator_rejects_file_clouds():
    params = {**_valid_cloud_params(), "particle_cloud(1)%cloud_geometry": 21, "particle_cloud(1)%filepath": "c.dat", "particle_cloud(1)%pdf": "exp(-x**2)"}
    with pytest.raises(CaseConstraintError, match="pdf requires a sampled cloud"):
        CaseValidator(params).validate("simulation")


def test_pdf_codegen_only_for_expressions():
    from mfc.case import Case

    def pdf_fpp(params):
        return Case({"n": 31, "p": 31, "num_particle_clouds": 2, **params})._Case__get_analytic_cloud_pdf_fpp(False)

    assert pdf_fpp({"particle_cloud(1)%pdf": 1.0}) == ""
    fpp = pdf_fpp({"particle_cloud(2)%pdf": "exp(-(x - xc)**2/lx)"})
    assert "#:def particle_cloud_pdf()" in fpp
    assert "case (2)" in fpp and "case (1)" not in fpp
    assert "pdf = exp((-(x - particle_cloud(2)%x_centroid) ** 2) / particle_cloud(2)%length_x)" in fpp


def test_pdf_codegen_rejects_unknown_names():
    from mfc.case import Case
    from mfc.common import MFCException

    with pytest.raises(MFCException, match="particle_cloud\\(1\\)%pdf"):
        Case({"n": 31, "p": 31, "num_particle_clouds": 1, "particle_cloud(1)%pdf": "exp(-q)"})._Case__get_analytic_cloud_pdf_fpp(False)


def _patch_cloud_params(**mods):
    params = {**_valid_cloud_params(), "particle_cloud(1)%cloud_geometry": 8, "particle_cloud(1)%region_radius": 0.5}
    for key in ("shell_inner_radius", "shell_outer_radius"):
        del params[f"particle_cloud(1)%{key}"]
    return {**params, **{f"particle_cloud(1)%{k}": v for k, v in mods.items()}}


@pytest.mark.parametrize("packing", [1, 2])
def test_patch_cloud_validator_accepts_sphere(packing):
    CaseValidator(_patch_cloud_params(packing_method=packing)).validate("simulation")


def test_patch_cloud_validator_requires_shape_fields():
    with pytest.raises(CaseConstraintError, match="region_radius must be > 0"):
        CaseValidator(_patch_cloud_params(region_radius=-1.0)).validate("simulation")
    with pytest.raises(CaseConstraintError, match="model_id must be between 1 and num_stl_models"):
        CaseValidator(_patch_cloud_params(cloud_geometry=12)).validate("simulation")


def test_patch_cloud_validator_rejects_wrong_dimension():
    with pytest.raises(CaseConstraintError, match="2D patch shape"):
        CaseValidator(_patch_cloud_params(cloud_geometry=2)).validate("simulation")


def test_lattice_accepts_pdf_in_every_sampled_geometry():
    CaseValidator(_patch_cloud_params(packing_method=2, pdf="x**2")).validate("simulation")
    CaseValidator({**_valid_cloud_params(), "particle_cloud(1)%packing_method": 2, "particle_cloud(1)%pdf": "x**2"}).validate("simulation")
    CaseValidator(_patch_cloud_params(cloud_geometry=1, packing_method=2, pdf="x**2", length_x=0.5, length_y=0.5, length_z=0.5)).validate("simulation")


def test_cloud_geometry_accepts_names():
    from mfc.case import Case

    assert Case({"particle_cloud(1)%cloud_geometry": "hemisphere_shell"}).params["particle_cloud(1)%cloud_geometry"] == 20
    assert Case({"particle_cloud(1)%cloud_geometry": "STL_3D"}).params["particle_cloud(1)%cloud_geometry"] == 12


def _ib_cloud(**mods):
    params = {k: v for k, v in _valid_cloud_params().items() if k != "particle_cloud(1)%mass"}
    return {**params, **{f"particle_cloud(1)%{k}": v for k, v in mods.items()}}


def test_size_distribution_accepts_lognormal_and_classes_with_density():
    CaseValidator(_ib_cloud(density=10.0, radius_sigma=0.4)).validate("simulation")
    CaseValidator(_ib_cloud(density=10.0, **{"radius_ratio(2)": 2.0, "radius_fraction(1)": 0.7, "radius_fraction(2)": 0.3})).validate("simulation")


@pytest.mark.parametrize(
    "mods, message",
    [
        ({}, "needs density"),
        ({"density": 10.0, "mass": 1.0}, "density or mass, not both"),
        ({"mass": 1.0, "radius_sigma": 0.4}, "with varying sizes give density"),
        ({"density": 10.0, "radius_sigma": 0.4, "radius_fraction(1)": 1.0}, "pick one"),
        ({"density": 10.0, "radius_sigma": -0.1}, "radius_sigma must be > 0"),
        ({"density": 10.0, "radius_min": 0.01}, "need radius_sigma"),
        ({"density": 10.0, "radius_sigma": 0.4, "cloud_geometry": 1, "packing_method": 2}, "need rejection packing"),
    ],
)
def test_size_distribution_validator_rejects(mods, message):
    with pytest.raises(CaseConstraintError, match=message):
        CaseValidator(_ib_cloud(**mods)).validate("simulation")


def test_file_cloud_needs_no_radius():
    params = {k: v for k, v in _ib_cloud(cloud_geometry=21, filepath="c.dat", density=1.0).items() if not k.endswith("%radius")}
    CaseValidator(params).validate("simulation")
    with pytest.raises(CaseConstraintError, match="radius must be specified"):
        CaseValidator({k: v for k, v in _ib_cloud(density=1.0).items() if not k.endswith("%radius")}).validate("simulation")
