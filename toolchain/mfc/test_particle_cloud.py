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
        "particle_cloud(1)%cloud_geometry": 2,
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
    params = {**_valid_cloud_params(), "particle_cloud(1)%cloud_geometry": 3}
    del params["particle_cloud(1)%packing_method"]
    with pytest.raises(CaseConstraintError, match="requires filepath"):
        CaseValidator(params).validate("simulation")
    params["particle_cloud(1)%filepath"] = "centres.dat"
    CaseValidator(params).validate("simulation")


@pytest.mark.parametrize("name, code", [("IB", 1), ("ib", 1), ("EL", 2), ("el", 2)])
def test_particle_type_accepts_names(name, code):
    from mfc.case import Case

    assert Case({"particle_cloud(3)%particle_type": name}).params["particle_cloud(3)%particle_type"] == code
