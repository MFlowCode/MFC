import json

import pytest

from mfc import amr_check
from mfc.common import MFCException

BASE = {"m": 49, "n": 39, "p": 0, "dt": 1e-3, "t_step_start": 0, "t_step_stop": 100, "t_step_save": 50, "weno_order": 5, "time_stepper": 3}


def test_default_block_is_the_central_half_of_each_active_dimension():
    assert amr_check.block(BASE) == [(12, 36), (10, 29)]


def test_a_block_given_in_the_case_is_kept():
    case = {**BASE, "amr_block_beg(1)": 6, "amr_block_end(1)": 30, "amr_block_beg(2)": 6, "amr_block_end(2)": 25}
    assert amr_check.block(case) == [(6, 30), (6, 25)]


def test_variants_differ_only_in_the_grid_and_the_amr_block():
    amr, fine = amr_check.variants({**BASE, "amr_tag_eps": 0.1}, steps=2)
    assert (amr["amr"], amr["amr_regrid_int"], amr["amr_block_beg(2)"], amr["amr_block_end(2)"]) == ("T", 0, 10, 29)
    assert (fine["m"], fine["n"], fine["p"]) == (99, 79, 0)
    assert not any(k.startswith("amr") for k in fine)
    for case in (amr, fine):
        assert (case["t_step_stop"], case["t_step_save"], case["parallel_io"], case["dt"]) == (2, 2, "F", 1e-3)


def test_subcycling_halves_the_fine_run_time_step():
    _, fine = amr_check.variants({**BASE, "amr_subcycle": "T"}, steps=2)
    assert (fine["dt"], fine["t_step_stop"], fine["t_step_save"]) == (5e-4, 4, 4)


def test_an_adaptive_time_step_is_rejected():
    with pytest.raises(MFCException, match="fixed dt"):
        amr_check.variants({**BASE, "cfl_adap_dt": "T"}, steps=1)


def test_a_restarted_case_is_rejected():
    with pytest.raises(MFCException, match="t_step_start"):
        amr_check.variants({**BASE, "t_step_start": 10}, steps=1)


def test_the_written_case_keeps_json_only_values(tmp_path):
    params = {**BASE, "flag": True, "off": False, "none": None}
    with open(amr_check.write_case(str(tmp_path / "amr"), params)) as f:
        assert json.load(f) == params


def test_files_named_relative_to_the_case_become_absolute(tmp_path):
    (tmp_path / "mech.yaml").write_text("")
    out = amr_check.absolute_paths({"cantera_file": "mech.yaml", "surface_cantera_file": "h2o2.yaml", "weno_order": 5, "amr": "T"}, str(tmp_path))
    assert out == {"cantera_file": str(tmp_path / "mech.yaml"), "surface_cantera_file": "h2o2.yaml", "weno_order": 5, "amr": "T"}


def test_reach_counts_the_stencil_of_every_stage_in_coarse_cells():
    assert amr_check.reach(BASE, steps=1) == 7  # 3 stages x 4 fine cells = 12 fine = 6 coarse, plus one
    assert amr_check.reach({**BASE, "amr_subcycle": "T"}, steps=1) == 13
    assert amr_check.reach({**BASE, "weno_order": 3, "time_stepper": 1}, steps=2) == 4


def _fine(value):
    return {1: {(i, j): value(i, j) for i in range(16) for j in range(16)}}


def test_depth_errors_compare_each_block_cell_with_the_mean_of_its_fine_cells():
    fine = _fine(lambda i, j: float(i + 10 * j))
    block = [(2, 5), (2, 5)]
    coarse = amr_check.restrict(fine, block)
    assert coarse[1][(2, 3)] == pytest.approx(0.25 * sum(fine[1][(4 + a, 6 + b)] for a in (0, 1) for b in (0, 1)))
    assert amr_check.depth_errors(coarse, coarse, block) == {0: 0.0, 1: 0.0}
    off = {1: {**coarse[1], (2, 4): coarse[1][(2, 4)] + 1.0, (3, 3): coarse[1][(3, 3)] + 0.5}}
    scale = max(abs(v) for v in coarse[1].values())
    assert amr_check.depth_errors(off, coarse, block) == pytest.approx({0: 1.0 / scale, 1: 0.5 / scale})


def test_verdict_needs_matching_initial_data_an_active_core_and_an_exact_core():
    ok = dict(initial={0: 0.0, 8: 0.0}, final={0: 1e-3, 8: 1e-14}, active={0: 0.2, 8: 0.1}, reach=7, tol=1e-8)
    assert amr_check.verdict(**ok)[0]
    assert not amr_check.verdict(**{**ok, "final": {0: 1e-3, 8: 1e-5}})[0]
    assert "initial" in amr_check.verdict(**{**ok, "initial": {0: 0.0, 8: 1e-4}})[1]
    assert "did not change" in amr_check.verdict(**{**ok, "active": {0: 0.2, 8: 0.0}})[1]
    assert "too small" in amr_check.verdict(**{**ok, "reach": 9})[1]
