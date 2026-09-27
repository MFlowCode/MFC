"""External analytic material files are resolved before MFC case validation."""

import pytest
import yaml

from mfc.common import MFCException
from mfc.materials import resolve_materials
from mfc.run.input import MFCInputFile

_JWL = {"A": 6.0, "B": 0.15, "R1": 4.0, "R2": 1.0, "omega": 0.3, "rho0": 1.0}


def _write_material(path, family="jwl", parameters=None, release_status="public"):
    parameters = {**_JWL, "Q": 2.0} if parameters is None else parameters
    path.write_text(
        yaml.safe_dump(
            {
                "material": {"name": "synthetic_material", "eos": {"family": family, "parameters": parameters}},
                "provenance": {"citation": "synthetic example", "release_status": release_status},
            }
        )
    )


def test_material_expands_into_runtime_case_parameters(tmp_path, capsys):
    path = tmp_path / "products.yaml"
    _write_material(path, parameters={**_JWL, "A": "3.712e11", "Q": 2.0, "qv": -2.0})
    case = MFCInputFile("case.py", str(tmp_path), {"fluid_pp(1)%material_file": str(path)})
    assert case.params["fluid_pp(1)%eos"] == 4
    assert case.params["fluid_pp(1)%jwl_a"] == 3.712e11
    assert case.params["fluid_pp(1)%qv"] == -2.0
    assert "fluid_pp(1)%jwl_Q" not in case.params
    assert "material Q is metadata" in capsys.readouterr().out
    assert "fluid_pp(1)%material_file" not in case.params
    inp = case.get_inp("simulation")
    assert "fluid_pp(1)%jwl_a = 371200000000.0" in inp
    assert "fluid_pp(1)%qv = -2.0" in inp


@pytest.mark.parametrize(
    ("reactant", "pi_inf"),
    [
        ({"fluid_pp(1)%eos": "stiffened_gas", "fluid_pp(1)%pi_inf": 0.5}, 0.5),
        ({"fluid_pp(1)%eos": "mie_gruneisen", "fluid_pp(1)%mg_rho0": 2.0, "fluid_pp(1)%mg_c0": 1.0, "fluid_pp(1)%mg_s": 1.5, "fluid_pp(1)%mg_gruneisen": 0.4}, 0.0),
    ],
)
def test_material_q_sets_reactive_qv(tmp_path, reactant, pi_inf):
    path = tmp_path / "products.yaml"
    _write_material(path, parameters={**_JWL, "rho0": 2.0, "Q": 3.0})
    loaded = resolve_materials({"reactive_burn": "T", **reactant, "fluid_pp(2)%material_file": str(path)}, str(tmp_path))
    assert loaded["fluid_pp(1)%qv"] == pytest.approx(3.0 - pi_inf / 2.0)
    assert loaded["fluid_pp(2)%qv"] == 0.0


def test_material_q_rejects_misplaced_or_duplicated_energy(tmp_path):
    path = tmp_path / "products.yaml"
    _write_material(path)
    with pytest.raises(MFCException, match="belongs to the products"):
        resolve_materials({"reactive_burn": "T", "fluid_pp(1)%material_file": str(path)}, str(tmp_path))
    with pytest.raises(MFCException, match="sets both reactant and product qv"):
        resolve_materials({"reactive_burn": "T", "fluid_pp(1)%qv": 1.0, "fluid_pp(2)%material_file": str(path)}, str(tmp_path))
    _write_material(path, parameters={**_JWL, "Q": -1.0})
    with pytest.raises(MFCException, match="positive JWL Q"):
        resolve_materials({"fluid_pp(2)%material_file": str(path)}, str(tmp_path))


def test_material_search_order(tmp_path, monkeypatch):
    case_dir, public_dir = tmp_path / "case", tmp_path / "public"
    case_dir.mkdir()
    public_dir.mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("MFC_PUBLIC_MATERIAL_DIR", str(public_dir))
    for directory, value in ((tmp_path, 6.0), (case_dir, 7.0), (public_dir, 8.0)):
        _write_material(directory / "products.yaml", parameters={**_JWL, "A": value})
    params = {"fluid_pp(1)%material_file": "products.yaml"}
    for directory, expected in ((tmp_path, 6.0), (case_dir, 7.0), (public_dir, 8.0)):
        assert resolve_materials(params, str(case_dir))["fluid_pp(1)%jwl_a"] == expected
        (directory / "products.yaml").unlink()


def test_material_rejects_conflicts_and_bad_provenance(tmp_path):
    path = tmp_path / "products.yaml"
    _write_material(path)
    directive = {"fluid_pp(2)%material_file": str(path)}
    with pytest.raises(MFCException, match="conflicts with fluid_pp\\(2\\)%jwl_a"):
        resolve_materials({**directive, "fluid_pp(2)%jwl_a": 9.0}, str(tmp_path))
    _write_material(path, parameters={**_JWL, "qv": 1.0})
    with pytest.raises(MFCException, match="conflicts with fluid_pp\\(2\\)%qv"):
        resolve_materials({**directive, "fluid_pp(2)%qv": 2.0}, str(tmp_path))
    with pytest.raises(MFCException, match="conflicts with fluid_pp\\(2\\)%eos"):
        resolve_materials({**directive, "fluid_pp(2)%eos": "ideal_gas"}, str(tmp_path))
    with pytest.raises(MFCException, match="beyond num_fluids"):
        resolve_materials({**directive, "num_fluids": 1}, str(tmp_path))
    _write_material(path, release_status="")
    with pytest.raises(MFCException, match="provenance.citation and provenance.release_status"):
        resolve_materials(directive, str(tmp_path))


def test_material_rejects_nonfinite_values_and_missing_file(tmp_path):
    path = tmp_path / "products.yaml"
    _write_material(path, parameters={**_JWL, "A": ".nan"})
    with pytest.raises(MFCException, match="finite numeric"):
        resolve_materials({"fluid_pp(1)%material_file": str(path)}, str(tmp_path))
    with pytest.raises(MFCException, match="Searched"):
        resolve_materials({"fluid_pp(1)%material_file": "missing.yaml"}, str(tmp_path))


@pytest.mark.parametrize(
    ("family", "prefix", "parameters"),
    [
        ("mie_gruneisen", "mg", {"rho0": 1.0, "c0": 1.0, "s": 1.5, "gruneisen": 0.4, "gruneisen_a": 0.2, "t0": 300.0, "s2": 0.1, "s3": 0.02}),
        ("vinet", "vinet", {"k0": 10.0, "k0p": 4.0, "rho0": 1.0, "gruneisen": 1.2, "gruneisen_a": 0.1, "t0": 300.0}),
    ],
)
def test_material_loads_family_parameters(tmp_path, family, prefix, parameters):
    path = tmp_path / "material.yaml"
    _write_material(path, family=family, parameters={**parameters, "cv": 2.0, "qv": 3.0})
    loaded = resolve_materials({"fluid_pp(1)%material_file": str(path)}, str(tmp_path))
    assert loaded["fluid_pp(1)%eos"] == family
    assert all(loaded[f"fluid_pp(1)%{prefix}_{name}"] == value for name, value in parameters.items())
    assert loaded["fluid_pp(1)%cv"] == 2.0
    assert loaded["fluid_pp(1)%qv"] == 3.0
