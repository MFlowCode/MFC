"""External analytic material files are resolved before MFC case validation."""

import pytest
import yaml

from mfc.common import MFCException
from mfc.materials import resolve_materials
from mfc.run.input import MFCInputFile

_JWL = {"A": 6.0, "B": 0.15, "R1": 4.0, "R2": 1.0, "omega": 0.3, "rho0": 1.0}


def _material(directory, parameters=None, family="jwl", release_status="public"):
    """Write directory/products.yaml and return its path."""
    path = directory / "products.yaml"
    material = {"name": "synthetic_material", "eos": {"family": family, "parameters": {**_JWL, "Q": 2.0} if parameters is None else parameters}}
    path.write_text(yaml.safe_dump({"material": material, "provenance": {"citation": "synthetic example", "release_status": release_status}}))
    return str(path)


def test_material_expands_into_runtime_case_parameters(tmp_path):
    path = _material(tmp_path, {**_JWL, "A": "3.712e11", "Q": 2.0, "qv": -2.0})
    case = MFCInputFile("case.py", str(tmp_path), {"fluid_pp(1)%material_file": path})
    assert case.params["fluid_pp(1)%eos"] == 4
    assert not {"fluid_pp(1)%q", "fluid_pp(1)%material_file"} & case.params.keys()
    inp = case.get_inp("simulation")
    assert "fluid_pp(1)%jwl_a = 371200000000.0" in inp
    assert "fluid_pp(1)%qv = -2.0" in inp


@pytest.mark.parametrize(
    ("reactant", "pi_inf"),
    [
        ({"fluid_pp(1)%eos": "stiffened_gas", "fluid_pp(1)%pi_inf": 0.5}, 0.5),
        ({"fluid_pp(1)%eos": 4, **{f"fluid_pp(1)%jwl_{k.lower()}": v for k, v in {**_JWL, "rho0": 2.0}.items()}}, -0.4675971238515872),
    ],
)
def test_material_q_sets_reactive_qv(tmp_path, reactant, pi_inf):
    path = _material(tmp_path, {**_JWL, "rho0": 2.0, "Q": 3.0})
    loaded = resolve_materials({"reactive_burn": "T", **reactant, "fluid_pp(2)%material_file": path}, str(tmp_path))
    assert loaded["fluid_pp(1)%qv"] == pytest.approx(3.0 - pi_inf / 2.0)
    assert loaded["fluid_pp(2)%qv"] == 0.0


def test_material_search_order(tmp_path, monkeypatch):
    case_dir, public_dir = tmp_path / "case", tmp_path / "public"
    case_dir.mkdir()
    public_dir.mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("MFC_PUBLIC_MATERIAL_DIR", str(public_dir))
    for directory, value in ((tmp_path, 6.0), (case_dir, 7.0), (public_dir, 8.0)):
        _material(directory, {**_JWL, "A": value})
    params = {"fluid_pp(1)%material_file": "products.yaml"}
    for directory, expected in ((tmp_path, 6.0), (case_dir, 7.0), (public_dir, 8.0)):
        assert resolve_materials(params, str(case_dir))["fluid_pp(1)%jwl_a"] == expected
        (directory / "products.yaml").unlink()


@pytest.mark.parametrize(
    ("case", "parameters", "release_status", "error"),
    [
        ({"fluid_pp(2)%jwl_a": 9.0}, _JWL, "public", "conflicts with fluid_pp\\(2\\)%jwl_a"),
        ({"fluid_pp(2)%qv": 2.0}, {**_JWL, "qv": 1.0}, "public", "conflicts with fluid_pp\\(2\\)%qv"),
        ({"fluid_pp(2)%eos": "ideal_gas"}, _JWL, "public", "conflicts with fluid_pp\\(2\\)%eos"),
        ({}, _JWL, "", "nonempty provenance citation"),
        ({}, {**_JWL, "A": ".nan"}, "public", "finite numeric"),
        ({}, {**_JWL, "a": 1.0}, "public", "names a parameter twice"),
        ({"reactive_burn": "T", "fluid_pp(1)%qv": 1.0}, {**_JWL, "Q": 2.0}, "public", "sets both reactant and product qv"),
    ],
)
def test_material_rejects_invalid_input(tmp_path, case, parameters, release_status, error):
    path = _material(tmp_path, parameters, release_status=release_status)
    with pytest.raises(MFCException, match=error):
        resolve_materials({**case, "fluid_pp(2)%material_file": path}, str(tmp_path))


def test_material_rejects_misplaced_q_and_missing_file(tmp_path):
    _material(tmp_path)
    with pytest.raises(MFCException, match="belongs to the products"):
        resolve_materials({"reactive_burn": "T", "fluid_pp(1)%material_file": "products.yaml"}, str(tmp_path))
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
    loaded = resolve_materials({"fluid_pp(1)%material_file": _material(tmp_path, {**parameters, "cv": 2.0, "qv": 3.0}, family)}, str(tmp_path))
    expected = {"eos": family, "cv": 2.0, "qv": 3.0, **{f"{prefix}_{name}": value for name, value in parameters.items()}}
    assert {key: loaded[f"fluid_pp(1)%{key}"] for key in expected} == expected
