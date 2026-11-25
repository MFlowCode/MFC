"""Compile the production PR routines and check them against independent thermodynamics."""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def pr_driver(tmp_path_factory):
    if os.environ.get("MFC_SKIP_COMPILER_TESTS") == "1":
        pytest.skip("Fortran compilation disabled for precheck")
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the PR thermodynamic checks")
    directory = tmp_path_factory.mktemp("pr")
    source = (ROOT / "src/common/m_eos.fpp").read_text()
    routines = source[source.index("        subroutine s_pr_check_density") : source.rindex("    #:endif")]
    driver = directory / "pr.f90"
    driver.write_text(
        """module production_pr
use iso_fortran_env, only: wp => real64
implicit none
real(wp) :: pr_a, pr_b, pr_r, pr_c, pr_d, pr_cv, qvs(1)
contains
"""
        + routines
        + """
end module
program check_pr
use production_pr
implicit none
real(wp) :: rho, t, p, e, cv, cp, c2, omega
pr_r = 8314.46261815324_wp/28.014_wp
pr_cv = 2.5_wp*pr_r
omega = 0.0372_wp
pr_c = 0.37464_wp + 1.54226_wp*omega - 0.26992_wp*omega**2
pr_d = pr_c/sqrt(126.192_wp)
pr_c = 1._wp + pr_c
pr_a = 0.45724_wp*(pr_r*126.192_wp)**2/3395800._wp
pr_b = 0.07780_wp*pr_r*126.192_wp/3395800._wp
qvs = 0._wp
read(*,*) rho,t
call s_pr_properties(rho,t,p,e,cv,cp,c2)
write(*,'(7ES25.16)') p,e,cv,cp,c2,f_pr_temperature_energy(rho,e),f_pr_temperature_pressure(rho,p)
end program
"""
    )
    executable = directory / "pr"
    subprocess.run([compiler, "-ffree-line-length-none", "-fcheck=all", "-J", str(directory), str(driver), "-o", str(executable)], check=True, cwd=directory)
    return executable


def evaluate(driver, rho, temperature):
    result = subprocess.run([str(driver)], input=f"{rho} {temperature}\n", text=True, capture_output=True, check=True)
    return np.fromstring(result.stdout, sep=" ")


@pytest.mark.parametrize("rho,temperature", [(0.01, 200), (20, 150), (150, 200), (500, 200), (150, 400)])
def test_pr_round_trip_and_derivatives(pr_driver, rho, temperature):
    p, e, cv, cp, c2, te, tp = evaluate(pr_driver, rho, temperature)
    assert te == pytest.approx(temperature, rel=2e-13)
    assert tp == pytest.approx(temperature, rel=2e-13)
    dt = temperature * 1e-4
    dr = rho * 1e-4
    tplus, tminus = evaluate(pr_driver, rho, temperature + dt), evaluate(pr_driver, rho, temperature - dt)
    rplus, rminus = evaluate(pr_driver, rho + dr, temperature), evaluate(pr_driver, rho - dr, temperature)
    dpdt = (tplus[0] - tminus[0]) / (2 * dt)
    dpdr = (rplus[0] - rminus[0]) / (2 * dr)
    dedr = (rplus[1] - rminus[1]) / (2 * dr)
    assert cv == pytest.approx((tplus[1] - tminus[1]) / (2 * dt), rel=2e-8)
    assert dedr == pytest.approx((p - temperature * dpdt) / rho**2, rel=2e-6)
    assert c2 == pytest.approx(dpdr + temperature * dpdt**2 / (rho**2 * cv), rel=2e-8)
    assert cp == pytest.approx(cv + temperature * dpdt**2 / (rho**2 * dpdr), rel=2e-8)


def test_pr_against_cantera(pr_driver):
    ct = pytest.importorskip("cantera")
    species = ct.Species("N2", {"N": 2})
    cp = 3.5 * ct.gas_constant
    species.thermo = ct.ConstantCp(1, 5000, 101325, [298.15, cp * 298.15, 0, cp])
    species.update_user_data(
        {
            "equation-of-state": {
                "model": "Peng-Robinson",
                "acentric-factor": 0.0372,
                "a": 0.45724 * (ct.gas_constant * 126.192) ** 2 / 3395800,
                "b": 0.07780 * ct.gas_constant * 126.192 / 3395800,
            }
        }
    )
    gas = ct.Solution(thermo="Peng-Robinson", species=[species])
    for rho, temperature in [(20, 150), (150, 200), (500, 200), (150, 400)]:
        gas.TD = temperature, rho
        actual = evaluate(pr_driver, rho, temperature)
        expected = [gas.P, gas.int_energy_mass, gas.cv_mass, gas.cp_mass, gas.sound_speed**2]
        # Cantera uses higher precision PR critical constants to infer Tc from a and b.
        np.testing.assert_allclose(actual[:5], expected, rtol=8e-5)


@pytest.mark.parametrize("rho,temperature", [(-1, 200), (2000, 200), (150, -1), (150, 80)])
def test_pr_rejects_invalid_states(pr_driver, rho, temperature):
    result = subprocess.run([str(pr_driver)], input=f"{rho} {temperature}\n", text=True, capture_output=True, check=False)
    assert result.returncode != 0
    assert "PR" in result.stderr


@pytest.mark.parametrize(
    "override",
    [
        {"num_fluids": 2},
        {"model_eqns": 3},
        {"wave_speeds": 2},
        {"riemann_solver": 1},
        {"fluid_pp(1)%cv": 0},
        {"fluid_pp(1)%pr_tc": None},
        {"fluid_pp(1)%pr_pc": -1},
        {"fluid_pp(1)%pr_omega": 0.8},
        {"chemistry": "T"},
        {"ib": "T"},
        {"bc_x%beg": -6},
    ],
)
def test_pr_rejects_unsupported_configuration(override):
    from .case_validator import CaseValidator

    params = {
        "num_fluids": 1,
        "model_eqns": 2,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fluid_pp(1)%eos": 6,
        "fluid_pp(1)%pr_tc": 126.192,
        "fluid_pp(1)%pr_pc": 3395800,
        "fluid_pp(1)%pr_r": 296.8,
        "fluid_pp(1)%pr_omega": 0.0372,
        "fluid_pp(1)%cv": 742.0,
    }
    validator = CaseValidator(params)
    validator.check_peng_robinson()
    assert not validator.errors
    validator = CaseValidator(params | override)
    validator.check_peng_robinson()
    assert validator.errors


@pytest.mark.parametrize("variant", ["normal_ideal", "normal_pr", "supercritical_pr"])
def test_isothermal_boundary_layer_configuration(variant):
    import importlib.util

    from .case import Case
    from .case_validator import CaseValidator

    path = ROOT / "examples/2D_peng_robinson_boundary_layer/setup_case.py"
    spec = importlib.util.spec_from_file_location("boundary_layer_setup", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    case = Case(module.case(variant))
    case.validate_params()
    validator = CaseValidator(case.params)
    validator.validate("pre_process")
    validator.validate("simulation")
    assert not validator.errors
    state = module.reference(variant)
    # The existing diffusion timestep uses ideal cv, a lower bound for supported PR states.
    for t in (state["wall_temperature"], state["temperature"]):
        assert module.properties(state["rho"], t, state["ideal"])[2] >= module.CV


@pytest.mark.parametrize(
    "folder",
    [
        "2D_boundary_layer_ideal_gas",
        "2D_boundary_layer_pr_normal",
        "2D_boundary_layer_pr_supercritical",
        "2D_kelvin_helmholtz_ideal_gas",
        "2D_kelvin_helmholtz_pr_normal",
        "2D_kelvin_helmholtz_pr_supercritical",
    ],
)
def test_standalone_eos_example_configuration(folder, monkeypatch, tmp_path):
    import dataclasses
    import importlib.util
    import json
    import sys

    from . import state
    from .case import Case
    from .case_validator import CaseValidator

    monkeypatch.setattr(state, "gCFG", dataclasses.replace(state.CFG(), mpi=True))
    path = ROOT / "examples" / folder / "case.py"
    spec = importlib.util.spec_from_file_location("standalone_boundary_layer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    case = Case(module.case())
    case.validate_params()
    validator = CaseValidator(case.params)
    validator.validate("pre_process")
    validator.validate("simulation")
    assert not validator.errors
    result = subprocess.run([sys.executable, str(path)], cwd=tmp_path, capture_output=True, text=True, check=True)
    assert json.loads(result.stdout) == module.case()


def test_kh_xenon_thermodynamics():
    """Independently check the new KH fluid and its dilute/dense initial states."""
    import importlib.util

    ct = pytest.importorskip("cantera")
    path = ROOT / "examples/2D_kelvin_helmholtz_pr_supercritical/case.py"
    spec = importlib.util.spec_from_file_location("kh_xenon", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    species = ct.Species("Xe", {"Xe": 1})
    cp = 2.5 * ct.gas_constant
    species.thermo = ct.ConstantCp(1, 5000, 101325, [298.15, cp * 298.15, 0, cp])
    species.update_user_data(
        {
            "equation-of-state": {
                "model": "Peng-Robinson",
                "acentric-factor": module.OMEGA,
                "a": 0.45724 * (ct.gas_constant * module.TC) ** 2 / module.PC,
                "b": 0.07780 * ct.gas_constant * module.TC / module.PC,
            }
        }
    )
    gas = ct.Solution(thermo="Peng-Robinson", species=[species])
    for rho, temperature in [(2.66669133, 600), (5.33338266, 300), (487.86108, 410.7), (975.72216, 330)]:
        gas.TD = temperature, rho
        p, e, cv, cp, sound = module.properties(rho, temperature)
        np.testing.assert_allclose([p, e, cv, cp, sound], [gas.P, gas.int_energy_mass, gas.cv_mass, gas.cp_mass, gas.sound_speed], rtol=2e-4)
        assert module.temperature_from_energy(rho, e) == pytest.approx(temperature, rel=1e-12)
        assert module.temperature_from_pressure(rho, p) == pytest.approx(temperature, rel=1e-12)
