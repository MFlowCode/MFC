"""Compile generated kernels and compare their numerical interface with Cantera."""

import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import cantera as ct
import numpy as np
import pytest

from mfc.thermochem import generate_fortran

ROOT = Path(__file__).resolve().parents[2]
MECHANISMS = [
    "h2o2.yaml",
    "gri30.yaml",
    str(ROOT / "examples/3D_reacting_mixing_layer/sandiego.yaml"),
    str(ROOT / "examples/2D_reactive_shock_bubble/h2o2_xe.yaml"),
]

DRIVER = """
program reference
    use m_thermochem
    implicit none
    integer :: ierr
    real(KIND) :: t, pressure, rho, cp, cv, h, e, mw, mu, lambda, te, th, pcheck
    real(KIND) :: y(num_species), omega(num_species), creation(num_species), destruction(num_species)
    real(KIND) :: diffusion(num_species), enthalpy(num_species)
    do
        read(*,*,iostat=ierr) t, pressure, y
        if (ierr /= 0) exit
        call get_density(pressure, t, y, rho)
        call get_pressure(rho, t, y, pcheck)
        call get_mixture_specific_heat_cp_mass(t, y, cp)
        call get_mixture_specific_heat_cv_mass(t, y, cv)
        call get_mixture_enthalpy_mass(t, y, h)
        call get_mixture_energy_mass(t, y, e)
        call get_mixture_molecular_weight(y, mw)
        call get_temperature(e, 0.9_KIND*t, y, .true., te)
        call get_temperature(h, 1.1_KIND*t, y, .false., th)
        call get_mixture_viscosity_mixavg(t, y, mu)
        call get_mixture_thermal_conductivity_mixavg(t, y, lambda)
        call get_species_mass_diffusivities_mixavg(pressure, t, y, diffusion)
        call get_species_enthalpies_rt(t, enthalpy)
        call get_net_production_rates(rho, t, y, omega)
        call get_creation_destruction_rates(rho, t, y, creation, destruction)
        write(*,'(*(ES25.16E3,1X))') rho, pcheck, cp, cv, h, e, mw, te, th, mu, lambda, &
            diffusion, enthalpy, omega, creation, destruction
    end do
end program
"""


def fypp(directory, name, source):
    """Preprocess generated Fypp source the way MFC's CMake build does."""
    executable = shutil.which("fypp") or str(Path(sys.executable).with_name("fypp"))
    fpp, f90 = directory / f"{name}.fpp", directory / f"{name}.f90"
    fpp.write_text(source)
    include = ["-I", str(ROOT / "src/common/include"), "-I", str(ROOT / "src/common")]
    defines = ["-D", 'MFC_COMPILER="GNU"', "-D", "MFC_CASE_OPTIMIZATION=False", "-D", "chemistry=False"]
    subprocess.run([executable, "-m", "re", *include, *defines, "--no-folding", "--line-length=999", str(fpp), str(f90)], check=True, capture_output=True, text=True)
    return f90


def compile_kernel(directory, gas, precision="dp", offload=None, *, source=None, driver_source=DRIVER, extra_flags=(), extra_sources=()):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required to validate generated Fortran")
    module = fypp(directory, "m_thermochem", source if source is not None else generate_fortran(gas))
    driver = directory / "driver.f90"
    driver.write_text(driver_source.replace("KIND", "wp").replace("use m_thermochem", "use m_precision_select, only: wp\n    use m_thermochem", 1))
    executable = directory / "reference"
    flags = {None: [], "acc": ["-fopenacc", "-DMFC_OpenACC"], "mp": ["-fopenmp", "-DMFC_OpenMP"]}[offload]
    flags += {"dp": [], "sp": ["-DMFC_SINGLE_PRECISION"]}[precision]
    sources = [ROOT / "src/common/m_precision_select.f90", *extra_sources, module, driver]
    subprocess.run(
        [compiler, "-cpp", "-O0", "-Wconversion", "-Werror=conversion", *flags, *extra_flags, *map(str, sources), "-o", str(executable)],
        cwd=directory,
        check=True,
        capture_output=True,
        text=True,
    )
    return executable


def reference_states(gas):
    rng = np.random.default_rng(23)
    for temperature, pressure in [(300, ct.one_atm), (999, 0.1 * ct.one_atm), (1001, ct.one_atm), (1800, 10 * ct.one_atm), (2800, ct.one_atm)]:
        y = rng.uniform(0.01, 1, gas.n_species)
        gas.TPY = temperature, pressure, y / y.sum()
        yield gas.T, gas.P, gas.Y
    # Zero concentrations and the pure-species diffusion fallback.
    gas.TPX = 800, ct.one_atm, "N2:1"
    yield gas.T, gas.P, gas.Y
    gas.TPX = 1200, ct.one_atm, "H2:2,O2:1,N2:3.76"
    yield gas.T, gas.P, gas.Y


def compare_kernel(executable, gas, precision="dp"):
    states = list(reference_states(gas))
    inputs = "\n".join(" ".join(map(str, [t, p, *y])) for t, p, y in states) + "\n"
    result = subprocess.run([str(executable)], input=inputs, capture_output=True, text=True, check=True)
    rows = np.array([np.fromstring(line, sep=" ") for line in result.stdout.splitlines()])
    assert len(rows) == len(states)
    rtol = 3e-5 if precision == "sp" else 2e-11
    for actual, (t, p, y) in zip(rows, states):
        gas.TPY = t, p, y
        diffusion = gas.mix_diff_coeffs.copy()
        # Preserve the old generator's pure-species self-diffusion convention.
        # Cantera's mix_diff_coeffs may instead return zero in this degenerate limit.
        pure = np.flatnonzero(y == 1)
        for k in pure:
            diffusion[k] = gas.binary_diff_coeffs[k, k]
        expected = np.concatenate(
            (
                [gas.density, gas.P, gas.cp_mass, gas.cv_mass, gas.enthalpy_mass, gas.int_energy_mass, gas.mean_molecular_weight, gas.T, gas.T, gas.viscosity, gas.thermal_conductivity],
                diffusion,
                gas.standard_enthalpies_RT,
                gas.net_production_rates,
                gas.creation_rates,
                gas.destruction_rates,
            )
        )
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=2e-8 if precision == "sp" else 1e-10)
        ns = gas.n_species
        net, creation, destruction = actual[11 + 2 * ns :].reshape(3, ns)
        np.testing.assert_allclose(net, creation - destruction, rtol=rtol, atol=rtol * max(1, np.max(creation)))
        elements = np.array([[gas.n_atoms(k, el) for k in range(ns)] for el in gas.element_names])
        residual = elements @ net
        assert np.max(np.abs(residual)) < rtol * max(1, np.max(np.abs(net)))


@pytest.mark.parametrize("mechanism", MECHANISMS)
def test_generated_mechanisms(tmp_path, mechanism):
    gas = ct.Solution(mechanism)
    compare_kernel(compile_kernel(tmp_path, gas), gas)


@pytest.mark.parametrize("precision,offload", [("sp", None), ("dp", "acc"), ("dp", "mp")])
def test_precision_and_directives(tmp_path, precision, offload):
    gas = ct.Solution("h2o2.yaml")
    compare_kernel(compile_kernel(tmp_path, gas, precision, offload), gas, precision)


def test_rejects_unsupported_thermo():
    species = ct.Species("N2", {"N": 2})
    species.thermo = ct.ConstantCp(200, 5000, ct.one_atm, [300, 0, 0, 30000])
    species.transport = ct.Solution("h2o2.yaml").species("N2").transport
    gas = ct.Solution(thermo="ideal-gas", kinetics="gas", species=[species], transport_model="mixture-averaged")
    with pytest.raises(ValueError, match="Species N2.*NASA7"):
        generate_fortran(gas)


def test_rejects_custom_orders():
    gas = ct.Solution("h2o2.yaml")
    reaction = ct.Reaction(equation="H2 + O => H + OH", rate=ct.ArrheniusRate(1e8, 0, 0))
    reaction.orders = {"H2": 0.5}
    custom = ct.Solution(thermo="ideal-gas", kinetics="gas", species=gas.species(), reactions=[reaction], transport_model="mixture-averaged")
    with pytest.raises(ValueError, match="Reaction 1.*custom reaction orders"):
        generate_fortran(custom)


@pytest.mark.parametrize("mode", ["double", "single", "mixed"])
def test_solver_working_precision(tmp_path, mode):
    """The module the toolchain writes compiles and agrees with Cantera in every precision mode."""
    from mfc.run import input as input_module

    case = input_module.MFCInputFile("case.py", str(tmp_path), {"chemistry": "T", "cantera_file": "h2o2.yaml"})
    case.get_fpp = lambda target: ""
    target = SimpleNamespace(name="simulation", isDependency=False, get_staging_dirpath=lambda case: str(tmp_path))
    case.generate_fpp(target)
    source = (tmp_path / "modules/simulation/m_thermochem.fpp").read_text()
    precision = "sp" if mode == "single" else "dp"
    flags = ["-DMFC_MIXED_PRECISION"] if mode == "mixed" else []
    gas = ct.Solution("h2o2.yaml")
    compare_kernel(compile_kernel(tmp_path, gas, precision, source=source, extra_flags=flags), gas, precision)


def test_long_species_names(tmp_path):
    names = ["nitrogen_reference_one", "nitrogen_reference_two"]
    nitrogen = ct.Solution("h2o2.yaml").species("N2")
    species = []
    for name in names:
        data = dict(nitrogen.input_data)
        data["name"] = name
        species.append(ct.Species.from_dict(data))
    gas = ct.Solution(thermo="ideal-gas", kinetics="gas", species=species, transport_model="mixture-averaged")
    driver = """
program names
    use m_thermochem
    implicit none
    character(len=50) :: name
    integer :: i, index
    do i = 1, num_species
        call get_species_name(i, name)
        call get_species_index(name, index)
        if (index /= i) stop 1
        print *, trim(name)
    end do
end program
"""
    executable = compile_kernel(tmp_path, gas, driver_source=driver)
    result = subprocess.run([str(executable)], capture_output=True, text=True, check=True)
    assert result.stdout.split() == names


@pytest.mark.parametrize("precision", ["sp", "dp"])
def test_zero_concentration_falloff(tmp_path, precision):
    gas = ct.Solution("h2o2.yaml")
    driver = """
program falloff
    use m_thermochem
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(KIND) :: concentrations(num_species), rates(num_reactions)
    concentrations = 0.0_KIND
    call get_fwd_rate_coefficients(1200.0_KIND, concentrations, rates)
    if (.not. all(ieee_is_finite(rates))) stop 1
    print *, rates
end program
"""
    executable = compile_kernel(tmp_path, gas, precision, driver_source=driver, extra_flags=["-ffpe-trap=invalid,zero,overflow"])
    result = subprocess.run([str(executable)], capture_output=True, text=True, check=True)
    rates = np.fromstring(result.stdout, sep=" ")
    for i, reaction in enumerate(gas.reactions()):
        if isinstance(reaction.rate, ct.FalloffRate):
            assert rates[i] == 0
