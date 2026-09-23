"""Compile the generated surface-chemistry module and compare it with Cantera's interface kinetics."""

import subprocess

import cantera as ct
import numpy as np
import pytest

from mfc.test_thermochem import ROOT, compile_kernel
from mfc.thermochem import generate_fortran, generate_surface_fortran

EXAMPLE = ROOT / "examples/2D_ibm_reacting_surface"
GAS = str(EXAMPLE / "carbon_gasphase_reduced_gri11.yaml")
SURFACE = EXAMPLE / "carbon_surface_bradley_11species.yaml"

DRIVER = """
program reference
    use m_thermochem
    use m_surface_thermochem
    implicit none
    integer :: ierr
    real(KIND) :: t, pressure, rho, q
    real(KIND) :: y(num_species), omega(num_species)
    do
        read(*,*,iostat=ierr) t, pressure, y
        if (ierr /= 0) exit
        call get_density(pressure, t, y, rho)
        call get_surface_net_production_rates(rho, t, y, omega)
        call get_surface_reaction_heat_flux(rho, t, y, q)
        write(*,'(*(ES25.16E3,1X))') omega, q
    end do
end program
"""


def carbon_surface(reactions=None):
    """The example mechanism, optionally with its reactions replaced by a YAML list."""
    text = SURFACE.read_text()
    if reactions is not None:
        text = text[: text.index("bradley-porous-reactions:")] + "bradley-porous-reactions:\n" + reactions
    return ct.Interface(yaml=text, name="carbon_surface")


def single_reaction(mechanism, phase, index):
    """One reaction of a Cantera surface mechanism, so each guard is exercised on its own."""
    full = ct.Interface(mechanism, phase)
    return ct.Interface(thermo="ideal-surface", species=full.species(), kinetics="surface", reactions=[full.reaction(index)], adjacent=[full.adjacent["gas"]])


def run_kernel(tmp_path, gas, surface, precision="dp", offload=None):
    """Compile m_thermochem + m_surface_thermochem and evaluate them at a set of gas states."""
    scalar_type = f"real({precision})"
    source = generate_fortran(gas, scalar_type=scalar_type, offload=offload)
    driver = generate_surface_fortran(gas, surface, scalar_type=scalar_type, offload=offload) + DRIVER
    executable = compile_kernel(tmp_path, gas, precision, offload, source=source, driver_source=driver)
    rng = np.random.default_rng(7)
    # Both sides of graphite's NASA7 midpoint (1000 K) and a state with no O/OH to consume.
    states = []
    for temperature in (300.0, 999.0, 1001.0, 1800.0, 2800.0):
        y = rng.uniform(0.01, 1, gas.n_species)
        states.append((temperature, ct.one_atm, y / y.sum()))
    gas.TPX = 1200.0, 2 * ct.one_atm, "O2:1,N2:3.76,H2O:0.5"
    states.append((gas.T, gas.P, gas.Y))
    inputs = "\n".join(" ".join(map(str, [t, pres, *y])) for t, pres, y in states) + "\n"
    result = subprocess.run([str(executable)], input=inputs, capture_output=True, text=True, check=True)
    return states, np.array([np.fromstring(line, sep=" ") for line in result.stdout.splitlines()])


def cantera_reference(surface, t, pressure, y):
    """Gas-species net production rates and the heat released, from Cantera."""
    for phase in surface.adjacent.values():
        phase.TP = t, pressure
    surface.adjacent["gas"].TPY = t, pressure, y
    surface.TP = t, pressure
    gas = surface.adjacent["gas"]
    omega = [surface.net_production_rates[surface.kinetics_species_index(name)] for name in gas.species_names]
    # NASA7 enthalpies for every species: the bulk's reference-pressure enthalpy, without the v*(P - P_ref)
    # term Cantera's condensed-phase standard state adds (~5e-6 of the reaction enthalpy for graphite).
    species = {sp.name: sp for phase in surface.adjacent.values() for sp in phase.species()}
    delta_h = [sum(nu * species[name].thermo.h(t) for name, nu in r.products.items()) - sum(nu * species[name].thermo.h(t) for name, nu in r.reactants.items()) for r in surface.reactions()]
    q = -np.dot(delta_h, surface.net_rates_of_progress)
    return np.array([*omega, q])


@pytest.mark.parametrize("precision,offload", [("dp", None), ("sp", None), ("dp", "acc"), ("dp", "mp")])
def test_carbon_surface_matches_cantera(tmp_path, precision, offload):
    surface = carbon_surface()
    gas = ct.Solution(GAS)
    states, rows = run_kernel(tmp_path, gas, surface, precision, offload)
    assert len(rows) == len(states)
    rtol = 3e-5 if precision == "sp" else 1e-12
    for actual, state in zip(rows, states):
        expected = cantera_reference(surface, *state)
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=rtol * np.max(np.abs(expected)))


def test_surface_mass_is_the_carbon_gasified(tmp_path):
    """Net gas-phase mass production equals the carbon mass leaving the solid: 12.011 per C(gr)."""
    surface = carbon_surface()
    gas = ct.Solution(GAS)
    states, rows = run_kernel(tmp_path, gas, surface)
    for actual, state in zip(rows, states):
        cantera_reference(surface, *state)
        carbon = -surface.net_production_rates[surface.kinetics_species_index("C(gr)")]
        np.testing.assert_allclose(actual[:-1] @ gas.molecular_weights, carbon * surface.adjacent["graphite"].molecular_weights[0], rtol=1e-12)


def test_without_a_surface_mechanism_the_module_returns_zero(tmp_path):
    """The simulation imports the module unconditionally, so chemistry without a surface still builds."""
    gas = ct.Solution("h2o2.yaml")
    _, rows = run_kernel(tmp_path, gas, None)
    assert not rows.any()


@pytest.mark.parametrize(
    "index,match",
    [
        (2, "rate type 'sticking-Arrhenius'"),
        (1, "coverage-dependent"),
        (11, "reversible"),
        (0, "surface-site species"),
    ],
)
def test_unsupported_ptcombust_reactions_are_refused(index, match):
    """Each would otherwise be emitted as a wrong rate law, silently."""
    surface = single_reaction("ptcombust.yaml", "Pt_surf", index)
    with pytest.raises(ValueError, match=match):
        generate_surface_fortran(ct.Solution("ptcombust.yaml", "gas"), surface)


def test_a_reversible_bulk_reaction_is_refused():
    """The reverse branch needs interface equilibrium constants, which are not generated."""
    surface = carbon_surface("- equation: C(gr) + OH <=> CO + H\n  rate-constant: {A: 1.65, b: 0.5, Ea: 0.0}\n")
    with pytest.raises(ValueError, match="reversible"):
        generate_surface_fortran(ct.Solution(GAS), surface)
