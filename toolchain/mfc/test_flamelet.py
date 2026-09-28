"""Physical invariants of the Cantera-only mixing-layer initialization."""

import cantera as ct
import numpy as np
import pytest

from mfc.flamelet import create_simulation_fields, density, mechanism_fingerprint, streams
from mfc.test_thermochem import ROOT


@pytest.mark.parametrize("cold", [True, False])
@pytest.mark.parametrize("mechanism,fuel_fraction", [("h2o2.yaml", 1.0), (str(ROOT / "examples/3D_reacting_mixing_layer/sandiego.yaml"), 0.5)])
def test_mixing_layer_profiles(mechanism, fuel_fraction, cold):
    gas = ct.Solution(mechanism)
    pressure = ct.one_atm
    oxidizer, fuel, _ = streams(gas, "H2", pressure, 500, 300, 0.21, fuel_fraction, 0.001, 0.3)
    coords = np.linspace(-0.01, 0.01, 401)
    fields = create_simulation_fields(gas, pressure, 500, 300, coords, 0.001, oxidizer, fuel, 100, cold)
    y = fields.mass_fractions
    assert np.all(np.isfinite(y)) and np.min(y) >= 0
    np.testing.assert_allclose(y.sum(axis=0), 1, atol=1e-12)
    np.testing.assert_allclose(y[:, 0], fuel[1], atol=1e-12)
    np.testing.assert_allclose(y[:, -1], oxidizer[1], atol=1e-12)
    np.testing.assert_allclose(fields.temperature[[0, -1]], [300, 500], atol=1e-6)
    np.testing.assert_allclose(fields.velocity[[0, -1]], [fuel[4], oxidizer[4]], atol=1e-8)
    rho = density(gas, fields.pressure, fields.temperature, y)
    enthalpy = oxidizer[0] + (fuel[0] - oxidizer[0]) * fields.mixture_fraction
    atoms = np.array([[gas.n_atoms(k, el) * gas.atomic_weight(el) / gas.molecular_weights[k] for k in range(gas.n_species)] for el in gas.element_names])
    mixed_y = oxidizer[1][:, None] + (fuel[1] - oxidizer[1])[:, None] * fields.mixture_fraction
    # Unity-Lewis counterflow conserves the stream elemental mixture fractions.
    np.testing.assert_allclose(atoms @ y, atoms @ mixed_y, atol=2e-5)
    for i in range(len(coords)):
        gas.TPY = fields.temperature[i], pressure, y[:, i]
        assert gas.density == pytest.approx(rho[i], rel=1e-12)
        # Cantera's HP inversion tolerance is scaled by cp*T, including near h=0.
        assert gas.enthalpy_mass == pytest.approx(enthalpy[i], abs=1e-8 * gas.cp_mass * gas.T)
    if cold:
        np.testing.assert_allclose(y, mixed_y, atol=1e-14)
        assert fields.temperature.max() <= 500 + 1e-6
    else:
        assert fields.temperature.max() > 1500
        assert y[gas.species_index("H2O")].max() > 0.05


def test_mechanism_fingerprint_ignores_state_but_tracks_rates():
    gas = ct.Solution("h2o2.yaml")
    original = mechanism_fingerprint(gas)
    gas.TPX = 1500, 2 * ct.one_atm, "H2:1"
    assert mechanism_fingerprint(gas) == original
    reaction = gas.reaction(2)
    rate = reaction.rate
    reaction.rate = ct.ArrheniusRate(2 * rate.pre_exponential_factor, rate.temperature_exponent, rate.activation_energy)
    gas.modify_reaction(2, reaction)
    assert mechanism_fingerprint(gas) != original
