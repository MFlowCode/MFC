"""R1: proves the EOS_FAMILIES registry changed nothing observable.

Expected values are hard-coded (captured from the tree before the registry existed), not derived
from the registry itself — asserting the registry against itself would be circular.
"""

from ..params import REGISTRY
from ..params.definitions import CONSTRAINTS

_EXPECTED_NAMES = {"stiffened_gas": 1, "ideal_gas": 2, "mie_gruneisen": 3, "jwl": 4, "vinet": 5}
_EXPECTED_LABELS = {1: "stiffened-gas", 2: "ideal-gas", 3: "Mie-Gruneisen", 4: "JWL", 5: "Vinet"}
_EXPECTED_CHOICES = [1, 2, 3, 4, 5]

# Captured from `sorted(n for n in REGISTRY.all_params if n.startswith("fluid_pp(1)%"))`
# on the unmodified tree (HEAD 949a131c), before eos_families.py existed.
_EXPECTED_FLUID_PP1_PARAMS = [
    "fluid_pp(1)%G",
    "fluid_pp(1)%K",
    "fluid_pp(1)%Re(1)",
    "fluid_pp(1)%Re(2)",
    "fluid_pp(1)%cv",
    "fluid_pp(1)%eos",
    "fluid_pp(1)%gamma",
    "fluid_pp(1)%hb_m",
    "fluid_pp(1)%jwl_a",
    "fluid_pp(1)%jwl_b",
    "fluid_pp(1)%jwl_omega",
    "fluid_pp(1)%jwl_r1",
    "fluid_pp(1)%jwl_r2",
    "fluid_pp(1)%jwl_rho0",
    "fluid_pp(1)%jwl_t0",
    "fluid_pp(1)%mg_c0",
    "fluid_pp(1)%mg_gruneisen",
    "fluid_pp(1)%mg_gruneisen_a",
    "fluid_pp(1)%mg_rho0",
    "fluid_pp(1)%mg_s",
    "fluid_pp(1)%mg_s2",
    "fluid_pp(1)%mg_s3",
    "fluid_pp(1)%mg_t0",
    "fluid_pp(1)%mu_bulk",
    "fluid_pp(1)%mu_max",
    "fluid_pp(1)%mu_min",
    "fluid_pp(1)%nn",
    "fluid_pp(1)%non_newtonian",
    "fluid_pp(1)%pi_inf",
    "fluid_pp(1)%qv",
    "fluid_pp(1)%qvp",
    "fluid_pp(1)%tau0",
    "fluid_pp(1)%vinet_gruneisen",
    "fluid_pp(1)%vinet_gruneisen_a",
    "fluid_pp(1)%vinet_k0",
    "fluid_pp(1)%vinet_k0p",
    "fluid_pp(1)%vinet_rho0",
    "fluid_pp(1)%vinet_t0",
]


def test_eos_names_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["names"] == _EXPECTED_NAMES


def test_eos_value_labels_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["value_labels"] == _EXPECTED_LABELS


def test_eos_choices_unchanged():
    assert CONSTRAINTS["fluid_pp(1)%eos"]["choices"] == _EXPECTED_CHOICES


def test_fluid_pp_parameter_set_unchanged():
    names = sorted(n for n in REGISTRY.all_params if n.startswith("fluid_pp(1)%"))
    assert names == _EXPECTED_FLUID_PP1_PARAMS
