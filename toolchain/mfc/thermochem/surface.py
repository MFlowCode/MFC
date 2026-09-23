"""Generate the heterogeneous surface-chemistry module for reacting immersed boundaries.

The interface mechanism is a bulk solid reacting with gas species: rates depend on gas
concentrations only, and the solid enters the reaction enthalpy through its NASA7 thermo.
Adsorbed-site balances are not modelled, so mechanisms that need them are refused.
"""

from functools import partial
from pathlib import Path

import cantera as ct
import pymbolic.primitives as p
from mako.template import Template

from . import expressions
from .fortran import OFFLOAD_DIRECTIVES, FortranExpressionMapper, check_options, float_to_fortran, wrap_code


def _net_stoich(reaction, name):
    return reaction.products.get(name, 0.0) - reaction.reactants.get(name, 0.0)


def _bulk_species(gas, surface):
    """Species of the adjacent non-gas phases, by name."""
    return {species.name: species for phase in surface.adjacent.values() for species in phase.species() if species.name not in gas.species_names}


def validate_surface_mechanism(gas, surface):
    """Reject rate laws the generator would otherwise emit silently wrong."""
    bulk = _bulk_species(gas, surface)
    for i, reaction in enumerate(surface.reactions()):
        label = f"Surface reaction {i + 1} ({reaction.equation})"
        rate = reaction.rate
        # An allowlist, not a hasattr() probe: sticking and Blowers-Masel rates also expose
        # pre_exponential_factor, but for a sticking rate it is a dimensionless probability and for
        # Blowers-Masel the activation energy is shifted by the reaction enthalpy at runtime.
        if not isinstance(rate, (ct.ArrheniusRate, ct.InterfaceArrheniusRate)):
            raise ValueError(f"{label}: rate type {rate.type!r} is not supported; only (interface-)Arrhenius rates are")
        if dict(getattr(rate, "coverage_dependencies", {}) or {}):
            raise ValueError(f"{label}: coverage-dependent rates are not supported")
        if reaction.reversible:
            raise ValueError(f"{label}: reversible surface reactions are not supported; write the reaction with '=>'")
        if rate.pre_exponential_factor <= 0:
            raise ValueError(f"{label}: Arrhenius pre-exponential factors must be positive")
        for name in sorted(set(reaction.reactants) | set(reaction.products) | set(reaction.orders)):
            if name in gas.species_names:
                continue
            if name in surface.species_names:
                raise ValueError(f"{label}: surface-site species {name!r} are not supported; the generator has no site balance")
            if name not in bulk:
                raise ValueError(f"{label}: species {name!r} is not in the gas or an adjacent phase")
            if _net_stoich(reaction, name) != 0 and not isinstance(bulk[name].thermo, ct.NasaPoly2):
                raise ValueError(f"{label}: bulk species {name!r} needs NASA7 thermo for the reaction enthalpy")


def _rate_of_progress_expr(gas, reaction, concentrations, temperature):
    """k(T) times gas concentrations raised to their orders; bulk species have unit activity."""
    orders = dict(reaction.reactants)
    orders.update(reaction.orders)
    progress = expressions.rate_coefficient_expr(reaction.rate, temperature)
    for name, order in orders.items():
        if name in gas.species_names and order != 0:
            progress = progress * concentrations[gas.species_index(name)] ** order
    return progress


def _production_rate_expr(surface, species, rates_of_progress):
    return sum(_net_stoich(r, species) * rates_of_progress[i] for i, r in enumerate(surface.reactions()) if _net_stoich(r, species) != 0)


def _reaction_enthalpy_rt_expr(gas, reaction, h_rt, bulk_h_rt):
    """Delta H / (R T) of one reaction from gas and bulk standard-state enthalpies."""
    terms = [(name, h_rt[k]) for k, name in enumerate(gas.species_names)] + list(bulk_h_rt.items())
    return sum(_net_stoich(reaction, name) * h for name, h in terms if _net_stoich(reaction, name) != 0)


def generate_surface_fortran(gas, surface=None, module_name="m_surface_thermochem", scalar_type="real(dp)", offload=None):
    """Emit get_surface_net_production_rates and get_surface_reaction_heat_flux.

    Without a surface mechanism the module still exists, with both routines returning zero, because
    the simulation imports it unconditionally.
    """
    check_options(module_name, scalar_type, offload)
    kind = "sp" if scalar_type == "real(sp)" else "dp"
    temperature = p.Variable("temperature")
    reactions = []
    bulk = []
    production_rates = []
    if surface is not None:
        validate_surface_mechanism(gas, surface)
        reactions = surface.reactions()
        species = _bulk_species(gas, surface)
        participating = sorted({name for r in reactions for name in set(r.reactants) | set(r.products) if name in species and _net_stoich(r, name) != 0})
        bulk = [(name, expressions.poly_to_enthalpy_expr(species[name].thermo, "temperature")) for name in participating]
        for k, name in enumerate(gas.species_names):
            expr = _production_rate_expr(surface, name, p.Variable("rates_of_progress"))
            if expr != 0:
                production_rates.append((k, expr))

    bulk_h_rt = {name: p.Variable("bulk_h_rt")[j] for j, (name, _) in enumerate(bulk)}
    template = Template(filename=str(Path(__file__).with_name("surface.f90.mako")))
    return wrap_code(
        template.render(
            module_name=module_name,
            real_type=scalar_type,
            kind=kind,
            gpu_routine=f"#define GPU_ROUTINE(name) {OFFLOAD_DIRECTIVES[offload]}",
            cgm=FortranExpressionMapper(kind),
            float_to_fortran=partial(float_to_fortran, kind=kind),
            gas=gas,
            reactions=reactions,
            bulk=bulk,
            rate_of_progress=lambda r: _rate_of_progress_expr(gas, r, p.Variable("concentrations"), temperature),
            production_rates=production_rates,
            reaction_enthalpy_rt=lambda r: _reaction_enthalpy_rt_expr(gas, r, p.Variable("h_rt"), bulk_h_rt),
        )
    )
