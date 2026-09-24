__copyright__ = """
Copyright (C) 2020 University of Illinois Board of Trustees
"""

__license__ = """
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from functools import singledispatch

import cantera as ct
import numpy as np
import pymbolic.primitives as p

# {{{ polynomial processing


def nasa7_conditional(t, poly, part_gen):
    # FIXME: Should check minTemp, maxTemp
    return p.If(
        p.Comparison(t, ">", poly.coeffs[0]),
        part_gen(poly.coeffs[1:8], t),
        part_gen(poly.coeffs[8:15], t),
    )


@singledispatch
def poly_to_expr(poly):
    raise TypeError(f"unexpected argument type in poly_to_expr: {type(poly)}")


@poly_to_expr.register
def _(poly: ct.NasaPoly2, arg_name):
    def gen(c, t):
        assert len(c) == 7
        return c[0] + c[1] * t + c[2] * t**2 + c[3] * t**3 + c[4] * t**4

    return nasa7_conditional(p.Variable(arg_name), poly, gen)


@singledispatch
def poly_to_enthalpy_expr(poly, arg_name):
    raise TypeError("unexpected argument type in poly_to_enthalpy_expr: " f"{type(poly)}")


@poly_to_enthalpy_expr.register
def _(poly: ct.NasaPoly2, arg_name):
    def gen(c, t):
        assert len(c) == 7
        return c[0] + c[1] / 2 * t + c[2] / 3 * t**2 + c[3] / 4 * t**3 + c[4] / 5 * t**4 + c[5] / t

    return nasa7_conditional(p.Variable(arg_name), poly, gen)


@singledispatch
def poly_to_entropy_expr(poly, arg_name):
    raise TypeError("unexpected argument type in poly_to_entropy_expr: " f"{type(poly)}")


@poly_to_entropy_expr.register
def _(poly: ct.NasaPoly2, arg_name):
    log = p.Variable("log")

    def gen(c, t):
        assert len(c) == 7
        return c[0] * log(t) + c[1] * t + c[2] / 2 * t**2 + c[3] / 3 * t**3 + c[4] / 4 * t**4 + c[6]

    return nasa7_conditional(p.Variable(arg_name), poly, gen)


# {{{ Data-handling helper


def _zeros_like(argument):
    # FIXME: This mishandles NaNs.
    return 0 * argument


# }}}


# {{{ Transport polynomials & mixture rules


def viscosity_polynomial_expr(c, t):
    """Generate code for viscosity polynomials

    :returns: Viscosity polynomial expression with coefficients c in terms of
    the temperature t as a :class:`pymbolic.primitives.Expression`.
    """
    assert len(c) == 5
    return p.Variable("sqrt")(t) * (c[0] + c[1] * p.Variable("log")(t) + c[2] * p.Variable("log")(t) ** 2 + c[3] * p.Variable("log")(t) ** 3 + c[4] * p.Variable("log")(t) ** 4) ** 2


def conductivity_polynomial_expr(c, t):
    """Generate code for conductivity polynomials

    :returns: Conductivity polynomial expression with coefficients c in terms
    of the temperature t as a :class:`pymbolic.primitives.Expression`.
    """
    assert len(c) == 5
    return p.Variable("sqrt")(t) * (c[0] + c[1] * p.Variable("log")(t) + c[2] * p.Variable("log")(t) ** 2 + c[3] * p.Variable("log")(t) ** 3 + c[4] * p.Variable("log")(t) ** 4)


def diffusivity_polynomial_expr(c, t):
    """Generate code for diffusivity polynomials

    :returns: Diffusivity polynomial expression with coefficients c in terms
    of the temperature t as a :class:`pymbolic.primitives.Expression`.
    """
    assert len(c) == 5
    return p.Variable("sqrt")(t) * t * (c[0] + c[1] * p.Variable("log")(t) + c[2] * p.Variable("log")(t) ** 2 + c[3] * p.Variable("log")(t) ** 3 + c[4] * p.Variable("log")(t) ** 4)


def viscosity_mixture_rule_wilke_expr(sol: ct.Solution, sp, x, mu):
    """Generate code for species mixture rule. See [Kee_2003]_, chapter 12.

    :returns: Expression for the Wilke viscosity mixture rule
        for species *sp* in terms of species mole fractions *w*
        and viscosities *mu* as a :class:`pymbolic.primitives.Expression`
    """
    w = sol.molecular_weights
    sqrt = p.Variable("sqrt")
    return sum([x[j] * (1 + sqrt((mu[sp] / mu[j]) * np.sqrt(w[j] / w[sp]))) ** 2 / np.sqrt(8 * (1 + (w[sp] / w[j]))) for j in range(sol.n_species)])


def diffusivity_mixture_rule_denom_expr(sol: ct.Solution, j_sp, x, bdiff):
    """See [Kee_2003]_, chapter 12 for details.
    :returns: The denominator expression to the mixture rule
    for mixture-averaged species diffusivities in terms
    of the species mole fractions *x* and binary diffusivities *bdiff* as a
    :class:`pymbolic.primitives.Expression`
    """
    return sum(x[i_sp] / bdiff[i_sp][j_sp] for i_sp in range(sol.n_species))


# }}}


# {{{ Equilibrium constants


def equilibrium_constants_expr(sol: ct.Solution, reaction_index, gibbs_rt):
    """Generate code for equilibrium constants.

    :returns: Equilibrium constant expression for reaction with
        index *reaction_index* in terms of the species Gibbs
        functions *gibbs_rt* as a :class:`pymbolic.primitives.Expression`
    """
    indices_reac = [sol.species_index(sp) for sp in sol.reaction(reaction_index).reactants]
    indices_prod = [sol.species_index(sp) for sp in sol.reaction(reaction_index).products]

    # Stoichiometric coefficients
    nu_reac = [sol.reactant_stoich_coeff(sol.species_index(sp), reaction_index) for sp in sol.reaction(reaction_index).reactants]
    nu_prod = [sol.product_stoich_coeff(sol.species_index(sp), reaction_index) for sp in sol.reaction(reaction_index).products]

    sum_r = sum(nu_reac_i * gibbs_rt[indices_reac_i] for indices_reac_i, nu_reac_i in zip(indices_reac, nu_reac))
    sum_p = sum(nu_prod_i * gibbs_rt[indices_prod_i] for indices_prod_i, nu_prod_i in zip(indices_prod, nu_prod))

    # Check if reaction is termolecular
    sum_nu_net = sum(nu_prod) - sum(nu_reac)
    if sum_nu_net != 0:
        return sum_p - sum_r - sum_nu_net * p.Variable("c0")
    else:
        return sum_p - sum_r


# }}}


# {{{ Rate coefficients


def rate_coefficient_expr(rate_coeff: ct.Arrhenius, t):
    """
    :returns: The rate coefficient expression for *rate_coeff* in terms
        of the temperature *t* as a :class:`pymbolic.primitives.Expression`
    """
    # Rate parameters
    a = rate_coeff.pre_exponential_factor
    b = rate_coeff.temperature_exponent
    t_a = rate_coeff.activation_energy / ct.gas_constant
    if t_a == 0:
        # Weakly temperature-dependent rate
        return a * t**b
    else:
        # Modified Arrhenius
        return p.Variable("exp")(np.log(a) + b * p.Variable("log")(t) - t_a / t)


def third_body_efficiencies_expr(sol: ct.Solution, react: ct.Reaction, c):
    """
    :returns: The third-body concentration expression for reaction *react* in
    terms of the species concentrations *c* as a
        :class:`pymbolic.primitives.Expression`
    """

    efficiencies = [react.third_body.efficiencies[sp] for sp in react.third_body.efficiencies]
    indices_nondef = [sol.species_index(sp) for sp in react.third_body.efficiencies]
    indices_default = [i for i in range(sol.n_species) if i not in indices_nondef]
    sum_nondef = sum(eff_i * c[index_i] for eff_i, index_i in zip(np.array(efficiencies), indices_nondef))
    sum_default = react.third_body.default_efficiency * sum(c[i] for i in indices_default)
    return sum_nondef + sum_default


def troe_falloff_center_expr(react: ct.Reaction, t):
    """
    :returns: The Troe falloff center expression for reaction *react* in
    terms of the temperature *t* as a
    :class:`pymbolic.primitives.Expression`
    """

    if isinstance(react.rate, ct.TroeRate):
        troe_params = react.rate.falloff_coeffs
    elif isinstance(react.rate, ct.LindemannRate):
        return 1
    else:
        raise ValueError("Unexpected value of 'rate.type': " f" '{react.rate.type}'")

    troe_1 = (1.0 - troe_params[0]) * p.Variable("exp")(-t / troe_params[1])
    troe_2 = troe_params[0] * p.Variable("exp")(-t / troe_params[2])
    if len(troe_params) == 3:
        return p.Variable("log10")(troe_1 + troe_2)
    elif len(troe_params) == 4:
        troe_3 = p.Variable("exp")(-troe_params[3] / t)
        return p.Variable("log10")(troe_1 + troe_2 + troe_3)
    else:
        raise ValueError("Unexpected length of 'troe_params': " f" '{len(troe_params)}'")
    return


def troe_falloff_factor_expr(react: ct.Reaction, i, red_pressure, falloff_center):
    """
    :returns: The Troe falloff factor expression for reaction
    *react* in terms of reduced pressure *red_pressure* and the
    falloff center *falloff_center* as a
    :class:`pymbolic.primitives.Expression`

    """
    if isinstance(react.rate, ct.TroeRate):
        # Fortran MERGE need not short-circuit: guard the logarithm's argument.
        positive = p.Comparison(red_pressure[i], ">", 0)
        log_rp = p.Variable("log10")(p.If(positive, red_pressure[i], 1))
        c = -0.4 - 0.67 * falloff_center[i]
        n = 0.75 - 1.27 * falloff_center[i]
        return p.If(positive, (log_rp + c) / (n - 0.14 * (log_rp + c)), -1 / 0.14)
    elif isinstance(react.rate, ct.LindemannRate):
        return 0
    else:
        raise ValueError("Unexpected value of 'rate.type': " f" '{react.rate.type}'")


def falloff_function_expr(react: ct.Reaction, i, falloff_factor, falloff_center):
    """
    :returns: Falloff function expression for reaction *react* in
    terms of the temperature *t*, falloff width factor
    *falloff_factor*, and falloff center *falloff_center* as a
    :class:`pymbolic.primitives.Expression`

    """

    falloff_type = react.reaction_type.split("-")[1]

    if falloff_type == "Troe":
        return 10 ** (falloff_center[i] / (1 + falloff_factor[i] ** 2))
    elif falloff_type == "Lindemann":
        return 1
    else:
        raise ValueError("Unexpected value of 'falloff_type': " f" '{falloff_type}'")


# }}}


# {{{ Rates of progress


def rate_of_progress_expr(sol: ct.Solution, reaction_index, c, k_fwd, log_k_eq):
    """
    :returns: Rate of progress expression for reaction with
    index *reaction_index* in terms of species concentrations *c*
    with rate coefficients *k_fwd* and equilbrium constants *k_eq*
    as a :class:`pymbolic.primitives.Expression`
    """
    indices_reac = [sol.species_index(sp) for sp in sol.reaction(reaction_index).reactants]
    indices_prod = [sol.species_index(sp) for sp in sol.reaction(reaction_index).products]

    if sol.reaction(reaction_index).orders:
        nu_reac = [sol.reaction(reaction_index).orders[sp] for sp in sol.reaction(reaction_index).orders]
    else:
        nu_reac = [sol.reaction(reaction_index).reactants[sp] for sp in sol.reaction(reaction_index).reactants]

    r_fwd = np.prod([c[index] ** nu for index, nu in zip(indices_reac, nu_reac)])

    if sol.reaction(reaction_index).reversible:
        nu_prod = [sol.reaction(reaction_index).products[sp] for sp in sol.reaction(reaction_index).products]
        r_rev = np.prod([c[index] ** nu for index, nu in zip(indices_prod, nu_prod)])
        return k_fwd[reaction_index] * (r_fwd - p.Variable("exp")(log_k_eq[reaction_index]) * r_rev)
    else:
        return k_fwd[reaction_index] * r_fwd


# }}}


# {{{ Species production rates


def production_rate_expr(sol: ct.Solution, species, r_net):
    """
    :returns: Species production rate for species *species* in terms of
        the net reaction rates of progress *r_net* as a
        :class:`pymbolic.primitives.Expression`
    """
    ones = _zeros_like(r_net[0]) + 1.0
    indices_fwd = [i for i, react in enumerate(sol.reactions()) if species in react.reactants]
    indices_rev = [i for i, react in enumerate(sol.reactions()) if species in react.products]
    nu_fwd = [sol.reactant_stoich_coeff(sol.species_index(species), react_index) for react_index in indices_fwd]
    nu_rev = [sol.product_stoich_coeff(sol.species_index(species), prod_index) for prod_index in indices_rev]
    sum_fwd = sum(nu * r_net[index] for nu, index in zip(nu_fwd, indices_fwd))
    sum_rev = sum(nu * r_net[index] for nu, index in zip(nu_rev, indices_rev))
    return (sum_rev - sum_fwd) * ones


def fwd_rate_of_progress_expr(sol: ct.Solution, reaction_index, c, k_fwd):
    """
    :returns: Forward rate of progress (non-negative) for reaction
        *reaction_index*, i.e. the forward part of
        :func:`rate_of_progress_expr`.
    """
    indices_reac = [sol.species_index(sp) for sp in sol.reaction(reaction_index).reactants]
    if sol.reaction(reaction_index).orders:
        nu_reac = [sol.reaction(reaction_index).orders[sp] for sp in sol.reaction(reaction_index).orders]
    else:
        nu_reac = [sol.reaction(reaction_index).reactants[sp] for sp in sol.reaction(reaction_index).reactants]
    r_fwd = np.prod([c[index] ** nu for index, nu in zip(indices_reac, nu_reac)])
    return k_fwd[reaction_index] * r_fwd


def rev_rate_of_progress_expr(sol: ct.Solution, reaction_index, c, k_fwd, log_k_eq):
    """
    :returns: Reverse rate of progress (non-negative) for reaction
        *reaction_index*; zero for irreversible reactions.
    """
    if not sol.reaction(reaction_index).reversible:
        return _zeros_like(c[0])
    indices_prod = [sol.species_index(sp) for sp in sol.reaction(reaction_index).products]
    nu_prod = [sol.reaction(reaction_index).products[sp] for sp in sol.reaction(reaction_index).products]
    r_rev = np.prod([c[index] ** nu for index, nu in zip(indices_prod, nu_prod)])
    return k_fwd[reaction_index] * p.Variable("exp")(log_k_eq[reaction_index]) * r_rev


def _species_reaction_stoich(sol: ct.Solution, species):
    """
    :returns: ``(idx_reactant, idx_product, nu_reactant, nu_product)`` -- the
        reaction indices in which *species* appears as a reactant and as a
        product, with the matching stoichiometric coefficients. Shared by the
        creation/destruction splits, which differ only in how these are summed.
    """
    si = sol.species_index(species)
    idx_reactant = [i for i, react in enumerate(sol.reactions()) if species in react.reactants]
    idx_product = [i for i, react in enumerate(sol.reactions()) if species in react.products]
    nu_reactant = [sol.reactant_stoich_coeff(si, i) for i in idx_reactant]
    nu_product = [sol.product_stoich_coeff(si, i) for i in idx_product]
    return idx_reactant, idx_product, nu_reactant, nu_product


def creation_rate_expr(sol: ct.Solution, species, r_fwd, r_rev):
    """
    :returns: Species creation rate for *species*: created as a product by
        forward reactions and as a reactant by reverse reactions. Mirrors
        Cantera's ``creation_rates``.
    """
    ones = _zeros_like(r_fwd[0]) + 1.0
    idx_reactant, idx_product, nu_reactant, nu_product = _species_reaction_stoich(sol, species)
    made = sum(nu * r_fwd[i] for nu, i in zip(nu_product, idx_product)) + sum(nu * r_rev[i] for nu, i in zip(nu_reactant, idx_reactant))
    return made * ones


def destruction_rate_expr(sol: ct.Solution, species, r_fwd, r_rev):
    """
    :returns: Species destruction rate for *species*: consumed as a reactant by
        forward reactions and as a product by reverse reactions. Mirrors
        Cantera's ``destruction_rates``. creation - destruction == net.
    """
    ones = _zeros_like(r_fwd[0]) + 1.0
    idx_reactant, idx_product, nu_reactant, nu_product = _species_reaction_stoich(sol, species)
    lost = sum(nu * r_fwd[i] for nu, i in zip(nu_reactant, idx_reactant)) + sum(nu * r_rev[i] for nu, i in zip(nu_product, idx_product))
    return lost * ones


# }}}

# vim
