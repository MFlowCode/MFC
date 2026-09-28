"""Generate the mechanism-specific Fortran module used by MFC. See LICENSE.

Adapted from the Pyrometheus 1.1.1 Fortran emitter.
"""

import shlex
from functools import partial
from numbers import Integral
from pathlib import Path

import cantera as ct
import pymbolic.primitives as p
from mako.template import Template
from pymbolic.mapper.stringifier import PREC_CALL, PREC_NONE, PREC_PRODUCT, StringifyMapper

from . import expressions

# {{{ code generation helpers


def pad_fortran(line, width):
    line += " " * (width - 1 - len(line))
    line += "&"
    return line


def wrap_line_base(line, level=0, width=80, indentation="    ", pad_func=lambda string, amount: string, lex_func=None):
    """
    The input is a line of code at the given indentation level. Return the list
    of lines that results from wrapping the line to the given width. Lines
    subsequent to the first line in the returned list are padded with extra
    indentation. The initial indentation level is not included in the input or
    output lines.

    The `pad_func` argument is a function that adds line continuations. The
    `lex_func` argument returns the list of tokens in the line.
    """
    if lex_func is None:
        lex_func = partial(shlex.split, posix=False)

    tokens = lex_func(line)
    resulting_lines = []
    at_line_start = True
    indentation_len = len(level * indentation)
    current_line = ""
    padding_width = width - indentation_len
    for index, word in enumerate(tokens):
        has_next_word = index < len(tokens) - 1
        word_len = len(word)
        if not at_line_start:
            next_len = indentation_len + len(current_line) + 1 + word_len
            if next_len < width or (not has_next_word and next_len == width):
                # The word goes on the same line.
                current_line += " " + word
            else:
                # The word goes on the next line.
                resulting_lines.append(pad_func(current_line, padding_width))
                at_line_start = True
                current_line = indentation
        if at_line_start:
            current_line += word
            at_line_start = False
    resulting_lines.append(current_line)
    return resulting_lines


def count_leading_spaces(s):
    n = 0
    while n < len(s) and s[n] == " ":
        n += 1
    return n


def wrap_code(s, indent=4):
    lines = s.split("\n")
    result_lines = []
    for ln in lines:
        # Fypp directives must stay on one line.
        if ln.lstrip().startswith(("$:", "#:", "@:")):
            result_lines.append(ln)
            continue
        nspaces = count_leading_spaces(ln)
        level, remainder = divmod(nspaces, indent)

        if remainder != 0:
            raise ValueError(f"indentation of '{ln}' is not a multiple of " f"{indent}")

        result_lines.extend((level * indent) * " " + subln for subln in wrap_line_base(ln, level=level, indentation=" " * indent, pad_func=pad_fortran))

    return "\n".join(result_lines)


def float_to_fortran(num, kind):
    result = f"{num}"
    if "." not in result and "e" not in result.lower():
        result += ".0"
    result += f"_{kind}"
    if num < 0:
        result = "(%s)" % result
    return result


def str_np(ary, kind):
    return ", ".join(float_to_fortran(entry, kind) for entry in ary)


# }}}


def validate_mechanism(sol):
    """Reject unsupported physics before emitting a partially valid module."""
    if sol.thermo_model != "ideal-gas":
        raise ValueError(f"MFC thermochemistry requires ideal-gas thermodynamics, got {sol.thermo_model!r}")
    if sol.transport_model not in ("mixture-averaged", "multicomponent", "unity-Lewis-number"):
        raise ValueError(f"MFC thermochemistry requires gas transport data, got {sol.transport_model!r}")
    for species in sol.species():
        if not isinstance(species.thermo, ct.NasaPoly2):
            raise ValueError(f"Species {species.name}: MFC thermochemistry supports NASA7 polynomials only")
    supported = {"Arrhenius", "three-body-Arrhenius", "falloff-Troe", "falloff-Lindemann"}
    for i, reaction in enumerate(sol.reactions()):
        label = f"Reaction {i + 1} ({reaction.equation})"
        if reaction.reaction_type not in supported:
            raise ValueError(f"{label}: unsupported rate type {reaction.reaction_type!r}")
        if reaction.orders:
            raise ValueError(f"{label}: custom reaction orders are not supported")
        rates = (reaction.rate.low_rate, reaction.rate.high_rate) if reaction.reaction_type.startswith("falloff") else (reaction.rate,)
        if any(rate.pre_exponential_factor <= 0 for rate in rates):
            raise ValueError(f"{label}: Arrhenius pre-exponential factors must be positive")


def generate_fortran(solution, module_name="m_thermochem"):
    """Emit MFC's thermodynamic, kinetics and transport interface from Cantera as Fypp source.

    Precision (wp) and offload directives ($:GPU_ROUTINE) are resolved by MFC's build, so one
    source serves every configuration.
    """
    import re

    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,62}", module_name):
        raise ValueError(f"Invalid Fortran module name: {module_name!r}")
    validate_mechanism(solution)
    falloff = [(i, r) for i, r in enumerate(solution.reactions()) if r.reaction_type.startswith("falloff")]
    three_body = [(i, r) for i, r in enumerate(solution.reactions()) if r.reaction_type == "three-body-Arrhenius"]
    template = Template(filename=str(Path(__file__).with_name("module.fpp.mako")))
    return wrap_code(
        template.render(
            ct=ct,
            sol=solution,
            str_np=partial(str_np, kind="wp"),
            cgm=FortranExpressionMapper("wp"),
            Variable=p.Variable,
            float_to_fortran=partial(float_to_fortran, kind="wp"),
            species_name_length=max(map(len, solution.species_names)),
            module_name=module_name,
            ce=expressions,
            falloff_reactions=falloff,
            falloff_indices={i for i, _ in falloff},
            three_body_reactions=three_body,
        )
    )


# {{{ fortran expression generation


class FortranExpressionMapper(StringifyMapper):
    """Converts expressions to Fortran code."""

    def __init__(self, kind):
        super().__init__()
        self.kind = kind

    def map_constant(self, expr, enclosing_prec):
        if isinstance(expr, bool):
            if expr:
                return ".true."
            else:
                return ".false."
        else:
            return float_to_fortran(expr, self.kind)

    def map_variable(self, expr, enclosing_prec):
        return expr.name

    def map_lookup(self, expr, enclosing_prec):
        return self.parenthesize_if_needed(self.format("%s%%%s", self.rec(expr.aggregate, PREC_CALL), expr.name), enclosing_prec, PREC_CALL)

    def map_subscript(self, expr, enclosing_prec):
        def get_base_and_indices(expr):
            if not hasattr(expr, "aggregate") or not hasattr(expr, "index"):
                return expr, []

            # Get current level indices
            if isinstance(expr.index, tuple):
                current_indices = list(expr.index)
            else:
                current_indices = [expr.index]

            # Only recurse if aggregate is another subscript
            if hasattr(expr.aggregate, "aggregate") and hasattr(expr.aggregate, "index"):
                base, prev_indices = get_base_and_indices(expr.aggregate)
                return base, prev_indices + current_indices
            else:
                return expr.aggregate, current_indices

        # Get base array and all indices
        base_array, all_indices = get_base_and_indices(expr)

        # Convert zero-based expression indices without applying real-literal kinds.
        def convert_index(idx):
            if isinstance(idx, Integral):
                return str(idx + 1)
            return f"({self.rec(idx, PREC_NONE)} + 1)"

        # Format indices, converting floats to integers and adding 1
        index_str = ", ".join(convert_index(idx) for idx in all_indices)

        # Format the final expression
        return self.parenthesize_if_needed(self.format("%s(%s)", self.rec(base_array, PREC_CALL), index_str), enclosing_prec, PREC_CALL)

    def map_product(self, expr, enclosing_prec, *args, **kwargs):
        # This differs from the superclass only by adding spaces
        # around the operator, which provide an opportunity for
        # line breaking.
        return self.parenthesize_if_needed(self.join_rec(" * ", expr.children, PREC_PRODUCT, *args, **kwargs), enclosing_prec, PREC_PRODUCT)

    def map_logical_not(self, expr, enclosing_prec):
        from pymbolic.mapper.stringifier import PREC_UNARY

        return self.parenthesize_if_needed(".not. " + self.rec(expr.child, PREC_UNARY), enclosing_prec, PREC_UNARY)

    def map_logical_or(self, expr, enclosing_prec):
        from pymbolic.mapper.stringifier import PREC_LOGICAL_OR

        return self.parenthesize_if_needed(self.join_rec(" .or. ", expr.children, PREC_LOGICAL_OR), enclosing_prec, PREC_LOGICAL_OR)

    def map_logical_and(self, expr, enclosing_prec):
        from pymbolic.mapper.stringifier import PREC_LOGICAL_AND

        return self.parenthesize_if_needed(self.join_rec(" .and. ", expr.children, PREC_LOGICAL_AND), enclosing_prec, PREC_LOGICAL_AND)

    def map_if(self, expr, enclosing_prec):
        return self.format("merge(%s)" % self.join_rec(", ", [expr.then, expr.else_, expr.condition], PREC_NONE))


# }}}
