"""AST-based translation of analytic case expressions (Python syntax) to Fortran.

Replaces the legacy regex token substitution: expressions are parsed, every
variable name is resolved through an explicit map, and anything unrecognized
is a load-time error instead of silently-wrong generated Fortran.
"""

import ast
from typing import Dict

# Fortran intrinsics permitted in analytic expressions; emitted unchanged.
INTRINSICS = {
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "atan2",
    "sinh",
    "cosh",
    "tanh",
    "exp",
    "log",
    "log10",
    "sqrt",
    "abs",
    "min",
    "max",
    "mod",
    "sign",
}

# Names with no mapping that resolve on the Fortran side (m_constants).
PASSTHROUGH = {"pi"}

_BINOPS = {ast.Add: "+", ast.Sub: "-", ast.Mult: "*", ast.Div: "/", ast.Pow: "**"}
# Precedence: additive < multiplicative < unary minus < power.
_PREC = {ast.Add: 1, ast.Sub: 1, ast.Mult: 2, ast.Div: 2, ast.Pow: 4}
_UNARY_PREC = 3


class AnalyticExprError(Exception):
    """Raised when an analytic expression cannot be translated to Fortran."""


def fortranize_expr(expr: str, var_map: Dict[str, str]) -> str:
    """Translate a Python-syntax analytic expression into a Fortran expression.

    var_map maps case-file variable names (x, lx, r, ...) to the Fortran
    entities they stand for. Intrinsics and PASSTHROUGH names are kept as-is.
    """
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as e:
        raise AnalyticExprError(f"invalid analytic expression {expr!r}: {e.msg}") from e
    return _emit(tree.body, expr, var_map, 0)


def _emit(node, expr: str, var_map: Dict[str, str], parent_prec: int) -> str:
    if isinstance(node, ast.BinOp):
        op = type(node.op)
        if op not in _BINOPS:
            raise AnalyticExprError(f"unsupported operator in {expr!r}")
        prec = _PREC[op]
        if op is ast.Pow:  # right-associative
            left = _emit(node.left, expr, var_map, prec + 1)
            right = _emit(node.right, expr, var_map, prec)
        else:
            left = _emit(node.left, expr, var_map, prec)
            right = _emit(node.right, expr, var_map, prec + 1)
        s = f"{left} {_BINOPS[op]} {right}"
        return f"({s})" if prec < parent_prec else s
    if isinstance(node, ast.UnaryOp):
        if not isinstance(node.op, (ast.USub, ast.UAdd)):
            raise AnalyticExprError(f"unsupported operator in {expr!r}")
        sign = "-" if isinstance(node.op, ast.USub) else "+"
        s = f"{sign}{_emit(node.operand, expr, var_map, _UNARY_PREC)}"
        # Fortran forbids an operator directly following another operator
        # (`y * -x` is invalid), so parenthesize inside any binary context.
        return f"({s})" if parent_prec > 1 else s
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in INTRINSICS or node.keywords:
            name = node.func.id if isinstance(node.func, ast.Name) else ast.unparse(node.func)
            raise AnalyticExprError(f"unknown function '{name}' in {expr!r}. Allowed: {', '.join(sorted(INTRINSICS))}")
        args = ", ".join(_emit(a, expr, var_map, 0) for a in node.args)
        return f"{node.func.id}({args})"
    if isinstance(node, ast.Name):
        if node.id in var_map:
            return var_map[node.id]
        if node.id in PASSTHROUGH:
            return node.id
        valid = ", ".join(sorted(var_map) + sorted(PASSTHROUGH))
        raise AnalyticExprError(f"unknown variable '{node.id}' in {expr!r}. Available: {valid} (plus intrinsics)")
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            s = repr(node.value)
            if isinstance(node.value, float) and ("e" in s or "E" in s):
                s += "_wp"  # exponent form defaults to single precision in Fortran
            return s
        raise AnalyticExprError(f"unsupported literal {node.value!r} in {expr!r}")
    raise AnalyticExprError(f"unsupported syntax in {expr!r}: {type(node).__name__}")
