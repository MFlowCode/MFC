"""
Shared AST analyzer for case_validator.py.

Extracts all `self.prohibit(...)` rules from CaseValidator, determines
which parameter "triggers" each rule, and provides convenience functions
for both doc generators (parameters.md and case_constraints.md).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set
from collections import defaultdict


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Rule:
    method: str                      # e.g. "check_igr_simulation"
    lineno: int                      # line number of the prohibit call
    params: List[str]                # case parameter names used in condition
    message: str                     # user-friendly error message
    stages: Set[str] = field(default_factory=set)  # e.g. {"simulation", "pre_process"}
    trigger: Optional[str] = None    # param that "owns" this rule


# ---------------------------------------------------------------------------
# F-string message extraction
# ---------------------------------------------------------------------------

def _extract_message(node: ast.AST) -> Optional[str]:
    """
    Extract the message string from a prohibit() call's second argument.

    Handles both plain strings (ast.Constant) and f-strings (ast.JoinedStr).
    For f-strings, FormattedValue expressions are replaced with their
    unparsed source representation, giving a readable approximation.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value

    if isinstance(node, ast.JoinedStr):
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant):
                parts.append(str(value.value))
            elif isinstance(value, ast.FormattedValue):
                # Unparse the expression to get a readable approximation
                try:
                    parts.append(ast.unparse(value.value))
                except Exception:  # pylint: disable=broad-except
                    parts.append("?")
            else:
                parts.append("?")
        return "".join(parts)

    return None


# ---------------------------------------------------------------------------
# AST analysis: methods, call graph, rules
# ---------------------------------------------------------------------------

class CaseValidatorAnalyzer(ast.NodeVisitor):
    """
    Analyzes the CaseValidator class:

    - collects all methods
    - builds a call graph between methods
    - extracts all self.prohibit(...) rules
    """

    def __init__(self):
        super().__init__()
        self.in_case_validator = False
        self.current_method: str | None = None

        self.methods: Dict[str, ast.FunctionDef] = {}
        self.call_graph: Dict[str, Set[str]] = defaultdict(set)
        self.rules: List[Rule] = []

        # Stack of {local_name -> param_name} maps, one per method
        self.local_param_stack: List[Dict[str, str]] = []

        # {method_name -> trigger_param} from guard detection
        self._method_guards: Dict[str, str] = {}

    # --- top-level entrypoint ---

    def visit_ClassDef(self, node: ast.ClassDef):
        if node.name == "CaseValidator":
            self.in_case_validator = True
            # collect methods
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    self.methods[item.name] = item
            # now analyze all methods
            for method in self.methods.values():
                self._analyze_method(method)
            self.in_case_validator = False
        else:
            self.generic_visit(node)

    # --- per-method analysis ---

    def _analyze_method(self, func: ast.FunctionDef):
        """Analyze a single method: local param mapping, call graph, rules."""
        self.current_method = func.name
        local_param_map = self._build_local_param_map(func)
        self.local_param_stack.append(local_param_map)

        # Detect method guard pattern
        guard = _extract_method_guard(func, local_param_map)
        if guard:
            self._method_guards[func.name] = guard

        self.generic_visit(func)
        self.local_param_stack.pop()
        self.current_method = None

    def _build_local_param_map(self, func: ast.FunctionDef) -> Dict[str, str]:  # pylint: disable=too-many-nested-blocks
        """
        Look for assignments like:
            igr = self.get('igr', 'F') == 'T'
            model_eqns = self.get('model_eqns')
        and record local_name -> 'param_name'.
        """
        m: Dict[str, str] = {}
        for stmt in func.body:  # pylint: disable=too-many-nested-blocks
            if isinstance(stmt, ast.Assign):
                # Handle both direct calls and comparisons
                value = stmt.value
                # Unwrap comparisons like "self.get('igr', 'F') == 'T'"
                if isinstance(value, ast.Compare):
                    value = value.left

                if isinstance(value, ast.Call):
                    call = value
                    if (  # pylint: disable=too-many-boolean-expressions
                        isinstance(call.func, ast.Attribute)
                        and isinstance(call.func.value, ast.Name)
                        and call.func.value.id == "self"
                        and call.func.attr == "get"
                        and call.args
                        and isinstance(call.args[0], ast.Constant)
                        and isinstance(call.args[0].value, str)
                    ):
                        param_name = call.args[0].value
                        for target in stmt.targets:
                            if isinstance(target, ast.Name):
                                m[target.id] = param_name
        return m

    # --- visit calls to build call graph + rules ---

    def visit_Call(self, node: ast.Call):
        # record method call edges: self.some_method(...)
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "self"
            and isinstance(node.func.attr, str)
        ):
            callee = node.func.attr
            if self.current_method is not None:
                # method call on self
                self.call_graph[self.current_method].add(callee)

        # detect self.prohibit(<condition>, "<message>")
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "self"
            and node.func.attr == "prohibit"
            and len(node.args) >= 2
        ):
            condition, msg_node = node.args[0], node.args[1]
            msg = _extract_message(msg_node)
            if msg is not None:
                local_map = self.local_param_stack[-1] if self.local_param_stack else {}
                params = sorted(self._extract_params(condition))
                trigger = self._determine_trigger(params, condition, local_map)
                rule = Rule(
                    method=self.current_method or "<unknown>",
                    lineno=node.lineno,
                    params=params,
                    message=msg,
                    trigger=trigger,
                )
                self.rules.append(rule)

        self.generic_visit(node)

    def _determine_trigger(self, _params: List[str], condition: ast.AST,
                           local_map: Dict[str, str]) -> Optional[str]:
        """Determine trigger param: method guard first, then condition fallback."""
        # 1. Method guard (high confidence)
        if self.current_method and self.current_method in self._method_guards:
            return self._method_guards[self.current_method]

        # 2. Condition first-param fallback
        return _extract_trigger_from_condition(condition, local_map)

    def _extract_params(self, condition: ast.AST) -> Set[str]:
        """
        Collect parameter names used in the condition via:
          - local variables mapped from self.get(...)
          - direct self.get('param_name', ...) calls
        """
        params: Set[str] = set()
        local_map = self.local_param_stack[-1] if self.local_param_stack else {}

        for node in ast.walk(condition):
            # local names
            if isinstance(node, ast.Name) and node.id in local_map:
                params.add(local_map[node.id])

            # direct self.get('param_name')
            if isinstance(node, ast.Call):
                if (  # pylint: disable=too-many-boolean-expressions
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "self"
                    and node.func.attr == "get"
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)
                ):
                    params.add(node.args[0].value)

        return params


# ---------------------------------------------------------------------------
# Trigger detection helpers
# ---------------------------------------------------------------------------

def _extract_method_guard(func: ast.FunctionDef, local_param_map: Dict[str, str]) -> Optional[str]:
    """
    Detect early-return guard patterns like:
        if not bubbles_euler:
            return
    The guarded variable's param is the trigger for all rules in that method.
    """
    for stmt in func.body:  # pylint: disable=too-many-nested-blocks
        if not isinstance(stmt, ast.If):
            continue

        # Check for "if not <var>: return" pattern
        test = stmt.test
        if isinstance(test, ast.UnaryOp) and isinstance(test.op, ast.Not):
            if isinstance(test.operand, ast.Name):
                var_name = test.operand.id
                # Check body is just "return"
                if (len(stmt.body) == 1
                        and isinstance(stmt.body[0], ast.Return)
                        and stmt.body[0].value is None):
                    if var_name in local_param_map:
                        return local_param_map[var_name]

        # Check for "if <var> != <value>: return" pattern
        # e.g. "if recon_type != 1: return" or "if model_eqns != 3: return"
        if isinstance(test, ast.Compare) and len(test.ops) == 1:
            if isinstance(test.ops[0], ast.NotEq):
                if isinstance(test.left, ast.Name):
                    var_name = test.left.id
                    if (len(stmt.body) == 1
                            and isinstance(stmt.body[0], ast.Return)
                            and stmt.body[0].value is None):
                        if var_name in local_param_map:
                            return local_param_map[var_name]

    return None


def _extract_trigger_from_condition(condition: ast.AST, local_param_map: Dict[str, str]) -> Optional[str]:
    """
    Fallback trigger detection: walk the condition AST left-to-right,
    return the first parameter name found.
    """
    # Walk in source order (left-to-right in boolean expressions)
    for node in ast.walk(condition):
        if isinstance(node, ast.Name) and node.id in local_param_map:
            return local_param_map[node.id]
        if isinstance(node, ast.Call):
            if (  # pylint: disable=too-many-boolean-expressions
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "self"
                and node.func.attr == "get"
                and node.args
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)
            ):
                return node.args[0].value
    return None


# ---------------------------------------------------------------------------
# Stage inference from validate_* roots and call graph
# ---------------------------------------------------------------------------

STAGE_ROOTS: Dict[str, List[str]] = {
    "common": ["validate_common"],
    "simulation": ["validate_simulation"],
    "pre_process": ["validate_pre_process"],
    "post_process": ["validate_post_process"],
}


def compute_method_stages(call_graph: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    """
    For each stage (simulation/pre_process/post_process/common), starting from
    validate_* roots, walk the call graph and record which methods belong to which stages.
    """
    method_stages: Dict[str, Set[str]] = defaultdict(set)

    def dfs(start: str, stage: str):
        stack = [start]
        visited: Set[str] = set()
        while stack:
            m = stack.pop()
            if m in visited:
                continue
            visited.add(m)
            method_stages[m].add(stage)
            for nxt in call_graph.get(m, ()):
                if nxt not in visited:
                    stack.append(nxt)

    for stage, roots in STAGE_ROOTS.items():
        for root in roots:
            dfs(root, stage)

    return method_stages


# ---------------------------------------------------------------------------
# Classification of messages for nicer grouping
# ---------------------------------------------------------------------------

def classify_message(msg: str) -> str:
    """
    Roughly classify rule messages for nicer grouping.

    Returns one of: "requirement", "incompatibility", "range", "other".
    """
    text = msg.lower()

    if (  # pylint: disable=too-many-boolean-expressions
        "not compatible" in text
        or "does not support" in text
        or "cannot be used" in text
        or "must not" in text
        or "is not supported" in text
        or "incompatible" in text
        or "untested" in text
        or "not available" in text
        or "is not compatible" in text
        or "activate only one" in text
    ):
        return "incompatibility"

    if (  # pylint: disable=too-many-boolean-expressions
        "requires" in text
        or "must be set if" in text
        or "must be specified" in text
        or "must be set with" in text
        or "can only be enabled if" in text
        or "must be set for" in text
        or "must be set when" in text
    ):
        return "requirement"

    if (  # pylint: disable=too-many-boolean-expressions
        "must be between" in text
        or "must be positive" in text
        or "must be non-negative" in text
        or "must be greater than" in text
        or "must be less than" in text
        or "must be at least" in text
        or "must be <=" in text
        or "must be >=" in text
        or "must be odd" in text
        or "divisible by" in text
        or "must be 1" in text
        or "must be 'T' or 'F'" in text
    ):
        return "range"

    return "other"


# Optional: nicer display names / categories (you can extend this)
FEATURE_META = {
    "igr": {"title": "Iterative Generalized Riemann (IGR)", "category": "solver"},
    "bubbles_euler": {"title": "Euler-Euler Bubble Model", "category": "bubbles"},
    "bubbles_lagrange": {"title": "Euler-Lagrange Bubble Model", "category": "bubbles"},
    "qbmm": {"title": "Quadrature-Based Moment Method (QBMM)", "category": "bubbles"},
    "polydisperse": {"title": "Polydisperse Bubble Dynamics", "category": "bubbles"},
    "mhd": {"title": "Magnetohydrodynamics (MHD)", "category": "physics"},
    "alt_soundspeed": {"title": "Alternative Sound Speed", "category": "physics"},
    "surface_tension": {"title": "Surface Tension Model", "category": "physics"},
    "hypoelasticity": {"title": "Hypoelasticity", "category": "physics"},
    "hyperelasticity": {"title": "Hyperelasticity", "category": "physics"},
    "relax": {"title": "Phase Change (Relaxation)", "category": "physics"},
    "viscous": {"title": "Viscosity", "category": "physics"},
    "acoustic_source": {"title": "Acoustic Sources", "category": "physics"},
    "ib": {"title": "Immersed Boundaries", "category": "geometry"},
    "cyl_coord": {"title": "Cylindrical Coordinates", "category": "geometry"},
    "weno_order": {"title": "WENO Order", "category": "numerics"},
    "muscl_order": {"title": "MUSCL Order", "category": "numerics"},
    "riemann_solver": {"title": "Riemann Solver", "category": "numerics"},
    "model_eqns": {"title": "Model Equations", "category": "fundamentals"},
    "num_fluids": {"title": "Number of Fluids", "category": "fundamentals"},
}


def feature_title(param: str) -> str:
    meta = FEATURE_META.get(param)
    if meta and "title" in meta:
        return meta["title"]
    return param


# ---------------------------------------------------------------------------
# Convenience: full analysis pipeline
# ---------------------------------------------------------------------------

_DEFAULT_VALIDATOR_PATH = Path(__file__).resolve().parent.parent / "case_validator.py"


def analyze_case_validator(path: Optional[Path] = None) -> Dict:
    """
    Parse case_validator.py and return extracted rules indexed multiple ways.

    Returns a dict with:
        rules: List[Rule]          - all rules
        by_trigger: Dict[str, List[Rule]]  - rules indexed by trigger param
        by_param: Dict[str, List[Rule]]    - rules indexed by all mentioned params
    """
    if path is None:
        path = _DEFAULT_VALIDATOR_PATH

    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(path))

    analyzer = CaseValidatorAnalyzer()
    analyzer.visit(tree)

    # Infer stages per method from call graph
    method_stages = compute_method_stages(analyzer.call_graph)

    # Attach stages to rules
    for r in analyzer.rules:
        r.stages = method_stages.get(r.method, set())

    # Build indices
    by_trigger: Dict[str, List[Rule]] = defaultdict(list)
    by_param: Dict[str, List[Rule]] = defaultdict(list)

    for r in analyzer.rules:
        if r.trigger:
            by_trigger[r.trigger].append(r)
        for p in r.params:
            by_param[p].append(r)

    return {
        "rules": analyzer.rules,
        "by_trigger": dict(by_trigger),
        "by_param": dict(by_param),
        "call_graph": analyzer.call_graph,
        "methods": analyzer.methods,
    }
