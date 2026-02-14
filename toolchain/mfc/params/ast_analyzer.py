"""
Shared AST analyzer for case_validator.py.

Extracts all `self.prohibit(...)` and `self.warn(...)` rules from
CaseValidator, determines which parameter "triggers" each rule, and
provides convenience functions for both doc generators (parameters.md
and case_constraints.md).
"""

import ast
import re
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
    lineno: int                      # line number of the prohibit/warn call
    params: List[str]                # case parameter names used in condition
    message: str                     # user-friendly error/warning message
    stages: Set[str] = field(default_factory=set)  # e.g. {"simulation", "pre_process"}
    trigger: Optional[str] = None    # param that "owns" this rule
    severity: str = "error"          # "error" (prohibit) or "warning" (warn)


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


def _resolve_fstring(node: ast.JoinedStr, subs: Dict[str, str]) -> Optional[str]:
    """Resolve a JoinedStr (f-string) by substituting known loop variables."""
    parts: List[str] = []
    for v in node.values:
        if isinstance(v, ast.Constant):
            parts.append(str(v.value))
        elif isinstance(v, ast.FormattedValue):
            if isinstance(v.value, ast.Name) and v.value.id in subs:
                parts.append(subs[v.value.id])
            else:
                try:
                    parts.append(ast.unparse(v.value))
                except Exception:  # pylint: disable=broad-except
                    parts.append("?")
        else:
            parts.append("?")
    return "".join(parts)


def _resolve_message(msg_node: ast.AST, subs: Dict[str, str]) -> Optional[str]:
    """Resolve a prohibit/warn message, substituting loop variables in f-strings."""
    if isinstance(msg_node, ast.Constant) and isinstance(msg_node.value, str):
        return msg_node.value
    if isinstance(msg_node, ast.JoinedStr):
        return _resolve_fstring(msg_node, subs)
    return None


def _is_self_get(call: ast.Call) -> bool:
    """Check if a Call node is self.get(...)."""
    return (isinstance(call.func, ast.Attribute)
            and isinstance(call.func.value, ast.Name)
            and call.func.value.id == "self"
            and call.func.attr == "get"
            and bool(call.args))


# ---------------------------------------------------------------------------
# AST analysis: methods, call graph, rules
# ---------------------------------------------------------------------------

class CaseValidatorAnalyzer(ast.NodeVisitor):  # pylint: disable=too-many-instance-attributes
    """
    Analyzes the CaseValidator class:

    - collects all methods
    - builds a call graph between methods
    - extracts all self.prohibit(...) and self.warn(...) rules
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

        # Stack of {alias_name -> [source_param, ...]} maps (parallel to local_param_stack)
        self.alias_map_stack: List[Dict[str, List[str]]] = []

        # {method_name -> trigger_param} from guard detection
        self._method_guards: Dict[str, str] = {}

        # Line numbers of prohibit/warn calls handled by loop expansion (skip in visit_Call)
        self._expanded_prohibit_lines: Set[int] = set()

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
        alias_map = self._build_alias_map(func, local_param_map)
        self.local_param_stack.append(local_param_map)
        self.alias_map_stack.append(alias_map)

        # Detect method guard pattern
        guard = _extract_method_guard(func, local_param_map)
        if guard:
            self._method_guards[func.name] = guard

        # Expand literal-list for-loops before generic_visit
        self._expand_literal_loops(func, local_param_map)

        self.generic_visit(func)

        # Enrich rules with params from if-guard conditions
        self._enrich_rules_with_if_guards(func, local_param_map, alias_map)

        self.alias_map_stack.pop()
        self.local_param_stack.pop()
        self.current_method = None

    def _enrich_rules_with_if_guards(self, func: ast.FunctionDef,
                                     local_param_map: Dict[str, str],
                                     alias_map: Dict[str, List[str]]):
        """
        After rules are extracted, walk the function body for ast.If nodes.
        For each if-block, extract guard params from the test condition and add
        them to every rule whose lineno falls within the block's line range.
        """
        for node in ast.walk(func):  # pylint: disable=too-many-nested-blocks
            if not isinstance(node, ast.If):
                continue
            # Extract params from the if-test condition
            guard_params = _extract_test_params(node.test, local_param_map, alias_map)
            if not guard_params:
                continue
            # Determine line ranges for body and orelse
            ranges = []
            if node.body:
                body_start = node.body[0].lineno
                body_end = node.body[-1].end_lineno or node.body[-1].lineno
                ranges.append((body_start, body_end))
            if node.orelse:
                else_start = node.orelse[0].lineno
                else_end = node.orelse[-1].end_lineno or node.orelse[-1].lineno
                ranges.append((else_start, else_end))
            # Enrich matching rules
            for rule in self.rules:
                if rule.method != func.name:
                    continue
                for rng_start, rng_end in ranges:
                    if rng_start <= rule.lineno <= rng_end:
                        for gp in guard_params:
                            if gp not in rule.params:
                                rule.params.append(gp)
                        break

    def _build_local_param_map(self, func: ast.FunctionDef) -> Dict[str, str]:  # pylint: disable=too-many-nested-blocks
        """
        Look for assignments like:
            igr = self.get('igr', 'F') == 'T'
            model_eqns = self.get('model_eqns')
        and record local_name -> 'param_name'.

        Uses ast.walk to find assignments at any nesting depth (inside if/for/with blocks).
        """
        m: Dict[str, str] = {}
        for node in ast.walk(func):  # pylint: disable=too-many-nested-blocks
            if isinstance(node, ast.Assign):
                # Handle both direct calls and comparisons
                value = node.value
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
                        for target in node.targets:
                            if isinstance(target, ast.Name):
                                m[target.id] = param_name
        return m

    @staticmethod
    def _build_alias_map(func: ast.FunctionDef,
                         local_param_map: Dict[str, str]) -> Dict[str, List[str]]:
        """
        Detect boolean alias assignments like:
            variable_dt = cfl_dt or cfl_adap_dt       # BoolOp(Or)
            has_output = rho_wrt or E_wrt or ...       # BoolOp(Or) with many operands
            skip_check = cyl_coord and dir in [...]    # BoolOp(And)

        Returns {alias_name -> [source_param, ...]}.
        """
        alias_map: Dict[str, List[str]] = {}
        for node in ast.walk(func):
            if not isinstance(node, ast.Assign):
                continue
            if not isinstance(node.value, ast.BoolOp):
                continue
            # Collect source params from BoolOp operands
            sources: List[str] = []
            for operand in node.value.values:
                if isinstance(operand, ast.Name) and operand.id in local_param_map:
                    sources.append(local_param_map[operand.id])
            if not sources:
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    alias_map[target.id] = sources
        return alias_map

    # --- literal-list for-loop expansion ---

    def _expand_literal_loops(self, func: ast.FunctionDef, local_param_map: Dict[str, str]):
        """Expand `for var in [x, y, z]:` loops into concrete Rules."""
        self._expand_loop_stmts(func.body, func.name, local_param_map, {})

    def _expand_loop_stmts(self, stmts: list, method_name: str,
                            parent_map: Dict[str, str], subs: Dict[str, str]):
        """Recursively find literal-list for-loops and create expanded Rules."""
        for stmt in stmts:
            if (isinstance(stmt, ast.For)
                    and isinstance(stmt.target, ast.Name)
                    and isinstance(stmt.iter, ast.List)
                    and all(isinstance(e, ast.Constant) for e in stmt.iter.elts)):
                var = stmt.target.id
                for elt in stmt.iter.elts:
                    new_subs = {**subs, var: str(elt.value)}
                    loop_map = self._resolve_loop_gets(stmt.body, new_subs)
                    merged = {**parent_map, **loop_map}
                    # Detect loop-body guard: `if not <var>: continue`
                    loop_guard = self._detect_loop_guard(stmt.body, merged)
                    # Recurse for nested literal-list loops
                    self._expand_loop_stmts(stmt.body, method_name, merged, new_subs)
                    # Create Rules for prohibit calls at this level
                    self._create_loop_rules(stmt.body, method_name, merged, new_subs, loop_guard)
            elif subs:
                # Inside an expanded loop: recurse into if/else blocks
                if isinstance(stmt, ast.If):
                    self._expand_loop_stmts(stmt.body, method_name, parent_map, subs)
                    if stmt.orelse:
                        self._expand_loop_stmts(stmt.orelse, method_name, parent_map, subs)

    @staticmethod
    def _detect_loop_guard(stmts: list, local_map: Dict[str, str]) -> Optional[str]:
        """Detect `if not <var>: continue` guard pattern in a loop body."""
        for stmt in stmts:
            if not isinstance(stmt, ast.If):
                continue
            test = stmt.test
            if (isinstance(test, ast.UnaryOp) and isinstance(test.op, ast.Not)
                    and isinstance(test.operand, ast.Name)):
                if (len(stmt.body) == 1 and isinstance(stmt.body[0], ast.Continue)):
                    var_name = test.operand.id
                    if var_name in local_map:
                        return local_map[var_name]
        return None

    @staticmethod
    def _resolve_loop_gets(stmts: list, subs: Dict[str, str]) -> Dict[str, str]:
        """Resolve self.get() assignments in loop body, substituting f-string loop vars."""
        m: Dict[str, str] = {}
        for stmt in stmts:
            for node in ast.walk(stmt):
                if not isinstance(node, ast.Assign):
                    continue
                value = node.value
                if isinstance(value, ast.Compare):
                    value = value.left
                if not isinstance(value, ast.Call) or not _is_self_get(value):
                    continue
                arg = value.args[0]
                if isinstance(arg, ast.JoinedStr):
                    resolved = _resolve_fstring(arg, subs)
                    if resolved is None:
                        continue
                    param_name = resolved
                elif isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    param_name = arg.value
                else:
                    continue
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        m[target.id] = param_name
        return m

    def _create_loop_rules(self, stmts: list, method_name: str,  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
                            local_map: Dict[str, str], subs: Dict[str, str],
                            loop_guard: Optional[str] = None):
        """Create Rules for self.prohibit()/self.warn() calls found in loop body statements."""
        for stmt in stmts:
            # Skip nested literal-list for-loops (handled by recursion)
            if (isinstance(stmt, ast.For)
                    and isinstance(stmt.target, ast.Name)
                    and isinstance(stmt.iter, ast.List)
                    and all(isinstance(e, ast.Constant) for e in stmt.iter.elts)):
                continue
            for node in ast.walk(stmt):
                if not isinstance(node, ast.Call):
                    continue
                if not (isinstance(node.func, ast.Attribute)
                        and isinstance(node.func.value, ast.Name)
                        and node.func.value.id == "self"
                        and node.func.attr in ("prohibit", "warn")
                        and len(node.args) >= 2):
                    continue
                severity = "warning" if node.func.attr == "warn" else "error"
                condition, msg_node = node.args[0], node.args[1]
                msg = _resolve_message(msg_node, subs)
                if msg is None:
                    continue
                param_set = self._extract_params_with_subs(condition, local_map, subs)
                # Use loop-body guard as trigger and include it in params
                if loop_guard:
                    trigger = loop_guard
                    param_set.add(loop_guard)
                else:
                    trigger = self._determine_trigger(
                        sorted(param_set), condition, local_map)
                params = sorted(param_set)
                rule = Rule(
                    method=method_name,
                    lineno=node.lineno,
                    params=params,
                    message=msg,
                    trigger=trigger,
                    severity=severity,
                )
                self.rules.append(rule)
                self._expanded_prohibit_lines.add(node.lineno)

    def _extract_params_with_subs(self, condition: ast.AST,
                                   local_map: Dict[str, str],
                                   subs: Dict[str, str]) -> Set[str]:
        """Like _extract_params but also resolves JoinedStr self.get() args."""
        params: Set[str] = set()
        for node in ast.walk(condition):
            if isinstance(node, ast.Name) and node.id in local_map:
                params.add(local_map[node.id])
            if isinstance(node, ast.Call) and _is_self_get(node):
                arg = node.args[0]
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    params.add(arg.value)
                elif isinstance(arg, ast.JoinedStr):
                    resolved = _resolve_fstring(arg, subs)
                    if resolved:
                        params.add(resolved)
        return params

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

        # detect self.prohibit(<condition>, "<message>") and self.warn(<condition>, "<message>")
        # Skip calls already handled by loop expansion
        if (  # pylint: disable=too-many-boolean-expressions
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "self"
            and node.func.attr in ("prohibit", "warn")
            and len(node.args) >= 2
            and node.lineno not in self._expanded_prohibit_lines
        ):
            severity = "warning" if node.func.attr == "warn" else "error"
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
                    severity=severity,
                )
                self.rules.append(rule)

        self.generic_visit(node)

    def _determine_trigger(self, _params: List[str], condition: ast.AST,
                           local_map: Dict[str, str]) -> Optional[str]:
        """Determine trigger param: method guard first, then condition fallback."""
        # 1. Method guard (high confidence)
        if self.current_method and self.current_method in self._method_guards:
            return self._method_guards[self.current_method]

        # 2. Condition first-param fallback (with alias resolution)
        alias_map = self.alias_map_stack[-1] if self.alias_map_stack else {}
        return _extract_trigger_from_condition(condition, local_map, alias_map)

    def _extract_params(self, condition: ast.AST) -> Set[str]:
        """
        Collect parameter names used in the condition via:
          - local variables mapped from self.get(...)
          - boolean aliases (variable_dt → [cfl_dt, cfl_adap_dt])
          - direct self.get('param_name', ...) calls
        """
        params: Set[str] = set()
        local_map = self.local_param_stack[-1] if self.local_param_stack else {}
        alias_map = self.alias_map_stack[-1] if self.alias_map_stack else {}

        for node in ast.walk(condition):
            # local names
            if isinstance(node, ast.Name):
                if node.id in local_map:
                    params.add(local_map[node.id])
                elif node.id in alias_map:
                    params.update(alias_map[node.id])

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


def _extract_test_params(test: ast.AST, local_param_map: Dict[str, str],
                         alias_map: Dict[str, List[str]]) -> Set[str]:
    """Extract parameter names from an if-test condition, resolving aliases."""
    params: Set[str] = set()
    for node in ast.walk(test):
        if isinstance(node, ast.Name):
            if node.id in local_param_map:
                params.add(local_param_map[node.id])
            elif node.id in alias_map:
                params.update(alias_map[node.id])
        if isinstance(node, ast.Call) and _is_self_get(node):
            arg = node.args[0]
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                params.add(arg.value)
    return params


def _extract_trigger_from_condition(condition: ast.AST, local_param_map: Dict[str, str],
                                    alias_map: Optional[Dict[str, List[str]]] = None) -> Optional[str]:
    """
    Fallback trigger detection: walk the condition AST left-to-right,
    return the first parameter name found. Resolves aliases to their first source param.
    """
    if alias_map is None:
        alias_map = {}
    # Walk in source order (left-to-right in boolean expressions)
    for node in ast.walk(condition):
        if isinstance(node, ast.Name):
            if node.id in local_param_map:
                return local_param_map[node.id]
            if node.id in alias_map and alias_map[node.id]:
                return alias_map[node.id][0]
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
        or re.search(r"must be 1\b", text) is not None
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
