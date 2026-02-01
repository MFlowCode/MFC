"""
Constraint Coverage Analysis Tool.

Analyzes which validation constraints are exercised by the test cases.
This helps identify gaps in test coverage before refactoring.
"""

import ast
import json
from pathlib import Path
from typing import Dict, List, Any
from dataclasses import dataclass


@dataclass
class ConstraintInfo:
    """
    Information about a single constraint.

    Attributes:
        method: Name of the check_* method containing this constraint.
        line_number: Line number of the prohibit() call start (1-indexed).
            For multi-line calls, this is the first line.
        message: Error message shown when constraint is violated.
        condition_code: Unparsed source code of the condition expression.
    """
    method: str
    line_number: int
    message: str
    condition_code: str


def _extract_message(msg_node: ast.expr) -> str:
    """Extract message string from AST node."""
    if isinstance(msg_node, ast.Constant):
        return msg_node.value
    if isinstance(msg_node, ast.JoinedStr):
        # f-string - extract the static parts
        return "".join(
            p.value if isinstance(p, ast.Constant) else "{...}"
            for p in msg_node.values
        )
    return "<dynamic>"


def _is_prohibit_call(node: ast.AST) -> bool:
    """Check if node is a self.prohibit() call with enough arguments."""
    if not isinstance(node, ast.Call):
        return False
    if not isinstance(node.func, ast.Attribute):
        return False
    return node.func.attr == 'prohibit' and len(node.args) >= 2


def _find_case_validator_class(tree: ast.Module) -> ast.ClassDef:
    """Find the CaseValidator class in the AST."""
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "CaseValidator":
            return node
    return None


def extract_constraints_from_validator() -> List[ConstraintInfo]:
    """Parse case_validator.py and extract all prohibit() calls."""
    validator_path = Path(__file__).parent.parent / "case_validator.py"

    with open(validator_path, 'r', encoding='utf-8') as f:
        source = f.read()

    tree = ast.parse(source)
    validator_class = _find_case_validator_class(tree)
    if validator_class is None:
        return []

    constraints: List[ConstraintInfo] = []

    # Iterate through each method in the class
    for item in validator_class.body:
        if not isinstance(item, ast.FunctionDef):
            continue

        # Walk only within this method to find prohibit() calls
        for node in ast.walk(item):
            if not _is_prohibit_call(node):
                continue

            message = _extract_message(node.args[1])
            try:
                condition_code = ast.unparse(node.args[0])
            except (ValueError, TypeError, AttributeError):
                # ast.unparse can fail on malformed AST or missing attributes
                condition_code = "<complex>"

            # Note: node.lineno points to the start of the prohibit() call.
            # For multi-line calls, this is the first line, not where the
            # condition or message appears.
            constraints.append(ConstraintInfo(
                method=item.name,
                line_number=node.lineno,
                message=message,
                condition_code=condition_code
            ))

    return constraints


def _count_prohibit_calls(func_node: ast.FunctionDef) -> int:
    """Count self.prohibit() calls in a function."""
    count = 0
    for subnode in ast.walk(func_node):
        if _is_prohibit_call(subnode):
            count += 1
    return count


def extract_check_methods() -> Dict[str, Dict[str, Any]]:
    """Extract all check_* methods from validator with their stage."""
    validator_path = Path(__file__).parent.parent / "case_validator.py"

    with open(validator_path, 'r', encoding='utf-8') as f:
        source = f.read()

    methods = {}
    tree = ast.parse(source)
    validator_class = _find_case_validator_class(tree)

    if validator_class is None:
        return methods

    for item in validator_class.body:
        if not isinstance(item, ast.FunctionDef):
            continue
        if not item.name.startswith("check_"):
            continue

        docstring = ast.get_docstring(item) or ""
        methods[item.name] = {
            "line_number": item.lineno,
            "docstring": docstring.split('\n')[0] if docstring else "",
            "prohibit_count": _count_prohibit_calls(item),
        }

    return methods


_VALIDATE_METHOD_TO_STAGE = {
    "validate_common": "common",
    "validate_pre_process": "pre_process",
    "validate_simulation": "simulation",
    "validate_post_process": "post_process",
}


def _find_check_calls(func_node: ast.FunctionDef) -> List[str]:
    """Find all self.check_* method calls in a function."""
    calls = []
    for subnode in ast.walk(func_node):
        if not isinstance(subnode, ast.Call):
            continue
        if not isinstance(subnode.func, ast.Attribute):
            continue
        if subnode.func.attr.startswith("check_"):
            calls.append(subnode.func.attr)
    return calls


def extract_validate_dispatch() -> Dict[str, List[str]]:
    """Extract which check methods are called for each stage."""
    validator_path = Path(__file__).parent.parent / "case_validator.py"

    with open(validator_path, 'r', encoding='utf-8') as f:
        source = f.read()

    dispatch = {stage: [] for stage in _VALIDATE_METHOD_TO_STAGE.values()}
    tree = ast.parse(source)
    validator_class = _find_case_validator_class(tree)

    if validator_class is None:
        return dispatch

    for item in validator_class.body:
        if not isinstance(item, ast.FunctionDef):
            continue
        stage = _VALIDATE_METHOD_TO_STAGE.get(item.name)
        if stage is None:
            continue
        dispatch[stage].extend(_find_check_calls(item))

    return dispatch


def generate_coverage_report() -> Dict[str, Any]:
    """Generate a comprehensive coverage report."""
    constraints = extract_constraints_from_validator()
    methods = extract_check_methods()
    dispatch = extract_validate_dispatch()

    # Group constraints by method
    by_method = {}
    for c in constraints:
        if c.method not in by_method:
            by_method[c.method] = []
        by_method[c.method].append({
            "line": c.line_number,
            "message": c.message,
            "condition": c.condition_code[:80] + "..." if len(c.condition_code) > 80 else c.condition_code
        })

    # Calculate coverage per stage
    stage_coverage = {}
    for stage, check_methods in dispatch.items():
        total_constraints = 0
        for method_name in check_methods:
            if method_name in methods:
                total_constraints += methods[method_name]["prohibit_count"]
        stage_coverage[stage] = {
            "methods": check_methods,
            "method_count": len(check_methods),
            "constraint_count": total_constraints,
        }

    # Add common constraints to all stages
    common_constraints = stage_coverage.get("common", {}).get("constraint_count", 0)
    for stage in ["pre_process", "simulation", "post_process"]:
        if stage in stage_coverage:
            stage_coverage[stage]["total_with_common"] = (
                stage_coverage[stage]["constraint_count"] + common_constraints
            )

    return {
        "summary": {
            "total_constraints": len(constraints),
            "total_check_methods": len(methods),
            "methods_with_most_constraints": sorted(
                [(name, info["prohibit_count"]) for name, info in methods.items()],
                key=lambda x: -x[1]
            )[:10],
        },
        "stage_coverage": stage_coverage,
        "methods": methods,
        "constraints_by_method": by_method,
    }


def print_coverage_report():
    """Print coverage report to console."""
    report = generate_coverage_report()

    print("=" * 70)
    print("MFC Case Validator Constraint Coverage Report")
    print("=" * 70)

    print(f"\nTotal constraints (self.prohibit calls): {report['summary']['total_constraints']}")
    print(f"Total check methods: {report['summary']['total_check_methods']}")

    print("\nMethods with most constraints:")
    for method, count in report['summary']['methods_with_most_constraints']:
        print(f"  {method}: {count} constraints")

    print("\nConstraints by stage:")
    for stage, info in report['stage_coverage'].items():
        total = info.get('total_with_common', info['constraint_count'])
        print(f"  {stage}:")
        print(f"    Methods: {info['method_count']}")
        print(f"    Constraints: {info['constraint_count']} (+ common = {total})")

    print("\n" + "=" * 70)
    print("Detailed constraint listing (top methods):")
    print("=" * 70)

    for method, count in report['summary']['methods_with_most_constraints'][:5]:
        print(f"\n{method} ({count} constraints):")
        if method in report['constraints_by_method']:
            for c in report['constraints_by_method'][method][:5]:
                print(f"  L{c['line']}: {c['message'][:60]}")
            if len(report['constraints_by_method'][method]) > 5:
                print(f"  ... and {len(report['constraints_by_method'][method]) - 5} more")


def save_coverage_report(output_path: Path = None):
    """Save coverage report to JSON file."""
    if output_path is None:
        output_path = Path(__file__).parent / "constraint_coverage.json"

    report = generate_coverage_report()

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2)

    return output_path


if __name__ == "__main__":
    print_coverage_report()
    path = save_coverage_report()
    print(f"\nReport saved to: {path}")
