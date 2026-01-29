"""
Mutation Testing for Validator Coverage.

Takes valid example cases and systematically mutates parameters
to verify the validator catches invalid configurations.
"""

import json
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple
from dataclasses import dataclass

from ..case_validator import CaseValidator


@dataclass
class MutationResult:
    """Result of a mutation test."""
    case_name: str
    param_name: str
    original_value: Any
    mutated_value: Any
    validator_caught: bool
    errors: List[str]


# Mutations to apply to parameters
MUTATIONS = {
    # Numeric parameters: try invalid values
    "m": [0, -1, None],
    "n": [-1, -10],
    "p": [-1, -5],
    "dt": [0, -1e-6, None],
    "t_step_start": [-1],
    "t_step_stop": [-1],
    "t_step_save": [0, -1],
    "num_fluids": [0, -1],
    "num_patches": [0, -1],
    "model_eqns": [0, 5, 10, -1],
    "weno_order": [0, 2, 4, 6, 8],
    "time_stepper": [0, 6, -1],
    "riemann_solver": [0, 10, -1],

    # Boolean parameters: try invalid strings
    "bubbles_euler": ["X", "yes", "1"],
    "mpp_lim": ["X", "yes"],
    "cyl_coord": ["X", "maybe"],

    # Boundary conditions
    "bc_x%beg": [None, 100, -100],
    "bc_x%end": [None, 100, -100],
    "bc_y%beg": [100, -100],
    "bc_y%end": [100, -100],

    # Remove required parameters
    "x_domain%beg": [None],
    "x_domain%end": [None],
}


def load_example_case(case_path: Path) -> Dict[str, Any]:
    """Load parameters from an example case file."""
    result = subprocess.run(
        ["python3", str(case_path)],
        capture_output=True,
        text=True,
        cwd=case_path.parent,
        timeout=30,
        check=False
    )
    if result.returncode != 0:
        return None
    return json.loads(result.stdout.strip())


def run_mutation(params: Dict[str, Any], param_name: str,
                 mutated_value: Any) -> Tuple[bool, List[str]]:
    """Apply mutation and check if validator catches it."""
    mutated_params = params.copy()

    if mutated_value is None:
        # Remove the parameter
        mutated_params.pop(param_name, None)
    else:
        mutated_params[param_name] = mutated_value

    validator = CaseValidator(mutated_params)

    try:
        validator.validate_pre_process()
    except Exception:
        pass

    try:
        validator.validate_simulation()
    except Exception:
        pass

    try:
        validator.validate_post_process()
    except Exception:
        pass

    return len(validator.errors) > 0, validator.errors


def run_mutations_on_case(case_name: str, params: Dict[str, Any]) -> List[MutationResult]:
    """Run all applicable mutations on a case."""
    results = []

    for param_name, mutations in MUTATIONS.items():
        if param_name not in params:
            continue

        original = params[param_name]

        for mutated_value in mutations:
            # Skip if mutation is same as original
            if mutated_value == original:
                continue

            caught, errors = run_mutation(params, param_name, mutated_value)

            results.append(MutationResult(
                case_name=case_name,
                param_name=param_name,
                original_value=original,
                mutated_value=mutated_value,
                validator_caught=caught,
                errors=errors[:3]  # Limit for memory
            ))

    return results


def run_mutation_tests(max_cases: int = 10) -> Dict[str, Any]:
    """Run mutation tests on example cases."""
    examples_dir = Path(__file__).parent.parent.parent.parent / "examples"
    case_files = sorted(examples_dir.glob("**/case.py"))[:max_cases]

    all_results = []
    cases_tested = 0

    for case_file in case_files:
        case_name = str(case_file.relative_to(examples_dir).parent)
        params = load_example_case(case_file)

        if params is None:
            continue

        cases_tested += 1
        results = run_mutations_on_case(case_name, params)
        all_results.extend(results)

    # Summarize
    total = len(all_results)
    caught = sum(1 for r in all_results if r.validator_caught)
    missed = sum(1 for r in all_results if not r.validator_caught)

    # Group by parameter
    by_param = {}
    for r in all_results:
        if r.param_name not in by_param:
            by_param[r.param_name] = {"caught": 0, "missed": 0}
        if r.validator_caught:
            by_param[r.param_name]["caught"] += 1
        else:
            by_param[r.param_name]["missed"] += 1

    return {
        "cases_tested": cases_tested,
        "total_mutations": total,
        "caught": caught,
        "missed": missed,
        "catch_rate": caught / total * 100 if total > 0 else 0,
        "by_param": by_param,
        "missed_details": [r for r in all_results if not r.validator_caught][:20],
    }


def print_mutation_report():
    """Print mutation test results."""
    print("Running mutation tests on example cases...")
    print("(This tests that the validator catches invalid parameter values)")
    print()

    results = run_mutation_tests(max_cases=20)

    print("=" * 70)
    print("MUTATION TEST RESULTS")
    print("=" * 70)
    print(f"\nCases tested: {results['cases_tested']}")
    print(f"Total mutations: {results['total_mutations']}")
    print(f"Caught by validator: {results['caught']}")
    print(f"Missed by validator: {results['missed']}")
    print(f"Catch rate: {results['catch_rate']:.1f}%")

    print("\n" + "-" * 70)
    print("BY PARAMETER:")
    print("-" * 70)
    for param, data in sorted(results["by_param"].items(),
                               key=lambda x: -x[1]["missed"]):
        total = data["caught"] + data["missed"]
        rate = data["caught"] / total * 100 if total > 0 else 0
        status = "OK" if data["missed"] == 0 else "GAPS"
        print(f"  {param}: {data['caught']}/{total} caught ({rate:.0f}%) [{status}]")

    if results["missed_details"]:
        print("\n" + "-" * 70)
        print("SAMPLE OF UNCAUGHT MUTATIONS (potential validator gaps):")
        print("-" * 70)
        for r in results["missed_details"][:10]:
            print(f"  {r.case_name}")
            print(f"    {r.param_name}: {r.original_value} -> {r.mutated_value}")
            print(f"    No validation error raised!")
            print()


if __name__ == "__main__":
    print_mutation_report()
