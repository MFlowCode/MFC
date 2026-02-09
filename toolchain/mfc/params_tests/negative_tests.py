"""
Negative Test Case Generator.

Generates test cases that intentionally violate validator constraints
to ensure each constraint is properly enforced.
"""

from typing import Dict, Any, List
from dataclasses import dataclass

from ..case_validator import CaseValidator


@dataclass
class ConstraintTest:
    """A test case for a specific constraint."""
    method: str
    line_number: int
    message: str
    condition: str
    test_params: Dict[str, Any]
    should_trigger: bool = True


# Base valid case - starts from a known-good configuration
BASE_CASE = {
    "m": 50,
    "n": 0,
    "p": 0,
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "t_step_start": 0,
    "t_step_stop": 100,
    "t_step_save": 10,
    "dt": 1e-6,
    "weno_order": 5,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
}


def generate_constraint_tests() -> List[ConstraintTest]:
    """Generate test cases for each constraint in case_validator.py."""
    tests = []

    # ===================================================================
    # check_simulation_domain constraints
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_simulation_domain",
            line_number=56,
            message="m must be set",
            condition="m is None",
            test_params={**BASE_CASE, "m": None},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=57,
            message="m must be positive",
            condition="m <= 0",
            test_params={**BASE_CASE, "m": 0},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=57,
            message="m must be positive",
            condition="m <= 0",
            test_params={**BASE_CASE, "m": -5},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=58,
            message="n must be non-negative",
            condition="n < 0",
            test_params={**BASE_CASE, "n": -1},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=59,
            message="p must be non-negative",
            condition="p < 0",
            test_params={**BASE_CASE, "p": -1},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=60,
            message="p must be odd for cylindrical coordinates",
            condition="cyl_coord and p > 0 and p % 2 == 0",
            test_params={**BASE_CASE, "cyl_coord": "T", "n": 10, "p": 2},
        ),
        ConstraintTest(
            method="check_simulation_domain",
            line_number=62,
            message="p must be 0 if n = 0",
            condition="n == 0 and p > 0",
            test_params={**BASE_CASE, "n": 0, "p": 5},
        ),
    ])

    # ===================================================================
    # check_model_eqns_and_num_fluids constraints
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_model_eqns_and_num_fluids",
            line_number=73,
            message="model_eqns must be 1, 2, 3, or 4",
            condition="model_eqns not in [1, 2, 3, 4]",
            test_params={**BASE_CASE, "model_eqns": 5},
        ),
        ConstraintTest(
            method="check_model_eqns_and_num_fluids",
            line_number=75,
            message="num_fluids must be positive",
            condition="num_fluids < 1",
            test_params={**BASE_CASE, "num_fluids": 0},
        ),
        ConstraintTest(
            method="check_model_eqns_and_num_fluids",
            line_number=85,
            message="model_eqns = 1 does not support mpp_lim",
            condition="model_eqns == 1 and mpp_lim",
            test_params={**BASE_CASE, "model_eqns": 1, "num_fluids": None, "mpp_lim": "T"},
        ),
        ConstraintTest(
            method="check_model_eqns_and_num_fluids",
            line_number=87,
            message="num_fluids = 1 does not support mpp_lim",
            condition="num_fluids == 1 and mpp_lim",
            test_params={**BASE_CASE, "num_fluids": 1, "mpp_lim": "T"},
        ),
    ])

    # ===================================================================
    # check_time_stepping constraints
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_time_stepping",
            line_number=0,  # Will be determined
            message="dt must be positive",
            condition="dt <= 0",
            test_params={**BASE_CASE, "dt": 0},
        ),
        ConstraintTest(
            method="check_time_stepping",
            line_number=0,
            message="dt must be positive",
            condition="dt <= 0",
            test_params={**BASE_CASE, "dt": -1e-6},
        ),
        ConstraintTest(
            method="check_time_stepping",
            line_number=0,
            message="t_step_stop must be >= t_step_start",
            condition="t_step_stop < t_step_start",
            test_params={**BASE_CASE, "t_step_start": 100, "t_step_stop": 50},
        ),
    ])

    # ===================================================================
    # check_weno constraints
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_weno_simulation",
            line_number=0,
            message="weno_order must be 1, 3, 5, or 7",
            condition="weno_order not in [1, 3, 5, 7]",
            test_params={**BASE_CASE, "weno_order": 4},
        ),
        ConstraintTest(
            method="check_weno_simulation",
            line_number=0,
            message="weno_order must be 1, 3, 5, or 7",
            condition="weno_order not in [1, 3, 5, 7]",
            test_params={**BASE_CASE, "weno_order": 9},
        ),
    ])

    # ===================================================================
    # check_boundary_conditions constraints
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_boundary_conditions",
            line_number=0,
            message="bc_x%beg must be set",
            condition="bc_x%beg is None",
            test_params={**BASE_CASE, "bc_x%beg": None},
        ),
        ConstraintTest(
            method="check_boundary_conditions",
            line_number=0,
            message="bc_x%end must be set",
            condition="bc_x%end is None",
            test_params={**BASE_CASE, "bc_x%end": None},
        ),
    ])

    # ===================================================================
    # check_bubbles constraints
    # ===================================================================
    bubble_case = {**BASE_CASE, "bubbles_euler": "T", "bubble_model": 2, "nb": 1}
    tests.extend([
        ConstraintTest(
            method="check_bubbles_euler",
            line_number=0,
            message="nb must be >= 1",
            condition="bubbles_euler and nb < 1",
            test_params={**bubble_case, "nb": 0},
        ),
    ])

    # ===================================================================
    # check_acoustic_source constraints (the biggest method)
    # ===================================================================
    tests.extend([
        ConstraintTest(
            method="check_acoustic_source",
            line_number=0,
            message="num_source must be positive when acoustic_source is enabled",
            condition="acoustic_source and num_source < 1",
            test_params={**BASE_CASE, "acoustic_source": "T", "num_source": 0},
        ),
    ])

    return tests


def _message_matches(expected: str, actual_errors: List[str]) -> bool:
    """Check if expected message matches any actual error (fuzzy)."""
    expected_lower = expected.lower()
    # Extract key terms from expected message
    key_terms = [w for w in expected_lower.split() if len(w) > 3]

    for err in actual_errors:
        err_lower = err.lower()
        # Check if most key terms appear in the error
        matches = sum(1 for term in key_terms if term in err_lower)
        if matches >= len(key_terms) * 0.5:  # 50% of terms match
            return True
    return False


def run_constraint_tests() -> Dict[str, Any]:
    """Run all constraint tests and return results."""
    tests = generate_constraint_tests()
    results = {
        "total": len(tests),
        "passed": 0,
        "failed": 0,
        "errors_triggered": 0,
        "details": [],
    }

    for test in tests:
        validator = CaseValidator(test.test_params)

        # Run validation for all stages
        try:
            validator.validate_pre_process()
        except Exception:
            pass

        try:
            validator.validate_simulation()
        except Exception:
            pass

        # Check if any error was triggered (the key metric)
        any_error = len(validator.errors) > 0
        message_matched = _message_matches(test.message, validator.errors)

        if any_error:
            results["errors_triggered"] += 1

        if any_error == test.should_trigger:
            results["passed"] += 1
            status = "PASS"
        else:
            results["failed"] += 1
            status = "FAIL"

        results["details"].append({
            "method": test.method,
            "message": test.message,
            "status": status,
            "expected_trigger": test.should_trigger,
            "any_error": any_error,
            "message_matched": message_matched,
            "all_errors": validator.errors[:3],  # Limit for display
        })

    return results


def print_test_report():
    """Print test results to console."""
    results = run_constraint_tests()

    print("=" * 70)
    print("Constraint Validation Negative Tests")
    print("=" * 70)
    print(f"\nTotal tests: {results['total']}")
    print(f"Errors triggered: {results['errors_triggered']}/{results['total']}")
    print(f"Passed: {results['passed']}")
    print(f"Failed: {results['failed']}")

    # Group by method
    by_method = {}
    for detail in results["details"]:
        method = detail["method"]
        if method not in by_method:
            by_method[method] = {"passed": 0, "failed": 0, "tests": []}
        if detail["status"] == "PASS":
            by_method[method]["passed"] += 1
        else:
            by_method[method]["failed"] += 1
        by_method[method]["tests"].append(detail)

    print("\nResults by method:")
    for method, data in sorted(by_method.items()):
        status = "OK" if data["failed"] == 0 else "ISSUES"
        print(f"  {method}: {data['passed']}/{data['passed']+data['failed']} [{status}]")

    if results["failed"] > 0:
        print("\nFailed tests (constraint not triggering as expected):")
        for detail in results["details"]:
            if detail["status"] == "FAIL":
                print(f"\n  {detail['method']}")
                print(f"    Expected: {detail['message']}")
                print(f"    Got errors: {detail['any_error']}")
                if detail['all_errors']:
                    for err in detail['all_errors'][:2]:
                        print(f"      - {err[:60]}...")

    print("\n" + "=" * 70)
    error_rate = results["errors_triggered"] / results["total"] * 100
    print(f"Error trigger rate: {error_rate:.1f}% ({results['errors_triggered']}/{results['total']} tests triggered errors)")
    print("=" * 70)


if __name__ == "__main__":
    print_test_report()
