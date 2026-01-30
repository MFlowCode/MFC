"""
Test Safety Net Runner.

Main entry point for building and verifying the parameter validation test suite.
"""
# pylint: disable=import-outside-toplevel

import sys
import json
import argparse
from pathlib import Path

from .inventory import (
    export_parameter_inventory,
    save_inventory,
    print_inventory_summary
)
from .snapshot import (
    capture_all_examples,
    save_snapshots,
    load_snapshots,
    compare_snapshots,
    print_comparison_report
)
from .coverage import (
    generate_coverage_report,
    print_coverage_report,
    save_coverage_report
)


def get_data_dir() -> Path:
    """Get the directory for storing test data."""
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)
    return data_dir


def build_safety_net(verbose: bool = True):
    """
    Build the complete test safety net.

    This captures:
    1. Parameter inventory
    2. Validation snapshots from all examples
    3. Constraint coverage analysis
    """
    data_dir = get_data_dir()

    if verbose:
        print("=" * 70)
        print("Building Parameter Validation Safety Net")
        print("=" * 70)

    # 1. Parameter inventory
    if verbose:
        print("\n[1/3] Exporting parameter inventory...")
    inventory_path = data_dir / "param_inventory.json"
    save_inventory(inventory_path)
    inventory = export_parameter_inventory()
    if verbose:
        print(f"  Total parameters: {inventory['metadata']['total_parameters']}")
        print(f"  Saved to: {inventory_path}")

    # 2. Validation snapshots
    if verbose:
        print("\n[2/3] Capturing validation snapshots from examples...")
    snapshots = capture_all_examples()
    snapshots_path = data_dir / "validation_snapshots.json"
    save_snapshots(snapshots, snapshots_path)
    load_errors = sum(1 for s in snapshots.values() if s.load_error)
    if verbose:
        print(f"  Total cases: {len(snapshots)}")
        print(f"  Load errors: {load_errors}")
        print(f"  Saved to: {snapshots_path}")

    # 3. Constraint coverage
    if verbose:
        print("\n[3/3] Analyzing constraint coverage...")
    coverage_path = data_dir / "constraint_coverage.json"
    save_coverage_report(coverage_path)
    coverage = generate_coverage_report()
    if verbose:
        print(f"  Total constraints: {coverage['summary']['total_constraints']}")
        print(f"  Check methods: {coverage['summary']['total_check_methods']}")
        print(f"  Saved to: {coverage_path}")

    if verbose:
        print("\n" + "=" * 70)
        print("Safety net built successfully!")
        print("=" * 70)
        print(f"\nData stored in: {data_dir}")
        print("\nFiles created:")
        print(f"  - param_inventory.json     ({inventory['metadata']['total_parameters']} params)")
        print(f"  - validation_snapshots.json ({len(snapshots)} cases)")
        print(f"  - constraint_coverage.json  ({coverage['summary']['total_constraints']} constraints)")

    return {
        "inventory": inventory,
        "snapshots": snapshots,
        "coverage": coverage,
    }


def _print_if(verbose: bool, *args, **kwargs):
    """Print only if verbose mode is enabled."""
    if verbose:
        print(*args, **kwargs)


def _print_changes_report(differences: dict, verbose: bool):
    """Print report when validation has changed."""
    if not verbose:
        return
    print("\n" + "=" * 70)
    print("VALIDATION CHANGED!")
    print("=" * 70)
    if differences['changed_validation']:
        print(f"  {len(differences['changed_validation'])} cases have different validation results")
    if differences['removed_cases']:
        print(f"  {len(differences['removed_cases'])} cases were removed")
    print("\nIf this is expected, run 'build' to update the safety net.")


def verify_safety_net(verbose: bool = True) -> bool:
    """
    Verify that current validation matches the captured safety net.

    Returns True if validation is unchanged, False if there are differences.
    """
    data_dir = get_data_dir()
    snapshots_path = data_dir / "validation_snapshots.json"

    if not snapshots_path.exists():
        _print_if(verbose, "ERROR: Safety net not found. Run 'build' first.")
        return False

    _print_if(verbose, "=" * 70)
    _print_if(verbose, "Verifying Parameter Validation Against Safety Net")
    _print_if(verbose, "=" * 70)

    _print_if(verbose, "\nLoading saved snapshots...")
    old_snapshots = load_snapshots(snapshots_path)
    _print_if(verbose, f"  Loaded {len(old_snapshots.get('snapshots', {}))} cases")

    _print_if(verbose, "\nCapturing current validation results...")
    new_snapshots = capture_all_examples()
    _print_if(verbose, f"  Captured {len(new_snapshots)} cases")

    _print_if(verbose, "\nComparing results...")
    differences = compare_snapshots(old_snapshots, new_snapshots)

    if verbose:
        print_comparison_report(differences)

    has_changes = bool(differences['changed_validation'] or differences['removed_cases'])
    if has_changes:
        _print_changes_report(differences, verbose)
        return False

    _print_if(verbose, "\n" + "=" * 70)
    _print_if(verbose, "VALIDATION UNCHANGED - All tests pass!")
    _print_if(verbose, "=" * 70)

    return True


def show_summary():
    """Show summary of captured safety net data."""
    data_dir = get_data_dir()

    print("=" * 70)
    print("Parameter Validation Safety Net Summary")
    print("=" * 70)

    # Inventory
    inventory_path = data_dir / "param_inventory.json"
    if inventory_path.exists():
        with open(inventory_path) as f:
            inventory = json.load(f)
        print("\nParameter Inventory:")
        print(f"  Total parameters: {inventory['metadata']['total_parameters']}")
        print(f"  By stage:")
        print(f"    Common: {inventory['metadata']['common_count']}")
        print(f"    Pre-process: {inventory['metadata']['pre_process_count']}")
        print(f"    Simulation: {inventory['metadata']['simulation_count']}")
        print(f"    Post-process: {inventory['metadata']['post_process_count']}")
    else:
        print("\nParameter Inventory: NOT FOUND")

    # Snapshots
    snapshots_path = data_dir / "validation_snapshots.json"
    if snapshots_path.exists():
        with open(snapshots_path) as f:
            snapshots = json.load(f)
        print("\nValidation Snapshots:")
        print(f"  Total cases: {snapshots['metadata']['total_cases']}")
        print(f"  Load errors: {snapshots['metadata']['load_errors']}")
        print(f"  Validation errors: {snapshots['metadata']['validation_errors']}")
    else:
        print("\nValidation Snapshots: NOT FOUND")

    # Coverage
    coverage_path = data_dir / "constraint_coverage.json"
    if coverage_path.exists():
        with open(coverage_path) as f:
            coverage = json.load(f)
        print("\nConstraint Coverage:")
        print(f"  Total constraints: {coverage['summary']['total_constraints']}")
        print(f"  Check methods: {coverage['summary']['total_check_methods']}")
        print("  Top methods by constraint count:")
        for method, count in coverage['summary']['methods_with_most_constraints'][:5]:
            print(f"    {method}: {count}")
    else:
        print("\nConstraint Coverage: NOT FOUND")


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Parameter Validation Test Safety Net"
    )
    parser.add_argument(
        "command",
        choices=["build", "verify", "summary", "inventory", "coverage",
                 "negative", "mutation"],
        help="Command to run"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Reduce output verbosity"
    )

    args = parser.parse_args()
    verbose = not args.quiet

    if args.command == "build":
        build_safety_net(verbose=verbose)
    elif args.command == "verify":
        success = verify_safety_net(verbose=verbose)
        sys.exit(0 if success else 1)
    elif args.command == "summary":
        show_summary()
    elif args.command == "inventory":
        print_inventory_summary()
    elif args.command == "coverage":
        print_coverage_report()
    elif args.command == "negative":
        from .negative_tests import print_test_report
        print_test_report()
    elif args.command == "mutation":
        from .mutation_tests import print_mutation_report
        print_mutation_report()


if __name__ == "__main__":
    main()
