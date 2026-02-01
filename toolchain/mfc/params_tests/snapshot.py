"""
Validation Snapshot Tool.

Captures validation results from case files for regression testing.
This allows us to verify that refactoring doesn't change validation behavior.
"""

import json
import hashlib
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, asdict

from ..case_validator import CaseValidator


@dataclass
class ValidationResult:
    """Result of validating a single case file for a single stage."""
    case_path: str
    stage: str
    success: bool
    errors: List[str]
    param_hash: str  # Hash of parameters for change detection
    error_count: int


@dataclass
class CaseSnapshot:
    """Complete validation snapshot for a case file."""
    case_path: str
    param_count: int
    param_hash: str
    stages: Dict[str, ValidationResult]
    load_error: Optional[str] = None


def hash_params(params: Dict[str, Any]) -> str:
    """Create a hash of parameters for change detection."""
    # Sort keys for consistent hashing
    sorted_items = sorted(params.items(), key=lambda x: x[0])
    param_str = json.dumps(sorted_items, sort_keys=True, default=str)
    return hashlib.md5(param_str.encode()).hexdigest()[:12]


def validate_case_for_stage(params: Dict[str, Any], stage: str) -> ValidationResult:
    """Run validation for a specific stage and capture results."""
    validator = CaseValidator(params)

    try:
        if stage == "pre_process":
            validator.validate_pre_process()
        elif stage == "simulation":
            validator.validate_simulation()
        elif stage == "post_process":
            validator.validate_post_process()
        else:
            raise ValueError(f"Unknown stage: {stage}")

        return ValidationResult(
            case_path="",  # Will be filled in by caller
            stage=stage,
            success=len(validator.errors) == 0,
            errors=validator.errors.copy(),
            param_hash=hash_params(params),
            error_count=len(validator.errors)
        )
    except (ValueError, KeyError, TypeError, AttributeError) as e:
        # Catch expected validation errors, not programming bugs like SystemExit
        return ValidationResult(
            case_path="",
            stage=stage,
            success=False,
            errors=[f"Exception during validation: {type(e).__name__}: {str(e)}"],
            param_hash=hash_params(params),
            error_count=1
        )


def load_case_params(case_path: Path) -> Dict[str, Any]:
    """Load parameters from a case file by running it and capturing JSON output."""
    # MFC case files print JSON to stdout when run
    result = subprocess.run(
        ["python3", str(case_path)],
        capture_output=True,
        text=True,
        cwd=case_path.parent,
        timeout=30,
        check=False
    )

    if result.returncode != 0:
        raise ValueError(f"Case file failed: {result.stderr[:200]}")

    # Parse the JSON output
    output = result.stdout.strip()
    if not output:
        raise ValueError("Case file produced no output")

    return json.loads(output)


def capture_case_snapshot(case_path: Path) -> CaseSnapshot:
    """Capture complete validation snapshot for a case file."""
    case_path = Path(case_path)

    try:
        params = load_case_params(case_path)
    except (ValueError, json.JSONDecodeError, subprocess.TimeoutExpired,
            subprocess.SubprocessError, OSError, FileNotFoundError) as e:
        # Catch expected case loading errors, not programming bugs
        return CaseSnapshot(
            case_path=str(case_path),
            param_count=0,
            param_hash="",
            stages={},
            load_error=f"{type(e).__name__}: {str(e)}"
        )

    stages = {}
    for stage in ["pre_process", "simulation", "post_process"]:
        result = validate_case_for_stage(params, stage)
        result.case_path = str(case_path)
        stages[stage] = result

    return CaseSnapshot(
        case_path=str(case_path),
        param_count=len(params),
        param_hash=hash_params(params),
        stages=stages
    )


def capture_all_examples(examples_dir: Path = None) -> Dict[str, CaseSnapshot]:
    """Capture validation snapshots for all example cases."""
    if examples_dir is None:
        examples_dir = Path(__file__).parent.parent.parent.parent / "examples"

    snapshots = {}

    # Find all case.py files
    case_files = sorted(examples_dir.glob("**/case.py"))

    for case_file in case_files:
        relative_path = case_file.relative_to(examples_dir)
        case_name = str(relative_path.parent)

        print(f"  Capturing: {case_name}...", end=" ", flush=True)
        try:
            snapshot = capture_case_snapshot(case_file)
            snapshots[case_name] = snapshot

            if snapshot.load_error:
                print(f"LOAD ERROR: {snapshot.load_error[:50]}")
            else:
                errors = sum(s.error_count for s in snapshot.stages.values())
                if errors > 0:
                    print(f"ERRORS: {errors}")
                else:
                    print("OK")
        except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as e:
            # Catch expected errors during capture, not programming bugs
            print(f"EXCEPTION: {e}")
            snapshots[case_name] = CaseSnapshot(
                case_path=str(case_file),
                param_count=0,
                param_hash="",
                stages={},
                load_error=f"Capture exception: {type(e).__name__}: {str(e)}"
            )

    return snapshots


def snapshot_to_dict(snapshot: CaseSnapshot) -> Dict[str, Any]:
    """Convert snapshot to JSON-serializable dict."""
    result = asdict(snapshot)
    # Convert ValidationResult objects in stages
    result["stages"] = {
        stage: asdict(vr) for stage, vr in snapshot.stages.items()
    }
    return result


def save_snapshots(snapshots: Dict[str, CaseSnapshot], output_path: Path = None):
    """Save snapshots to JSON file."""
    if output_path is None:
        output_path = Path(__file__).parent / "validation_snapshots.json"

    data = {
        "metadata": {
            "total_cases": len(snapshots),
            "load_errors": sum(1 for s in snapshots.values() if s.load_error),
            "validation_errors": sum(
                sum(stage.error_count for stage in s.stages.values())
                for s in snapshots.values() if not s.load_error
            ),
        },
        "snapshots": {
            name: snapshot_to_dict(snapshot)
            for name, snapshot in snapshots.items()
        }
    }

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)

    return output_path


def load_snapshots(input_path: Path) -> Dict[str, Any]:
    """Load snapshots from JSON file."""
    with open(input_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def compare_snapshots(
    old_snapshots: Dict[str, Any],
    new_snapshots: Dict[str, CaseSnapshot]
) -> Dict[str, Any]:
    """Compare old and new snapshots, report differences."""
    differences = {
        "new_cases": [],
        "removed_cases": [],
        "changed_validation": [],
        "unchanged": [],
    }

    old_cases = set(old_snapshots.get("snapshots", {}).keys())
    new_cases = set(new_snapshots.keys())

    differences["new_cases"] = sorted(new_cases - old_cases)
    differences["removed_cases"] = sorted(old_cases - new_cases)

    for case_name in sorted(old_cases & new_cases):
        old_snap = old_snapshots["snapshots"][case_name]
        new_snap = snapshot_to_dict(new_snapshots[case_name])

        # Compare validation results
        changed = False
        changes = []

        for stage in ["pre_process", "simulation", "post_process"]:
            old_stage = old_snap.get("stages", {}).get(stage, {})
            new_stage = new_snap.get("stages", {}).get(stage, {})

            old_errors = set(old_stage.get("errors", []))
            new_errors = set(new_stage.get("errors", []))

            if old_errors != new_errors:
                changed = True
                changes.append({
                    "stage": stage,
                    "old_error_count": len(old_errors),
                    "new_error_count": len(new_errors),
                    "added_errors": sorted(new_errors - old_errors),
                    "removed_errors": sorted(old_errors - new_errors),
                })

        if changed:
            differences["changed_validation"].append({
                "case": case_name,
                "changes": changes
            })
        else:
            differences["unchanged"].append(case_name)

    return differences


def print_comparison_report(differences: Dict[str, Any]):
    """Print a human-readable comparison report."""
    print("=" * 60)
    print("Validation Comparison Report")
    print("=" * 60)

    print(f"\nNew cases: {len(differences['new_cases'])}")
    for case in differences['new_cases'][:5]:
        print(f"  + {case}")
    if len(differences['new_cases']) > 5:
        print(f"  ... and {len(differences['new_cases']) - 5} more")

    print(f"\nRemoved cases: {len(differences['removed_cases'])}")
    for case in differences['removed_cases'][:5]:
        print(f"  - {case}")

    print(f"\nChanged validation: {len(differences['changed_validation'])}")
    for item in differences['changed_validation'][:10]:
        print(f"\n  {item['case']}:")
        for change in item['changes']:
            print(f"    [{change['stage']}] {change['old_error_count']} -> {change['new_error_count']} errors")
            for err in change['added_errors'][:2]:
                print(f"      + {err[:60]}...")
            for err in change['removed_errors'][:2]:
                print(f"      - {err[:60]}...")

    print(f"\nUnchanged: {len(differences['unchanged'])}")


if __name__ == "__main__":
    print("Capturing validation snapshots for all examples...")
    all_snapshots = capture_all_examples()
    path = save_snapshots(all_snapshots)
    print(f"\nSnapshots saved to: {path}")
