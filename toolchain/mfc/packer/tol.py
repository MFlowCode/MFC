import math, typing

from .pack   import Pack
from .errors import compute_error, AverageError, Error

Tolerance = Error

def is_close(error: Error, tolerance: Tolerance) -> bool:
    if error.absolute <= tolerance.absolute:
        return True

    if math.isnan(error.relative):
        return True

    if error.relative <= tolerance.relative:
        return True

    return False


def _format_error_diagnostics(max_abs_info, max_rel_info) -> str:
    """Format diagnostic information for maximum errors among failing variables."""
    diagnostic_msg = ""

    if max_abs_info:
        filepath, val_idx, g_val, c_val, abs_err, rel_err = max_abs_info
        rel_error_str = f"{rel_err:.2E}" if not math.isnan(rel_err) else "NaN"
        diagnostic_msg += f"\n\nDiagnostics - Maximum absolute error among FAILING variables:\n" \
                         f" - File: {filepath}\n" \
                         f" - Variable n°{val_idx+1}\n" \
                         f" - Candidate: {c_val}\n" \
                         f" - Golden: {g_val}\n" \
                         f" - Absolute Error: {abs_err:.2E}\n" \
                         f" - Relative Error: {rel_error_str}"

    if max_rel_info:
        filepath, val_idx, g_val, c_val, rel_err, abs_err = max_rel_info
        diagnostic_msg += f"\n\nDiagnostics - Maximum relative error among FAILING variables:\n" \
                         f" - File: {filepath}\n" \
                         f" - Variable n°{val_idx+1}\n" \
                         f" - Candidate: {c_val}\n" \
                         f" - Golden: {g_val}\n" \
                         f" - Relative Error: {rel_err:.2E}\n" \
                         f" - Absolute Error: {abs_err:.2E}"

    return diagnostic_msg


# pylint: disable=too-many-return-statements
def compare(candidate: Pack, golden: Pack, tol: Tolerance) -> typing.Tuple[Error, str]:
    # Keep track of the average error
    avg_err = AverageError()

    # Compare entry-count
    if len(candidate.entries) != len(golden.entries):
        return None, "Line count does not match."

    # For every entry in the golden's pack
    for gFilepath, gEntry in golden.entries.items():
        # Find the corresponding entry in the candidate's pack
        cEntry = candidate.find(gFilepath)

        if cEntry is None:
            return None, f"No reference to {gFilepath} in the candidate's pack."

        # Compare variable-count
        if len(gEntry.doubles) != len(cEntry.doubles):
            return None, f"Variable count didn't match for {gFilepath}."

        # Check if each variable is within tolerance
        for valIndex, (gVal, cVal) in enumerate(zip(gEntry.doubles, cEntry.doubles)):
            # Keep track of the error and average errors
            error = compute_error(cVal, gVal)
            avg_err.push(error)

            def raise_err_with_failing_diagnostics(msg: str):
                # Find maximum errors among FAILING variables only
                max_abs_info, max_rel_info = find_maximum_errors_among_failing(candidate, golden, tol)
                diagnostic_msg = _format_error_diagnostics(max_abs_info, max_rel_info)

                return None, f"""\
Variable n°{valIndex+1} (1-indexed) in {gFilepath} {msg}:
  - Candidate:   {cVal}
  - Golden:      {gVal}
  - Error:       {error}
  - Tolerance:   {tol}{diagnostic_msg}
"""

            if math.isnan(gVal):
                return raise_err_with_failing_diagnostics("is NaN in the golden file")
            if math.isnan(cVal):
                return raise_err_with_failing_diagnostics("is NaN in the pack file")
            if not is_close(error, tol):
                return raise_err_with_failing_diagnostics("is not within tolerance")

    # Return the average relative error
    return avg_err.get(), None


def find_maximum_errors_among_failing(candidate: Pack, golden: Pack, tol: Tolerance) -> typing.Tuple[typing.Optional[typing.Tuple[str, int, float, float, float, float]], typing.Optional[typing.Tuple[str, int, float, float, float, float]]]:
    """
    Scan all files to find the maximum absolute and relative errors among FAILING variables only.
    A variable fails if is_close(error, tol) returns False.

    Returns tuple of:
    - max_abs_info: (filepath, var_index, golden_val, candidate_val, absolute_error, relative_error)
    - max_rel_info: (filepath, var_index, golden_val, candidate_val, relative_error, absolute_error)
    """
    max_abs_error = -1.0
    max_abs_info = None

    max_rel_error = -1.0
    max_rel_info = None

    for gFilepath, gEntry in golden.entries.items():
        cEntry = candidate.find(gFilepath)
        if cEntry is None:
            continue

        for valIndex, (gVal, cVal) in enumerate(zip(gEntry.doubles, cEntry.doubles)):
            # Skip NaN values in golden or candidate
            if math.isnan(gVal) or math.isnan(cVal):
                continue

            error = compute_error(cVal, gVal)

            # Only consider variables that FAIL tolerance
            if is_close(error, tol):
                continue  # Skip variables that pass tolerance

            # Track maximum absolute error among failing variables
            if error.absolute > max_abs_error:
                max_abs_error = error.absolute
                max_abs_info = (gFilepath, valIndex, gVal, cVal, error.absolute, error.relative)

            # Track maximum relative error among failing variables (only if it's not NaN)
            if not math.isnan(error.relative) and error.relative > max_rel_error:
                max_rel_error = error.relative
                max_rel_info = (gFilepath, valIndex, gVal, cVal, error.relative, error.absolute)

    return max_abs_info, max_rel_info
