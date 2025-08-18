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

            def raise_err_with_max_diagnostics(msg: str):
                # Find maximum errors across ALL files for diagnostics
                max_errors = find_maximum_errors(candidate, golden)
                diagnostic_msg = format_diagnostic_message(max_errors)

                return None, f"""\
Variable n°{valIndex+1} (1-indexed) in {gFilepath} {msg}:
  - Candidate:   {cVal}
  - Golden:      {gVal}
  - Error:       {error}
  - Tolerance:   {tol}{diagnostic_msg}
"""

            if math.isnan(gVal):
                return raise_err_with_max_diagnostics("is NaN in the golden file")
            if math.isnan(cVal):
                return raise_err_with_max_diagnostics("is NaN in the pack file")
            if not is_close(error, tol):
                return raise_err_with_max_diagnostics("is not within tolerance")

    # Return the average relative error
    return avg_err.get(), None


def format_diagnostic_message(max_errors: typing.Tuple[typing.Optional[typing.Tuple[str, int, float, float, float, float]], typing.Optional[typing.Tuple[str, int, float, float, float, float]]]) -> str:
    """Format the diagnostic message showing maximum errors."""
    max_abs_info, max_rel_info = max_errors
    diagnostic_msg = ""
    
    if max_abs_info:
        filepath, var_idx, golden_val, candidate_val, abs_error, rel_error = max_abs_info
        rel_error_str = f"{rel_error:.2E}" if not math.isnan(rel_error) else "NaN"
        diagnostic_msg += f"\n\nDiagnostics - Maximum absolute error across ALL files:\n" \
                         f" - File: {filepath}\n" \
                         f" - Variable n°{var_idx+1}\n" \
                         f" - Candidate: {candidate_val}\n" \
                         f" - Golden: {golden_val}\n" \
                         f" - Absolute Error: {abs_error:.2E}\n" \
                         f" - Relative Error: {rel_error_str}"

    if max_rel_info:
        filepath, var_idx, golden_val, candidate_val, rel_error, abs_error = max_rel_info
        diagnostic_msg += f"\n\nDiagnostics - Maximum relative error across ALL files:\n" \
                         f" - File: {filepath}\n" \
                         f" - Variable n°{var_idx+1}\n" \
                         f" - Candidate: {candidate_val}\n" \
                         f" - Golden: {golden_val}\n" \
                         f" - Relative Error: {rel_error:.2E}\n" \
                         f" - Absolute Error: {abs_error:.2E}"

    return diagnostic_msg


def find_maximum_errors(candidate: Pack, golden: Pack) -> typing.Tuple[typing.Optional[typing.Tuple[str, int, float, float, float, float]], typing.Optional[typing.Tuple[str, int, float, float, float, float]]]:
    """
    Scan all files to find the maximum absolute and relative errors.
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

            # Track maximum absolute error
            if error.absolute > max_abs_error:
                max_abs_error = error.absolute
                max_abs_info = (gFilepath, valIndex, gVal, cVal, error.absolute, error.relative)

            # Track maximum relative error (only if it's not NaN)
            if not math.isnan(error.relative) and error.relative > max_rel_error:
                max_rel_error = error.relative
                max_rel_info = (gFilepath, valIndex, gVal, cVal, error.relative, error.absolute)

    return max_abs_info, max_rel_info
