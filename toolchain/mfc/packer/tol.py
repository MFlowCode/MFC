import math, typing

from ..common import MFCException

from .pack   import Pack
from .errors import compute_error, AverageError, Error

class Tolerance(Error):
    pass


def is_close(error: Error, tolerance: Tolerance) -> bool:
    if error.absolute <= tolerance.absolute:
        return True

    if math.isnan(error.relative):
        return True

    if error.relative <= tolerance.relative:
        return True

    return False


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

        if cEntry == None:
            return None, f"No reference to {gFilepath} in the candidate's pack."

        # Compare variable-count
        if len(gEntry.doubles) != len(cEntry.doubles):
            return None, f"Variable count didn't match for {gFilepath}."

        # Check if each variable is within tolerance
        for valIndex, (gVal, cVal) in enumerate(zip(gEntry.doubles, cEntry.doubles)):
            # Keep track of the error and average errors
            error = compute_error(cVal, gVal)
            avg_err.push(error)

            def raise_err(msg: str):
                return None, f"""\
Variable n°{valIndex+1} (1-indexed) in {gFilepath} {msg}:
  - Candidate:   {cVal}
  - Golden:      {gVal}
  - Error:       {error}
  - Tolerance:   {tol}
"""

            if math.isnan(gVal):
                return raise_err("is NaN in the golden file")

            if math.isnan(cVal):
                return raise_err("is NaN in the pack file")

            if not is_close(error, tol):
                return raise_err("is not within tolerance")

    # Return the average relative error
    return avg_err.get(), None
