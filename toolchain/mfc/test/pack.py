import os, re, math, dataclasses

from pathlib import Path

from .  import case
from .. import common

from ..common import MFCException


@dataclasses.dataclass(repr=False)
class Error:
    absolute: float
    relative: float

    def __repr__(self) -> str:
        return f"abs: {self.absolute:.2E}, rel: {self.relative:.2E}"


def compute_error(measured: float, expected: float) -> Error:
    absolute = abs(measured - expected)

    if expected != 0:
        relative = absolute / abs(expected)
    elif measured == expected:
        relative = 0
    else:
        relative = float("NaN")

    return Error(absolute, relative)


class AverageError:
    accumulated: Error
    count:       int

    def __init__(self) -> None:
        self.accumulated = Error(0, 0)
        self.count       = 0

    def get(self) -> Error:
        if self.count == 0:
            return Error(0, 0)

        return Error(self.accumulated.absolute / self.count,
                     self.accumulated.relative / self.count)

    def push(self, error: Error) -> None:
        # Do not include nans in the result
        # See: compute_error()
        if math.isnan(error.relative):
            return

        self.accumulated.absolute += error.absolute
        self.accumulated.relative += error.relative

        self.count += 1

    def __repr__(self) -> str:
        return self.get().__repr__()


@dataclasses.dataclass(repr=False)
class PackEntry:
    filepath: str
    doubles:  list

    def __repr__(self) -> str:
        return f"{self.filepath} {' '.join([ str(d) for d in self.doubles ])}"


@dataclasses.dataclass
class Pack:
    entries: list

    def save(self, filepath: str):
        common.file_write(filepath, '\n'.join([ str(e) for e in self.entries ]))


def load(filepath: str) -> Pack:
    entries: list = []

    for line in common.file_read(filepath).splitlines():
        if common.isspace(line):
            continue

        arr = line.split(' ')

        entries.append(PackEntry(
            filepath=arr[0],
            doubles=[ float(d) for d in arr[1:] ]
        ))

    return Pack(entries)


def generate(case: case.TestCase) -> Pack:
    entries = []

    case_dir = case.get_dirpath()
    D_dir    = os.path.join(case_dir, "D")

    for filepath in list(Path(D_dir).rglob("*.dat")):
        short_filepath = str(filepath).replace(f'{case_dir}', '')[1:].replace("\\", "/")

        data_content = common.file_read(filepath)

        # 2 or more (contiguous) spaces
        pattern = r"([ ]{2,})"

        numbers_str = re.sub(pattern, " ", data_content.replace('\n', '')).strip()

        doubles: list = [ float(e) for e in numbers_str.split(' ') ]

        for double in doubles:
            if math.isnan(double):
                raise MFCException(f"Test {case}: A NaN was found while generating a pack file.")

        entries.append(PackEntry(short_filepath,doubles))

    return Pack(entries)


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


def check_tolerance(case: case.TestCase, candidate: Pack, golden: Pack, tol: float) -> Error:
    # Keep track of the average error
    avg_err = AverageError()

    # Compare entry-count
    if len(candidate.entries) != len(golden.entries):
        raise MFCException(f"Test {case}: Line count didn't match.")

    # For every entry in the golden's pack
    for gIndex, gEntry in enumerate(golden.entries):
        # Find the corresponding entry in the candidate's pack
        cIndex, cEntry = common.find(lambda i, e: e.filepath == gEntry.filepath, candidate.entries)

        if cIndex == None:
            raise MFCException(f"Test {case}: No reference of {gEntry.filepath} in the candidate's pack.")

        filepath: str = gEntry.filepath

        # Compare variable-count
        if len(gEntry.doubles) != len(cEntry.doubles):
            raise MFCException(f"Test {case}: Variable count didn't match for {filepath}.")

        # Check if each variable is within tolerance
        for valIndex, (gVal, cVal) in enumerate(zip(gEntry.doubles, cEntry.doubles)):
            # Keep track of the error and average errors
            error = compute_error(cVal, gVal)
            avg_err.push(error)

            def raise_err(msg: str):
                raise MFCException(f"""\
Test {case}: Variable n°{valIndex+1} (1-indexed) in {filepath} {msg}:
  - Description: {case.trace}
  - Candidate:   {cVal}
  - Golden:      {gVal}
  - Error:       {error}
  - Tolerance:   {tol}
""")

            if math.isnan(gVal):
                raise_err("is NaN in the golden file")

            if math.isnan(cVal):
                raise_err("is NaN in the pack file")

            if not is_close(error, Tolerance(absolute=tol, relative=tol)):
                raise_err("is not within tolerance")

    # Return the average relative error
    return avg_err.get()
