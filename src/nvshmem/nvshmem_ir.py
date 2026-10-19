"""Parse developer-authored MPI/Fypp regions into the NVSHMEM transport IR."""

from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from pathlib import Path

POINT_TO_POINT_CALLS = {
    "BSEND",
    "BSEND_INIT",
    "CANCEL",
    "IBSEND",
    "IMPROBE",
    "IMRECV",
    "IPROBE",
    "IRECV",
    "IRSEND",
    "ISEND",
    "ISSEND",
    "MPROBE",
    "MRECV",
    "PACK",
    "PROBE",
    "RECV",
    "RECV_INIT",
    "REQUEST_FREE",
    "RSEND",
    "RSEND_INIT",
    "SEND",
    "SENDRECV",
    "SENDRECV_REPLACE",
    "SEND_INIT",
    "SSEND",
    "SSEND_INIT",
    "START",
    "STARTALL",
    "TEST",
    "TESTALL",
    "TESTANY",
    "TESTSOME",
    "UNPACK",
    "WAIT",
    "WAITALL",
    "WAITANY",
    "WAITSOME",
}
_BEGIN = re.compile(r"^\s*!\s*MFC_COMM_BEGIN\s+(.*)$")
_END = re.compile(r"^\s*!\s*MFC_COMM_END\s*$")
_CALL = re.compile(rf"\bcall\s+MPI_({'|'.join(sorted(POINT_TO_POINT_CALLS, key=len, reverse=True))})\s*\(", re.IGNORECASE)


@dataclass(frozen=True)
class MpiCall:
    operation: str
    arguments: tuple[str, ...]
    line: int


@dataclass(frozen=True)
class TransportRegion:
    name: str
    mode: str
    pattern: str
    path: Path
    begin_line: int
    end_line: int
    calls: tuple[MpiCall, ...]


def _without_comments(source: str) -> str:
    """Blank Fortran comments while preserving offsets and line numbers."""
    output = list(source)
    quote: str | None = None
    index = 0
    while index < len(source):
        char = source[index]
        if char == "\n":
            quote = None
        elif quote is None and char in {"'", '"'}:
            quote = char
        elif quote is not None and char == quote:
            if index + 1 < len(source) and source[index + 1] == quote:
                index += 1
            else:
                quote = None
        elif quote is None and char == "!":
            while index < len(source) and source[index] != "\n":
                output[index] = " "
                index += 1
            continue
        index += 1
    return "".join(output)


def _split_arguments(arguments: str) -> tuple[str, ...]:
    values: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    index = 0
    while index < len(arguments):
        char = arguments[index]
        if quote is None and char in {"'", '"'}:
            quote = char
        elif quote is not None and char == quote:
            if index + 1 < len(arguments) and arguments[index + 1] == quote:
                index += 1
            else:
                quote = None
        elif quote is None:
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
            elif char == "," and depth == 0:
                values.append(arguments[start:index])
                start = index + 1
        index += 1
    values.append(arguments[start:])
    return tuple(re.sub(r"\s*&\s*", " ", value).strip() for value in values)


def _parse_calls(source: str) -> list[MpiCall]:
    clean = _without_comments(source)
    calls: list[MpiCall] = []
    for match in _CALL.finditer(clean):
        depth = 1
        index = match.end()
        quote: str | None = None
        while index < len(clean) and depth > 0:
            char = clean[index]
            if quote is None and char in {"'", '"'}:
                quote = char
            elif quote is not None and char == quote:
                quote = None
            elif quote is None:
                if char == "(":
                    depth += 1
                elif char == ")":
                    depth -= 1
            index += 1
        if depth != 0:
            line = clean.count("\n", 0, match.start()) + 1
            raise ValueError(f"Unterminated MPI_{match.group(1).upper()} call at line {line}")
        line = clean.count("\n", 0, match.start()) + 1
        calls.append(MpiCall(match.group(1).upper(), _split_arguments(clean[match.end() : index - 1]), line))
    return calls


def _metadata(text: str, path: Path, line: int) -> dict[str, str]:
    fields: dict[str, str] = {}
    for token in shlex.split(text):
        if "=" not in token:
            raise ValueError(f"Invalid communication-region field '{token}' at {path}:{line}")
        key, value = token.split("=", 1)
        fields[key] = value
    missing = {"name", "mode", "pattern"} - fields.keys()
    if missing:
        raise ValueError(f"Missing communication-region fields {sorted(missing)} at {path}:{line}")
    if fields["mode"] not in {"device", "mpi_only"}:
        raise ValueError(f"Invalid communication-region mode at {path}:{line}")
    return fields


def parse_source(path: Path) -> list[TransportRegion]:
    source = path.read_text()
    lines = source.splitlines()
    ranges: list[tuple[int, int, dict[str, str]]] = []
    active: tuple[int, dict[str, str]] | None = None
    for line_number, line in enumerate(lines, 1):
        begin = _BEGIN.match(line)
        if begin:
            if active is not None:
                raise ValueError(f"Nested communication region at {path}:{line_number}")
            active = (line_number, _metadata(begin.group(1), path, line_number))
        elif _END.match(line):
            if active is None:
                raise ValueError(f"Communication-region end without begin at {path}:{line_number}")
            ranges.append((active[0], line_number, active[1]))
            active = None
    if active is not None:
        raise ValueError(f"Unterminated communication region at {path}:{active[0]}")

    calls = _parse_calls(source)
    regions: list[TransportRegion] = []
    covered: set[int] = set()
    for begin, end, fields in ranges:
        selected = tuple(call for call in calls if begin < call.line < end)
        covered.update(call.line for call in selected)
        regions.append(TransportRegion(fields["name"], fields["mode"], fields["pattern"], path, begin, end, selected))
    uncovered = [call for call in calls if call.line not in covered]
    if uncovered:
        locations = ", ".join(f"MPI_{call.operation}@{call.line}" for call in uncovered)
        raise ValueError(f"Unclassified point-to-point MPI in {path}: {locations}")
    return regions


def parse_sources(paths: list[Path]) -> list[TransportRegion]:
    regions = [region for path in paths for region in parse_source(path)]
    names = [region.name for region in regions]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate communication-region names: {duplicates}")
    return regions
