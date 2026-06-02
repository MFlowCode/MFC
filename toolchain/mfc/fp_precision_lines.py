"""FP-stability precision-lines transform (Tier 2).

A fypp #:for/#:def expansion emits many generated computations that all carry
the same cpp line marker (`# N "file.fpp"`), so DWARF — and therefore Verrou —
collapse every expanded instance onto one .fpp line.  This transform removes the
fypp line markers from a generated .f90 so the compiler attributes each statement
to the generated file's own physical line (which *is* distinct per expanded
instance), and records a sidecar mapping each surviving physical line back to
(file, original .fpp line, instance index).  Genuine cpp directives
(#if/#define/#endif/...) are preserved so conditional compilation is unchanged.

When the stripped .f90 is compiled, Verrou attributes — and fp-stability ranks
and isolates via --source — per expanded instance rather than per source line.
Used only by a dedicated precision build (MFC_FP_PRECISION_LINES); the normal
build is unaffected.  The mechanism (stripped markers -> instance-distinct
physical-line attribution -> per-instance Verrou --source isolation, surviving
the cpp #if layer) is validated against gfortran + Verrou.
"""

import json
import os
import re

# A fypp line marker: "# <number> "<file>"" possibly with trailing flags.  A cpp
# conditional/define directive (#if, #define, #endif, ...) has a word, not a
# number, after the '#', so the two are unambiguous.
_FYPP_MARKER = re.compile(r'^#\s+(\d+)\s+"([^"]+)"')
# Any other preprocessor directive line (kept, but it is not a .fpp source line,
# so it neither consumes a source-line increment nor gets a sidecar entry).
_CPP_DIRECTIVE = re.compile(r"^\s*#")


def strip_markers(lines: list) -> tuple:
    """Strip fypp line markers; return (output_lines, sidecar).

    sidecar maps each 1-based physical output line number to
    {"file", "line", "instance"}: the .fpp file, the .fpp line that physical
    line came from (auto-incremented within a marker region), and how many times
    that marker's (file, line) had been seen before (0 = first/real occurrence,
    >=1 = an expanded instance).
    """
    seen = {}
    out = []
    sidecar = {}
    cur_file = None
    cur_line = None
    cur_instance = None
    for raw in lines:
        m = _FYPP_MARKER.match(raw)
        if m:
            cur_file = m.group(2)
            cur_line = int(m.group(1))
            cur_instance = seen.get((cur_file, cur_line), 0)
            seen[(cur_file, cur_line)] = cur_instance + 1
            continue  # drop the marker line
        out.append(raw)
        if cur_file is None or _CPP_DIRECTIVE.match(raw):
            # cpp directives are kept verbatim but are not .fpp source lines
            continue
        sidecar[len(out)] = {"file": cur_file, "line": cur_line, "instance": cur_instance}
        cur_line += 1  # subsequent physical source lines map to the next .fpp line
    return out, sidecar


def transform_file(in_path: str, out_path: str, sidecar_path: str) -> int:
    """Strip a generated .f90 to its precision-lines variant.

    Reads in_path, writes the marker-stripped source to out_path and the sidecar
    JSON to sidecar_path.  Returns the number of mapped physical lines.
    """
    with open(in_path) as fh:
        lines = fh.readlines()
    out, sidecar = strip_markers(lines)
    with open(out_path, "w") as fh:
        fh.writelines(out)
    with open(sidecar_path, "w") as fh:
        json.dump({str(k): v for k, v in sidecar.items()}, fh)
    return len(sidecar)


# --- consumption side (Tier 2): locating and querying the sidecars ---


def sidecar_dir_for_binary(sim_bin: str) -> str:
    """Map a precision simulation binary path to its sidecar directory.

    .../build/install/<hash>/bin/simulation -> .../build/staging/<hash>/fypp/simulation
    """
    bin_dir = os.path.dirname(os.path.abspath(sim_bin))  # .../install/<hash>/bin
    hash_dir = os.path.dirname(bin_dir)  # .../install/<hash>
    cfg_hash = os.path.basename(hash_dir)
    build_root = os.path.dirname(os.path.dirname(hash_dir))  # .../build
    return os.path.join(build_root, "staging", cfg_hash, "fypp", "simulation")


def sidecar_path(sidecar_dir: str, fpp_file: str) -> str:
    """Sidecar JSON path for a .fpp file: <dir>/<basename>.linemap.json."""
    return os.path.join(sidecar_dir, os.path.basename(fpp_file) + ".linemap.json")


def load_sidecar(path: str) -> dict:
    """Load a sidecar JSON into {physical_line:int -> {file, line, instance}}."""
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        raw = json.load(fh)
    return {int(k): v for k, v in raw.items()}


def instances_of(sidecar: dict, fpp_file: str, fpp_line: int) -> list:
    """Return [(physical_line, instance), ...] (sorted by physical line) for every
    expanded instance of fpp_file:fpp_line, matched by basename."""
    base = os.path.basename(fpp_file)
    hits = [(physline, entry["instance"]) for physline, entry in sidecar.items() if os.path.basename(entry["file"]) == base and entry["line"] == fpp_line]
    return sorted(hits)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 4:
        sys.exit("usage: fp_precision_lines.py <in.f90> <out.f90> <sidecar.json>")
    transform_file(sys.argv[1], sys.argv[2], sys.argv[3])
