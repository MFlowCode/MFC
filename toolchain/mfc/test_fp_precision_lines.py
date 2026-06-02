"""Unit tests for the fp-stability precision-lines transform (Tier 2, P1).

A fypp #:for/#:def expansion re-marks many generated computations with the same
cpp line marker (`# N "file.fpp"`), so DWARF — and Verrou — collapse every
expanded instance onto one .fpp line.  strip_markers removes the fypp line
markers so the compiler attributes to the generated .f90's own (instance-
distinct) physical lines, and emits a sidecar mapping each surviving physical
line back to (file, original .fpp line, instance index).  Genuine cpp directives
(#if/#define/...) are kept so conditional compilation still works.
"""

import os

from mfc.fp_precision_lines import (
    instances_of,
    sidecar_dir_for_binary,
    sidecar_path,
    strip_markers,
)


def test_strips_fypp_markers_and_keeps_code():
    out, sidecar = strip_markers(['# 700 "real.fpp"\n', "  x = a - b\n"])
    assert out == ["  x = a - b\n"]
    assert sidecar == {1: {"file": "real.fpp", "line": 700, "instance": 0}}


def test_keeps_cpp_conditional_directives():
    lines = ['# 700 "real.fpp"\n', "#if defined(FOO)\n", "  x = 1\n", "#endif\n"]
    out, _ = strip_markers(lines)
    assert out == ["#if defined(FOO)\n", "  x = 1\n", "#endif\n"]


def test_repeated_marker_increments_instance():
    lines = ['# 700 "real.fpp"\n', "  s1 = x\n", '# 700 "real.fpp"\n', "  s2 = y\n"]
    out, sidecar = strip_markers(lines)
    assert out == ["  s1 = x\n", "  s2 = y\n"]
    assert sidecar[1] == {"file": "real.fpp", "line": 700, "instance": 0}
    assert sidecar[2] == {"file": "real.fpp", "line": 700, "instance": 1}


def test_distinguishes_fypp_marker_from_cpp_directive():
    # no fypp line markers here -> nothing stripped, no origin recorded
    lines = ["#define X 1\n", "#if X\n", "  a = 1\n", "#endif\n"]
    out, sidecar = strip_markers(lines)
    assert out == lines
    assert sidecar == {}


def test_source_line_auto_increments_within_a_region():
    lines = ['# 700 "real.fpp"\n', "  a = 1\n", "  b = 2\n"]
    _, sidecar = strip_markers(lines)
    assert sidecar[1]["line"] == 700
    assert sidecar[2]["line"] == 701


# --- Tier 2 consumption: locating + querying sidecars ---


def test_instances_of_returns_physical_lines_for_a_source_line():
    sidecar = {
        7: {"file": "/abs/src/simulation/m_weno.fpp", "line": 241, "instance": 0},
        11: {"file": "/abs/src/simulation/m_weno.fpp", "line": 241, "instance": 1},
        20: {"file": "/abs/src/simulation/m_weno.fpp", "line": 999, "instance": 0},
    }
    # matched by basename; the repo-relative path from a dd_line hotspot still matches
    assert instances_of(sidecar, "src/simulation/m_weno.fpp", 241) == [(7, 0), (11, 1)]


def test_instances_of_empty_when_no_match():
    sidecar = {7: {"file": "m_weno.fpp", "line": 241, "instance": 0}}
    assert instances_of(sidecar, "m_weno.fpp", 999) == []
    assert instances_of(sidecar, "m_other.fpp", 241) == []


def test_instances_of_sorted_by_physical_line():
    sidecar = {
        30: {"file": "f.fpp", "line": 5, "instance": 2},
        10: {"file": "f.fpp", "line": 5, "instance": 0},
        20: {"file": "f.fpp", "line": 5, "instance": 1},
    }
    assert instances_of(sidecar, "f.fpp", 5) == [(10, 0), (20, 1), (30, 2)]


def test_sidecar_dir_for_binary_maps_install_to_staging():
    got = sidecar_dir_for_binary("/x/build/install/HASH/bin/simulation")
    assert got == os.path.join("/x/build/staging/HASH/fypp/simulation")


def test_sidecar_path_uses_fpp_basename_and_linemap_suffix():
    got = sidecar_path("/x/staging/HASH/fypp/simulation", "src/simulation/m_weno.fpp")
    assert got == os.path.join("/x/staging/HASH/fypp/simulation", "m_weno.fpp.linemap.json")


def test_cpp_directives_do_not_consume_a_source_line_increment():
    # the #else line must not advance the .fpp source line nor get a sidecar entry
    lines = ['# 700 "real.fpp"\n', "  a = 1\n", "#else\n", "  b = 2\n"]
    out, sidecar = strip_markers(lines)
    assert out == ["  a = 1\n", "#else\n", "  b = 2\n"]
    assert sidecar[1]["line"] == 700  # a = 1
    assert 2 not in sidecar  # #else: kept, but not a source line
    assert sidecar[3]["line"] == 701  # b = 2 (not 702)


def test_sidecar_line_numbers_are_physical_output_lines():
    # output physical line numbers (1-based, after stripping) are the keys
    lines = ['# 10 "f"\n', "  a = 1\n", '# 20 "f"\n', "  b = 2\n", "  c = 3\n"]
    out, sidecar = strip_markers(lines)
    assert out == ["  a = 1\n", "  b = 2\n", "  c = 3\n"]
    assert sidecar[1] == {"file": "f", "line": 10, "instance": 0}
    assert sidecar[2] == {"file": "f", "line": 20, "instance": 0}
    assert sidecar[3] == {"file": "f", "line": 21, "instance": 0}
