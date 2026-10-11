"""Tests for ./mfc.sh context (context_expand and context)."""

import re
from pathlib import Path

from mfc import context_expand

REPO = Path(__file__).resolve().parents[2]


def test_fypp_flags_match_cmake():
    text = (REPO / "cmake" / "Fypp.cmake").read_text()
    start = text.index("COMMAND  ${FYPP_EXE}")
    block = text[start : text.index("DEPENDS", start)]
    cmake_flags = sorted(re.sub(r"\s+", " ", f) for f in re.findall(r"(--[\w-]+(?:=\S+)?|-m\s+\w+)", block))
    flags = context_expand.FYPP_FLAGS
    ours = sorted([f for f in flags if f.startswith("--")] + [f"{flags[0]} {flags[1]}"])
    assert cmake_flags == ours


def test_expand_one_target_produces_markers(tmp_path):
    files, errors = context_expand.expand(REPO, ["post_process"], tmp_path)
    assert errors == []
    names = {p.name for p in files["post_process"]}
    assert "m_derived_variables.fpp.f90" in names and "m_helper_basic.fpp.f90" in names
    text = next(p for p in files["post_process"] if p.name == "m_helper_basic.fpp.f90").read_text()
    assert re.search(r'^# \d+ ".*src/common/m_helper_basic\.fpp"', text, re.M)


def test_expand_cache_is_reused(tmp_path):
    files, _ = context_expand.expand(REPO, ["post_process"], tmp_path)
    p = files["post_process"][0].resolve()
    mtime = p.stat().st_mtime_ns
    files2, _ = context_expand.expand(REPO, ["post_process"], tmp_path)
    assert files2["post_process"][0].resolve() == p and p.stat().st_mtime_ns == mtime


def test_expand_reports_fypp_error_and_continues(tmp_path):
    fake = tmp_path / "repo"
    (fake / "src" / "common" / "include").mkdir(parents=True)
    (fake / "src" / "post_process").mkdir(parents=True)
    (fake / "src" / "common" / "m_ok.fpp").write_text("module m_ok\nend module m_ok\n")
    (fake / "src" / "post_process" / "m_bad.fpp").write_text("#:if\nmodule m_bad\nend module m_bad\n")
    files, errors = context_expand.expand(fake, ["post_process"], tmp_path / "cache", write_includes=False)
    assert any("m_bad.fpp" in e for e in errors)
    assert any(p.name == "m_ok.fpp.f90" for p in files["post_process"])
