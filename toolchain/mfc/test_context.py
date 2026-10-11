"""Tests for ./mfc.sh context (context_expand and context)."""

import re
from pathlib import Path

from mfc import context, context_expand

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


# --- context: index of expanded files ---------------------------------------------------------------------------

EXPANDED = """# 1 "/r/src/common/m_g.fpp"
module m_g
    implicit none
    integer :: dir_idx, nsteps
    type t_box
        real :: a, b
    end type t_box
contains
    subroutine s_set_dir(d)
        integer, intent(in) :: d
        dir_idx = d
    end subroutine s_set_dir
end module m_g
# 1 "/r/src/simulation/m_k.fpp"
module m_k
    use m_g, only: dir_idx
    implicit none
contains
    subroutine s_kernel(q)
        real, intent(inout) :: q(:)
        integer :: i
        real :: a
        do i = 1, 10
            a = q(i)
            q(i) = a + dir_idx
        end do
        call s_helper(q)
    contains
        subroutine s_helper(x)
            real, intent(inout) :: x(:)
            x = 0.
        end subroutine s_helper
    end subroutine s_kernel
end module m_k
"""


def _write(tmp_path, text, target="simulation"):
    p = tmp_path / target / "x.f90"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def _index(tmp_path, text=EXPANDED, target="simulation"):
    return context.build_index({target: [_write(tmp_path, text, target)]})


def test_markers_map_lines(tmp_path):
    rows = context.read_expanded(_write(tmp_path, EXPANDED))
    assert rows[0] == ("module m_g", "/r/src/common/m_g.fpp", 1)
    assert ("module m_k", "/r/src/simulation/m_k.fpp", 1) in rows


def test_type_components_are_not_module_vars(tmp_path):
    idx = _index(tmp_path)
    assert {"dir_idx", "nsteps"} <= set(idx.module_vars)
    assert "a" not in idx.module_vars and "b" not in idx.module_vars


def test_routine_locals_calls_and_use_only(tmp_path):
    idx = _index(tmp_path)
    k = next(r for r in idx.routines if r.name == "s_kernel")
    assert {"i", "a", "q"} <= k.locals
    assert "s_helper" in k.calls
    assert k.module == "m_k" and k.src_file == "/r/src/simulation/m_k.fpp"
    assert idx.module_uses[("simulation", "m_k")] == {"m_g": {"dir_idx"}}
    assert "dir_idx" in next(r for r in idx.routines if r.name == "s_set_dir").assigns
