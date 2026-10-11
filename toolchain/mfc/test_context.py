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


# --- context: diff mapping and report ---------------------------------------------------------------------------

EXPANDED_SRC = [ln for ln in EXPANDED.splitlines() if not ln.startswith("#")]

DIFF_KERNEL = """diff --git a/src/simulation/m_k.fpp b/src/simulation/m_k.fpp
--- a/src/simulation/m_k.fpp
+++ b/src/simulation/m_k.fpp
@@ -10 +10 @@
-            a = q(i)
+            a = 2*q(i)
"""

COPIES = EXPANDED.replace(
    "        do i = 1, 10\n            a = q(i)\n",
    '        do i = 1, 10\n# 10 "/r/src/simulation/m_k.fpp"\n            a = q(i)\n            q(i) = a + dir_idx\n# 10 "/r/src/simulation/m_k.fpp"\n            a = q(i)\n',
)


def _render(idx, diff, side="new", **kw):
    changed, text, new = context.parse_diff(diff, Path("/r"), side)
    return context.render(idx, changed, text, new, Path("/r"), read_source=lambda f: EXPANDED_SRC, **kw)


def test_parse_diff_sides():
    changed, text, new = context.parse_diff(DIFF_KERNEL, Path("/r"), "new")
    assert changed == {"/r/src/simulation/m_k.fpp": {10}} and new == []
    assert any("2*q(i)" in t for t in text["/r/src/simulation/m_k.fpp"])
    assert context.parse_diff(DIFF_KERNEL, Path("/r"), "old")[0] == {"/r/src/simulation/m_k.fpp": {10}}


def test_report_changed_routine_callee_and_global(tmp_path):
    out = _render(_index(tmp_path), DIFF_KERNEL)
    assert "`s_kernel`" in out
    assert "`s_helper`" in out.split("## Callees")[1]
    globals_ = out.split("## Module variables")[1]
    assert "`dir_idx`" in globals_ and "`s_set_dir`" in globals_
    assert "`a`" not in globals_ and "`i`" not in globals_


def test_report_flags_generated_copies(tmp_path):
    assert "emitted 2 times" in _render(_index(tmp_path, COPIES), DIFF_KERNEL)


def test_no_src_changes():
    diff = "diff --git a/docs/x.md b/docs/x.md\n--- a/docs/x.md\n+++ b/docs/x.md\n@@ -1 +1 @@\n-a\n+b\n"
    assert _render(context.Index([], {}, {}), diff).startswith("No changes under src/.")


def test_new_file_is_listed_not_dropped():
    diff = "diff --git a/src/simulation/m_new.fpp b/src/simulation/m_new.fpp\nnew file mode 100644\n--- /dev/null\n+++ b/src/simulation/m_new.fpp\n@@ -0,0 +1,2 @@\n+module m_new\n+end module m_new\n"
    out = _render(context.Index([], {}, {}), diff, side="old")
    assert "New files" in out and "src/simulation/m_new.fpp" in out


def test_inline_is_capped(tmp_path):
    many = "".join(f'# 1 "/r/src/simulation/m_{i}.fpp"\nmodule m_{i}\ncontains\n    subroutine s_{i}()\n        x = 1\n    end subroutine s_{i}\nend module m_{i}\n' for i in range(8))
    idx = _index(tmp_path, many)
    changed = {f"/r/src/simulation/m_{i}.fpp": {4} for i in range(8)}
    out = context.render(idx, changed, {}, [], Path("/r"), read_source=lambda f: ["module m", "contains", "subroutine s", "x = 1", "end subroutine", "end module"])
    assert out.count("```fortran") == 5 and "not shown inline" in out


# --- context: changed macros ------------------------------------------------------------------------------------


def test_macro_sites_through_wrapper_filtered_by_argument(tmp_path):
    inc = tmp_path / "src" / "common" / "include"
    inc.mkdir(parents=True)
    (tmp_path / "src" / "simulation").mkdir(parents=True)
    (inc / "omp.fpp").write_text("#:def OMP_LOOP(collapse=None, reduction=None, &\n    & copy=None)\n    #:set red = reduction\n#:enddef\n")
    (inc / "par.fpp").write_text("#:def GPU_LOOP(**kw)\n    $:OMP_LOOP(**kw)\n#:enddef\n")
    (tmp_path / "src" / "simulation" / "m_a.fpp").write_text("$:GPU_LOOP(collapse=3)\n$:GPU_LOOP(reduction='[[s]]', reductionOp='[+]')\n")
    (tmp_path / "src" / "simulation" / "m_b.fpp").write_text("$:GPU_LOOP(collapse=2, &\n    & reduction='[[t]]')\n")
    f = str(inc / "omp.fpp")
    out = "\n".join(context.macro_sections({f: {3}}, {f: ["    #:set red = reduction"]}, tmp_path))
    assert "OMP_LOOP" in out and "GPU_LOOP" in out
    assert "2 call sites passing reduction" in out
    assert "src/simulation/m_a.fpp: 1" in out and "src/simulation/m_b.fpp: 1" in out


# --- context: command ---------------------------------------------------------------------------------------------


def test_default_base_falls_back_without_origin(tmp_path):
    import subprocess  # pylint: disable=import-outside-toplevel

    env = {"GIT_AUTHOR_NAME": "a", "GIT_AUTHOR_EMAIL": "a@b", "GIT_COMMITTER_NAME": "a", "GIT_COMMITTER_EMAIL": "a@b", "PATH": "/usr/bin:/bin"}
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "commit", "-q", "--allow-empty", "-m", "x"], check=True, env=env)
    base, note = context.default_base(tmp_path)
    assert base == "HEAD" and "origin/master" in note


def test_parse_diff_ignores_context_lines():
    diff = """--- a/src/simulation/m_k.fpp
+++ b/src/simulation/m_k.fpp
@@ -7,7 +7,8 @@
         integer :: i
         real :: a
         do i = 1, 10
-            a = q(i)
+            a = 2*q(i)
+            a = a + 1
             q(i) = a + dir_idx
         end do
         call s_helper(q)
"""
    assert context.parse_diff(diff, Path("/r"), "new")[0] == {"/r/src/simulation/m_k.fpp": {10, 11}}
    assert context.parse_diff(diff, Path("/r"), "old")[0] == {"/r/src/simulation/m_k.fpp": {10}}


def test_parse_diff_old_side_maps_insertions_to_preceding_line():
    diff = """--- a/src/simulation/m_k.fpp
+++ b/src/simulation/m_k.fpp
@@ -8,4 +8,5 @@
         real :: a
         do i = 1, 10
             a = q(i)
+            a = a + 1
             q(i) = a + dir_idx
"""
    assert context.parse_diff(diff, Path("/r"), "old")[0] == {"/r/src/simulation/m_k.fpp": {10}}


def test_parse_diff_new_side_maps_deletions_to_following_line():
    diff = """--- a/src/simulation/m_k.fpp
+++ b/src/simulation/m_k.fpp
@@ -9,3 +9,2 @@
         do i = 1, 10
-            a = q(i)
             q(i) = a + dir_idx
"""
    assert context.parse_diff(diff, Path("/r"), "new")[0] == {"/r/src/simulation/m_k.fpp": {10}}


def test_changed_continuation_line_maps_to_its_statement(tmp_path):
    expanded = '# 1 "/r/src/simulation/m_c.fpp"\nmodule m_c\ncontains\n    subroutine s_c()\n# 5 "/r/src/simulation/m_c.fpp"\n        call s_x(1, 2)\n# 7 "/r/src/simulation/m_c.fpp"\n    end subroutine s_c\nend module m_c\n'
    source = ["module m_c", "contains", "    subroutine s_c()", "", "        call s_x(1, &", "            & 2)", "    end subroutine s_c", "end module m_c"]
    idx = _index(tmp_path, expanded)
    out = context.render(idx, {"/r/src/simulation/m_c.fpp": {6}}, {}, [], Path("/r"), read_source=lambda f: source)
    assert "`s_c`" in out


def test_one_macro_expansion_is_not_reported_as_copies(tmp_path):
    # A macro call expands to several lines that all map to the call's line: one copy, not several.
    expanded = EXPANDED.replace(
        "        do i = 1, 10\n            a = q(i)\n",
        '        do i = 1, 10\n# 10 "/r/src/simulation/m_k.fpp"\n!$acc parallel loop\n# 10 "/r/src/simulation/m_k.fpp"\n!$omp target teams loop\n# 10 "/r/src/simulation/m_k.fpp"\n            a = q(i)\n',
    )
    assert "emitted" not in _render(_index(tmp_path, expanded), DIFF_KERNEL)


def test_same_named_routines_in_different_files_are_kept_apart(tmp_path):
    pre = '# 1 "/r/src/pre_process/m_gp.fpp"\nmodule m_gp\ncontains\n    subroutine s_defaults()\n        x = 1\n    end subroutine s_defaults\nend module m_gp\n'
    sim = '# 1 "/r/src/simulation/m_gp.fpp"\nmodule m_gp\ncontains\n    subroutine s_defaults()\n        y = 2\n    end subroutine s_defaults\nend module m_gp\n'
    idx = context.build_index({"pre_process": [_write(tmp_path, pre, "pre_process")], "simulation": [_write(tmp_path, sim, "simulation")]})
    out = context.render(idx, {"/r/src/simulation/m_gp.fpp": {4}}, {}, [], Path("/r"), read_source=lambda f: ["module m_gp", "contains", "subroutine s_defaults()", "y = 2", "end subroutine", "end module"])
    assert "src/simulation/m_gp.fpp" in out and "src/pre_process/m_gp.fpp" not in out
    assert "(in simulation)" in out
