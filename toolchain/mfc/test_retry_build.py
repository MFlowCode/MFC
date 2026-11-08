"""Tests for .github/scripts/retry-build.sh.

Frontier's concurrent build variants share one build/ directory, and a compute
node cannot reach PyPI or the dependency sources. A retry must therefore clear
only the slugs its own failed attempt touched: never build/venv, the
build/pyproject.toml install marker, prebuilt dependencies, or another job's slugs.
"""

import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[2]

KEPT = ["venv/bin/python", "pyproject.toml", "staging/fftw/CMakeCache.txt", "install/fftw/lib/libfftw3.a", "staging/gpu-mp-other/x.o", "install/gpu-mp-other/bin/simulation"]
OWN = ["staging/gpu-mp-own/x.o", "install/gpu-mp-own/bin/simulation"]


def _run(tmp_path, attempt_script):
    shutil.copytree(REPO / ".github", tmp_path / ".github")
    for rel in KEPT + OWN:
        (tmp_path / "build" / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / "build" / rel).write_text("")
    (tmp_path / "attempt.sh").write_text(attempt_script)
    return subprocess.run(
        ["bash", "-c", "source .github/scripts/retry-build.sh && retry_build bash attempt.sh"],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "MFC_BUILD_RETRY_DELAY": "0"},
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )


# Attempt 1 logs its slug and fails; attempt 2 checks what survived and succeeds.
ATTEMPT = """
n=$(cat count 2>/dev/null || echo 0); echo $((n + 1)) > count
if [ "$n" = 0 ]; then echo gpu-mp-own >> "$MFC_SLUG_LOG"; exit 1; fi
ls -R build > survived.txt
"""


def test_retry_clears_only_the_failed_attempts_slugs(tmp_path):
    result = _run(tmp_path, ATTEMPT)
    assert result.returncode == 0, result.stdout + result.stderr
    for rel in KEPT:
        assert (tmp_path / "build" / rel).exists(), f"retry removed {rel}"
    for rel in OWN:
        assert not (tmp_path / "build" / rel).exists(), f"retry kept the failed slug's {rel}"


def test_a_failure_before_any_slug_removes_nothing(tmp_path):
    result = _run(tmp_path, 'n=$(cat count 2>/dev/null || echo 0); echo $((n + 1)) > count; [ "$n" != 0 ]\n')
    assert result.returncode == 0, result.stdout + result.stderr
    for rel in KEPT + OWN:
        assert (tmp_path / "build" / rel).exists(), f"retry removed {rel}"


def test_the_toolchain_logs_each_slug_it_builds_but_not_installed_dependencies(tmp_path, monkeypatch):
    from mfc import build

    log = tmp_path / "slugs"
    monkeypatch.setenv("MFC_SLUG_LOG", str(log))
    for target, slug in ((build.SYSCHECK, "gpu-mp-own"), (build.HIPFORT, "hipfort")):
        monkeypatch.setattr(target, "get_slug", lambda case, s=slug: s)
        monkeypatch.setattr(target, "is_buildable", lambda: True)
        monkeypatch.setattr(target, "is_installed", lambda case: True)
        monkeypatch.setattr(target, "is_configured", lambda case: True)
        monkeypatch.setattr(target, "build", lambda case: None)
        monkeypatch.setattr(target, "install", lambda case: None)
    monkeypatch.setattr(build.SYSCHECK, "requires", SimpleNamespace(compute=lambda: [build.HIPFORT]))
    monkeypatch.setattr(build.HIPFORT, "requires", SimpleNamespace(compute=lambda: []))
    monkeypatch.setattr(build, "_stale_dependencies", lambda *a, **k: [])

    vars(build)["__build_target"](build.SYSCHECK, None)

    assert log.read_text().split() == ["gpu-mp-own"]
