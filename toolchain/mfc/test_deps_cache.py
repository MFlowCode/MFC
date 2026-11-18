"""Tests for .github/scripts/deps-cache.sh.

Every Frontier job rebuilt LAPACK, hipfort, HDF5 and the other dependencies from
source on the login node (14-35 minutes), although they change only with the
compiler, the loaded modules or the pinned versions. The cache stores the
installed libraries and each dependency's install record under a key made of
those inputs, so a job restores them in seconds.

A wrong entry is worse than no entry: a compute node cannot download sources, so
the build tool's own "rebuild a stale dependency" path fails there. Hence the key
covers everything the dependency build reads, an entry is published atomically
and only after a successful build, and an entry that is incomplete or damaged is
a miss.
"""

import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

DEP_FILES = {
    "staging/lapack/CMakeCache.txt": "CMAKE_Fortran_COMPILER:FILEPATH=/opt/ftn\n",
    "staging/lapack/install_manifest.txt": "liblapack.a\n",
    "install/lapack/lib64/liblapack.a": "lapack-bytes",
    # satisfied by a system library: an install record and nothing installed
    "staging/hdf5/CMakeCache.txt": "CMAKE_C_COMPILER:FILEPATH=/opt/cc\n",
    "staging/hdf5/install_manifest.txt": "",
}
NOT_DEPS = {
    "staging/lapack/CMakeFiles/big.o": "scratch",  # build scratch is not cached
    "staging/gpu-acc-0123456789/CMakeCache.txt": "solver",  # MFC targets are not dependencies
    "install/gpu-acc-0123456789/bin/syscheck": "solver",
    "venv/bin/python": "venv",
}


def _workspace(tmp_path, name):
    ws = tmp_path / name
    shutil.copytree(REPO / ".github", ws / ".github")
    (ws / "toolchain" / "dependencies").mkdir(parents=True)
    (ws / "toolchain" / "dependencies" / "CMakeLists.txt").write_text("GIT_TAG v3.12.0\n")
    return ws


def _populate(ws, files):
    for rel, text in files.items():
        path = ws / "build" / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)


def _sh(ws, cache, script, key_extra="--gpu acc"):
    env = {"PATH": "/usr/bin:/bin", "MFC_DEPS_CACHE_DIR": str(cache), "HOME": str(ws)}
    return subprocess.run(
        ["bash", "-c", f"source .github/scripts/deps-cache.sh && deps_cache_init frontier {key_extra!r} && {script}"],
        cwd=ws,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )


def _entries(cache):
    return [p for p in cache.iterdir() if not p.name.startswith(".")] if cache.exists() else []


def test_a_saved_entry_is_restored_into_another_checkout(tmp_path):
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, {**DEP_FILES, **NOT_DEPS})

    assert _sh(a, cache, "deps_cache_restore").stdout.count("miss") == 1
    saved = _sh(a, cache, "deps_cache_save")
    assert saved.returncode == 0, saved.stdout + saved.stderr
    assert len(_entries(cache)) == 1

    restored = _sh(b, cache, "deps_cache_restore")
    assert restored.returncode == 0 and "hit" in restored.stdout, restored.stdout + restored.stderr
    for rel, text in DEP_FILES.items():
        assert (b / "build" / rel).read_text() == text, rel
    for rel in NOT_DEPS:
        assert not (b / "build" / rel).exists(), f"{rel} must not be cached"


def test_nothing_is_saved_without_an_installed_dependency(tmp_path):
    cache = tmp_path / "cache"
    a = _workspace(tmp_path, "a")
    _populate(a, {"staging/lapack/CMakeCache.txt": "configured, never installed\n", **NOT_DEPS})
    assert _sh(a, cache, "deps_cache_save").returncode == 0
    assert not _entries(cache)


def test_a_damaged_entry_is_a_miss(tmp_path):
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, DEP_FILES)
    _sh(a, cache, "deps_cache_save")
    (_entries(cache)[0] / "install" / "lapack" / "lib64" / "liblapack.a").unlink()

    restored = _sh(b, cache, "deps_cache_restore")
    assert restored.returncode == 0 and "miss" in restored.stdout, restored.stdout + restored.stderr
    assert not (b / "build" / "staging").exists() and not (b / "build" / "install").exists()


def test_an_unfinished_entry_is_a_miss(tmp_path):
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, DEP_FILES)
    _sh(a, cache, "deps_cache_save")
    (_entries(cache)[0] / ".complete").unlink()
    assert "miss" in _sh(b, cache, "deps_cache_restore").stdout
    assert not (b / "build" / "staging").exists()


def test_the_key_follows_everything_the_dependency_build_reads(tmp_path):
    cache = tmp_path / "cache"
    a = _workspace(tmp_path, "a")

    def key(**kw):
        return _sh(a, cache, 'echo "$deps_cache_key"', **kw).stdout.strip()

    base = key()
    assert len(base) == 64 and key() == base
    assert key(key_extra="--gpu mp") != base
    (a / "toolchain" / "dependencies" / "CMakeLists.txt").write_text("GIT_TAG v3.13.0\n")
    assert key() != base


def test_an_existing_entry_is_never_replaced(tmp_path):
    """Two jobs can finish the same dependency build at once; the first publish wins and the loser cleans up."""
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, DEP_FILES)
    _populate(b, {**DEP_FILES, "install/lapack/lib64/liblapack.a": "a different build"})
    _sh(a, cache, "deps_cache_save")
    second = _sh(b, cache, "deps_cache_save")
    assert second.returncode == 0, second.stdout + second.stderr
    entry = _entries(cache)
    assert len(entry) == 1 and (entry[0] / "install" / "lapack" / "lib64" / "liblapack.a").read_text() == "lapack-bytes"
    assert not [p for p in cache.iterdir() if p.name.startswith(".tmp")], "the losing job left its temporary folder behind"


def test_an_unusable_cache_folder_disables_the_cache(tmp_path):
    a = _workspace(tmp_path, "a")
    _populate(a, DEP_FILES)
    blocked = tmp_path / "file"
    blocked.write_text("")
    for step in ("deps_cache_restore", "deps_cache_save"):
        result = _sh(a, blocked / "cache", step)
        assert result.returncode == 0, result.stdout + result.stderr
    assert (a / "build" / "staging" / "lapack" / "install_manifest.txt").exists()


def test_off_disables_the_cache(tmp_path):
    a = _workspace(tmp_path, "a")
    _populate(a, DEP_FILES)
    assert _sh(a, "off", "deps_cache_save").returncode == 0
    assert not (a / "off").exists()


def test_a_damaged_entry_is_replaced_by_the_next_save(tmp_path):
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, DEP_FILES)
    _sh(a, cache, "deps_cache_save")
    lib = _entries(cache)[0] / "install" / "lapack" / "lib64" / "liblapack.a"
    lib.unlink()

    assert "saved" in _sh(a, cache, "deps_cache_save").stdout
    assert lib.read_text() == "lapack-bytes" and "hit" in _sh(b, cache, "deps_cache_restore").stdout
    assert [p.name for p in cache.iterdir()] == [_entries(cache)[0].name], "leftovers in the cache folder"


def test_a_copy_that_fails_publishes_nothing(tmp_path):
    """An entry with a library missing would pass its own completeness check: the census is taken from what was copied."""
    cache = tmp_path / "cache"
    a, b = _workspace(tmp_path, "a"), _workspace(tmp_path, "b")
    _populate(a, {**DEP_FILES, "install/silo/lib/libsilo.a": "silo", "staging/silo/install_manifest.txt": "libsilo.a\n"})
    unreadable = a / "build" / "install" / "lapack" / "lib64" / "liblapack.a"
    unreadable.chmod(0)
    try:
        saved = _sh(a, cache, "deps_cache_save")
    finally:
        unreadable.chmod(0o644)
    assert saved.returncode == 0 and "not saved" in saved.stdout, saved.stdout + saved.stderr
    assert not _entries(cache), "a partial entry was published"
    assert not list(cache.iterdir()), "the failed save left files behind"
    assert "miss" in _sh(b, cache, "deps_cache_restore").stdout
