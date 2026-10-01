"""Tests for detecting installed dependencies that cannot serve the current build.

Dependencies are keyed by name alone, so a LAPACK built under one compiler module
(amdflang, say) was silently reused after switching to another (gfortran), and
post_process then failed to link. Likewise, a dependency the superbuild satisfied
with a system library installs nothing; in an environment where that library is
not visible, the dependent target cannot configure, and nothing ever rebuilt it.
"""

import os

from mfc.build import _stale_dependencies, cmake_cache_compilers, compiler_mismatches

GFORTRAN = "/opt/gcc/bin/gfortran"
AMDFLANG = "/opt/rocm/bin/amdflang"
GCC = "/opt/gcc/bin/cc"


def write_cache(dirpath, fortran=None, c=GCC):
    os.makedirs(dirpath, exist_ok=True)
    lines = ["# This is the CMakeCache file.", "CMAKE_BUILD_TYPE:STRING=Release"]
    if fortran:
        lines.append(f"CMAKE_Fortran_COMPILER:FILEPATH={fortran}")
    if c:
        lines.append(f"CMAKE_C_COMPILER:FILEPATH={c}")
    lines.append("CMAKE_CXX_COMPILER:STRING=CMAKE_CXX_COMPILER-NOTFOUND")
    with open(os.path.join(dirpath, "CMakeCache.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


class Deps:
    def __init__(self, deps):
        self.deps = deps

    def compute(self):
        return self.deps


class FakeTarget:
    def __init__(self, root, name, is_dependency, deps=(), installed=True, buildable=True, install_files=1):
        self.root = root
        self.name = name
        self.isDependency = is_dependency
        self.requires = Deps(list(deps))
        self.installed = installed
        self.buildable = buildable
        install = self.get_install_dirpath(None)
        os.makedirs(install, exist_ok=True)
        for i in range(install_files):
            with open(os.path.join(install, f"lib{i}.a"), "w") as f:
                f.write("x")

    def get_staging_dirpath(self, _case):
        return os.path.join(self.root, "staging", self.name)

    def get_install_dirpath(self, _case):
        return os.path.join(self.root, "install", self.name)

    def is_installed(self, _case):
        return self.installed

    def is_buildable(self):
        return self.buildable


def test_compilers_are_read_from_the_cache_and_notfound_entries_ignored(tmp_path):
    write_cache(str(tmp_path), fortran=GFORTRAN)
    assert cmake_cache_compilers(str(tmp_path)) == {"Fortran": GFORTRAN, "C": GCC}


def test_a_missing_cache_has_no_compilers(tmp_path):
    assert cmake_cache_compilers(str(tmp_path / "nope")) == {}


def test_only_languages_both_sides_recorded_are_compared():
    assert compiler_mismatches({"Fortran": GFORTRAN, "C": GCC}, {"Fortran": AMDFLANG}) == ["Fortran"]
    assert compiler_mismatches({"Fortran": GFORTRAN}, {"C": GCC}) == []
    assert compiler_mismatches({}, {"Fortran": AMDFLANG}) == []


def test_symlinked_paths_to_one_compiler_match(tmp_path):
    real = tmp_path / "amdllvm"
    real.write_text("")
    link = tmp_path / "amdflang"
    link.symlink_to(real)
    assert compiler_mismatches({"Fortran": str(real)}, {"Fortran": str(link)}) == []


def test_a_dependency_built_by_another_fortran_compiler_is_stale(tmp_path):
    root = str(tmp_path)
    lapack = FakeTarget(root, "lapack", True)
    fftw = FakeTarget(root, "fftw", True)
    post = FakeTarget(root, "post_process", False, deps=[lapack, fftw])
    write_cache(lapack.get_staging_dirpath(None), fortran=AMDFLANG)
    write_cache(fftw.get_staging_dirpath(None), fortran=GFORTRAN)
    write_cache(post.get_staging_dirpath(None), fortran=GFORTRAN)

    stale = _stale_dependencies(post, None, include_system_found=False)

    assert [dep.name for dep, _ in stale] == ["lapack"]
    assert AMDFLANG in stale[0][1] and GFORTRAN in stale[0][1]


def test_transitive_dependencies_are_checked(tmp_path):
    root = str(tmp_path)
    hdf5 = FakeTarget(root, "hdf5", True)
    silo = FakeTarget(root, "silo", True, deps=[hdf5])
    post = FakeTarget(root, "post_process", False, deps=[silo])
    write_cache(hdf5.get_staging_dirpath(None), fortran=AMDFLANG)
    write_cache(silo.get_staging_dirpath(None), fortran=GFORTRAN)
    write_cache(post.get_staging_dirpath(None), fortran=GFORTRAN)

    assert [dep.name for dep, _ in _stale_dependencies(post, None, include_system_found=False)] == ["hdf5"]


def test_a_dependency_cycle_terminates(tmp_path):
    root = str(tmp_path)
    a = FakeTarget(root, "a", True)
    b = FakeTarget(root, "b", True, deps=[a])
    a.requires = Deps([b])
    post = FakeTarget(root, "post_process", False, deps=[a])
    write_cache(a.get_staging_dirpath(None), fortran=AMDFLANG)
    write_cache(b.get_staging_dirpath(None), fortran=GFORTRAN)
    write_cache(post.get_staging_dirpath(None), fortran=GFORTRAN)

    assert [dep.name for dep, _ in _stale_dependencies(post, None, include_system_found=False)] == ["a"]


def test_system_found_dependencies_are_flagged_only_when_asked(tmp_path):
    # Same compiler, but the superbuild found FFTW on the system and installed nothing.
    root = str(tmp_path)
    fftw = FakeTarget(root, "fftw", True, install_files=0)
    sim = FakeTarget(root, "simulation", False, deps=[fftw])
    write_cache(fftw.get_staging_dirpath(None), fortran=GFORTRAN)
    write_cache(sim.get_staging_dirpath(None), fortran=GFORTRAN)

    assert _stale_dependencies(sim, None, include_system_found=False) == []
    assert [dep.name for dep, _ in _stale_dependencies(sim, None, include_system_found=True)] == ["fftw"]


def test_system_provided_and_uninstalled_dependencies_are_left_alone(tmp_path):
    root = str(tmp_path)
    sys_fftw = FakeTarget(root, "fftw", True, buildable=False, install_files=0)  # --sys-fftw
    lapack = FakeTarget(root, "lapack", True, installed=False)
    post = FakeTarget(root, "post_process", False, deps=[sys_fftw, lapack])
    write_cache(sys_fftw.get_staging_dirpath(None), fortran=AMDFLANG)
    write_cache(lapack.get_staging_dirpath(None), fortran=AMDFLANG)
    write_cache(post.get_staging_dirpath(None), fortran=GFORTRAN)

    assert _stale_dependencies(post, None, include_system_found=True) == []


def test_matching_toolchains_are_not_stale(tmp_path):
    root = str(tmp_path)
    lapack = FakeTarget(root, "lapack", True)
    post = FakeTarget(root, "post_process", False, deps=[lapack])
    write_cache(lapack.get_staging_dirpath(None), fortran=GFORTRAN)
    write_cache(post.get_staging_dirpath(None), fortran=GFORTRAN)

    assert _stale_dependencies(post, None, include_system_found=True) == []
