import os, typing, dataclasses

from .state     import ARG, CFG
from .printer   import cons
from .          import common
from .run.input import MFCInputFile


@dataclasses.dataclass
class MFCTarget:
    @dataclasses.dataclass
    class Dependencies:
        all: typing.List[str]
        cpu: typing.List[str]
        gpu: typing.List[str]

        def compute(self) -> typing.List[str]:
            r  = self.all[:]
            r += self.gpu[:] if ARG("gpu") else self.cpu[:]

            return r

    name:         str              # Name of the target
    flags:        typing.List[str] # Extra flags to pass to CMake
    isDependency: bool             # Is it a dependency of an MFC target?
    isDefault:    bool             # Should it be built by default? (unspecified -t | --targets)
    isRequired:   bool             # Should it always be built? (no matter what -t | --targets is)
    requires:     Dependencies     # Build dependencies of the target


FFTW          = MFCTarget('fftw',          ['-DMFC_FFTW=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []))
HDF5          = MFCTarget('hdf5',          ['-DMFC_HDF5=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []))
SILO          = MFCTarget('silo',          ['-DMFC_SILO=ON'],          True,  False, False, MFCTarget.Dependencies(["hdf5"], [], []))
PRE_PROCESS   = MFCTarget('pre_process',   ['-DMFC_PRE_PROCESS=ON'],   False, True,  False, MFCTarget.Dependencies([], [], []))
SIMULATION    = MFCTarget('simulation',    ['-DMFC_SIMULATION=ON'],    False, True,  False, MFCTarget.Dependencies([], ["fftw"], []))
POST_PROCESS  = MFCTarget('post_process',  ['-DMFC_POST_PROCESS=ON'],  False, True,  False, MFCTarget.Dependencies(['fftw', 'silo'], [], []))
SYSCHECK      = MFCTarget('syscheck',      ['-DMFC_SYSCHECK=ON'],      False, False, True,  MFCTarget.Dependencies([], [], []))
DOCUMENTATION = MFCTarget('documentation', ['-DMFC_DOCUMENTATION=ON'], False, False, False, MFCTarget.Dependencies([], [], []))

TARGETS: typing.List[MFCTarget] = [ FFTW, HDF5, SILO, PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK, DOCUMENTATION ]

def get_mfc_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isDefault ]


def get_dependencies_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isDependency ]


def get_required_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isRequired ]


def get_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS ]


def get_target(name: str) -> MFCTarget:
    for target in TARGETS:
        if target.name == name:
            return target

    raise common.MFCException(f"Target '{name}' does not exist.")


# Get path to directory that will store the build files
def get_build_dirpath(target: MFCTarget) -> str:
    return os.sep.join([
        os.getcwd(),
        "build",
        [CFG().make_slug(), 'dependencies'][int(target.isDependency)],
        target.name
    ])


# Get the directory that contains the target's CMakeLists.txt
def get_cmake_dirpath(target: MFCTarget) -> str:
    # The CMakeLists.txt file is located:
    #  * Regular:    <root>/CMakelists.txt
    #  * Dependency: <root>/toolchain/dependencies/CMakelists.txt
    return os.sep.join([
        os.getcwd(),
        os.sep.join(["toolchain", "dependencies"]) if target.isDependency else "",
    ])


def get_install_dirpath(target: MFCTarget) -> str:
    # The install directory is located:
    # Regular:    <root>/build/install/<configuration_slug>
    # Dependency: <root>/build/install/dependencies (shared)
    return os.sep.join([
        os.getcwd(),
        "build",
        "install",
        'dependencies' if target.isDependency else CFG().make_slug()
    ])


def get_dependency_install_dirpath() -> str:
    # Since dependencies share the same install directory, we can just return
    # the install directory of the first dependency we find.
    for target in TARGETS:
        if target.isDependency:
            return get_install_dirpath(target)

    raise common.MFCException("No dependency target found.")


def is_target_configured(target: MFCTarget) -> bool:
    # We assume that if the CMakeCache.txt file exists, then the target is
    # configured. (this isn't perfect, but it's good enough for now)
    return os.path.isfile(
        os.sep.join([get_build_dirpath(target), "CMakeCache.txt"])
    )


def clean_target(name: str):
    cons.print(f"[bold]Cleaning [magenta]{name}[/magenta]:[/bold]")
    cons.indent()

    target = get_target(name)

    build_dirpath = get_build_dirpath(target)

    if not os.path.isdir(build_dirpath):
        cons.print("Target not configured. Nothing to clean.")
        cons.unindent()
        return

    clean = ["cmake", "--build",  build_dirpath, "--target", "clean",
                      "--config", "Debug" if ARG("debug") else "Release" ]

    if ARG("verbose"):
        clean.append("--verbose")

    common.system(clean, exception_text=f"Failed to clean the [bold magenta]{name}[/bold magenta] target.")

    cons.unindent()


def clean_targets(targets: typing.List[str]):
    for target in targets:
        clean_target(target)


def clean():
    clean_targets(ARG("targets"))


def build_target(name: str, history: typing.List[str] = None):
    cons.print(f"[bold]Building [magenta]{name}[/magenta]:[/bold]")
    cons.indent()

    if history is None:
        history = []

    if ARG("no_build"):
        cons.print("--no-build specified, skipping...")
        cons.unindent()
        return

    if name in history:
        cons.print("Already built, skipping...")
        cons.unindent()
        return

    history.append(name)

    target = get_target(name)

    if target.isDependency and ARG(f"no_{target.name}"):
        cons.print(f"--no-{target.name} given, skipping...")
        cons.unindent()
        return

    build_dirpath   = get_build_dirpath(target)
    cmake_dirpath   = get_cmake_dirpath(target)
    install_dirpath = get_install_dirpath(target)
    
    install_prefixes = ';'.join([install_dirpath, get_dependency_install_dirpath()])

    flags: list = target.flags.copy() + [
        # Disable CMake warnings intended for developers (us).
        # See: https://cmake.org/cmake/help/latest/manual/cmake.1.html.
        f"-Wno-dev",
        # Save a compile_commands.json file with the compile commands used to
        # build the configured targets. This is mostly useful for debugging.
        # See: https://cmake.org/cmake/help/latest/variable/CMAKE_EXPORT_COMPILE_COMMANDS.html.
        f"-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        # Set build type (e.g Debug, Release, etc.).
        # See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
        f"-DCMAKE_BUILD_TYPE={'Debug' if ARG('debug') else 'Release'}",
        # Used by FIND_PACKAGE (/FindXXX) to search for packages, with the
        # second heighest level of priority, still letting users manually
        # specify <PackageName>_ROOT, which has precedence over CMAKE_PREFIX_PATH.
        # See: https://cmake.org/cmake/help/latest/command/find_package.html.
        f"-DCMAKE_PREFIX_PATH={install_prefixes}",
        # First directory that FIND_LIBRARY searches.
        # See: https://cmake.org/cmake/help/latest/command/find_library.html.
        f"-DCMAKE_FIND_ROOT_PATH={install_prefixes}",
        # Location prefix to install bin/, lib/, include/, etc.
        # See: https://cmake.org/cmake/help/latest/command/install.html.
        f"-DCMAKE_INSTALL_PREFIX={install_dirpath}",
    ]

    if not target.isDependency:
        flags.append(f"-DMFC_MPI={    'ON' if ARG('mpi') else 'OFF'}")
        flags.append(f"-DMFC_OpenACC={'ON' if ARG('gpu') else 'OFF'}")

    configure = ["cmake"] + flags + ["-S", cmake_dirpath, "-B", build_dirpath]
    build     = ["cmake", "--build",  build_dirpath,
                          "--target", name,
                          "-j",       ARG("jobs"),
                          "--config", 'Debug' if ARG('debug') else 'Release']
    if ARG('verbose'):
        build.append("--verbose")

    install   = ["cmake", "--install", build_dirpath]

    if not is_target_configured(target):
        for dependency_name in target.requires.compute():
            build_target(dependency_name, history)

        common.delete_directory(build_dirpath)
        common.create_directory(build_dirpath)

        if common.system(configure, no_exception=True) != 0:
            raise common.MFCException(f"Failed to configure the [bold magenta]{name}[/bold magenta] target.")

    if not target.isDependency and ARG("command") == "build":
        MFCInputFile("", "", {}).generate(name, bOnlyFPPs = True)

    common.system(build,   exception_text=f"Failed to build the [bold magenta]{name}[/bold magenta] target.")
    common.system(install, exception_text=f"Failed to install the [bold magenta]{name}[/bold magenta] target.")
    
    cons.print(no_indent=True)
    cons.unindent()


def build_targets(targets):
    for target in set(targets).union(set(get_required_target_names())):
        build_target(target)


def build():
    build_targets(ARG("targets"))
