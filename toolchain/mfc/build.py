import re, os, typing, dataclasses

from .common    import MFCException, system, delete_directory, create_directory
from .state     import ARG, CFG
from .printer   import cons
from .run.input import MFCInputFile


@dataclasses.dataclass
class MFCTarget:
    @dataclasses.dataclass
    class Dependencies:
        all: typing.List
        cpu: typing.List
        gpu: typing.List

        def compute(self) -> typing.Set:
            r  = self.all[:]
            r += self.gpu[:] if ARG("gpu") else self.cpu[:]

            return r

    name:         str              # Name of the target
    flags:        typing.List[str] # Extra flags to pass to CMake
    isDependency: bool             # Is it a dependency of an MFC target?
    isDefault:    bool             # Should it be built by default? (unspecified -t | --targets)
    isRequired:   bool             # Should it always be built? (no matter what -t | --targets is)
    requires:     Dependencies     # Build dependencies of the target
    
    def __hash__(self) -> int:
        return hash(self.name)
    
    # Get path to directory that will store the build files
    def get_build_dirpath(self) -> str:
        return os.sep.join([
            os.getcwd(),
            "build",
            [CFG().make_slug(), 'dependencies'][int(self.isDependency)],
            self.name
        ])
    
    # Get the directory that contains the target's CMakeLists.txt
    def get_cmake_dirpath(self) -> str:
        # The CMakeLists.txt file is located:
        #  * Regular:    <root>/CMakelists.txt
        #  * Dependency: <root>/toolchain/dependencies/CMakelists.txt
        return os.sep.join([
            os.getcwd(),
            os.sep.join(["toolchain", "dependencies"]) if self.isDependency else "",
        ])

    def get_install_dirpath(self) -> str:
        # The install directory is located:
        # Regular:    <root>/build/install/<configuration_slug>
        # Dependency: <root>/build/install/dependencies (shared)
        return os.sep.join([
            os.getcwd(),
            "build",
            "install",
            'dependencies' if self.isDependency else CFG().make_slug()
        ])
    
    def get_install_binpath(self) -> str:
        # <root>/install/<slug>/bin/<target>
        return os.sep.join([self.get_install_dirpath(), "bin", self.name])

    def is_configured(self) -> bool:
        # We assume that if the CMakeCache.txt file exists, then the target is
        # configured. (this isn't perfect, but it's good enough for now)
        return os.path.isfile(
            os.sep.join([self.get_build_dirpath(), "CMakeCache.txt"])
        )

    def get_configuration_txt(self) -> typing.Optional[dict]:
        if not self.is_configured():
            return None

        configpath = os.path.join(self.get_build_dirpath(), "configuration.txt")
        if not os.path.exists(configpath):
            return None

        with open(configpath) as f:
            return f.read()

        return None


    def build(self, history: typing.Set[str] = None):
        if history is None:
            history = set()

        if self.name in history:
            return

        history.add(self.name)
        
        build_targets(REQUIRED_TARGETS, history)
        
        cons.print(f"[bold]Building [magenta]{self.name}[/magenta]:[/bold]")
        cons.indent()

        if ARG("no_build"):
            cons.print("--no-build specified, skipping...")
            cons.unindent()
            return

        if self.isDependency and ARG(f"no_{self.name}"):
            cons.print(f"--no-{self.name} given, skipping...")
            cons.unindent()
            return

        build_dirpath   = self.get_build_dirpath()
        cmake_dirpath   = self.get_cmake_dirpath()
        install_dirpath = self.get_install_dirpath()
        
        install_prefixes = ';'.join([install_dirpath, get_dependency_install_dirpath()])

        flags: list = self.flags.copy() + [
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

        if not self.isDependency:
            flags.append(f"-DMFC_MPI={    'ON' if ARG('mpi') else 'OFF'}")
            flags.append(f"-DMFC_OpenACC={'ON' if ARG('gpu') else 'OFF'}")

        configure = ["cmake"] + flags + ["-S", cmake_dirpath, "-B", build_dirpath]
        build     = ["cmake", "--build",  build_dirpath,
                            "--target", self.name,
                            "-j",       ARG("jobs"),
                            "--config", 'Debug' if ARG('debug') else 'Release']
        if ARG('verbose'):
            build.append("--verbose")

        install   = ["cmake", "--install", build_dirpath]

        if not self.is_configured():
            build_targets(self.requires.compute(), history)

            delete_directory(build_dirpath)
            create_directory(build_dirpath)

            if system(configure, no_exception=True) != 0:
                raise MFCException(f"Failed to configure the [bold magenta]{self.name}[/bold magenta] target.")

        if not self.isDependency and ARG("command") == "build":
            MFCInputFile("", "", {}).generate(self, bOnlyFPPs = True)

        system(build,   exception_text=f"Failed to build the [bold magenta]{self.name}[/bold magenta] target.")
        system(install, exception_text=f"Failed to install the [bold magenta]{self.name}[/bold magenta] target.")
        
        cons.print(no_indent=True)
        cons.unindent()

    def clean(self):
        cons.print(f"[bold]Cleaning [magenta]{self.name}[/magenta]:[/bold]")
        cons.indent()

        build_dirpath = self.get_build_dirpath()

        if not os.path.isdir(build_dirpath):
            cons.print("Target not configured. Nothing to clean.")
            cons.unindent()
            return

        clean = ["cmake", "--build",  build_dirpath, "--target", "clean",
                          "--config", "Debug" if ARG("debug") else "Release" ]

        if ARG("verbose"):
            clean.append("--verbose")

        system(clean, exception_text=f"Failed to clean the [bold magenta]{self.name}[/bold magenta] target.")

        cons.unindent()


FFTW          = MFCTarget('fftw',          ['-DMFC_FFTW=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []))
HDF5          = MFCTarget('hdf5',          ['-DMFC_HDF5=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []))
SILO          = MFCTarget('silo',          ['-DMFC_SILO=ON'],          True,  False, False, MFCTarget.Dependencies([HDF5], [], []))
PRE_PROCESS   = MFCTarget('pre_process',   ['-DMFC_PRE_PROCESS=ON'],   False, True,  False, MFCTarget.Dependencies([], [], []))
SIMULATION    = MFCTarget('simulation',    ['-DMFC_SIMULATION=ON'],    False, True,  False, MFCTarget.Dependencies([], [FFTW], []))
POST_PROCESS  = MFCTarget('post_process',  ['-DMFC_POST_PROCESS=ON'],  False, True,  False, MFCTarget.Dependencies([FFTW, SILO], [], []))
SYSCHECK      = MFCTarget('syscheck',      ['-DMFC_SYSCHECK=ON'],      False, False, True,  MFCTarget.Dependencies([], [], []))
DOCUMENTATION = MFCTarget('documentation', ['-DMFC_DOCUMENTATION=ON'], False, False, False, MFCTarget.Dependencies([], [], []))

TARGETS = { FFTW, HDF5, SILO, PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK, DOCUMENTATION }

DEFAULT_TARGETS    = { target for target in TARGETS if target.isDefault }
REQUIRED_TARGETS   = { target for target in TARGETS if target.isRequired }
DEPENDENCY_TARGETS = { target for target in TARGETS if target.isDependency }

TARGET_MAP = { target.name: target for target in TARGETS }

def get_target(target: typing.Union[str, MFCTarget]) -> MFCTarget:
    if isinstance(target, MFCTarget):
        return target
    
    if target in TARGET_MAP:
        return TARGET_MAP[target]

    raise MFCException(f"Target '{target}' does not exist.")


def get_dependency_install_dirpath() -> str:
    # Since dependencies share the same install directory, we can just return
    # the install directory of the first dependency we find.
    for target in TARGETS:
        if target.isDependency:
            return target.get_install_dirpath()

    raise MFCException("No dependency target found.")


def build_targets(targets: typing.Iterable[typing.Union[MFCTarget, str]], history: typing.Set[str] = None):
    if history is None:
        history = set()
    
    for target in targets:
        get_target(target).build(history)


def clean_targets(targets: typing.Iterable[typing.Union[MFCTarget, str]]):
    for target in targets:
        get_target(target).clean()


def get_configured_targets() -> typing.List[MFCTarget]:
    return [ target for target in TARGETS if target.is_configured() ]


def build():
    build_targets(ARG("targets"))


def clean():
    clean_targets(ARG("targets"))

