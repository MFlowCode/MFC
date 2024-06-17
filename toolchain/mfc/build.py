import os, typing, hashlib, dataclasses

from .printer import cons
from .common  import MFCException, system, delete_directory, create_directory, \
                     format_list_to_string
from .state   import ARG, CFG
from .run     import input

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
    runOrder:     int              # For MFC Targets: Order in which targets should logically run

    def __hash__(self) -> int:
        return hash(self.name)

    def get_slug(self, case: input.MFCInputFile) -> str:
        if self.isDependency:
            return self.name

        m = hashlib.sha256()
        m.update(self.name.encode())
        m.update(CFG().make_slug().encode())
        m.update(case.get_fpp(self, False).encode())

        return m.hexdigest()[:10]

    # Get path to directory that will store the build files
    def get_staging_dirpath(self, case: input.MFCInputFile) -> str:
        return os.sep.join([os.getcwd(), "build", "staging", self.get_slug(case) ])

    # Get the directory that contains the target's CMakeLists.txt
    def get_cmake_dirpath(self) -> str:
        # The CMakeLists.txt file is located:
        #  * Regular:    <root>/CMakelists.txt
        #  * Dependency: <root>/toolchain/dependencies/CMakelists.txt
        return os.sep.join([
            os.getcwd(),
            os.sep.join(["toolchain", "dependencies"]) if self.isDependency else "",
        ])

    def get_install_dirpath(self, case: input.MFCInputFile) -> str:
        # The install directory is located:
        # Regular:    <root>/build/install/<slug>
        # Dependency: <root>/build/install/dependencies (shared)
        return os.sep.join([
            os.getcwd(),
            "build",
            "install",
            'dependencies' if self.isDependency else self.get_slug(case),
        ])

    def get_install_binpath(self, case: input.MFCInputFile) -> str:
        # <root>/install/<slug>/bin/<target>
        return os.sep.join([self.get_install_dirpath(case), "bin", self.name])

    def is_configured(self, case: input.MFCInputFile) -> bool:
        # We assume that if the CMakeCache.txt file exists, then the target is
        # configured. (this isn't perfect, but it's good enough for now)
        return os.path.isfile(
            os.sep.join([self.get_staging_dirpath(case), "CMakeCache.txt"])
        )

    def get_configuration_txt(self, case: input.MFCInputFile) -> typing.Optional[dict]:
        if not self.is_configured(case):
            return None

        configpath = os.path.join(self.get_staging_dirpath(case), "configuration.txt")
        if not os.path.exists(configpath):
            return None

        with open(configpath) as f:
            return f.read()

    def is_buildable(self) -> bool:
        if ARG("no_build"):
            return False

        if self.isDependency and ARG(f"sys_{self.name}", False):
            return False

        return True

    def configure(self, case: input.MFCInputFile):
        build_dirpath   = self.get_staging_dirpath(case)
        cmake_dirpath   = self.get_cmake_dirpath()
        install_dirpath = self.get_install_dirpath(case)

        install_prefixes = ';'.join([install_dirpath, get_dependency_install_dirpath(case)])

        mod_dirs = ';'.join(['build/install/dependencies/include/hipfort/amdgcn'])

        flags: list = self.flags.copy() + [
            # Disable CMake warnings intended for developers (us).
            # See: https://cmake.org/cmake/help/latest/manual/cmake.1.html.
            "-Wno-dev",
            # Disable warnings about unused command line arguments. This is
            # useful for passing arguments to CMake that are not used by the
            # current target.
            "--no-warn-unused-cli",
            # Save a compile_commands.json file with the compile commands used to
            # build the configured targets. This is mostly useful for debugging.
            # See: https://cmake.org/cmake/help/latest/variable/CMAKE_EXPORT_COMPILE_COMMANDS.html.
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            # Set build type (e.g Debug, Release, etc.).
            # See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
            f"-DCMAKE_BUILD_TYPE={'Debug' if ARG('debug') else 'Release'}",
            # Used by FIND_PACKAGE (/FindXXX) to search for packages, with the
            # second highest level of priority, still letting users manually
            # specify <PackageName>_ROOT, which has precedence over CMAKE_PREFIX_PATH.
            # See: https://cmake.org/cmake/help/latest/command/find_package.html.
            f"-DCMAKE_PREFIX_PATH={install_prefixes}",
            # First directory that FIND_LIBRARY searches.
            # See: https://cmake.org/cmake/help/latest/command/find_library.html.
            f"-DCMAKE_FIND_ROOT_PATH={install_prefixes}",
            # Location prefix to install bin/, lib/, include/, etc.
            # See: https://cmake.org/cmake/help/latest/command/install.html.
            f"-DCMAKE_INSTALL_PREFIX={install_dirpath}",
            # Fortran .mod include directories. Currently used for the HIPFORT
            # dependency that has this missing from its config files.
            f"-DCMAKE_Fortran_MODULE_DIRECTORY={mod_dirs}",
        ]

        if ARG("verbose"):
            flags.append('--debug-find')

        if not self.isDependency:
            flags.append(f"-DMFC_MPI={    'ON' if ARG('mpi') else 'OFF'}")
            flags.append(f"-DMFC_OpenACC={'ON' if ARG('gpu') else 'OFF'}")
            flags.append(f"-DMFC_GCov={   'ON' if ARG('gcov') else 'OFF'}")
            flags.append(f"-DMFC_Unified={'ON' if ARG('unified') else 'OFF'}")

        command = ["cmake"] + flags + ["-S", cmake_dirpath, "-B", build_dirpath]

        delete_directory(build_dirpath)
        create_directory(build_dirpath)

        case.generate_fpp(self)

        if system(command).returncode != 0:
            raise MFCException(f"Failed to configure the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(no_indent=True)

    def build(self, case: input.MFCInputFile):
        case.generate_fpp(self)

        command = ["cmake", "--build",    self.get_staging_dirpath(case),
                            "--target",   self.name,
                            "--parallel", ARG("jobs"),
                            "--config",   'Debug' if ARG('debug') else 'Release']
        if ARG('verbose'):
            command.append("--verbose")

        if system(command).returncode != 0:
            raise MFCException(f"Failed to build the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(no_indent=True)

    def install(self, case: input.MFCInputFile):
        command = ["cmake", "--install", self.get_staging_dirpath(case)]

        if system(command).returncode != 0:
            raise MFCException(f"Failed to install the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(no_indent=True)

    def clean(self, case: input.MFCInputFile):
        build_dirpath = self.get_staging_dirpath(case)

        if not os.path.isdir(build_dirpath):
            return

        command = ["cmake", "--build",  build_dirpath, "--target", "clean",
                            "--config", "Debug" if ARG("debug") else "Release" ]

        if ARG("verbose"):
            command.append("--verbose")

        if system(command).returncode != 0:
            raise MFCException(f"Failed to clean the [bold magenta]{self.name}[/bold magenta] target.")

FFTW          = MFCTarget('fftw',          ['-DMFC_FFTW=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []), -1)
HDF5          = MFCTarget('hdf5',          ['-DMFC_HDF5=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []), -1)
SILO          = MFCTarget('silo',          ['-DMFC_SILO=ON'],          True,  False, False, MFCTarget.Dependencies([HDF5], [], []), -1)
HIPFORT       = MFCTarget('hipfort',       ['-DMFC_HIPFORT=ON'],       True,  False, False, MFCTarget.Dependencies([], [], []), -1)
PRE_PROCESS   = MFCTarget('pre_process',   ['-DMFC_PRE_PROCESS=ON'],   False, True,  False, MFCTarget.Dependencies([], [], []), 0)
SIMULATION    = MFCTarget('simulation',    ['-DMFC_SIMULATION=ON'],    False, True,  False, MFCTarget.Dependencies([], [FFTW], [HIPFORT]), 1)
POST_PROCESS  = MFCTarget('post_process',  ['-DMFC_POST_PROCESS=ON'],  False, True,  False, MFCTarget.Dependencies([FFTW, SILO], [], []), 2)
SYSCHECK      = MFCTarget('syscheck',      ['-DMFC_SYSCHECK=ON'],      False, False, True,  MFCTarget.Dependencies([], [], [HIPFORT]), -1)
DOCUMENTATION = MFCTarget('documentation', ['-DMFC_DOCUMENTATION=ON'], False, False, False, MFCTarget.Dependencies([], [], []), -1)

TARGETS = { FFTW, HDF5, SILO, HIPFORT, PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK, DOCUMENTATION }

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


def get_targets(targets: typing.List[typing.Union[str, MFCTarget]]) -> typing.List[MFCTarget]:
    return [ get_target(t) for t in targets ]


def get_dependency_install_dirpath(case: input.MFCInputFile) -> str:
    # Since dependencies share the same install directory, we can just return
    # the install directory of the first dependency we find.
    for target in TARGETS:
        if target.isDependency:
            return target.get_install_dirpath(case)

    raise MFCException("No dependency target found.")


def __build_target(target: typing.Union[MFCTarget, str], case: input.MFCInputFile, history: typing.Set[str] = None):
    if history is None:
        history = set()

    target = get_target(target)

    if target.name in history or not target.is_buildable():
        return

    history.add(target.name)

    for dep in target.requires.compute():
        # If we have already built and installed this target,
        # do not do so again. This can be inferred by whether
        # the target requesting this dependency is already configured.
        if dep.isDependency and target.is_configured(case):
            continue

        build([dep], case, history)

    if not target.is_configured(case):
        target.configure(case)

    target.build(case)
    target.install(case)


def get_configured_targets(case: input.MFCInputFile) -> typing.List[MFCTarget]:
    return [ target for target in TARGETS if target.is_configured(case) ]


def __generate_header(step_name: str, targets: typing.List):
    caseopt_info = "Generic Build"
    if ARG("case_optimization"):
        caseopt_info = f"Case Optimized for [magenta]{ARG('input')}[/magenta]"

    targets     = get_targets(targets)
    target_list = format_list_to_string([ t.name for t in targets ], 'magenta')

    return f"[bold]{step_name} | {target_list} | {caseopt_info}[/bold]"


def build(targets = None, case: input.MFCInputFile = None, history: typing.Set[str] = None):
    if history is None:
        history = set()
    if targets is None:
        targets = ARG("targets")

    targets = get_targets(list(REQUIRED_TARGETS) + targets)
    case    = case or input.load(ARG("input"), ARG("--"), {})
    case.validate_params()

    if len(history) == 0:
        cons.print(__generate_header("Build", targets))
        cons.print(no_indent=True)

    for target in targets:
        __build_target(target, case, history)

    if len(history) == 0:
        cons.print(no_indent=True)


def clean(targets = None, case: input.MFCInputFile = None):
    targets = get_targets(list(REQUIRED_TARGETS) + (targets or ARG("targets")))
    case    = case or input.load(ARG("input"), ARG("arguments"), {})

    cons.print(__generate_header("Clean", targets))
    cons.print(no_indent=True)

    for target in targets:
        if target.is_configured():
            target.clean()

    cons.print(no_indent=True)
    cons.unindent()
