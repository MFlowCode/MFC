import os, typing, hashlib, dataclasses, subprocess, re, time, sys, threading, queue

from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn, TaskProgressColumn

from .case    import Case
from .printer import cons
from .common  import MFCException, system, delete_directory, create_directory, \
                     format_list_to_string, debug
from .state   import ARG, CFG
from .run     import input
from .state   import gpuConfigOptions
from .user_guide import Tips


# Regex to parse build progress
# Ninja format: [42/156] Building Fortran object ...
_NINJA_PROGRESS_RE = re.compile(r'^\[(\d+)/(\d+)\]\s+(.*)$')
# Make format: [ 16%] Building Fortran object ... or [100%] Linking ...
_MAKE_PROGRESS_RE = re.compile(r'^\[\s*(\d+)%\]\s+(.*)$')


# pylint: disable=too-many-locals,too-many-branches,too-many-statements,too-many-nested-blocks
def _run_build_with_progress(command: typing.List[str], target_name: str, streaming: bool = False) -> subprocess.CompletedProcess:
    """
    Run a build command with a progress bar that parses ninja output.

    Args:
        command: The cmake build command to run
        target_name: Name of the target being built
        streaming: If True, print [X/Y] lines as they happen instead of progress bar (-v mode)

    Shows:
    - Progress bar with file count (e.g., 42/156)
    - Current file being compiled
    - Elapsed time

    Falls back to spinner with elapsed time if ninja progress can't be parsed.
    """
    cmd = [str(x) for x in command]

    # Check if we're in a TTY (interactive terminal)
    is_tty = sys.stdout.isatty()

    # Collect all output for error reporting
    all_stdout = []
    all_stderr = []

    # Start the process (can't use 'with' since process is used in multiple branches)
    process = subprocess.Popen(  # pylint: disable=consider-using-with
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1  # Line buffered
    )

    if streaming:
        # Streaming mode (-v): print build progress lines as they happen
        cons.print(f"  [bold blue]Building[/bold blue] [magenta]{target_name}[/magenta] [dim](-v)[/dim]...")
        start_time = time.time()

        # Read stdout and print matching lines
        for line in iter(process.stdout.readline, ''):
            all_stdout.append(line)
            stripped = line.strip()

            # Try ninja format first: [42/156] Action
            ninja_match = _NINJA_PROGRESS_RE.match(stripped)
            if ninja_match:
                completed = ninja_match.group(1)
                total = ninja_match.group(2)
                action = ninja_match.group(3)
                # Extract filename from action
                parts = action.split()
                if len(parts) >= 3:
                    filename = os.path.basename(parts[-1]).replace('.o', '').replace('.obj', '')
                    if len(filename) > 40:
                        filename = filename[:37] + "..."
                    cons.print(f"  [dim][{completed}/{total}][/dim] {filename}")
                continue

            # Try make format: [ 16%] Action
            make_match = _MAKE_PROGRESS_RE.match(stripped)
            if make_match:
                percent = make_match.group(1)
                action = make_match.group(2)
                # Extract filename from action
                parts = action.split()
                if len(parts) >= 3:
                    # Get the last part which is usually the file path
                    obj_path = parts[-1]
                    filename = os.path.basename(obj_path).replace('.o', '').replace('.obj', '')
                    if len(filename) > 40:
                        filename = filename[:37] + "..."
                    cons.print(f"  [dim][{percent:>3}%][/dim] {filename}")

        # Read any remaining stderr
        stderr = process.stderr.read()
        all_stderr.append(stderr)
        process.wait()

        elapsed = time.time() - start_time
        if elapsed > 5:
            cons.print(f"  [dim](build took {elapsed:.1f}s)[/dim]")

        return subprocess.CompletedProcess(cmd, process.returncode, ''.join(all_stdout), ''.join(all_stderr))

    if not is_tty:
        # Non-interactive, non-streaming: show message with elapsed time
        cons.print(f"  [bold blue]Building[/bold blue] [magenta]{target_name}[/magenta]...")
        start_time = time.time()
        stdout, stderr = process.communicate()
        elapsed = time.time() - start_time
        if elapsed > 5:  # Only show time for longer builds
            cons.print(f"  [dim](build took {elapsed:.1f}s)[/dim]")
        return subprocess.CompletedProcess(cmd, process.returncode, stdout, stderr)

    # Interactive: show progress bar
    current_file = ""
    total_files = 0
    completed_files = 0
    progress_detected = False

    # Create a custom progress display
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]Building[/bold blue] [magenta]{task.fields[target]}[/magenta]"),
        BarColumn(bar_width=30),
        TaskProgressColumn(),
        TextColumn("•"),
        TimeElapsedColumn(),
        TextColumn("[dim]{task.fields[current_file]}[/dim]"),
        console=cons.raw,
        transient=True,  # Remove progress bar when done
        refresh_per_second=4,
    ) as progress:
        # Start with indeterminate progress (total=None shows spinner behavior)
        task = progress.add_task(
            "build",
            total=None,
            target=target_name,
            current_file=""
        )

        # Use threads to read stdout and stderr concurrently
        stdout_queue = queue.Queue()
        stderr_queue = queue.Queue()

        def read_stdout():
            for line in iter(process.stdout.readline, ''):
                stdout_queue.put(line)
            stdout_queue.put(None)  # Signal EOF

        def read_stderr():
            for line in iter(process.stderr.readline, ''):
                stderr_queue.put(line)
            stderr_queue.put(None)  # Signal EOF

        stdout_thread = threading.Thread(target=read_stdout, daemon=True)
        stderr_thread = threading.Thread(target=read_stderr, daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        stdout_done = False
        stderr_done = False

        while not (stdout_done and stderr_done):
            # Check stdout
            try:
                line = stdout_queue.get_nowait()
                if line is None:
                    stdout_done = True
                else:
                    all_stdout.append(line)
                    stripped = line.strip()

                    # Try ninja format first: [42/156] Action
                    ninja_match = _NINJA_PROGRESS_RE.match(stripped)
                    if ninja_match:
                        completed_files = int(ninja_match.group(1))
                        total_files = int(ninja_match.group(2))
                        action = ninja_match.group(3)

                        # Extract just the filename from the action
                        if action:
                            parts = action.split()
                            if len(parts) >= 3:
                                obj_path = parts[-1]
                                current_file = os.path.basename(obj_path).replace('.o', '').replace('.obj', '')
                                if len(current_file) > 30:
                                    current_file = current_file[:27] + "..."

                        if not progress_detected:
                            progress_detected = True
                            progress.update(task, total=total_files)

                        progress.update(
                            task,
                            completed=completed_files,
                            current_file=current_file
                        )
                    else:
                        # Try make format: [ 16%] Action
                        make_match = _MAKE_PROGRESS_RE.match(stripped)
                        if make_match:
                            percent = int(make_match.group(1))
                            action = make_match.group(2)

                            # Extract filename from action
                            if action:
                                parts = action.split()
                                if len(parts) >= 3:
                                    obj_path = parts[-1]
                                    current_file = os.path.basename(obj_path).replace('.o', '').replace('.obj', '')
                                    if len(current_file) > 30:
                                        current_file = current_file[:27] + "..."

                            if not progress_detected:
                                progress_detected = True
                                # Make uses percentage, so set total to 100
                                progress.update(task, total=100)

                            progress.update(
                                task,
                                completed=percent,
                                current_file=current_file
                            )
            except queue.Empty:
                pass

            # Check stderr
            try:
                line = stderr_queue.get_nowait()
                if line is None:
                    stderr_done = True
                else:
                    all_stderr.append(line)
            except queue.Empty:
                pass

            # Small sleep to avoid busy waiting
            if not stdout_done or not stderr_done:
                time.sleep(0.01)

        # Wait for process to complete
        process.wait()

        # Ensure threads are done
        stdout_thread.join(timeout=1)
        stderr_thread.join(timeout=1)

    return subprocess.CompletedProcess(
        cmd,
        process.returncode,
        ''.join(all_stdout),
        ''.join(all_stderr)
    )


def _show_build_error(result: subprocess.CompletedProcess, stage: str):
    """Display build error details from captured subprocess output."""
    cons.print()
    cons.print(f"[bold red]{stage} Failed - Error Details:[/bold red]")

    # Show stdout if available (often contains the actual error for CMake)
    if result.stdout:
        stdout_text = result.stdout if isinstance(result.stdout, str) else result.stdout.decode('utf-8', errors='replace')
        stdout_lines = stdout_text.strip().split('\n')
        # Show last 40 lines to capture the relevant error
        if len(stdout_lines) > 40:
            stdout_lines = ['... (truncated) ...'] + stdout_lines[-40:]
        if stdout_lines and stdout_lines != ['']:
            cons.raw.print(Panel('\n'.join(stdout_lines), title="Output", border_style="yellow"))

    # Show stderr if available
    if result.stderr:
        stderr_text = result.stderr if isinstance(result.stderr, str) else result.stderr.decode('utf-8', errors='replace')
        stderr_lines = stderr_text.strip().split('\n')
        if len(stderr_lines) > 40:
            stderr_lines = ['... (truncated) ...'] + stderr_lines[-40:]
        if stderr_lines and stderr_lines != ['']:
            cons.raw.print(Panel('\n'.join(stderr_lines), title="Errors", border_style="red"))

    cons.print()

@dataclasses.dataclass
class MFCTarget:
    @dataclasses.dataclass
    class Dependencies:
        all: typing.List
        cpu: typing.List
        gpu: typing.List

        def compute(self) -> typing.Set:
            r  = self.all[:]
            r += self.gpu[:] if (ARG("gpu") != gpuConfigOptions.NONE.value) else self.cpu[:]

            return r

    name:         str              # Name of the target
    flags:        typing.List[str] # Extra flags to pass to CMakeMFCTarget
    isDependency: bool             # Is it a dependency of an MFC target?
    isDefault:    bool             # Should it be built by default? (unspecified -t | --targets)
    isRequired:   bool             # Should it always be built? (no matter what -t | --targets is)
    requires:     Dependencies     # Build dependencies of the target
    runOrder:     int              # For MFC Targets: Order in which targets should logically run

    def __hash__(self) -> int:
        return hash(self.name)

    def get_slug(self, case: Case ) -> str:
        if self.isDependency:
            return self.name

        m = hashlib.sha256()
        m.update(self.name.encode())
        m.update(CFG().make_slug().encode())
        m.update(case.get_fpp(self, False).encode())

        if case.params.get('chemistry', 'F') == 'T':
            m.update(case.get_cantera_solution().name.encode())

        return m.hexdigest()[:10]

    # Get path to directory that will store the build files
    def get_staging_dirpath(self, case: Case ) -> str:
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

    def get_install_dirpath(self, case: Case ) -> str:
        # The install directory is located <root>/build/install/<slug>
        return os.sep.join([os.getcwd(), "build", "install", self.get_slug(case)])

    def get_home_dirpath(self) -> str:
        return os.sep.join([os.getcwd()])

    def get_install_binpath(self, case: Case ) -> str:
        # <root>/install/<slug>/bin/<target>
        return os.sep.join([self.get_install_dirpath(case), "bin", self.name])

    def is_configured(self, case: Case ) -> bool:
        # We assume that if the CMakeCache.txt file exists, then the target is
        # configured. (this isn't perfect, but it's good enough for now)
        return os.path.isfile(
            os.sep.join([self.get_staging_dirpath(case), "CMakeCache.txt"])
        )

    def get_configuration_txt(self, case: Case ) -> typing.Optional[dict]:
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

    def configure(self, case: Case):
        build_dirpath   = self.get_staging_dirpath(case)
        cmake_dirpath   = self.get_cmake_dirpath()
        install_dirpath = self.get_install_dirpath(case)

        install_prefixes = ';'.join([
            t.get_install_dirpath(case) for t in self.requires.compute()
        ])

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
            # First directory that FIND_PACKAGE searches.
            # See: https://cmake.org/cmake/help/latest/variable/CMAKE_FIND_PACKAGE_REDIRECTS_DIR.html.
            f"-DCMAKE_FIND_PACKAGE_REDIRECTS_DIR={install_prefixes}",
            # Location prefix to install bin/, lib/, include/, etc.
            # See: https://cmake.org/cmake/help/latest/command/install.html.
            f"-DCMAKE_INSTALL_PREFIX={install_dirpath}",
            f"-DMFC_SINGLE_PRECISION={'ON' if (ARG('single') or ARG('mixed')) else 'OFF'}",
            f"-DMFC_MIXED_PRECISION={'ON' if ARG('mixed') else 'OFF'}"
        ]

        # Verbosity level 3 (-vvv): add cmake debug flags
        if ARG("verbose") >= 3:
            flags.append('--debug-find')

        if not self.isDependency:
            flags.append(f"-DMFC_MPI={    'ON' if ARG('mpi') else 'OFF'}")
            # flags.append(f"-DMFC_OpenACC={'ON' if ARG('acc') else 'OFF'}")
            # flags.append(f"-DMFC_OpenMP={'ON' if ARG('mp') else 'OFF'}")
            flags.append(f"-DMFC_OpenACC={'ON' if (ARG('gpu') == gpuConfigOptions.ACC.value) else 'OFF'}")
            flags.append(f"-DMFC_OpenMP={'ON' if (ARG('gpu') == gpuConfigOptions.MP.value) else 'OFF'}")
            flags.append(f"-DMFC_GCov={   'ON' if ARG('gcov') else 'OFF'}")
            flags.append(f"-DMFC_Unified={'ON' if ARG('unified') else 'OFF'}")
            flags.append(f"-DMFC_Fastmath={'ON' if ARG('fastmath') else 'OFF'}")

        command = ["cmake"] + flags + ["-S", cmake_dirpath, "-B", build_dirpath]

        delete_directory(build_dirpath)
        create_directory(build_dirpath)

        case.generate_fpp(self)

        debug(f"Configuring {self.name} in {build_dirpath}")
        debug(f"CMake flags: {' '.join(flags)}")

        verbosity = ARG('verbose')
        if verbosity >= 2:
            # -vv or higher: show raw cmake output
            level_str = "vv" + "v" * (verbosity - 2) if verbosity > 2 else "vv"
            cons.print(f"  [bold blue]Configuring[/bold blue] [magenta]{self.name}[/magenta] [dim](-{level_str})[/dim]...")
            if verbosity >= 2:
                cons.print(f"  [dim]$ {' '.join(str(c) for c in command)}[/dim]")
            cons.print()
            result = system(command, print_cmd=False)
        else:
            # Normal mode: capture output, show on error
            cons.print(f"  [bold blue]Configuring[/bold blue] [magenta]{self.name}[/magenta]...")
            result = system(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, print_cmd=False)

        if result.returncode != 0:
            cons.print(f"  [bold red]✗[/bold red] Configuration failed for [magenta]{self.name}[/magenta]")
            if verbosity < 2:
                _show_build_error(result, "Configuration")
            Tips.after_build_failure()
            raise MFCException(f"Failed to configure the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(f"  [bold green]✓[/bold green] Configured [magenta]{self.name}[/magenta]")
        cons.print(no_indent=True)

    def build(self, case: input.MFCInputFile):
        case.generate_fpp(self)

        command = ["cmake", "--build",    self.get_staging_dirpath(case),
                            "--target",   self.name,
                            "--parallel", ARG("jobs"),
                            "--config",   'Debug' if ARG('debug') else 'Release']

        verbosity = ARG('verbose')
        # -vv or higher: add cmake --verbose flag for full compiler commands
        if verbosity >= 2:
            command.append("--verbose")

        debug(f"Building {self.name} with {ARG('jobs')} parallel jobs")
        debug(f"Build command: {' '.join(str(c) for c in command)}")

        if verbosity >= 2:
            # -vv or higher: show raw compiler output (full verbose)
            level_str = "vv" + "v" * (verbosity - 2) if verbosity > 2 else "vv"
            cons.print(f"  [bold blue]Building[/bold blue] [magenta]{self.name}[/magenta] [dim](-{level_str})[/dim]...")
            cons.print(f"  [dim]$ {' '.join(str(c) for c in command)}[/dim]")
            cons.print()
            result = system(command, print_cmd=False)
        elif verbosity == 1:
            # -v: show ninja [X/Y] lines as they compile (streaming, no progress bar)
            result = _run_build_with_progress(command, self.name, streaming=True)
        else:
            # Default: show progress bar
            result = _run_build_with_progress(command, self.name, streaming=False)

        if result.returncode != 0:
            cons.print(f"  [bold red]✗[/bold red] Build failed for [magenta]{self.name}[/magenta]")
            if verbosity < 2:
                _show_build_error(result, "Build")
            Tips.after_build_failure()
            raise MFCException(f"Failed to build the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(f"  [bold green]✓[/bold green] Built [magenta]{self.name}[/magenta]")
        cons.print(no_indent=True)

    def install(self, case: input.MFCInputFile):
        command = ["cmake", "--install", self.get_staging_dirpath(case)]

        # Show progress indicator during install
        cons.print(f"  [bold blue]Installing[/bold blue] [magenta]{self.name}[/magenta]...")

        # Capture output to show detailed errors on failure
        result = system(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, print_cmd=False)
        if result.returncode != 0:
            cons.print(f"  [bold red]✗[/bold red] Install failed for [magenta]{self.name}[/magenta]")
            _show_build_error(result, "Install")
            raise MFCException(f"Failed to install the [bold magenta]{self.name}[/bold magenta] target.")

        cons.print(f"  [bold green]✓[/bold green] Installed [magenta]{self.name}[/magenta]")
        cons.print(no_indent=True)

#                         name             flags                       isDep  isDef  isReq  dependencies                        run order
FFTW          = MFCTarget('fftw',          ['-DMFC_FFTW=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []), -1)
HDF5          = MFCTarget('hdf5',          ['-DMFC_HDF5=ON'],          True,  False, False, MFCTarget.Dependencies([], [], []), -1)
SILO          = MFCTarget('silo',          ['-DMFC_SILO=ON'],          True,  False, False, MFCTarget.Dependencies([HDF5], [], []), -1)
LAPACK        = MFCTarget('lapack',        ['-DMFC_LAPACK=ON'],        True,  False, False, MFCTarget.Dependencies([],[],[]), -1)
HIPFORT       = MFCTarget('hipfort',       ['-DMFC_HIPFORT=ON'],       True,  False, False, MFCTarget.Dependencies([], [], []), -1)
PRE_PROCESS   = MFCTarget('pre_process',   ['-DMFC_PRE_PROCESS=ON'],   False, True,  False, MFCTarget.Dependencies([], [], []), 0)
SIMULATION    = MFCTarget('simulation',    ['-DMFC_SIMULATION=ON'],    False, True,  False, MFCTarget.Dependencies([], [FFTW], [HIPFORT]), 1)
POST_PROCESS  = MFCTarget('post_process',  ['-DMFC_POST_PROCESS=ON'],  False, True,  False, MFCTarget.Dependencies([FFTW, HDF5, SILO, LAPACK], [], []), 2)
SYSCHECK      = MFCTarget('syscheck',      ['-DMFC_SYSCHECK=ON'],      False, False, True,  MFCTarget.Dependencies([], [], [HIPFORT]), -1)
DOCUMENTATION = MFCTarget('documentation', ['-DMFC_DOCUMENTATION=ON'], False, False, False, MFCTarget.Dependencies([], [], []), -1)

TARGETS = { FFTW, HDF5, SILO, LAPACK, HIPFORT, PRE_PROCESS, SIMULATION, POST_PROCESS, SYSCHECK, DOCUMENTATION }

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


def __generate_header(case: input.MFCInputFile, targets: typing.List):
    feature_flags = [
        'Build',
        format_list_to_string([ t.name for t in get_targets(targets) ], 'magenta')
    ]
    if ARG("case_optimization"):
        feature_flags.append(f"Case Optimized: [magenta]{ARG('input')}[/magenta]")
    if case.params.get('chemistry', 'F') == 'T':
        feature_flags.append(f"Chemistry: [magenta]{case.get_cantera_solution().source}[/magenta]")

    return f"[bold]{' | '.join(feature_flags or ['Generic'])}[/bold]"


def build(targets = None, case: input.MFCInputFile = None, history: typing.Set[str] = None):
    if history is None:
        history = set()
    if isinstance(targets, (MFCTarget, str)):
        targets = [ targets ]
    if targets is None:
        targets = ARG("targets")

    targets = get_targets(list(REQUIRED_TARGETS) + targets)
    case    = case or input.load(ARG("input"), ARG("--"), {})
    case.validate_params()

    if len(history) == 0:
        cons.print(__generate_header(case, targets))
        cons.print(no_indent=True)

    for target in targets:
        __build_target(target, case, history)

    if len(history) == 0:
        cons.print(no_indent=True)
