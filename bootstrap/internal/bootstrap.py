import io
import os
import re
import sys
import copy
import shutil
import subprocess
import dataclasses
import urllib.request


import rich
import rich.tree
import rich.text
import rich.live
import rich.progress


import internal.run    as run
import internal.args   as args
import internal.conf   as conf
import internal.lock   as lock
import internal.user   as user
import internal.test   as test
import internal.common as common

class Bootstrap:
    def get_configuration_base_path(self, cc: str = None):
        if cc is None:
            cc = self.args["compiler_configuration"]

        return f'{common.MFC_SUBDIR}/{cc}'

    def get_target_base_path(self, name: str):
        default_cfg_name = self.conf.get_target_configuration_name(name, self.args["compiler_configuration"])
        cc = self.conf.get_target_configuration_folder_name(name, default_cfg_name)

        return f'{common.MFC_SUBDIR}/{cc}'

    def get_source_path(self, name: str):
        return f'{self.get_target_base_path(name)}/src/{name}'

    def get_build_path(self, name: str):
        return f'{self.get_target_base_path(name)}/build'

    def get_log_filepath(self, name: str):
        return f'{self.get_target_base_path(name)}/log/{name}.log'

    def get_temp_path(self, name: str):
        return f'{self.get_target_base_path(name)}/temp/{name}'

    def get_target_configuration(self, name: str, default: str) -> user.Configuration:
        return self.user.get_configuration(self.conf.get_target_configuration_name(name, default))

    def setup_directories(self):
        common.create_directory(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for cc in [ cc.name for cc in self.user.configurations ] + ["common"]:
                common.create_directory(f"{common.MFC_SUBDIR}/{cc}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory(f"{common.MFC_SUBDIR}/{cc}/{d}/{build_subdir}")

    def check_environment(self):
        self.console.print("[bold][u]Check Environment:[/u][/bold]")

        utilities = ["python3", "python3-config", "make", "git"]
        utilities += [ compiler for compiler in list(vars(self.user.build.compilers).values()) ]

        for utility in utilities:
            if shutil.which(utility) is None:
                raise common.MFCException(
                    f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

        self.console.print(f"> Found [bold blue]{len(utilities)}/{len(utilities)}[/bold blue] Command-line utilities [bold green]✓[/bold green]")

        # Run checks on the user's current compilers
        def compiler_str_replace(s: str):
            s = s.replace("${C}",       self.user.build.compilers.c)
            s = s.replace("${CPP}",     self.user.build.compilers.cpp)
            s = s.replace("${FORTRAN}", self.user.build.compilers.fortran)

            return s

        for check in self.conf.compiler_verions:
            check: conf.CompilerVersion

            # Check if used
            is_used_cmd = compiler_str_replace(check.is_used)
            if 0 != common.execute_shell_command(is_used_cmd, no_exception=True):
                continue

            version_fetch_cmd     = compiler_str_replace(check.fetch)
            version_fetch_cmd_out = subprocess.check_output(version_fetch_cmd, shell=True, encoding='UTF-8').split()[0]

            def get_ver_from_str(s: str) -> int:
                return int("".join([ n.zfill(4) for n in re.findall("[0-9]+", s) ]))

            version_num_fetched = get_ver_from_str(version_fetch_cmd_out)
            version_num_minimum = get_ver_from_str(check.minimum)

            if version_num_fetched >= version_num_minimum:
                self.console.print(
f"""> [bold blue]{check.name}[/bold blue] [bold magenta]v{version_fetch_cmd_out}[/bold magenta] \
>= [bold magenta]v{check.minimum}[/bold magenta] [bold green]✓[/bold green]""")
            else:
                raise common.MFCException(f"Compiler check {check.name} failed. Version v{check.minimum} minimum.")

        # TODO: MacOS Checks
        if sys.platform == "darwin": # MacOS
            pass


    def string_replace(self, dependency_name: str, string: str, recursive=True):
        dep       = self.conf.get_target(dependency_name)
        compilers = self.user.build.compilers

        configuration = self.get_target_configuration(dependency_name, self.args["compiler_configuration"])

        install_path = self.get_build_path (dependency_name)
        source_path  = self.get_source_path(dependency_name)

        flags = vars(copy.deepcopy(configuration))
        for lang in flags.keys():
            lang: str
            if "${CUDA:INSTALL_PATH}" in flags[lang]:
                matches = list(filter(lambda test_key: test_key in [ "CUDA_HOME", "CUDA_DIR" ], os.environ))

                if len(matches) == 0:
                    raise common.MFCException(f'''\
Failed to find where CUDA was installed for {dependency_name} with {configuration.name}/{lang}.
Please follow the instructions bellow:
- Make sure CUDA is installed and properly configured.
- Open mfc.conf.py.
- Locate section compilers -> configurations -> {configuration.name} -> {lang}:
- Replace $(CUDA:INSTALL_PATH) with the root path to your CUDA installation.
  "include" and "lib" should be folders directly accessible from this folder.

If you think MFC could (or should) be able to find it automatically for you system, you are welcome to file an issue on GitHub or a pull request with your changes to mfc.py at https://github.com/MFlowCode/MFC.
''')

                cuda_install_path = os.environ[matches[0]]

                flags[lang] = flags[lang].replace("${CUDA:INSTALL_PATH}", cuda_install_path)

        replace_list = [
            ("${MFC_ROOT_PATH}",     common.MFC_ROOTDIR),
            ("${CONFIGURE_OPTIONS}", f'--prefix="{install_path}"'),
            ("${SOURCE_PATH}",       source_path),
            ("${INSTALL_PATH}",      install_path),
            ("${INSTALL_PATH}",      install_path),
            ("${MAKE_OPTIONS}",      f'-j {self.args["jobs"]}'),
            ("${COMPILER_FLAGS}",    f'CFLAGS="{flags.get("c")}" CPPFLAGS="{flags.get("cpp")}" FFLAGS="{flags.get("fortran")}"'),
            ("${COMPILERS}",         f'CC="{compilers.c}" CXX="{compilers.cpp}" FC="{compilers.fortran}"')
        ]

        for e in replace_list:
            string = string.replace(*e)

        # Combine different assignments to flag variables (CFLAGS, FFLAGS, ...)
        for FLAG_NAME in [ "CFLAGS", "CPPFLAGS", "FFLAGS" ]:
            FLAG_PATTERN = f' {FLAG_NAME}=".*?"'

            matches = re.findall(FLAG_PATTERN, string)

            if len(matches) <= 1:
                continue

            for i in range(len(matches)):
                matches[i] = re.sub(f'^ {FLAG_NAME}="', ' ', matches[i])
                matches[i] = re.sub(r'"$', ' ', matches[i])

            string = re.sub(FLAG_PATTERN, ' ', string, len(matches) - 1)
            string = re.sub(FLAG_PATTERN, f' {FLAG_NAME}="{" ".join(matches)}"', string)

        # Fetch
        if recursive:
            for dep2_info in self.conf.targets:
                string = string.replace("${" + dep2_info.name         + ":", "${")
                string = string.replace("${" + dep2_info.name.upper() + ":", "${")
                string = self.string_replace(dep2_info.name, string, recursive=False)

        return string

    def is_build_satisfied(self, name: str):
        # Check if it hasn't been built before
        compiler_cfg = self.get_target_configuration(name, self.args.tree_get("compiler_configuration"))

        if not self.lock.does_target_exist(name, compiler_cfg.name):
            return False

        # Retrive CONF & LOCK descriptors
        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, compiler_cfg.name)

        # Check if any source file is newer than the previously built executable
        if conf_desc.fetch.method == "source":
            check_filepath = self.string_replace(conf_desc.name, conf_desc.fetch.params.check)
            if not os.path.isfile(check_filepath):
                return False
            
            last_check_date = os.path.getmtime(check_filepath)
            for subdir, dirs, files in os.walk(self.string_replace(conf_desc.name, conf_desc.fetch.params.source)):
                for file in files:
                    if os.path.getmtime(os.path.join(subdir, file)) > last_check_date:
                        return False
            
        # Check if it needs updating (LOCK & CONFIG descriptions don't match)
        if conf_desc.fetch.method != lock_desc.target.fetch.method    or \
           lock_desc.metadata.bCleaned                                or \
           conf_desc.fetch.params != lock_desc.target.fetch.params:
            return False

        # Check if build commands have changed
        if conf_desc.build != lock_desc.target.build:
            return False

        # Check if any of its dependencies needs updating
        for dependency_name in self.conf.get_dependency_names(name, recursive=True):
            if not self.is_build_satisfied(dependency_name):
                return False

        # Check for "scratch" flag
        if self.args.tree_get("scratch"):
            return False

        return True

    def build_target__clean_previous(self, name: str, depth: str):
        compiler_cfg = self.get_target_configuration(name, self.args.tree_get("compiler_configuration"))
        if not self.lock.does_unique_target_exist(name, compiler_cfg.name):
            return

        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, compiler_cfg.name)

        if ((    conf_desc.fetch.method != lock_desc.target.fetch.method
             and lock_desc.target.fetch.method in ["clone", "download"]
            ) or (self.args.tree_get("scratch"))):
            common.delete_directory_recursive(f'{common.MFC_SUBDIR}/{lock_desc.metadata.compiler_configuration}/src/{name}')

    def build_target__fetch(self, name: str, logfile: io.IOBase, depth: str):
        compiler_cfg = self.get_target_configuration(name, self.args.tree_get("compiler_configuration"))
        conf = self.conf.get_target(name)

        if conf.fetch.method in ["clone", "download"]:
            if conf.fetch.method == "clone":
                lock_matches = self.lock.get_target_matches(name, compiler_cfg.name)

                if ((   len(lock_matches)    == 1
                    and conf.fetch.params.git != self.lock.get_target(name, compiler_cfg.name).target.fetch.params.git)
                    or (self.args.tree_get("scratch"))):
                    self.console.print(f'{depth}GIT repository changed. Updating...')

                    common.delete_directory_recursive(self.get_source_path(name))

                if not os.path.isdir(self.get_source_path(name)):
                    self.console.print(f'{depth}Cloning repository...')

                    common.execute_shell_command(
                        f'git clone --recursive "{conf.fetch.params.git}" "{self.get_source_path(name)}" >> "{logfile.name}" 2>&1')

                self.console.print(f"{depth}Checking out [yellow]{conf.fetch.params.hash}[/yellow]...")

                common.execute_shell_command(
                    f'cd "{self.get_source_path(name)}" && git checkout "{conf.fetch.params.hash}" >> "{logfile.name}" 2>&1')
            elif conf.fetch.method == "download":
                self.console.print(f'{depth}Removing previously downloaded version...')

                common.delete_directory_recursive(self.get_source_path(name))

                download_link = conf.fetch.params.link.replace("${VERSION}", conf.fetch.params.version)
                filename = download_link.split("/")[-1]

                self.console.print(f'{depth}Downloading source...')

                common.create_directory(self.get_temp_path(name))

                download_path = f'{self.get_temp_path(name)}/{filename}'
                urllib.request.urlretrieve(download_link, download_path)

                self.console.print(f'{depth}Uncompressing archive...')

                common.uncompress_archive_to(download_path,
                                      f'{self.get_source_path(name)}')

                os.remove(download_path)
        elif conf.fetch.method == "source":
            dest_src: str = self.get_source_path(name)
            if os.path.isdir(dest_src):
                common.delete_directory_recursive(dest_src)

            # Copy files over
            shutil.copytree(self.string_replace(name, conf.fetch.params.source),
                            dest_src)
        elif conf.fetch.method == "collection":
            common.create_directory(self.get_source_path(name))
        else:
            raise common.MFCException(f'Dependency type "{conf.fetch.method}" is unsupported.')

    def build_target__build(self, name: str, logfile: io.IOBase, depth: str):
        conf = self.conf.get_target(name)

        if conf.fetch.method in ["clone", "download", "source"]:
            for cmd_idx, command in enumerate(rich.progress.track(conf.build, f"{depth}Building...")):
                command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
stdbuf -oL bash -c '{command}' >> "{logfile.name}" 2>&1""")

                logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                logfile.flush()

                def cmd_on_error():
                    print(logfile.read())

                cmd_exception_text=f"Above is the output of {name}'s build command that failed. (#{cmd_idx+1} in mfc.conf.yaml)"
                cmd_exception_text=cmd_exception_text+f"You can also view it by running:\n\ncat \"{logfile.name}\"\n"

                common.execute_shell_command(command, exception_text=cmd_exception_text, on_error=cmd_on_error)
        elif conf.fetch.method == "collection":
            pass
        else:
            raise common.MFCException(f'Unknown target type "{conf.fetch.method}".')


    def build_target__update_lock(self, name: str, depth: str):
        compiler_cfg = self.get_target_configuration(name, self.args["compiler_configuration"])
        conf = self.conf.get_target(name)

        self.console.print(f'{depth}Updating lock file...')

        new_entry = lock.LockTargetHolder({
            "target": dataclasses.asdict(conf),
            "metadata": {
                "compiler_configuration": compiler_cfg.name,
                "bCleaned": False
            }
        })

        # If the target - in the selected configuration - isn't already
        # in the lock file, we add a new target/metdata entry into it.
        if len(self.lock.get_target_matches(name, compiler_cfg.name)) == 0:
            self.lock.add_target(new_entry)
            self.lock.save()
            return

        # Otherwise, we simply need to update the existing entry.
        for index, dep in enumerate(self.lock.targets):
            if dep.target.name == name and dep.metadata.compiler_configuration == compiler_cfg.name:
                self.lock.targets[index] = new_entry
                self.lock.flush()
                self.lock.save()
                return

        # If for some reason we can't find the target, throw.
        raise common.MFCException(f"Failed to update the lock file for {name} in the {compiler_cfg.name} configuration.")

    def build_target(self, name: str, depth=""):
        prepend  =f"{depth}Package [bold blue]{name}[/bold blue]"
        prepend_u=f"{depth}[u]Package [bold blue]{name}[/bold blue][/u]"
        # Check if it needs to be (re)built
        if self.is_build_satisfied(name):
            self.console.print(f"{prepend_u} - satisfied [bold green]✓[/bold green]")
            return False

        dependencies = self.conf.get_dependency_names(name, recursive=False)
        if len(dependencies) > 0:
            self.console.print(f"{prepend} requires {', '.join([ f'[bold blue]{x}[/bold blue]' for x in dependencies ])}.")
            
            # Build its dependencies
            for dependency in dependencies:
                self.build_target(dependency, depth+"> ")

        self.console.print(f"{prepend}")

        common.create_file(self.get_log_filepath(name))

        with open(self.get_log_filepath(name), "r+") as logfile:
            self.build_target__clean_previous(name,          depth+"> ") # Clean any old build artifacts
            self.build_target__fetch         (name, logfile, depth+"> ") # Fetch Source Code
            self.build_target__build         (name, logfile, depth+"> ") # Build
            self.build_target__update_lock   (name,          depth+"> ") # Update LOCK

        self.console.print(f"{prepend_u} - Built [bold green]✓[/bold green]")

        return True

    def print_header(self):
        self.console.print(f"""[bold blue]
     ___            ___          ___
    /__/\          /  /\        /  /\\
   |  |::\        /  /:/_      /  /:/
   |  |:|:\      /  /:/ /\    /  /:/
 __|__|:|\:\    /  /:/ /:/   /  /:/  ___
/__/::::| \:\  /__/:/ /:/   /__/:/  /  /\\
\  \:\~~\__\/  \  \:\/:/    \  \:\ /  /:/
 \  \:\         \  \::/      \  \:\  /:/
  \  \:\         \  \:\       \  \:\/:/
   \  \:\         \  \:\       \  \::/
    \__\/          \__\/        \__\/
[/bold blue]\
""")

    def clean_target(self, name: str):
        if not self.is_build_satisfied(name):
            raise common.MFCException(f"Can't clean {name} because its build isn't satisfied.")

        self.console.print(f"Cleaning Package {name}")

        for dependency_name in self.conf.get_dependency_names(name, recursive=False):
            if not self.conf.is_target_common(dependency_name):
                self.clean_target(dependency_name)

        target = self.lock.get_target(name, self.args["compiler_configuration"])

        if not target.metadata.bCleaned:
            if os.path.isdir(self.string_replace(name, "${SOURCE_PATH}")):
                with open(self.get_log_filepath(name), "a") as log_file:
                    for cmd_idx, command in enumerate(target.target.clean):
                        self.console.print(f'> Cleaning [{cmd_idx+1}/{len(target.target.clean)}] (Logging to {log_file.name})...')

                        command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
stdbuf -oL bash -c '{command}' >> "{log_file.name}" 2>&1""")

                        log_file.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                        log_file.flush()

                        common.execute_shell_command(command)

            target.metadata.bCleaned = True

        common.delete_file(self.get_log_filepath(name))

        self.lock.flush()
        self.lock.save()

        self.console.print(f"Cleaning done [bold green]✓[/bold green]")

    def __init__(self):
        self.console = rich.console.Console()

        self.conf = conf.MFCConf()
        self.user = user.MFCUser()
        self.setup_directories()
        self.lock = lock.MFCLock()
        self.args = args.MFCArgs(self.conf, self.user)
        self.test = test.MFCTest(self)
        self.run  = run.MFCRun(self)

        self.print_header()
        self.check_environment()

        # Update symlink to current build
        if self.args["command"] == "build":
            common.update_symlink(f"{common.MFC_SUBDIR}/___current___", self.get_configuration_base_path())

        if self.args["command"] == "test":
            self.console.print("[bold][u]Test:[/u][/bold]")
            self.test.test()
        
        if self.args["command"] == "run":
            self.run.run()

        for target_name in [ x.name for x in self.conf.targets ]:
            if target_name in self.args["targets"]:
                if self.args["command"] == "build":
                    self.console.print("[bold][u]Build:[/u][/bold]")
                    self.build_target(target_name)
                if self.args["command"] == "clean":
                    self.console.print("[bold][u]Clean:[/u][/bold]")
                    self.clean_target(target_name)

        self.lock.save()
