import conf, lock, user, common, mfc

import rich, rich.progress

import io
import os
import re
import sys
import copy
import shutil
import subprocess
import dataclasses
import urllib.request

class MFCBuild:
    def __init__(self, mfc: mfc.MFCState) -> None:
        self.mfc = mfc
        self.check_environment()

    def get_mode_base_path(self, cc: str = None):
        if cc is None:
            cc = self.mfc.args["mode"]

        return f'{common.MFC_SUBDIR}/{cc}'

    def get_target_base_path(self, name: str):
        default_cfg_name = self.mfc.conf.get_desired_target_mode_name(name)
        cc = self.mfc.conf.get_target_mode_folder_name(name, default_cfg_name)

        return f'{common.MFC_SUBDIR}/{cc}'

    def get_source_path(self, name: str):
        return f'{self.get_target_base_path(name)}/src/{name}'

    def get_build_path(self, name: str):
        return f'{self.get_target_base_path(name)}/build'

    def get_log_filepath(self, name: str):
        return f'{self.get_target_base_path(name)}/log/{name}.log'

    def get_temp_path(self, name: str):
        return f'{self.get_target_base_path(name)}/temp/{name}'

    def get_desired_target_mode(self, name: str) -> user.Mode:
        return self.mfc.user.get_mode(self.mfc.conf.get_desired_target_mode_name(name))

    def setup_directories(self):
        common.create_directory(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for mode in [ mode.name for mode in self.mfc.user.modes ] + ["common"]:
                common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}/{build_subdir}")

    def check_environment(self):
        rich.print("[bold][u]Check Environment:[/u][/bold]")

        utilities = ["python3", "python3-config", "make", "git"]
        utilities += [ compiler for compiler in list(vars(self.mfc.user.build.compilers).values()) ]

        for utility in utilities:
            if shutil.which(utility) is None:
                raise common.MFCException(
                    f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

        rich.print(f"> Found [bold blue]{len(utilities)}/{len(utilities)}[/bold blue] Command-line utilities [bold green]✓[/bold green]")

        # Run checks on the user's current compilers
        def compiler_str_replace(s: str):
            s = s.replace("${CC}",  self.mfc.user.build.compilers.c)
            s = s.replace("${CXX}", self.mfc.user.build.compilers.cpp)
            s = s.replace("${FC}",  self.mfc.user.build.compilers.fortran)

            return s

        for compiler in self.mfc.conf.compilers:
            compiler: conf.Compiler

            # Check if used
            is_used_cmd = compiler_str_replace(compiler.is_used_cmd)
            if 0 != common.execute_shell_command(is_used_cmd, no_exception=True):
                continue

            if not compiler.supported:
                raise common.MFCException(f"{compiler.name} is unsupported by MFC. Please check the documentation for a list of supported compilers.")

            version_fetch_cmd     = compiler_str_replace(compiler.get_version)
            version_fetch_cmd_out = subprocess.check_output(version_fetch_cmd, shell=True, encoding='UTF-8').split()[0]

            def get_ver_from_str(s: str) -> int:
                return int("".join([ n.zfill(4) for n in re.findall("[0-9]+", s) ]))

            version_num_fetched = get_ver_from_str(version_fetch_cmd_out)
            version_num_minimum = get_ver_from_str(compiler.min_version)

            if version_num_fetched >= version_num_minimum:
                rich.print(f"""\
> [bold blue]{compiler.name}[/bold blue] [bold magenta]v{version_fetch_cmd_out}[/bold magenta] \
>= [bold magenta]v{compiler.min_version}[/bold magenta] [bold green]✓[/bold green]""")
            else:
                raise common.MFCException(f"{compiler.name} requires version v{compiler.min_version} or above.")

        # TODO: MacOS Checks
        if sys.platform == "darwin": # MacOS
            pass


    def get_cuda_libdirpath(self) -> str:
        matches = list(filter(lambda test_key: test_key in [ "CUDA_HOME", "CUDA_DIR" ], os.environ))

        if len(matches) == 0:
            raise common.MFCException(f'Failed to locate Cuda. If you have Cuda installed, please define $CUDA_HOME and/or $CUDA_DIR.')

        return os.environ[matches[0]].rstrip('/')


    def string_replace(self, dependency_name: str, string: str, recursive=True):
        dep       = self.mfc.conf.get_target(dependency_name)
        compilers = self.mfc.user.build.compilers

        mode = self.get_desired_target_mode(dependency_name)

        install_path = self.get_build_path (dependency_name)
        source_path  = self.get_source_path(dependency_name)

        flags = vars(copy.deepcopy(mode))
        for lang in flags.keys():
            lang: str
            if "${CUDA:INSTALL_PATH}" in flags[lang]:
                flags[lang] = flags[lang].replace("${CUDA:INSTALL_PATH}", self.get_cuda_libdirpath())

        replace_list = [
            ("${MFC_ROOT_PATH}",     common.MFC_ROOTDIR),
            ("${CONFIGURE_OPTIONS}", f'--prefix="{install_path}"'),
            ("${SOURCE_PATH}",       source_path),
            ("${INSTALL_PATH}",      install_path),
            ("${INSTALL_PATH}",      install_path),
            ("${MAKE_OPTIONS}",      f'-j {self.mfc.args["jobs"]}'),
            ("${COMPILER_FLAGS}",    f'CFLAGS="{flags.get("c")}" CPPFLAGS="{flags.get("cpp")}" FFLAGS="{flags.get("fortran")}"'),
            ("${COMPILERS}",         f'CC="{compilers.c}" CXX="{compilers.cpp}" FC="{compilers.fortran}"')
        ]

        for e in replace_list:
            string = string.replace(*e)
        
        if "${CUDA:INSTALL_PATH}" in string:
            string = string.replace("${CUDA:INSTALL_PATH}", self.get_cuda_libdirpath())

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
            for dep2_info in self.mfc.conf.targets:
                string = string.replace("${" + dep2_info.name         + ":", "${")
                string = string.replace("${" + dep2_info.name.upper() + ":", "${")
                string = self.string_replace(dep2_info.name, string, recursive=False)

        return string

    def check_build_status(self, name: str, bIgnoreCleans: bool):
        # Check if it hasn't been built before
        target_mode = self.get_desired_target_mode(name)

        if not self.mfc.lock.does_target_exist(name, target_mode.name):
            return False

        # Retrive CONF & LOCK descriptors
        conf_desc = self.mfc.conf.get_target(name)
        lock_desc = self.mfc.lock.get_target(name, target_mode.name)

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
           conf_desc.fetch.params != lock_desc.target.fetch.params:
            return False

        # Check if build commands have changed
        if conf_desc.build != lock_desc.target.build:
            return False

        # Check if any of its dependencies needs updating
        for dependency_name in self.mfc.conf.get_dependency_names(name, recursive=True):
            if not self.check_build_status(dependency_name, bIgnoreCleans):
                return False

        # Check if target was cleaned
        if not bIgnoreCleans:
            if self.mfc.lock.get_target(name, target_mode.name).metadata.bCleaned:
                return False

        # Check for "scratch" flag
        if self.mfc.args["scratch"]:
            return False

        return True

    def is_built(self, name: str):
        return self.check_build_status(name, bIgnoreCleans=True)

    def build_should_rebuild(self, name: str):
        return not self.check_build_status(name, bIgnoreCleans=False)

    def build_target__clean_previous(self, name: str, depth: str):
        target_mode = self.get_desired_target_mode(name)
        if not self.mfc.lock.does_unique_target_exist(name, target_mode.name):
            return

        conf_desc = self.mfc.conf.get_target(name)
        lock_desc = self.mfc.lock.get_target(name, target_mode.name)

        if ((    conf_desc.fetch.method != lock_desc.target.fetch.method
             and lock_desc.target.fetch.method in ["clone", "download"]
            ) or (self.mfc.args["scratch"])):
            common.delete_directory_recursive(f'{common.MFC_SUBDIR}/{lock_desc.metadata.mode}/src/{name}')

    def build_target__fetch(self, name: str, logfile: io.IOBase, depth: str):
        target_mode = self.get_desired_target_mode(name)
        conf = self.mfc.conf.get_target(name)

        if conf.fetch.method in ["clone", "download"]:
            if conf.fetch.method == "clone":
                lock_matches = self.mfc.lock.get_target_matches(name, target_mode.name)

                if ((   len(lock_matches)    == 1
                    and conf.fetch.params.git != self.mfc.lock.get_target(name, target_mode.name).target.fetch.params.git)
                    or (self.mfc.args["scratch"])):
                    rich.print(f'{depth}GIT repository changed. Updating...')

                    common.delete_directory_recursive(self.get_source_path(name))

                if not os.path.isdir(self.get_source_path(name)):
                    rich.print(f'{depth}Cloning repository...')

                    common.execute_shell_command(
                        f'git clone --recursive "{conf.fetch.params.git}" "{self.get_source_path(name)}" >> "{logfile.name}" 2>&1')

                rich.print(f"{depth}Checking out [yellow]{conf.fetch.params.hash}[/yellow]...")

                common.execute_shell_command(
                    f'cd "{self.get_source_path(name)}" && git checkout "{conf.fetch.params.hash}" >> "{logfile.name}" 2>&1')
            elif conf.fetch.method == "download":
                rich.print(f'{depth}Removing previously downloaded version...')

                common.delete_directory_recursive(self.get_source_path(name))

                download_link = conf.fetch.params.link.replace("${VERSION}", conf.fetch.params.version)
                filename = download_link.split("/")[-1]

                rich.print(f'{depth}Downloading source...')

                common.create_directory(self.get_temp_path(name))

                download_path = f'{self.get_temp_path(name)}/{filename}'
                urllib.request.urlretrieve(download_link, download_path)

                rich.print(f'{depth}Uncompressing archive...')

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
        conf = self.mfc.conf.get_target(name)

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

                cmd_exception_text=f"""\
Above is the output of {name}'s build command that failed. (#{cmd_idx+1} in mfc.conf.yaml)
  You can also view it by running:\n\ncat \"{logfile.name}\"\n\
"""

                common.execute_shell_command(command, exception_text=cmd_exception_text, on_error=cmd_on_error)
        elif conf.fetch.method == "collection":
            pass
        else:
            raise common.MFCException(f'Unknown target type "{conf.fetch.method}".')


    def build_target__update_lock(self, name: str, depth: str):
        target_mode = self.get_desired_target_mode(name)
        conf = self.mfc.conf.get_target(name)

        rich.print(f'{depth}Updating lock file...')

        new_entry = lock.LockTargetHolder({
            "target": dataclasses.asdict(conf),
            "metadata": {
                "mode":     target_mode.name,
                "bCleaned": False
            }
        })

        # If the target - in the selected mode - isn't already
        # in the lock file, we add a new target/metdata entry into it.
        if len(self.mfc.lock.get_target_matches(name, target_mode.name)) == 0:
            self.mfc.lock.add_target(new_entry)
            self.mfc.lock.save()
            return

        # Otherwise, we simply need to update the existing entry.
        for index, dep in enumerate(self.mfc.lock.targets):
            if dep.target.name == name and dep.metadata.mode == target_mode.name:
                self.mfc.lock.targets[index] = new_entry
                self.mfc.lock.save()
                return

        # If for some reason we can't find the target, throw.
        raise common.MFCException(f"Failed to update the lock file for {name} in the {target_mode.name} mode.")


    def build_target(self, name: str, depth="> "):
        # Proceed to build target
        prepend=f"{depth}Package [bold blue]{name}[/bold blue]"
        # Check if it needs to be (re)built
        if not self.build_should_rebuild(name):
            rich.print(f"{prepend} -> Satisfied [bold green]✓[/bold green]")
            return False

        dependencies = self.mfc.conf.get_dependency_names(name, recursive=False)
        if len(dependencies) > 0:
            format_list = common.format_list_to_string([ f'[bold blue]{x}[/bold blue]' for x in dependencies ])
            rich.print(f"{prepend} requires {format_list}.")

            # Build its dependencies
            for dependency in dependencies:
                self.build_target(dependency, depth+"> ")

        rich.print(f"{prepend}:")

        common.create_file(self.get_log_filepath(name))

        with open(self.get_log_filepath(name), "r+") as logfile:
            self.build_target__clean_previous(name,          depth+"> ") # Clean any old build artifacts
            self.build_target__fetch         (name, logfile, depth+"> ") # Fetch Source Code
            self.build_target__build         (name, logfile, depth+"> ") # Build
            self.build_target__update_lock   (name,          depth+"> ") # Update LOCK

        rich.print(f"{prepend} -> Built [bold green]✓[/bold green]")

        return True
