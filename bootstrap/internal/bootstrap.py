import io
import os
import re
import sys
import copy
import shutil
import colorama
import subprocess
import dataclasses
import urllib.request

import internal.args      as args
import internal.conf      as conf
import internal.lock      as lock
import internal.common    as common
import internal.user      as user
import internal.treeprint as treeprint
import internal.test      as test

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
        common.create_directory_safe(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for cc in [ cc.name for cc in self.user.configurations ] + ["common"]:
                common.create_directory_safe(f"{common.MFC_SUBDIR}/{cc}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory_safe(f"{common.MFC_SUBDIR}/{cc}/{d}/{build_subdir}")

    def check_environment(self):
        self.tree.print("Environment Checks")
        self.tree.indent()

        self.tree.print("Required command-line utilities...")
        self.tree.indent()

        required = ["python3", "python3-config", "make", "git"]
        required += [ compiler for compiler in list(vars(self.user.compilers).values()) ]

        for index, utility in enumerate(required):
            common.clear_line()
            self.tree.print(f"{index+1}/{len(required)} Checking for {utility}...", end='\r')

            if shutil.which(utility) is None:
                raise common.MFCException(
                    f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

        common.clear_line()
        self.tree.print(f"Found {len(required)}/{len(required)}. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})")
        self.tree.unindent()

        # Run checks on the user's current compilers
        def compiler_str_replace(s: str):
            s = s.replace("${C}",       self.user.compilers.c)
            s = s.replace("${CPP}",     self.user.compilers.cpp)
            s = s.replace("${FORTRAN}", self.user.compilers.fortran)

            return s

        self.tree.print("Compiler Checks")
        self.tree.indent()

        for check in self.conf.compiler_verions:
            check: conf.CompilerVersion

            # Check if used
            is_used_cmd = compiler_str_replace(check.is_used)
            if 0 != common.execute_shell_command_safe(is_used_cmd, no_exception=True):
                continue

            version_fetch_cmd     = compiler_str_replace(check.fetch)
            version_fetch_cmd_out = subprocess.check_output(version_fetch_cmd, shell=True, encoding='UTF-8').split()[0]

            def get_ver_from_str(s: str) -> int:
                return int("".join([ n.zfill(4) for n in re.findall("[0-9]+", s) ]))

            version_num_fetched = get_ver_from_str(version_fetch_cmd_out)
            version_num_minimum = get_ver_from_str(check.minimum)

            if version_num_fetched >= version_num_minimum:
                self.tree.print(f'{check.name}: v{version_fetch_cmd_out} >= v{check.minimum}. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})')

        self.tree.unindent()

        # TODO: MacOS Checks
        if sys.platform == "darwin": # MacOS
            pass

        self.tree.print(f'Passing. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})')
        self.tree.unindent()

    def string_replace(self, dependency_name: str, string: str, recursive=True):
        dep       = self.conf.get_target(dependency_name)
        compilers = self.user.compilers

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

        # Check if it needs updating (LOCK & CONFIG descriptions don't match)
        if conf_desc.fetch.method != lock_desc.target.fetch.method    or \
           lock_desc.metadata.bCleaned                                or \
           conf_desc.fetch.params != lock_desc.target.fetch.params:
            return False

        if conf_desc.fetch.method == "source":
            # TODO:
            # FIXME:
            pass

        # Check if any of its dependencies needs updating
        for dependency_name in self.conf.get_dependency_names(name, recursive=True):
            if not self.is_build_satisfied(dependency_name):
                return False

        # Check for "scratch" flag
        if self.args.tree_get("scratch"):
            return False

        return True

    def build_target__clean_previous(self, name: str):
        compiler_cfg = self.get_target_configuration(name, self.args.tree_get("compiler_configuration"))
        if not self.lock.does_unique_target_exist(name, compiler_cfg.name):
            return

        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, compiler_cfg.name)

        if ((    conf_desc.fetch.method != lock_desc.target.fetch.method
             and lock_desc.target.fetch.method in ["clone", "download"]
            ) or (self.args.tree_get("scratch"))):
            common.delete_directory_recursive_safe(f'{common.MFC_SUBDIR}/{lock_desc.metadata.compiler_configuration}/src/{name}')

    def build_target__fetch(self, name: str, logfile: io.IOBase):
        compiler_cfg = self.get_target_configuration(name, self.args.tree_get("compiler_configuration"))
        conf = self.conf.get_target(name)

        if conf.fetch.method in ["clone", "download"]:
            if conf.fetch.method == "clone":
                lock_matches = self.lock.get_target_matches(name, compiler_cfg.name)

                if ((   len(lock_matches)    == 1
                    and conf.fetch.params.git != self.lock.get_target(name, compiler_cfg.name).target.fetch.params.git)
                    or (self.args.tree_get("scratch"))):
                    self.tree.print(f'GIT repository changed. Updating...')

                    common.delete_directory_recursive_safe(self.get_source_path(name))

                if not os.path.isdir(self.get_source_path(name)):
                    self.tree.print(f'Cloning repository...')

                    common.execute_shell_command_safe(
                        f'git clone --recursive "{conf.fetch.params.git}" "{self.get_source_path(name)}" >> "{logfile.name}" 2>&1')

                self.tree.print(f'Checking out {conf.fetch.params.hash}...')

                common.execute_shell_command_safe(
                    f'cd "{self.get_source_path(name)}" && git checkout "{conf.fetch.params.hash}" >> "{logfile.name}" 2>&1')

                self.tree.print(f'Deleting .git/ to save space...')

                common.delete_directory_recursive_safe(f"{self.get_source_path(name)}/.git")
            elif conf.fetch.method == "download":
                self.tree.print(f'Removing previously downloaded version...')

                common.delete_directory_recursive_safe(self.get_source_path(name))

                download_link = conf.fetch.params.link.replace("${VERSION}", conf.fetch.params.version)
                filename = download_link.split("/")[-1]

                self.tree.print(f'Downloading source...')

                common.create_directory_safe(self.get_temp_path(name))

                download_path = f'{self.get_temp_path(name)}/{filename}'
                urllib.request.urlretrieve(download_link, download_path)

                self.tree.print(f'Uncompressing archive...')

                common.uncompress_archive_to(download_path,
                                      f'{self.get_source_path(name)}')

                os.remove(download_path)
        elif conf.fetch.method == "source":
            if os.path.isdir(self.get_source_path(name)):
                common.delete_directory_recursive_safe(self.get_source_path(name))

            shutil.copytree(self.string_replace(name, conf.fetch.params.source),
                            self.get_source_path(name))
        elif conf.fetch.method == "collection":
            common.create_directory_safe(self.get_source_path(name))
        else:
            raise common.MFCException(f'Dependency type "{conf.fetch.method}" is unsupported.')

    def build_target__build(self, name: str, logfile: io.IOBase):
        conf = self.conf.get_target(name)

        if conf.fetch.method in ["clone", "download", "source"]:
            for cmd_idx, command in enumerate(conf.build):
                self.tree.print_progress(f'Building', cmd_idx+1, len(conf.build))

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

                common.execute_shell_command_safe(command, exception_text=cmd_exception_text, on_error=cmd_on_error)
        elif conf.fetch.method == "collection":
            pass
        else:
            raise common.MFCException(f'Unknown target type "{conf.fetch.method}".')


    def build_target__update_lock(self, name: str):
        compiler_cfg = self.get_target_configuration(name, self.args["compiler_configuration"])
        conf = self.conf.get_target(name)

        self.tree.print(f'Updating lock file...')

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

    def build_target(self, name: str):
        possible_colors = [colorama.Fore.BLUE, colorama.Fore.CYAN, colorama.Fore.GREEN, colorama.Fore.YELLOW]

        color = possible_colors[self.tree.get_depth() % len(possible_colors)]

        self.tree.print(f"{color}Building Package {name}{colorama.Style.RESET_ALL}")
        self.tree.indent(color=color)

        # Check if it needs to be (re)built
        if self.is_build_satisfied(name):
            self.tree.print(f'Nothing to do ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})')
            self.tree.unindent()
            return False

        # Build its dependencies
        for dependency_name in self.conf.get_dependency_names(name, recursive=False):
            self.build_target(dependency_name)

        self.tree.print(f'Preparing build...')

        common.create_file_safe(self.get_log_filepath(name))

        with open(self.get_log_filepath(name), "r+") as logfile:
            self.build_target__clean_previous(name)          # Clean any old build artifacts
            self.build_target__fetch         (name, logfile) # Fetch Source Code
            self.build_target__build         (name, logfile) # Build
            self.build_target__update_lock   (name)          # Update LOCK

        self.tree.print(f'Done. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})')
        self.tree.unindent()

        return True

    def print_header(self):
        print(f"""{colorama.Fore.BLUE}
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
{colorama.Style.RESET_ALL}\

   +----------------------------------+
   |    Multi-component Flow Code     |
   |                                  |
   | https://github.com/MFlowCode/MFC |
   +----------------------------------+
""")

    def test_target(self, name: str): 
        # TODO: ERROR if name is not MFC
        self.test.test()

    def clean_target(self, name: str):
        if not self.is_build_satisfied(name):
            raise common.MFCException(f"Can't clean {name} because its build isn't satisfied.")

        self.tree.print(f"Cleaning Package {name}")
        self.tree.indent()

        for dependency_name in self.conf.get_dependency_names(name, recursive=False):
            if not self.conf.is_target_common(dependency_name):
                self.clean_target(dependency_name)

        target = self.lock.get_target(name, self.args["compiler_configuration"])

        if not target.metadata.bCleaned:
            if os.path.isdir(self.string_replace(name, "${SOURCE_PATH}")):
                with open(self.get_log_filepath(name), "a") as log_file:
                    for cmd_idx, command in enumerate(target.target.clean):
                        self.tree.print(f'Cleaning [{cmd_idx+1}/{len(target.target.clean)}] (Logging to {log_file.name})...')

                        command = self.string_replace(name, f"""\
        cd "${{SOURCE_PATH}}" && \
        stdbuf -oL bash -c '{command}' >> "{log_file.name}" 2>&1""")

                        log_file.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                        log_file.flush()

                        common.execute_shell_command_safe(command)

            target.metadata.bCleaned = True

        common.delete_file_safe(self.get_log_filepath(name))

        self.lock.flush()
        self.lock.save()

        self.tree.print(f"Cleaning done. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})")
        self.tree.unindent()

    def __init__(self):
        self.tree = treeprint.TreePrinter()

        self.conf = conf.MFCConf()
        self.user = user.MFCUser()
        self.setup_directories()
        self.lock = lock.MFCLock()
        self.args = args.MFCArgs(self.conf, self.user)
        self.test = test.MFCTest(self)

        self.print_header()
        self.check_environment()

        if self.args["set_current"] is not None:
            common.update_symlink(f"{common.MFC_SUBDIR}/___current___", self.get_configuration_base_path(self.args["set_current"]))

        # Update symlink to current build
        if self.args["build"]:
            common.update_symlink(f"{common.MFC_SUBDIR}/___current___", self.get_configuration_base_path())

        if self.args["test"]:
            self.test.test()

        for target_name in [ x.name for x in self.conf.targets ]:
            if target_name in self.args["targets"]:
                if self.args["build"]: self.build_target(target_name)
                if self.args["clean"]: self.clean_target(target_name)

        self.lock.save()
