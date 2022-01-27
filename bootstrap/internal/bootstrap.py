import io
import os
import re
import sys
import glob
import copy
import shutil
import colorama
import urllib.request

import internal.args   as args
import internal.conf   as conf
import internal.lock   as lock
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

    def setup_directories(self):
        common.create_directory_safe(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for cc in [ cc["name"] for cc in self.conf["configurations"] ] + ["common"]:
                common.create_directory_safe(f"{common.MFC_SUBDIR}/{cc}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory_safe(f"{common.MFC_SUBDIR}/{cc}/{d}/{build_subdir}")

    def check_environment(self):
        print("|--> Checking for the presence of required command-line utilities...", end='\r')

        required = ["python3", "python3-config", "make", "git"]
        required += self.conf["compilers"].values()

        for index, utility in enumerate(required):
            common.clear_print(f"|--> {index+1}/{len(required)} Checking for {utility}...", end='\r')

            if shutil.which(utility) is None:
                raise common.MFCException(
                    f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

        # TODO: MacOS Checks
        if sys.platform == "darwin": # MacOS
            pass
            #mfc_state.conf["compilers"]["mpi"]["fortran"]

        common.clear_print(f"|--> Build environment: Passing. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})")

    def string_replace(self, dependency_name: str, string: str, recursive=True):
        dep       = self.conf.get_target(dependency_name)
        compilers = self.conf["compilers"]

        compiler_cfg = self.conf.get_target_configuration(dependency_name, self.args["compiler_configuration"])

        install_path = self.get_build_path (dependency_name)
        source_path  = self.get_source_path(dependency_name)

        flags = copy.deepcopy(compiler_cfg["flags"])
        for lang in flags.keys():
            lang: str
            if "${CUDA:INSTALL_PATH}" in flags[lang]:
                matches = list(filter(lambda test_key: test_key in [ "CUDA_HOME", "CUDA_DIR" ], os.environ))

                if len(matches) == 0:
                    raise common.MFCException(f'''\
Failed to find where CUDA was installed for {dependency_name} with {compiler_cfg["name"]}/{lang}.
Please follow the instructions bellow:
- Make sure CUDA is installed and properly configured.
- Open mfc.conf.py.
- Locate section compilers -> configurations -> {compiler_cfg["name"]} -> {lang}:
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
            ("${COMPILER_FLAGS}",    f'CFLAGS="{flags.get("c", "")}" CPPFLAGS="{flags.get("c++", "")}" FFLAGS="{flags.get("fortran", "")}"'),
            ("${COMPILERS}",         f'CC="{compilers.get("c", "")}" CXX="{compilers.get ("c++", "")}" FC="{compilers.get("fortran", "")}"')
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
            for dep2_info in self.conf["targets"]:
                string = string.replace("${" + dep2_info["name"] + ":", "${")
                string = self.string_replace(dep2_info["name"], string, recursive=False)

        return string

    def is_build_satisfied(self, name: str, ignoreIfSource: bool = False):
        # Check if it hasn't been built before
        compiler_cfg = self.conf.get_target_configuration(name, self.args["compiler_configuration"])
        if not self.lock.does_target_exist(name, compiler_cfg["name"]):
            return False

        # Retrive CONF & LOCK descriptors
        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, compiler_cfg["name"])

        # Check if it needs updating (LOCK & CONFIG descriptions don't match)
        if conf_desc["type"] != lock_desc["type"]                or \
           lock_desc["bCleaned"]                                 or \
           conf_desc["type"] == "source" and not(ignoreIfSource) or \
           conf_desc[conf_desc["type"]] != conf_desc[conf_desc["type"]]:
            return False

        # Check if any of its dependencies needs updating
        for dependency_name in self.conf.get_dependency_names(name, recursive=True):
            if not self.is_build_satisfied(dependency_name, ignoreIfSource=ignoreIfSource):
                return False

        # Check for "scratch" flag
        if self.args["scratch"]:
            return False

        return True

    def build_target__clean_previous(self, name: str):
        compiler_cfg = self.conf.get_target_configuration(name, self.args["compiler_configuration"])
        if not self.lock.does_unique_target_exist(name, compiler_cfg["name"]):
            return

        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, compiler_cfg["name"])

        if ((    conf_desc["type"] != lock_desc["type"]
             and lock_desc["type"] in ["clone", "download"]
            ) or (self.args["scratch"])):
            common.delete_directory_recursive_safe(f'{common.MFC_SUBDIR}/{lock_desc["compiler_configuration"]}/src/{name}')

    def build_target__fetch(self, name: str, logfile: io.IOBase):
        compiler_cfg = self.conf.get_target_configuration(name, self.args["compiler_configuration"])
        conf = self.conf.get_target(name)

        if conf["type"] in ["clone", "download"]:
            if conf["type"] == "clone":
                lock_matches = self.lock.get_target_matches(name, compiler_cfg["name"])

                if ((   len(lock_matches)    == 1
                    and conf["clone"]["git"] != self.lock.get_target(name, compiler_cfg["name"])["clone"]["git"])
                    or (self.args["scratch"])):
                    common.clear_print(f'|--> Package {name}: GIT repository changed. Updating...', end='\r')

                    common.delete_directory_recursive_safe(self.get_source_path(name))

                if not os.path.isdir(self.get_source_path(name)):
                    common.clear_print(f'|--> Package {name}: Cloning repository...', end='\r')

                    common.execute_shell_command_safe(
                        f'git clone --recursive "{conf["clone"]["git"]}" "{self.get_source_path(name)}" >> "{logfile.name}" 2>&1')

                common.clear_print(f'|--> Package {name}: Checking out {conf["clone"]["hash"]}...', end='\r')

                common.execute_shell_command_safe(
                    f'cd "{self.get_source_path(name)}" && git checkout "{conf["clone"]["hash"]}" >> "{logfile.name}" 2>&1')
            elif conf["type"] == "download":
                common.clear_print(f'|--> Package {name}: Removing previously downloaded version...', end='\r')

                common.delete_directory_recursive_safe(self.get_source_path(name))

                download_link = conf[conf["type"]]["link"].replace("${VERSION}", conf[conf["type"]]["version"])
                filename = download_link.split("/")[-1]

                common.clear_print(f'|--> Package {name}: Downloading source...', end='\r')

                common.create_directory_safe(self.get_temp_path(name))

                download_path = f'{self.get_temp_path(name)}/{filename}'
                urllib.request.urlretrieve(download_link, download_path)

                common.clear_print(f'|--> Package {name}: Uncompressing archive...', end='\r')

                common.uncompress_archive_to(download_path,
                                      f'{self.get_source_path(name)}')

                os.remove(download_path)
        elif conf["type"] == "source":
            if os.path.isdir(self.get_source_path(name)):
                common.delete_directory_recursive_safe(self.get_source_path(name))

            shutil.copytree(self.string_replace(name, conf["source"]["source"]),
                            self.get_source_path(name))
        elif conf["type"] == "collection":
            common.create_directory_safe(self.get_source_path(name))
        else:
            raise common.MFCException(f'Dependency type "{conf["type"]}" is unsupported.')

    def build_target__build(self, name: str, logfile: io.IOBase):
        conf = self.conf.get_target(name)

        if conf["type"] in ["clone", "download", "source"]:
            for cmd_idx, command in enumerate(conf["build"]):
                command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
stdbuf -oL bash -c '{command}' >> "{logfile.name}" 2>&1""")

                logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                logfile.flush()

                common.clear_print(f'|--> Package {name}: Building (Logging to {logfile.name})...', end='\r')

                def cmd_on_error():
                    print(logfile.read())

                cmd_exception_text=f"Above is the output of {name}'s build command that failed. (#{cmd_idx+1} in mfc.conf.yaml)"
                cmd_exception_text=cmd_exception_text+f"You can also view it by running:\n\ncat \"{logfile.name}\"\n"

                common.execute_shell_command_safe(command, exception_text=cmd_exception_text, on_error=cmd_on_error)
        elif conf["type"] == "collection":
            pass
        else:
            raise common.MFCException(f'Unknown target type "{conf["type"]}".')


    def build_target__update_lock(self, name: str):
        compiler_cfg = self.conf.get_target_configuration(name, self.args["compiler_configuration"])
        conf = self.conf.get_target(name)

        common.clear_print(f'|--> Package {name}: Updating lock file...', end='\r')

        new_entry = dict(conf, **{"compiler_configuration": compiler_cfg["name"], "bCleaned": False})

        if len(self.lock.get_target_matches(name, compiler_cfg["name"])) == 0:
            self.lock["targets"].append(new_entry)
        else:
            for index, dep in enumerate(self.lock["targets"]):
                if dep["name"] == name:
                    self.lock["targets"][index] = new_entry

        self.lock.save()

    def build_target(self, name: str):
        # Check if it needs to be (re)built
        if self.is_build_satisfied(name):
            common.clear_print(f'|--> Package {name}: Nothing to do ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')
            return False

        # Build its dependencies
        for dependency_names in self.conf.get_dependency_names(name, recursive=False):
            self.build_target(dependency_names)

        common.clear_print(f'|--> Package {name}: Preparing build...', end='\r')

        common.create_file_safe(self.get_log_filepath(name))

        with open(self.get_log_filepath(name), "r+") as logfile:
            self.build_target__clean_previous(name)          # Clean any old build artifacts
            self.build_target__fetch         (name, logfile) # Fetch Source Code
            self.build_target__build         (name, logfile) # Build
            self.build_target__update_lock   (name)          # Update LOCK

        common.clear_print(f'|--> Package {name}: Done. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')

        return True

    def print_header(self):
        print(common.center_ansi_escaped_text(f"""{colorama.Fore.BLUE}
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
{colorama.Style.RESET_ALL}"""))

        print(common.center_ansi_escaped_text("""
+----------------------------------+
|    Multi-component Flow Code     |
|                                  |
| https://github.com/MFlowCode/MFC |
+----------------------------------+

"""))

    def test_target(self, name: str):
        if not self.is_build_satisfied(name, ignoreIfSource=True):
            raise common.MFCException(f"Can't test {name} because its build isn't satisfied.")

        with open(self.get_log_filepath(name), "w") as logfile:
            for command in self.conf.get_target(name)["test"]:
                command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
bash -c '{command}' >> "{logfile.name}" 2>&1""")

                logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                logfile.flush()

                common.execute_shell_command_safe(command)

                common.clear_print(f'|--> Package {name}: Testing (Logging to {logfile.name})...', end='\r')

        common.clear_print(f'|--> Package {name}: Tested. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')

    def clean_target(self, name: str):
        for target in self.lock["targets"]:
            if "clean" in target:
                with open(self.get_log_filepath(name), "a") as log_file:
                    for command in target["clean"]:
                        command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
stdbuf -oL bash -c '{command}' >> "{log_file.name}" 2>&1""")

                        log_file.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')
                        log_file.flush()

                target["bCleaned"] = True

        common.delete_directory_recursive_safe(self.get_source_path(name))
        common.delete_file_safe(self.get_log_filepath(name))

        self.lock["targets"] = list(filter(lambda x: x["name"] != name, self.lock["targets"]))
        self.lock.save()

    def __init__(self):
        self.conf = conf.MFCConf()
        self.setup_directories()
        self.lock = lock.MFCLock()
        self.args = args.MFCArgs(self.conf)

        self.print_header()
        self.check_environment()

        if self.args["set_current"] is not None:
            common.update_symlink(f"{common.MFC_SUBDIR}/___current___", self.get_configuration_base_path(self.args["set_current"]))

        # Update symlink to current build
        if self.args["build"]:
            common.update_symlink(f"{common.MFC_SUBDIR}/___current___", self.get_configuration_base_path())

        for target_name in [ x["name"] for x in self.conf["targets"] ]:
            if target_name in self.args["targets"]:
                if self.args["build"]: self.build_target(target_name)
                if self.args["test"]:  self.test_target(target_name)
                if self.args["clean"]: self.clean_target(target_name)

        self.lock.save()
