import os
import json
import colorama

import internal.common      as common
import internal.input_dicts as input_dicts

class MFCRun:
    def __init__(self, bootstrap):
        self.bootstrap = bootstrap

        # Aliases
        self.tree = self.bootstrap.tree
        self.args = self.bootstrap.args

    def get_case_dict(self):
        case:  dict = {}
        input: str  = self.args["input"].strip()

        self.tree.print(f"Fetching case dir from {input}...")

        if input.endswith(".py"):
            (output, err) = common.get_py_program_output(input)

            if err != 0:
                self.tree.print(f"Input file {input} terminated with a non-zero exit code. View the output bellow: ({colorama.Fore.RED}ERROR{colorama.Style.RESET_ALL})")
                for line in output.splitlines():
                    self.tree.print(line)

                raise common.MFCException(f"Input file {input} terminated with a non-zero exit code. View above.")

            case = json.loads(output)
        else:
            self.tree.print(f"Unrecognized input file format for '{input}'. Please check the extension. ({colorama.Fore.RED}ERROR{colorama.Style.RESET_ALL})")
            raise common.MFCException("Unrecognized input file format.")
        
        return case

    def create_input_file(self, target_name: str, case_dict: dict):
        MASTER_KEYS: list = input_dicts.get_keys(target_name)

        # Create Fortran-style input file content string
        dict_str = ""
        for key,val in case_dict.items():
            if key in MASTER_KEYS:
                dict_str += f"{key} = {val}\n"
        
        contents = f"&user_inputs\n{dict_str}&end/"

        # Save .inp input file
        dirpath  = os.path.abspath(os.path.dirname(self.args["input"]))
        filename = f"{target_name}.inp"
        common.file_write(f"{dirpath}/{filename}", contents)

    def detect_queue_system(self) -> str:
        SYSTEMS = {"PBS":   ["echo"], #FIXME: #qsub
                   "SGE":   ["qsub"], #FIXME: hmm
                   "SLURM": ["sbatch"]}
        for system,cmds in SYSTEMS.items():
            for cmd in cmds:
                if 0 == os.system(f"{cmd} -h > /dev/null 2>&1"):
                    self.tree.print(f"Detected the {colorama.Fore.MAGENTA}{system}{colorama.Style.RESET_ALL} queueing system.")
                    return system
        
        raise common.MFCException(f"Failed to detect a queueing system.")

    def get_binpath(self, target: str) -> str:
        return f'{self.bootstrap.get_build_path(target)}/bin/{target}'

    def get_ld(self) -> str:
        return f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{common.MFC_SUBDIR}/common/build/lib"'

    def get_case_dirpath(self) -> str:
        return os.path.abspath(os.path.dirname(self.args["input"]))

    def create_batch_file(self, system: str, target: str):
        case_dirpath = self.get_case_dirpath()

        if system == "PBS":
            BATCH_CONTENT: str = f"""\
#!/bin/sh -l
#PBS -l nodes={self.args["nodes"]}:ppn={self.args["tasks_per_node"]}
#PBS -l walltime={self.args["walltime"]}
#PBS -q {self.args["partition"]}
#PBS -N {target}

echo "================================================="
echo "| Starting job #{target}"
echo "| - Start-date: `date +%D`"
echo "| - Start-time: `date +%T`"
echo "================================================="

t_start=$(date +%s)

{self.get_ld()} \\
    mpiexec "{self.get_binpath(target)}"

code=$?

status_msg="{colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL}"
if [ "$code" -ne "0" ]; then
    status_msg="{colorama.Fore.RED}FAILED{colorama.Style.RESET_ALL}"
fi

t_stop=$(date +%s)

echo "================================================="
echo "| Finished job {target}: $status_msg"
echo "| - End-date: `date +%D`"
echo "| - End-time: `date +%T`"
echo "| - Total-time: $(expr $t_stop - $t_start)s"
echo "================================================="

exit $code
"""

            common.file_write(f"{case_dirpath}/{target}.sh", BATCH_CONTENT)
        else:
            raise common.MFCException(f"Can't create batch file for {system}.")

    def execute_batch_file(self, system: str):
        if system == "PBS":
            # handle
            pass
        else:
            raise common.MFCException(f"Running batch file for {system} is not supported.")

    def run(self):
        cc       = self.args["compiler_configuration"]
        input    = self.args["input"].strip()
        engine   = self.args["engine"]
        targets  = self.args["targets"]
        tasks_pn = self.args["tasks_per_node"]

        if targets[0] == "mfc":
            targets = ["pre_process", "simulation", "post_process"]

        self.tree.print(f"Running MFC:")
        self.tree.indent()

        self.tree.print(f"Target(s)     (-t)  {', '.join(targets)}")
        self.tree.print(f"Engine        (-e)  {engine}")
        self.tree.print(f"Config        (-cc) {cc}")
        self.tree.print(f"Input         (-i)  {input}")
        self.tree.print(f"Tasks (/node) (-n)  {tasks_pn}")

        for target_name in targets:
            self.tree.print(f"Running {colorama.Fore.MAGENTA}{target_name}{colorama.Style.RESET_ALL}:")
            self.tree.indent()
            
            if not self.bootstrap.is_build_satisfied(target_name):
                self.tree.print(f"Target {target_name} needs (re)building...")
                self.build_target(target_name)

            self.create_input_file(target_name, self.get_case_dict())

            if engine == 'serial':
                date = f"{colorama.Fore.CYAN}[{common.get_datetime_str()}]{colorama.Style.RESET_ALL}"
                bin  = self.get_binpath(target_name)
                
                cd   = f'cd {self.get_case_dirpath()}'
                ld   = self.get_ld()
                exec = f'mpiexec -np {tasks_pn} "{bin}"'

                self.tree.print(f"{date} Running...")
                common.execute_shell_command(f"{cd} && {ld} {exec}")

                self.tree.print(f"Done. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})")
            elif engine == 'parallel':
                queue_sys = self.detect_queue_system()

                self.create_batch_file(queue_sys, target_name)

                self.execute_batch_file(queue_sys)
            else:
                raise common.MFCException(f"Unsupported engine {engine}.")

            self.tree.unindent()

        self.tree.unindent()
