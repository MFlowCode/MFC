import os
import json
import rich

import internal.common      as common
import internal.input_dicts as input_dicts

class MFCRun:
    def __init__(self, bootstrap):
        self.bootstrap = bootstrap

        # Aliases
        self.console = self.bootstrap.console
        self.args    = self.bootstrap.args

    def get_case_dict(self):
        case:  dict = {}
        input: str  = self.args["i"].strip()

        self.console.print(f"> > Fetching case dir from {input}...")

        if input.endswith(".py"):
            (output, err) = common.get_py_program_output(input)

            if err != 0:
                self.console.print(f"> > Input file {input} terminated with a non-zero exit code. View the output bellow: [bold red]❌[/bold red]")
                for line in output.splitlines():
                    self.console.print(line)

                raise common.MFCException(f"> > Input file {input} terminated with a non-zero exit code. View above.")

            case = json.loads(output)
        else:
            self.console.print(f"> > Unrecognized input file format for '{input}'. Please check the extension. [bold red]✘[/bold red]")
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
        dirpath  = os.path.abspath(os.path.dirname(self.args["i"]))
        filename = f"{target_name}.inp"
        common.file_write(f"{dirpath}/{filename}", contents)

    def detect_queue_system(self) -> str:
        SYSTEMS = {"PBS":   ["echo"], #FIXME: #qsub
                   "SGE":   ["qsub"], #FIXME: hmm
                   "SLURM": ["sbatch"]}
        for system,cmds in SYSTEMS.items():
            for cmd in cmds:
                if 0 == os.system(f"{cmd} -h > /dev/null 2>&1"):
                    self.console.print(f"Detected the [bold magenta]{system}[/bold magenta] queueing system.")
                    return system
        
        raise common.MFCException(f"Failed to detect a queueing system.")

    def get_binpath(self, target: str) -> str:
        return f'{self.bootstrap.get_build_path(target)}/bin/{target}'

    def get_ld(self) -> str:
        return f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{common.MFC_SUBDIR}/common/build/lib"'

    def get_case_dirpath(self) -> str:
        return os.path.abspath(os.path.dirname(self.args["i"]))

    def create_batch_file(self, system: str, target: str):
        case_dirpath = self.get_case_dirpath()

        if system == "PBS":
            BATCH_CONTENT: str = f"""\
#!/bin/sh -l
#PBS -l nodes={self.args["N"]}:ppn={self.args["n"]}
#PBS -l walltime={self.args["w"]}
#PBS -q {self.args["p"]}
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

status_msg="{rich.pri}{colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL}"
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
        cc       = self.args["c"]
        input    = self.args["i"].strip()
        engine   = self.args["e"]
        targets  = self.args["t"]
        nodes    = self.args["N"]
        tasks_pn = self.args["n"]

        if targets[0] == "mfc":
            targets = ["pre_process", "simulation", "post_process"]

            self.console.print("[bold][u]Run:[/u][/bold]")

        self.console.print(f"> Targets       (-t)  {', '.join(targets)}")
        self.console.print(f"> Engine        (-e)  {engine}")
        self.console.print(f"> Config        (-c)  {cc}")
        self.console.print(f"> Input         (-i)  {input}")
        self.console.print(f"> Nodes         (-N)  {nodes}")
        self.console.print(f"> Tasks (/node) (-n)  {tasks_pn}")

        for target_name in targets:
            self.console.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")
            
            if not self.bootstrap.is_build_satisfied(target_name):
                self.console.print(f"> > Target {target_name} needs (re)building...")
                self.console.print(f"> > Please (re)build target {target_name}.")
                raise common.MFCException(f"Can't run mfc. {target_name} needs rebuilding.")

            self.create_input_file(target_name, self.get_case_dict())

            if engine == 'serial':
                date = f"> > [bold cyan][{common.get_datetime_str()}][/bold cyan]"
                bin  = self.get_binpath(target_name)
                
                cd   = f'cd {self.get_case_dirpath()}'
                ld   = self.get_ld()
                exec = f'mpiexec -np {tasks_pn} "{bin}"'

                self.console.print(f"{date} Running...")
                common.execute_shell_command(f"{cd} && {ld} {exec}")

                self.console.print(f"> > Done [bold green]✓[/bold green]")
            elif engine == 'parallel':
                queue_sys = self.detect_queue_system()

                self.create_batch_file(queue_sys, target_name)

                self.execute_batch_file(queue_sys)
            else:
                raise common.MFCException(f"Unsupported engine {engine}.")

