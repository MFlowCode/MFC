import os, json, rich, typing, time, dataclasses

from queue    import Queue
from datetime import timedelta

import common, case_dicts

@dataclasses.dataclass
class QueueSystem:
    name: str

    def is_active(self) -> bool:
        raise common.MFCException("QueueSystem::is_active: not implemented.")

    def gen_batch_header(self, args: dict, target_name: str) -> str:
        raise common.MFCException("QueueSystem::gen_batch_header: not implemented.")

    def gen_submit_cmd(self, filename: str) -> None:
        raise common.MFCException("QueueSystem::gen_submit_cmd: not implemented.")


class PBSSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("PBS")

    def is_active(self) -> bool:
        return 0 == os.system(f"qsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#PBS -A {args["account"]}
#PBS -l nodes={args["nodes"]}:ppn={args["cpus_per_node"]}
#PBS -l walltime={args["walltime"]}
#PBS -q {args["partition"]}
#PBS -N {job_name}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"qsub {filename}"


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")

    def is_active(self) -> bool:
        return 0 == os.system(f"bsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#BSUB -P {args["account"]}
#BSUB -W {args["walltime"][:-3]}
#BSUB -J {job_name}
#BSUB -nnodes {args["nodes"]}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"bsub {filename}"


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")

    def is_active(self) -> bool:
        return 0 == os.system(f"sbatch -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#SBATCH --job-name="{job_name}"
#SBATCH --time={args["walltime"]}
#SBATCH --nodes={args["nodes"]}
#SBATCH --ntasks-per-node={args["cpus_per_node"]}
#SBATCH --cpus-per-task={1}
#SBATCH --partition={args["partition"]}
#SBATCH --account={args["account"]}
#SBATCH --mail-user={args["email"]}
#SBATCH --mail-type=BEGIN, END, FAIL
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"sbatch {filename}"


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]


@dataclasses.dataclass
class Engine:
    name: str
    slug: str

    def run(self, mfc, target_name: str) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")


class SerialEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Serial", "serial")

    def run(self, mfc, target_name: str) -> None:
        self.mfc = mfc

        date = f"> > [bold cyan][{common.get_datetime_str()}][/bold cyan]"
        rich.print(f"{date} Running...")

        start_time = time.monotonic()
        common.execute_shell_command(self.mfc.run.get_exec_cmd(target_name))
        end_time   = time.monotonic()

        rich.print(f"> > Done [bold green]✓[/bold green] (in {timedelta(seconds=end_time - start_time)})")


class ParallelEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Parallel", "parallel")

    def run(self, mfc, target_name: str) -> None:
        self.mfc = mfc

        queue_sys = self.mfc.run.detect_queue_system()

        self.create_batch_file(queue_sys, target_name)

        self.execute_batch_file(queue_sys, target_name)

        self.remove_batch_file(queue_sys, target_name)

    def get_batch_filepath(self, system: QueueSystem, target_name: str):
        case_dirpath = self.mfc.run.get_case_dirpath()

        return os.path.abspath(f"{case_dirpath}/{target_name}.sh")

    def create_batch_file(self, system: QueueSystem, target_name: str):
        job_name = self.mfc.run.get_job_name(target_name)

        BATCH_CONTENT: str = f"""\
#!/usr/bin/env bash
{system.gen_batch_header(self.mfc.args, job_name)}

RED="\\u001b[31m";   CYAN="\\u001b[36m";   BLUE="\\u001b[34m";    WHITE="\\u001b[37m"
GREEN="\\u001b[32m"; YELLOW="\\u001b[33m"; MAGENTA="\\u001b[35m"; COLOR_RESET="\\033[m"

TABLE_FORMAT_LINE="| - %-14s %-25s - %-14s %-25s |\\n"
TABLE_HEADER="/‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\ \\n"
TABLE_FOOTER="\_______________________________________________________________________________________/\\n"
TABLE_TITLE_FORMAT="| %8s $MAGENTA%-51s$COLOR_RESET                          |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time:"    "$(date +%T)"                       "Start-date:"    "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{self.mfc.args["partition"]}"      "Walltime:"      "{self.mfc.args["walltime"]}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{self.mfc.args["account"]}"        "Nodes:"         "{self.mfc.args["nodes"]}")
$(printf "$TABLE_FORMAT_LINE" "CPUs (/node):"  "{self.mfc.args["cpus_per_node"]}"  "GPUs (/node):"  "{self.mfc.args["gpus_per_node"]}")
$(printf "$TABLE_FORMAT_LINE" "Input File:"    "{self.mfc.args["input"]}"          "Engine"         "{self.mfc.args["engine"]}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{system.name}"                     "Mode:"          "{self.mfc.args["mode"]}")
$(printf "$TABLE_FORMAT_LINE" "Email:"         "{self.mfc.args["email"]}"          "Job Name:"      "{job_name}")\\n
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Starting" "{job_name}"
printf "$TABLE_CONTENT"
printf "$TABLE_FOOTER"

t_start=$(date +%s)

echo ""

{self.mfc.run.get_exec_cmd(target_name)}

echo ""

code=$?

status_msg="SUCCESS"
if [ "$code" -ne "0" ]; then
    status_msg="FAILURE"
fi

t_stop="$(date +%s)"

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished" "{job_name}"
printf "$TABLE_FORMAT_LINE" "Total-time:"  "$(expr $t_stop - $t_start)s"  "Status:"    "$status_msg"
printf "$TABLE_CONTENT"
printf "$TABLE_FORMAT_LINE" "End-time:"    "$(date +%T)"                  "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

printf "\\nThank you for using MFC !\\n\\n"

exit $code
"""

        filepath = self.get_batch_filepath(system, target_name)

        common.file_write(filepath, BATCH_CONTENT)

    def remove_batch_file(self, system: QueueSystem, target_name: str):
        os.remove(self.get_batch_filepath(system, target_name))

    def execute_batch_file(self, system: QueueSystem, target_name: str):
        if 0 != os.system(system.gen_submit_cmd(self.get_batch_filepath(system, target_name))):
            raise common.MFCException(f"Running batch file for {system.name} failed.")


ENGINES = [ SerialEngine(), ParallelEngine() ]


class MFCRun:
    def __init__(self, mfc):
        self.mfc = mfc

    def get_case_dict(self):
        case:  dict = {}
        input: str  = self.mfc.args["input"].strip()

        rich.print(f"> > Fetching case dictionary from {input}...")

        if input.endswith(".py"):
            (output, err) = common.get_py_program_output(input)

            if err != 0:
                rich.print(f"> > Input file {input} terminated with a non-zero exit code. View the output bellow: [bold red]❌[/bold red]")
                for line in output.splitlines():
                    rich.print(line)

                raise common.MFCException(f"> > Input file {input} terminated with a non-zero exit code. View above.")

            case = json.loads(output)
        else:
            rich.print(f"> > Unrecognized input file format for '{input}'. Please check the extension. [bold red]✘[/bold red]")
            raise common.MFCException("Unrecognized input file format.")

        return case

    def get_job_name(self, target_name: str):
        return f'MFC-{str(self.mfc.args["name"])}-{target_name}'

    def get_input_filepath(self, target_name: str):
        dirpath  = os.path.abspath(os.path.dirname(self.mfc.args["input"]))
        filename = f"{target_name}.inp"

        return f"{dirpath}/{filename}"

    def create_input_file(self, target_name: str, case_dict: dict):
        MASTER_KEYS: list = case_dicts.get_input_dict_keys(target_name)

        # Create Fortran-style input file content string
        dict_str = ""
        for key,val in case_dict.items():
            if key in MASTER_KEYS:
                dict_str += f"{key} = {val}\n"

        contents = f"&user_inputs\n{dict_str}&end/\n"

        # Save .inp input file
        common.file_write(self.get_input_filepath(target_name), contents)

    def remove_input_file(self, target_name: str):
        os.remove(self.get_input_filepath(target_name))

    def detect_queue_system(self) -> str:
        for system in QUEUE_SYSTEMS:
            if system.is_active():
                rich.print(f"> > Detected the [bold magenta]{system.name}[/bold magenta] queue system.")
                return system

        raise common.MFCException(f"Failed to detect a queue system.")

    def get_binpath(self, target: str) -> str:
        return f'{self.mfc.build.get_build_path(target)}/bin/{target}'

    def get_ld(self) -> str:
        return f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{common.MFC_SUBDIR}/common/build/lib"'

    def get_case_dirpath(self) -> str:
        return os.path.abspath(os.path.dirname(self.mfc.args["input"]))

    def get_exec_cmd(self, target_name: str):
        bin = self.get_binpath(target_name)

        cd = f'cd "{self.get_case_dirpath()}"'
        ld = self.get_ld()

        def cmd_exists(cmd: str):
            return os.system(f"{cmd} -h > /dev/null 2>&1") == 0

        np = self.mfc.args["cpus_per_node"]*self.mfc.args["nodes"]

        if cmd_exists("jsrun"):
            # ORNL Summit: https://docs.olcf.ornl.gov/systems/summit_user_guide.html?highlight=lsf#launching-a-job-with-jsrun
                        
            if self.mfc.args["cpus_per_node"] != self.mfc.args["gpus_per_node"] \
               and self.mfc.args["gpus_per_node"] != 0:
               raise common.MFCException("JSRUN: Conflicting job execution parameters. If using GPUs, CPUs per node and GPUs per node must match.")

            # One resource set per CPU(Core)/GPU pair.
            rs=np
            cpus_per_rs=1
            gpus_per_rs=min(self.mfc.args["gpus_per_node"], 1)
            tasks_per_rs=1

            options = f'--smpiargs="-gpu" --nrs{rs} --cpu_per_rs{cpus_per_rs} --gpu_per_rs{gpus_per_rs} --tasks_per_rs{tasks_per_rs}'

            return f'{cd} && {ld} jsrun {options} "{bin}"'
        elif cmd_exists("srun"):
            options = f'-N {self.mfc.args["nodes"]} -n {self.mfc.args["cpus_per_node"]} -G {self.mfc.args["gpus_per_node"]}'
            
            if not self.mfc.args["account"].isspace():
                options += f' -A "{self.mfc.args["account"]}"'

            if not self.mfc.args["partition"].isspace():
                options += f' -p "{self.mfc.args["partition"]}"'

            return f'{cd} && {ld} srun {options} "{bin}"'
        elif cmd_exists("mpiexec"):            
            return f'{cd} && {ld} mpiexec -np {np} "{bin}"'
        elif cmd_exists("mpirun"):            
            return f'{cd} && {ld} mpirun -np {np} "{bin}"'
        else:
            raise common.MFCException("Not program capable of running an MPI program could be located.")

    def run(self):
        targets = self.mfc.args["targets"]
        if targets[0] == "mfc":
            targets = ["pre_process", "simulation", "post_process"]

        rich.print(f"""\
[bold][u]Run:[/u][/bold]
> Input               {self.mfc.args['input']}
> Job Name      (-#)  {self.mfc.args['name']}
> Engine        (-e)  {self.mfc.args['engine']}
> Mode          (-m)  {self.mfc.args['mode']}
> Targets       (-t)  {self.mfc.args['targets']}
> Nodes         (-N)  {self.mfc.args['nodes']}
> CPUs (/node)  (-n)  {self.mfc.args['cpus_per_node']}
> GPUs (/node)  (-g)  {self.mfc.args["gpus_per_node"]}
> Walltime      (-w)  {self.mfc.args["walltime"]}
> Partition     (-p)  {self.mfc.args["partition"]}
> Account       (-a)  {self.mfc.args["account"]}
> Email         (-@)  {self.mfc.args["email"]}
""")

        for target_name in targets:
            rich.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")

            if not self.mfc.build.is_built(target_name):
                rich.print(f"> > Target {target_name} needs (re)building...")
                self.mfc.build.build_target(target_name)

            self.create_input_file(target_name, self.get_case_dict())

            engine_slug = self.mfc.args["engine"]
            engine: Engine = None
            for candidate in ENGINES:
                candidate: Engine

                if candidate.slug == engine_slug:
                    engine = candidate
                    break

            if engine == None:
                raise common.MFCException(f"Unsupported engine {engine_slug}.")

            engine.run(self.mfc, target_name)

            self.remove_input_file(target_name)
