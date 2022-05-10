import os, re, json, rich, typing, time, dataclasses

import run.engines as engines

import common, run.case_dicts as case_dicts

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
            if os.system(f"{cmd} --help > /dev/null 2>&1") == 0:
                return True

            if os.system(f"{cmd} --h > /dev/null 2>&1") == 0:
                return True

            return False

        np = self.mfc.args["cpus_per_node"]*self.mfc.args["nodes"]

        options = ""
        if self.mfc.args["engine"] == "serial":
            for flag in self.mfc.args["flags"]:
                options += f"\"{flag}\" "

        if cmd_exists("jsrun"):
            # ORNL Summit: https://docs.olcf.ornl.gov/systems/summit_user_guide.html?highlight=lsf#launching-a-job-with-jsrun
            # We create one resource-set per CPU(Core)/GPU pair.
            rs=np
            cpus_per_rs=1
            gpus_per_rs=min(self.mfc.args["gpus_per_node"], 1)
            tasks_per_rs=1

            options += f'--smpiargs="-gpu" --nrs {rs} --cpu_per_rs {cpus_per_rs} --gpu_per_rs {gpus_per_rs} --tasks_per_rs {tasks_per_rs}'

            return f'{cd} && {ld} jsrun {options} "{bin}"'
        elif cmd_exists("srun"):
            options += f' -n {self.mfc.args["cpus_per_node"]}'

            if self.mfc.args["nodes"] != 1:
                options += f' -N {self.mfc.args["nodes"]}'

            # MFC binds its GPUs on its own, as long as they have been allocated
            # by the system's scheduler, or are present on your local machine,
            # if running in serial mode.
            #
            # if self.mfc.args["gpus_per_node"] != 0:
            #    options += f' -G {self.mfc.args["gpus_per_node"]}'


            if not common.isspace(self.mfc.args["account"]):
                options += f' -A "{self.mfc.args["account"]}"'

            if not common.isspace(self.mfc.args["partition"]):
                options += f' -p "{self.mfc.args["partition"]}"'

            return f'{cd} && {ld} srun {options} "{bin}"'
        elif cmd_exists("mpiexec"):
            options += f" -np {np}"

            return f'{cd} && {ld} mpiexec {options} "{bin}"'
        elif cmd_exists("mpirun"):
            options += f" -np {np}"
            
            return f'{cd} && {ld} mpirun {options} "{bin}"'
        else:
            raise common.MFCException("Not program capable of running an MPI program could be located.")


    def validate_job_options(self) -> None:
        if self.mfc.args["cpus_per_node"] != self.mfc.args["gpus_per_node"] \
            and self.mfc.args["gpus_per_node"] != 0:
            raise common.MFCException("RUN: Conflicting job execution parameters. If using GPUs, CPUs per node and GPUs per node must match.")

        if self.mfc.args["nodes"] <= 0:
            raise common.MFCException("RUN: At least one node must be requested.")

        if self.mfc.args["cpus_per_node"] <= 0:
            raise common.MFCException("RUN: At least one CPU per node must be requested.")

        if not common.isspace(self.mfc.args["email"]):
            # https://stackoverflow.com/questions/8022530/how-to-check-for-valid-email-address
            if not re.match(r"\"?([-a-zA-Z0-9.`?{}]+@\w+\.\w+)\"?", self.mfc.args["email"]):
                raise common.MFCException(f'RUN: {self.mfc.args["email"]} is not a valid e-mail address.')

        engines.get_engine(self.mfc.args["engine"]).validate_job_options(self.mfc)


    def run(self) -> None:
        if len(self.mfc.args["targets"]) == 0:
            rich.print(f"> No target selected.")
            return

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

        self.validate_job_options()

        engine = engines.get_engine(self.mfc.args["engine"])
        
        for target_name in engine.get_targets(self.mfc.args["targets"]):
            rich.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")

            if not self.mfc.build.is_built(target_name):
                rich.print(f"> > Target {target_name} needs (re)building...")
                self.mfc.build.build_target(target_name, "> > > ")

            self.create_input_file(target_name, self.get_case_dict())

            engine.run(self.mfc, target_name)

            self.remove_input_file(target_name)
