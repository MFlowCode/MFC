import re

import rich

import run.engines  as engines
import run.input    as input

import common


class MFCRun:
    def __init__(self, mfc):
        self.mfc = mfc

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

        # Load input file
        input_file = input.load(self.mfc.args["input"].strip())

        engine = engines.get_engine(self.mfc.args["engine"])
        engine.init(self.mfc, input_file)

        rich.print(f"""\
[bold][u]Run:[/u][/bold]
> Input               {self.mfc.args['input']}
> Job Name      (-#)  {self.mfc.args['name']}
> Engine        (-e)  {self.mfc.args['engine']}
> Mode          (-m)  {self.mfc.args['mode']}
> Targets       (-t)  {self.mfc.args['targets']}
{engine.get_args()}\
""")

        self.validate_job_options()

        for target_name in engine.get_targets(self.mfc.args["targets"]):
            rich.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")

            if not self.mfc.build.is_built(target_name):
                rich.print(f"> > Target {target_name} needs (re)building...")
                self.mfc.build.build_target(target_name, "> > > ")
                        
            # Create input file
            input_file.dump(target_name)

            engine.run(target_name)
