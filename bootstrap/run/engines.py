import re
import os
import time
import datetime
import dataclasses

import rich

import common

import run.queues   as queues
import run.mpi_bins as mpi_bins

from run.input import MFCInputFile

@dataclasses.dataclass
class Engine:
    name: str
    slug: str

    def init(self, mfc, input: MFCInputFile) -> None:
        self.mfc   = mfc
        self.input = input

        self._init()

    def _init(self) -> None:
        pass 

    def get_args(self) -> str:
        raise common.MFCException(f"MFCEngine::get_args: not implemented for {self.name}.")

    def get_targets(self, targets: list) -> list:
        raise common.MFCException(f"MFCEngine::get_targets: not implemented for {self.name}.")

    def run(self, target_name: str) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")

    def validate_job_options(self) -> None:
        raise common.MFCException(f"MFCEngine::validate_job_options: not implemented for {self.name}.")

    def get_binpath(self, target: str) -> str:
        return f'{self.mfc.build.get_build_path(target)}/bin/{target}'


class InteractiveEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Interactive", "interactive")

    def _init(self) -> None:
        self.mpibin = mpi_bins.get_binary(self.mfc.args)

    def get_targets(self, targets: list) -> list:
        if targets[0] == "mfc":
            return ["pre_process", "simulation", "post_process"]

        return targets

    def get_args(self) -> str:
        return f"""\
> MPI Binary    (-b)  {self.mpibin.bin} {f"[green](autodetect: {self.mpibin.name})[/green]" if self.mfc.args["binary"] == None else f"[yellow](override: {self.mpibin.name})[/yellow]"}
"""

    def get_exec_cmd(self, target_name: str):
        binpath = self.get_binpath(target_name)

        cd = f'cd "{self.input.case_dirpath}"'

        flags = ""
        for flag in self.mfc.args["flags"]:
            flags += f"\"{flag}\" "

        exec_params = self.mpibin.gen_params(self.mfc.args)

        return f'{cd} && {self.mpibin.bin} {exec_params} {flags} "{binpath}"'

    def run(self, target_name: str) -> None:
        date = f"> > [bold cyan][{common.get_datetime_str()}][/bold cyan]"
        rich.print(f"{date} Running...")

        start_time = time.monotonic()
        common.execute_shell_command(self.get_exec_cmd(target_name))
        end_time   = time.monotonic()

        rich.print(f"> > Done [bold green]✓[/bold green] (in {datetime.timedelta(seconds=end_time - start_time)})")

    def validate_job_options(self, mfc) -> None:
        if mfc.args["nodes"] != 1:
            raise common.MFCException("InteractiveEngine: Only node can be used with the interactive engine.")


class BatchEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Batch", "batch")

    def get_targets(self, targets: list) -> list:
        if targets[0] == "mfc" or len(targets) != 1:
            raise common.MFCException("The batch engine requires a unique target to run. Please specify -t <target_name>.")

        return targets

    def get_args(self) -> str:
        return f"""\
> Nodes         (-N)  {self.mfc.args['nodes']}
> CPUs (/node)  (-n)  {self.mfc.args['cpus_per_node']}
> GPUs (/node)  (-g)  {self.mfc.args["gpus_per_node"]}
> Walltime      (-w)  {self.mfc.args["walltime"]}
> Partition     (-p)  {self.mfc.args["partition"]}
> Account       (-a)  {self.mfc.args["account"]}
> Email         (-@)  {self.mfc.args["email"]}
"""

    def run(self, target_name: str) -> None:
        system = queues.get_system()
        rich.print(f"> > Detected the [bold magenta]{system.name}[/bold magenta] queue system.")

        self.create_batch_file(system, target_name)

        self.execute_batch_file(system, target_name)

    def get_batch_filepath(self, target_name: str):
        case_dirpath = self.input.case_dirpath

        return os.path.abspath(f"{case_dirpath}/{target_name}.sh")

    def generate_prologue(self, system: queues.QueueSystem,) -> str:
        return f"""\
TABLE_FORMAT_LINE="| - %-14s %-35s - %-14s %-35s |\\n"
TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_TITLE_FORMAT="| %8s %-96s |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time:"    "$(date +%T)"                       "Start-date:"    "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{self.mfc.args["partition"]}"      "Walltime:"      "{self.mfc.args["walltime"]}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{self.mfc.args["account"]}"        "Nodes:"         "{self.mfc.args["nodes"]}")
$(printf "$TABLE_FORMAT_LINE" "CPUs (/node):"  "{self.mfc.args["cpus_per_node"]}"  "GPUs (/node):"  "{self.mfc.args["gpus_per_node"]}")
$(printf "$TABLE_FORMAT_LINE" "Job Name:"      "{self.mfc.args["name"]}"           "Engine"         "{self.mfc.args["engine"]}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{system.name}"                     "Mode:"          "{self.mfc.args["mode"]}")
$(printf "$TABLE_FORMAT_LINE" "Email:"         "{self.mfc.args["email"]}"          ""               "")
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Starting" "{self.mfc.args["name"]} from {self.mfc.args["input"]}:"
printf "$TABLE_CONTENT\\n"
printf "$TABLE_FOOTER\\n"

cd "{self.input.case_dirpath}"

echo -e ":) Running MFC..."

t_start=$(date +%s)
"""

    def generate_epilogue(self) -> str:
        return f"""\
code=$?

t_stop="$(date +%s)"

printf "\\n$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished" "{self.mfc.args["name"]}:"
printf "$TABLE_FORMAT_LINE" "Total-time:"  "$(expr $t_stop - $t_start)s"  "Exit Code:" "$code"
printf "$TABLE_FORMAT_LINE" "End-time:"    "$(date +%T)"                  "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

exit $code
"""

    def evaluate_variable(self, var: str) -> str:
        v: str = var.strip()
        if v in self.mfc.args:
            return str(self.mfc.args[v])
        
        return None

    def evaluate_expression(self, expr: str) -> str:
        expr_original = expr[:]

        evaluated = self.evaluate_variable(expr)
        if evaluated is not None:
            if not common.isspace(evaluated):
                return evaluated
            else:
                # The expression is valid but it is empty
                return None

        # It may me a calculation. Try and parse it
        for var_candidate in re.split(r"[\*,\ ,\+,\-,\/,\\,\%,\,,\.,\^,\',\",\[,\],\(,\),\=]", expr):
            evaluated = self.evaluate_variable(var_candidate)

            if evaluated is not None and not common.isspace(evaluated):                
                expr = expr.replace(var_candidate, evaluated)
        
        # See if it computable
        try:
            # We assume eval is safe because we control the expression.
            return str(eval(expr))
        except Exception as exc:
            raise common.MFCException(f"BatchEngine: {expr_original} (interpreted as {expr}) is not a valid expression in the template file. Please check your spelling.")

    def batch_evaluate(self, s: str, system: queues.QueueSystem, target_name: str):
        replace_list = [
            ("{MFC::PROLOGUE}", self.generate_prologue(system)),
            ("{MFC::EPILOGUE}", self.generate_epilogue()),
            ("{MFC::BIN}",      self.get_binpath(target_name))
        ]

        for (key, value) in replace_list:
            s = s.replace(key, value)

        # Remove "#>" comments & redundant newlines
        s = re.sub(r"^#>.*\n", "",   s, flags=re.MULTILINE)
        s = re.sub(r"^\n{2,}", "\n", s, flags=re.MULTILINE)

        # Evaluate expressions of the form "{expression}"
        for match in re.findall(r"{[^\{]+}", s, flags=re.MULTILINE):
            repl = self.evaluate_expression(match[1:-1])

            if repl is not None:
                s = s.replace(match, repl)
            else:
                # If not specified, then remove the line it appears on
                s = re.sub(f"^.*\{match}.*$\n", "", s, flags=re.MULTILINE)

                rich.print(f"[bold yellow]Warning:[/bold yellow] [magenta]{match[1:-1]}[/magenta] was not specified. Thus, any line it figures on will be discarded.")

        return s

    def create_batch_file(self, system: queues.QueueSystem, target_name: str):
        common.file_write(self.get_batch_filepath(target_name),
                          self.batch_evaluate(system.template, system, target_name))

    def execute_batch_file(self, system: queues.QueueSystem, target_name: str):
        if 0 != os.system(system.gen_submit_cmd(self.get_batch_filepath(target_name))):
            raise common.MFCException(f"Running batch file for {system.name} failed. It can be found here: {self.get_batch_filepath(target_name)}. Please check the file for errors.")

    def validate_job_options(self, mfc) -> None:
        pass


ENGINES = [ InteractiveEngine(), BatchEngine() ]

def get_engine(slug: str) -> Engine:    
    engine: Engine = None
    for candidate in ENGINES:
        candidate: Engine

        if candidate.slug == slug:
            engine = candidate
            break

    if engine == None:
        raise common.MFCException(f"Unsupported engine {slug}.")

    return engine
