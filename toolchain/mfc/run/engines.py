import re
import os
import time
import copy
import queue
import typing
import datetime
import subprocess
import dataclasses
import multiprocessing

from ..util.printer import cons

from .. import build
from ..util import common

from ..run import queues
from ..run import mpi_bins

from ..run.input import MFCInputFile


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

    def get_args(self) -> typing.List[str]:
        raise common.MFCException(f"MFCEngine::get_args: not implemented for {self.name}.")

    def run(self, target_name: str) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")

    def validate_job_options(self) -> None:
        raise common.MFCException(f"MFCEngine::validate_job_options: not implemented for {self.name}.")

    def get_binpath(self, target: str) -> str:
        return os.sep.join([build.get_install_dirpath(), "bin", target])


def _interactive_working_worker(cmd, q):
	q.put(common.system(cmd, no_exception=True, stdout=subprocess.DEVNULL,
                                                stderr=subprocess.DEVNULL) == 0)

class InteractiveEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Interactive", "interactive")

    def _init(self) -> None:
        self.mpibin = mpi_bins.get_binary(self.mfc.args)
        self.bWorks = False # We don't know yet whether this engine works

    def get_args(self) -> str:
        return f"""\
Nodes         (-N)  {self.mfc.args['nodes']}
CPUs (/node)  (-n)  {self.mfc.args['cpus_per_node']}
GPUs (/node)  (-g)  {self.mfc.args["gpus_per_node"]}
MPI Binary    (-b)  {self.mpibin.bin}\
"""

    def get_exec_cmd(self, target_name: str) -> typing.List[str]:
        binpath = self.get_binpath(target_name)

        if not self.mfc.args["mpi"]:
            return [binpath]

        flags = self.mfc.args["flags"][:]

        return [self.mpibin.bin] + self.mpibin.gen_params(self.mfc.args) + flags + [binpath]


    def run(self, target_name: str) -> None:
        if not self.bWorks:
            # Fix MFlowCode/MFC#21: Check whether attempting to run a job will hang
            # forever. This can happen when using the wrong queue system.

            work_timeout = 10

            cons.print(f"Ensuring the [bold magenta]Interactive Engine[/bold magenta] works ({work_timeout}s timeout):")

            q = multiprocessing.Queue()

            p = multiprocessing.Process(
                target=_interactive_working_worker,
                args=([self.mpibin.bin] +
                      self.mpibin.gen_params(self.mfc.args) +
                      [ "hostname" ], q, ))
            p.start()
            p.join(work_timeout)

            try:
                self.bWorks = q.get(block=False)
            except queue.Empty as e:
                self.bWorks = False

            if p.is_alive() or not self.bWorks:
                raise common.MFCException(
                      "The [bold magenta]Interactive Engine[/bold magenta] appears to hang or exit with a non-zero status code. "
                    + "This may indicate that the wrong MPI binary is being used to launch parallel jobs. You can specify the correct one for your system "
                    + "using the <-b,--binary> option. For example:\n"
                    + " - ./mfc.sh run <myfile.py> -b mpirun\n"
                    + " - ./mfc.sh run <myfile.py> -b srun\n"
                    + f"[bold magenta]Reason[/bold magenta]: {'Time limit.' if p.is_alive() else 'Exit code.'}"
                )


        cons.print(f"Running [bold magenta]{target_name}[/bold magenta]:")
        cons.indent()

        if not self.mfc.args["dry_run"]:
            start_time = time.monotonic()
            common.system(self.get_exec_cmd(target_name), cwd=self.input.case_dirpath)
            end_time   = time.monotonic()
            cons.print(no_indent=True)

            cons.print(f"[bold green]Done[/bold green] (in {datetime.timedelta(seconds=end_time - start_time)})")

        cons.unindent()

    def validate_job_options(self, mfc) -> None:
        if mfc.args["nodes"] != 1:
            raise common.MFCException("InteractiveEngine: Only node can be used with the interactive engine.")


class BatchEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Batch", "batch")

    def get_args(self) -> str:
        return f"""\
Nodes         (-N)  {self.mfc.args['nodes']}
CPUs (/node)  (-n)  {self.mfc.args['cpus_per_node']}
GPUs (/node)  (-g)  {self.mfc.args["gpus_per_node"]}
Walltime      (-w)  {self.mfc.args["walltime"]}
Partition     (-p)  {self.mfc.args["partition"]}
Account       (-a)  {self.mfc.args["account"]}
Email         (-@)  {self.mfc.args["email"]}
"""

    def run(self, target_name: str) -> None:
        cons.print(f"Running [bold magenta]{target_name}[/bold magenta]:")
        cons.indent()

        system = queues.get_system()
        cons.print(f"Detected the [bold magenta]{system.name}[/bold magenta] queue system.")

        self.__create_batch_file(system, target_name)

        if not self.mfc.args["dry_run"]:
            cons.print(no_indent=True)
            self.__execute_batch_file(system, target_name)
            cons.print(no_indent=True)

            cons.print("[bold yellow]INFO:[/bold yellow] Batch file submitted! Please check your queue system for the job status.")
            cons.print("[bold yellow]INFO:[/bold yellow] If an error occurs, please check the generated batch file and error logs for more information.")
            cons.print("[bold yellow]INFO:[/bold yellow] You can modify the template batch file to your needs.")

        cons.unindent()

    def __get_batch_dirpath(self) -> str:
        return copy.copy(self.input.case_dirpath)

    def __get_batch_filename(self, target_name: str) -> str:
        return f"{target_name}.sh"

    def __get_batch_filepath(self, target_name: str):
        return os.path.abspath(os.sep.join([
            self.__get_batch_dirpath(),
            self.__get_batch_filename(target_name)
        ]))

    def __generate_prologue(self, system: queues.QueueSystem,) -> str:
        modules = f""

        if common.does_system_use_modules():
            modules = f"""\
printf ":) Loading modules...\\n"

module purge
module load {' '.join(common.get_loaded_modules())}
"""

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

{modules}

cd "{self.input.case_dirpath}"

echo -e ":) Running MFC..."

t_start=$(date +%s)
"""

    def __generate_epilogue(self) -> str:
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

    def __evaluate_variable(self, var: str) -> str:
        v: str = var.strip()
        if v in self.mfc.args:
            return str(self.mfc.args[v])

        return None

    def __evaluate_expression(self, expr: str) -> str:
        expr_original = expr[:]

        evaluated = self.__evaluate_variable(expr)
        if evaluated is not None:
            if not common.isspace(evaluated):
                return evaluated
            else:
                # The expression is valid but it is empty
                return None

        # It may be an expression. Try and parse it
        for var_candidate in re.split(r"[\*,\ ,\+,\-,\/,\\,\%,\,,\.,\^,\',\",\[,\],\(,\),\=]", expr):
            evaluated = self.__evaluate_variable(var_candidate)

            if evaluated is not None and not common.isspace(evaluated):
                expr = expr.replace(var_candidate, evaluated)

        # See if it computable
        try:
            # We assume eval is safe because we control the expression.
            return str(eval(expr))
        except Exception as exc:
            raise common.MFCException(f"BatchEngine: {expr_original} (interpreted as {expr}) is not a valid expression in the template file. Please check your spelling.")

    def __batch_evaluate(self, s: str, system: queues.QueueSystem, target_name: str):
        replace_list = [
            ("{MFC::PROLOGUE}", self.__generate_prologue(system)),
            ("{MFC::EPILOGUE}", self.__generate_epilogue()),
            ("{MFC::BIN}",      self.get_binpath(target_name))
        ]

        for (key, value) in replace_list:
            s = s.replace(key, value)

        # Remove "#>" comments & redundant newlines
        s = re.sub(r"^#>.*\n", "",   s, flags=re.MULTILINE)
        s = re.sub(r"^\n{2,}", "\n", s, flags=re.MULTILINE)

        # Evaluate expressions of the form "{expression}"
        for match in re.findall(r"{[^\{]+}", s, flags=re.MULTILINE):
            repl = self.__evaluate_expression(match[1:-1])

            if repl is not None:
                s = s.replace(match, repl)
            else:
                # If not specified, then remove the line it appears on
                s = re.sub(f"^.*\{match}.*$\n", "", s, flags=re.MULTILINE)

                cons.print(f"> > [bold yellow]Warning:[/bold yellow] [magenta]{match[1:-1]}[/magenta] was not specified. Thus, any line it figures on will be discarded.")

        return s

    def __create_batch_file(self, system: queues.QueueSystem, target_name: str):
        cons.print("> > Generating batch file...")
        filepath = self.__get_batch_filepath(target_name)
        cons.print("> > Evaluating template file...")
        content = self.__batch_evaluate(system.template, system, target_name)

        cons.print("> > Writing batch file...")
        common.file_write(filepath, content)

    def __execute_batch_file(self, system: queues.QueueSystem, target_name: str):
        # We CD to the case directory before executing the batch file so that
        # any files the queue system generates (like .err and .out) are created
        # in the correct directory.
        cmd = system.gen_submit_cmd(self.__get_batch_filename(target_name))

        if common.system(cmd, cwd=self.__get_batch_dirpath()) != 0:
            raise common.MFCException(f"Submitting batch file for {system.name} failed. It can be found here: {self.__get_batch_filepath(target_name)}. Please check the file for errors.")

    def validate_job_options(self, mfc) -> None:
        if len(mfc.args["targets"]) != 1:
            raise common.MFCException(f"The Batch engine requires a unique target (-t) to run.")
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

