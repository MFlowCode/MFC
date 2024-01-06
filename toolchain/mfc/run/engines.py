import re, os, time, copy, typing, datetime, subprocess, dataclasses, multiprocessing

from ..state     import ARG, ARGS
from ..printer   import cons
from ..build     import MFCTarget, SYSCHECK, get_targets
from ..common    import MFCException, does_command_exist, isspace, system
from ..common    import format_list_to_string, does_system_use_modules
from ..common    import get_loaded_modules, file_write
from ..run       import queues, mpi_bins
from ..run.input import MFCInputFile


def profiler_prepend():
    if ARG("ncu") is not None:
        if not does_command_exist("ncu"):
            raise MFCException("Failed to locate [bold green]NVIDIA Nsight Compute[/bold green] (ncu).")

        return ["ncu", "--nvtx", "--mode=launch-and-attach",
                       "--cache-control=none", "--clock-control=none"] + ARG("ncu")

    if ARG("nsys") is not None:
        if not does_command_exist("nsys"):
            raise MFCException("Failed to locate [bold green]NVIDIA Nsight Systems[/bold green] (nsys).")

        return ["nsys", "profile", "--stats=true", "--trace=mpi,nvtx,openacc"] + ARG("nsys")

    return []


@dataclasses.dataclass(init=False)
class Engine:
    name:  str
    slug:  str
    input: MFCInputFile

    def __init__(self, name: str, slug: str) -> None:
        self.name = name
        self.slug = slug

    def init(self, _input: MFCInputFile) -> None:
        self.input = _input

        self._init()

    def _init(self) -> None:
        pass

    def get_args(self) -> typing.List[str]:
        raise MFCException(f"MFCEngine::get_args: not implemented for {self.name}.")

    def run(self, targets: typing.List[MFCTarget]) -> None:
        raise MFCException(f"MFCEngine::run: not implemented for {self.name}.")


def _interactive_working_worker(cmd: typing.List[str], q: multiprocessing.Queue):
    """ Runs a command and puts the result in a queue. """
    cmd = [ str(_) for _ in cmd ]
    cons.print(f"$ {' '.join(cmd)}")
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    q.put(result)

class InteractiveEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Interactive", "interactive")

    # pylint: disable=attribute-defined-outside-init
    def _init(self) -> None:
        self.mpibin = mpi_bins.get_binary()

        # If using MPI, we don't know yet whether this engine works
        self.bKnowWorks = not ARG("mpi")

    def get_args(self) -> str:
        return f"""\
Nodes         (-N)  {ARG('nodes')}
Tasks (/node) (-n)  {ARG('tasks_per_node')}
MPI Binary    (-b)  {self.mpibin.bin}\
"""


    def __get_exec_cmd(self, target: MFCTarget) -> typing.List[str]:
        cmd = []
        if ARG("mpi"):
            cmd += [self.mpibin.bin] + self.mpibin.gen_params() + ARG("--")[:]

        cmd += profiler_prepend()

        cmd.append(target.get_install_binpath())

        return cmd

    def run(self, targets) -> None:
        targets = get_targets(targets)

        if not self.bKnowWorks:
            # Fix MFlowCode/MFC#21: Check whether attempting to run a job will hang
            # forever. This can happen when using the wrong queue system.

            work_timeout = 30

            cons.print(f"Ensuring the [bold magenta]Interactive Engine[/bold magenta] works ({work_timeout}s timeout) via [bold magenta]syscheck[/bold magenta]:")
            cons.print()
            cons.indent()

            q = multiprocessing.Queue()
            p = multiprocessing.Process(
                target=_interactive_working_worker,
                args=(
                    [self.mpibin.bin] + self.mpibin.gen_params() + [os.sep.join([SYSCHECK.get_install_dirpath(), "bin", "syscheck"])],
                    q,
                ))

            p.start()
            p.join(work_timeout)

            if p.is_alive():
                raise MFCException("""\
The [bold magenta]Interactive Engine[/bold magenta] appears to hang.
This may indicate that the wrong MPI binary is being used to launch parallel jobs. You can specify the correct one for your system
using the <-b,--binary> option. For example:
* ./mfc.sh run <myfile.py> -b mpirun
* ./mfc.sh run <myfile.py> -b srun
""")

            result = q.get(block=False)
            self.bKnowWorks = result.returncode == 0

            if not self.bKnowWorks:
                error_txt = """\
MFC's [bold magenta]syscheck[/bold magenta] (system check) failed to run successfully.
Please review the output bellow and ensure that your system is configured correctly:

"""

                if result is not None:
                    error_txt += f"""\
STDOUT:
{result.stdout}

STDERR:
{result.stderr}
"""
                else:
                    error_txt += f"Evaluation timed out after {work_timeout}s."

                raise MFCException(error_txt)

            cons.print()
            cons.unindent()

        for target in targets:
            cons.print(f"[bold]Running [magenta]{target.name}[/magenta][/bold]:")
            cons.indent()

            if not ARG("dry_run"):
                start_time = time.monotonic()
                env = os.environ.copy()
                if ARG('gpus') is not None:
                    env['CUDA_VISIBLE_DEVICES'] = ','.join([str(_) for _ in ARG('gpus')])

                system(
                    self.__get_exec_cmd(target), cwd=self.input.case_dirpath,
                    env=env
                )
                end_time = time.monotonic()
                cons.print(no_indent=True)

                cons.print(f"[bold green]Done[/bold green] (in {datetime.timedelta(seconds=end_time - start_time)})")

            cons.print()
            cons.unindent()


class BatchEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Batch", "batch")

    def get_args(self) -> str:
        return f"""\
Nodes         (-N)  {ARG('nodes')}
Tasks (/node) (-n)  {ARG('tasks_per_node')}
Walltime      (-w)  {ARG("walltime")}
Partition     (-p)  {ARG("partition")}
Account       (-a)  {ARG("account")}
Email         (-@)  {ARG("email")}
"""

    def run(self, targets) -> None:
        qsystem = queues.get_system()
        cons.print(f"Detected the [bold magenta]{qsystem.name}[/bold magenta] queue system.")

        targets = get_targets([SYSCHECK] + targets)

        cons.print(f"Running {format_list_to_string([_.name for _ in targets], 'bold magenta')}:")
        cons.indent()

        self.__create_batch_file(qsystem, targets)

        if not ARG("dry_run"):
            self.__execute_batch_file(qsystem)

            cons.print("[bold yellow]INFO:[/bold yellow] Batch file submitted! Please check your queue system for the job status.")
            cons.print("[bold yellow]INFO:[/bold yellow] If an error occurs, please check the generated batch file and error logs for more information.")
            cons.print("[bold yellow]INFO:[/bold yellow] You can modify the template batch file to your needs.")

        cons.unindent()

    def __get_batch_dirpath(self) -> str:
        return copy.copy(self.input.case_dirpath)

    def __get_batch_filename(self) -> str:
        return f"{ARG('name')}.sh"

    def __get_batch_filepath(self):
        return os.path.abspath(os.sep.join([
            self.__get_batch_dirpath(),
            self.__get_batch_filename()
        ]))

    def __generate_prologue(self, qsystem: queues.QueueSystem) -> str:
        modules = f""

        if does_system_use_modules():
            modules = f"""\
printf ":) Loading modules...\\n"

module purge
module load {' '.join(get_loaded_modules())}
"""

        return f"""\
TABLE_FORMAT_LINE="| - %-14s %-35s - %-14s %-35s |\\n"
TABLE_HEADER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_FOOTER="+-----------------------------------------------------------------------------------------------------------+ \\n"
TABLE_TITLE_FORMAT="| %8s %-96s |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time:"    "$(date +%T)"                       "Start-date:"    "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{ARG("partition")}"      "Walltime:"      "{ARG("walltime")}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{ARG("account")}"        "Nodes:"         "{ARG("nodes")}")
$(printf "$TABLE_FORMAT_LINE" "Job Name:"      "{ARG("name")}"           "Engine"         "{ARG("engine")}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{qsystem.name}"                     "Email:"         "{ARG("email")}")
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Starting" "{ARG("name")} from {ARG("input")}:"
printf "$TABLE_CONTENT\\n"
printf "$TABLE_FOOTER\\n"

{modules}

cd "{self.input.case_dirpath}"

t_start=$(date +%s)
"""

    def __generate_epilogue(self) -> str:
        return f"""\
code=$?

t_stop="$(date +%s)"

printf "\\n$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished" "{ARG("name")}:"
printf "$TABLE_FORMAT_LINE" "Total-time:"  "$(expr $t_stop - $t_start)s"  "Exit Code:" "$code"
printf "$TABLE_FORMAT_LINE" "End-time:"    "$(date +%T)"                  "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

exit $code
"""

    def __evaluate_expression(self, expr: str) -> str:
        # See if it computable
        try:
            # We assume eval is safe because we control the expression.
            # pylint: disable=eval-used
            r = str(eval(expr, ARGS()))
            return r if not isspace(r) else None
        except Exception as exc:
            raise MFCException(f"BatchEngine: '{expr}' is not a valid expression in the template file. Please check your spelling.") from exc

    def __batch_evaluate(self, s: str, qsystem: queues.QueueSystem, targets):
        targets = get_targets(targets)

        replace_list = [
            ("{MFC::PROLOGUE}", self.__generate_prologue(qsystem)),
            ("{MFC::PROFILER}", ' '.join(profiler_prepend())),
            ("{MFC::EPILOGUE}", self.__generate_epilogue()),
            ("{MFC::BINARIES}", ' '.join([f"'{target.get_install_binpath()}'" for target in targets])),
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
                s = re.sub(rf"^.*{match}.*$\n", "", s, flags=re.MULTILINE)

                cons.print(f"> > [bold yellow]Warning:[/bold yellow] [magenta]{match[1:-1]}[/magenta] was not specified. Thus, any line it appears on will be discarded.")

        return s

    def __create_batch_file(self, qsystem: queues.QueueSystem, targets: typing.List[MFCTarget]):
        cons.print("> Generating batch file...")
        filepath = self.__get_batch_filepath()
        cons.print("> Evaluating template file...")
        content = self.__batch_evaluate(qsystem.template, qsystem, targets)

        cons.print("> Writing batch file...")
        file_write(filepath, content)

    def __execute_batch_file(self, qsystem: queues.QueueSystem):
        # We CD to the case directory before executing the batch file so that
        # any files the queue system generates (like .err and .out) are created
        # in the correct directory.
        cmd = qsystem.gen_submit_cmd(self.__get_batch_filename())

        if system(cmd, cwd=self.__get_batch_dirpath()) != 0:
            raise MFCException(f"Submitting batch file for {qsystem.name} failed. It can be found here: {self.__get_batch_filepath()}. Please check the file for errors.")


ENGINES = [ InteractiveEngine(), BatchEngine() ]

def get_engine(slug: str) -> Engine:
    engine: Engine = None
    for candidate in ENGINES:
        candidate: Engine

        if candidate.slug == slug:
            engine = candidate
            break

    if engine is None:
        raise MFCException(f"Unsupported engine {slug}.")

    return engine
