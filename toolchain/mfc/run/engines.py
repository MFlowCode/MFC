import re, os, time, copy, queue, typing, datetime, subprocess, dataclasses, multiprocessing

from ..state     import ARG, ARGS
from ..printer   import cons
from ..          import build, common
from ..run       import queues, mpi_bins
from ..run.input import MFCInputFile


def profiler_prepend():
    if ARG("ncu") is not None:
        if not common.does_command_exist("ncu"):
            raise common.MFCException("Failed to locate [bold green]NVIDIA Nsight Compute[/bold green] (ncu).")

        return ["ncu", "--nvtx", "--mode=launch-and-attach",
                       "--cache-control=none", "--clock-control=none"] + ARG("ncu")

    if ARG("nsys") is not None:
        if not common.does_command_exist("nsys"):
            raise common.MFCException("Failed to locate [bold green]NVIDIA Nsight Systems[/bold green] (nsys).")

        return ["nsys", "profile", "--stats=true", "--trace=mpi,nvtx,openacc"] + ARG("nsys")

    return []


@dataclasses.dataclass
class Engine:
    name: str
    slug: str

    def init(self, input: MFCInputFile) -> None:
        self.input = input

        self._init()

    def _init(self) -> None:
        pass

    def get_args(self) -> typing.List[str]:
        raise common.MFCException(f"MFCEngine::get_args: not implemented for {self.name}.")

    def run(self, names: typing.List[str]) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")

    def get_binpath(self, target: str) -> str:
        # <root>/install/<slug>/bin/<target>
        prefix = build.get_install_dirpath(build.get_target(target))
        return os.sep.join([prefix, "bin", target])


def _interactive_working_worker(cmd: typing.List[str], q: multiprocessing.Queue):
    """ Runs a command and puts the result in a queue. """
    cmd = [ str(_) for _ in cmd ]
    cons.print(f"$ {' '.join(cmd)}")
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    q.put(result)

class InteractiveEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Interactive", "interactive")

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

    def get_exec_cmd(self, target_name: str) -> typing.List[str]:
        cmd = []

        if ARG("mpi"):
            cmd += [self.mpibin.bin] + self.mpibin.gen_params() + ARG("flags")[:]

        cmd += profiler_prepend()

        cmd.append(self.get_binpath(target_name))

        return cmd


    def run(self, names: typing.List[str]) -> None:
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
                    [self.mpibin.bin] + self.mpibin.gen_params() + [os.sep.join([build.get_install_dirpath(build.SYSCHECK), "bin", "syscheck"])],
                    q,
                ))

            p.start()
            p.join(work_timeout)

            if p.is_alive():
                raise common.MFCException("""\
The [bold magenta]Interactive Engine[/bold magenta] appears to hang.
This may indicate that the wrong MPI binary is being used to launch parallel jobs. You can specify the correct one for your system
using the <-b,--binary> option. For example:
* ./mfc.sh run <myfile.py> -b mpirun
* ./mfc.sh run <myfile.py> -b srun
""")

            result = q.get(block=False)
            self.bKnowWorks = result.returncode == 0

            if not self.bKnowWorks:
                error_txt = f"""\
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

                raise common.MFCException(error_txt)

            cons.print()
            cons.unindent()

        for name in names:
            cons.print(f"[bold]Running [magenta]{name}[/magenta][/bold]:")
            cons.indent()

            if not ARG("dry_run"):
                start_time = time.monotonic()
                common.system(self.get_exec_cmd(name), cwd=self.input.case_dirpath)
                end_time   = time.monotonic()
                cons.print(no_indent=True)

                cons.print(f"[bold green]Done[/bold green] (in {datetime.timedelta(seconds=end_time - start_time)})")

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

    def run(self, names: typing.List[str]) -> None:
        system = queues.get_system()
        cons.print(f"Detected the [bold magenta]{system.name}[/bold magenta] queue system.")

        cons.print(f"Running {common.format_list_to_string(names, 'bold magenta')}:")
        cons.indent()

        self.__create_batch_file(system, names)

        if not ARG("dry_run"):
            self.__execute_batch_file(system, names)

            cons.print("[bold yellow]INFO:[/bold yellow] Batch file submitted! Please check your queue system for the job status.")
            cons.print("[bold yellow]INFO:[/bold yellow] If an error occurs, please check the generated batch file and error logs for more information.")
            cons.print("[bold yellow]INFO:[/bold yellow] You can modify the template batch file to your needs.")

        cons.unindent()

    def __get_batch_dirpath(self) -> str:
        return copy.copy(self.input.case_dirpath)

    def __get_batch_filename(self, names: typing.List[str]) -> str:
        return f"{ARG('name')}.sh"

    def __get_batch_filepath(self, names: typing.List[str]):
        return os.path.abspath(os.sep.join([
            self.__get_batch_dirpath(),
            self.__get_batch_filename(names)
        ]))

    def __generate_prologue(self, system: queues.QueueSystem, names: typing.List[str]) -> str:
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
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{ARG("partition")}"      "Walltime:"      "{ARG("walltime")}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{ARG("account")}"        "Nodes:"         "{ARG("nodes")}")
$(printf "$TABLE_FORMAT_LINE" "Job Name:"      "{ARG("name")}"           "Engine"         "{ARG("engine")}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{system.name}"                     "Email:"         "{ARG("email")}")
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
            r = str(eval(expr, ARGS()))
            return r if not common.isspace(r) else None
        except Exception as exc:
            raise common.MFCException(f"BatchEngine: '{expr}' is not a valid expression in the template file. Please check your spelling.")

    def __batch_evaluate(self, s: str, system: queues.QueueSystem, names: typing.List[str]):
        replace_list = [
            ("{MFC::PROLOGUE}", self.__generate_prologue(system, names)),
            ("{MFC::PROFILER}", ' '.join(profiler_prepend())),
            ("{MFC::EPILOGUE}", self.__generate_epilogue()),
            ("{MFC::BINARIES}", ' '.join([f"'{self.get_binpath(x)}'" for x in ["syscheck"] + names])),
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

                cons.print(f"> > [bold yellow]Warning:[/bold yellow] [magenta]{match[1:-1]}[/magenta] was not specified. Thus, any line it appears on will be discarded.")

        return s

    def __create_batch_file(self, system: queues.QueueSystem, names: typing.List[str]):
        cons.print("> Generating batch file...")
        filepath = self.__get_batch_filepath(names)
        cons.print("> Evaluating template file...")
        content = self.__batch_evaluate(system.template, system, names)

        cons.print("> Writing batch file...")
        common.file_write(filepath, content)

    def __execute_batch_file(self, system: queues.QueueSystem, names: typing.List[str]):
        # We CD to the case directory before executing the batch file so that
        # any files the queue system generates (like .err and .out) are created
        # in the correct directory.
        cmd = system.gen_submit_cmd(self.__get_batch_filename(names))

        if common.system(cmd, cwd=self.__get_batch_dirpath()) != 0:
            raise common.MFCException(f"Submitting batch file for {system.name} failed. It can be found here: {self.__get_batch_filepath(target_name)}. Please check the file for errors.")


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
