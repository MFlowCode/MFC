import re, typing

from ..printer import cons
from ..state   import ARG

from .  import engines, input
from .. import common,  build


def validate_job_options() -> None:
    if not ARG("mpi") and any({ARG("nodes") > 1, ARG("tasks_per_node") > 1}):
        raise common.MFCException("RUN: Cannot run on more than one rank with --no-mpi.")

    if ARG("nodes") <= 0:
        raise common.MFCException("RUN: At least one node must be requested.")

    if ARG("tasks_per_node") <= 0:
        raise common.MFCException("RUN: At least one task per node must be requested.")

    if not common.isspace(ARG("email")):
        # https://stackoverflow.com/questions/8022530/how-to-check-for-valid-email-address
        if not re.match(r"\"?([-a-zA-Z0-9.`?{}]+@\w+\.\w+)\"?", ARG("email")):
            raise common.MFCException(f'RUN: {ARG("email")} is not a valid e-mail address.')


def run_targets(targets: typing.List[str]):    
    cons.print("[bold]Run[/bold]")
    cons.indent()

    if len(ARG("targets")) == 0:
        cons.print(f"> No target selected.")
        return

    input_file = input.load()

    engine = engines.get_engine(ARG("engine"))
    engine.init(input_file)

    cons.print(f"Configuration:")
    cons.indent()
    cons.print(f"""\
Input               {ARG('input')}
Job Name      (-#)  {ARG('name')}
Engine        (-e)  {ARG('engine')}
{engine.get_args()}\
""")
    cons.unindent()

    validate_job_options()

    for name in ARG("targets"):
        cons.print(f"Generating input files for [magenta]{name}[/magenta]...")
        cons.indent()
        cons.print()
        input_file.generate(name)
        cons.print()
        cons.unindent()

    build.build_targets(targets)

    engine.run(ARG("targets"))
    

def run_target(target: str):
    run_targets([target])


def run() -> None:
    run_targets(ARG("targets"))
