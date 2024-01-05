import re

from ..build   import get_targets, build
from ..printer import cons
from ..state   import ARG
from ..common  import MFCException, isspace

from . import engines, input


def validate_job_options() -> None:
    if not ARG("mpi") and any({ARG("nodes") > 1, ARG("tasks_per_node") > 1}):
        raise MFCException("RUN: Cannot run on more than one rank with --no-mpi.")

    if ARG("nodes") <= 0:
        raise MFCException("RUN: At least one node must be requested.")

    if ARG("tasks_per_node") <= 0:
        raise MFCException("RUN: At least one task per node must be requested.")

    if not isspace(ARG("email")):
        # https://stackoverflow.com/questions/8022530/how-to-check-for-valid-email-address
        if not re.match(r"\"?([-a-zA-Z0-9.`?{}]+@\w+\.\w+)\"?", ARG("email")):
            raise MFCException(f'RUN: {ARG("email")} is not a valid e-mail address.')


def run(targets = None):
    targets = get_targets(targets or ARG("targets"))

    build(targets)

    cons.print("[bold]Run[/bold]")
    cons.indent()

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

    for target in targets:
        cons.print(f"Generating input files for [magenta]{target.name}[/magenta]...")
        cons.indent()
        cons.print()
        input_file.generate_inp(target)
        cons.print()
        cons.unindent()

    engine.run(targets)
