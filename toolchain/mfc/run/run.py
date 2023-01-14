import re, typing

from ..printer import cons
from ..state   import ARG

from .  import engines, input
from .. import common,  build


def validate_job_options() -> None:
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

    cons.print("Generating input files...")
    for name in ARG("targets"):
        input_file.generate(name)

    build.build_targets(targets)

    engine.run(ARG("targets"))
    

def run_target(target: str):
    run_targets([target])


def run() -> None:
    run_targets(ARG("targets"))
