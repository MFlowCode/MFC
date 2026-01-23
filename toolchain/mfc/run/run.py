import re, os, sys, typing, dataclasses, shlex

from glob import glob

from mako.lookup   import TemplateLookup
from mako.template import Template

from ..build   import get_targets, build, REQUIRED_TARGETS, SIMULATION
from ..printer import cons
from ..state   import ARG, ARGS, CFG, gpuConfigOptions
from ..common  import MFCException, isspace, file_read, does_command_exist
from ..common  import MFC_TEMPLATE_DIR, file_write, system, MFC_ROOT_DIR
from ..common  import format_list_to_string, file_dump_yaml

from . import queues, input


def __validate_job_options() -> None:
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


def __profiler_prepend() -> typing.List[str]:
    if ARG("ncu") is not None:
        if not does_command_exist("ncu"):
            raise MFCException("Failed to locate [bold green]NVIDIA Nsight Compute[/bold green] (ncu).")

        return ["ncu", "--nvtx", "--mode=launch-and-attach",
                       "--cache-control=none", "--clock-control=none"] + ARG("ncu")

    if ARG("nsys") is not None:
        if not does_command_exist("nsys"):
            raise MFCException("Failed to locate [bold green]NVIDIA Nsight Systems[/bold green] (nsys).")

        return ["nsys", "profile", "--stats=true", "--trace=mpi,nvtx,openacc"] + ARG("nsys")

    if ARG("rcu") is not None:
        if not does_command_exist("rocprof-compute"):
            raise MFCException("Failed to locate [bold red]ROCM rocprof-compute[/bold red] (rocprof-compute).")

        return ["rocprof-compute", "profile", "-n", ARG("name").replace('-', '_').replace('.', '_')] + ARG("rcu") + ["--"]

    if ARG("rsys") is not None:
        if not does_command_exist("rocprof"):
            raise MFCException("Failed to locate [bold red]ROCM rocprof-systems[/bold red] (rocprof-systems).")

        return ["rocprof"] + ARG("rsys")

    return []


def get_baked_templates() -> dict:
    return {
        os.path.splitext(os.path.basename(f))[0] : file_read(f)
        for f in glob(os.path.join(MFC_TEMPLATE_DIR, "*.mako"))
    }


def __job_script_filepath() -> str:
    return os.path.abspath(os.sep.join([
        os.path.dirname(ARG("input")),
        f"{ARG('name')}.{'bat' if os.name == 'nt' else 'sh'}"
    ]))


def __get_template() -> Template:
    computer = ARG("computer")
    lookup   = TemplateLookup(directories=[MFC_TEMPLATE_DIR, os.path.join(MFC_TEMPLATE_DIR, "include")])
    baked    = get_baked_templates()

    if (content := baked.get(computer)) is not None:
        cons.print(f"Using baked-in template for [magenta]{computer}[/magenta].")
        return Template(content, lookup=lookup)

    if os.path.isfile(computer):
        cons.print(f"Using template from [magenta]{computer}[/magenta].")
        return Template(file_read(computer), lookup=lookup)

    raise MFCException(f"Failed to find a template for --computer '{computer}'. Baked-in templates are: {format_list_to_string(list(baked.keys()), 'magenta')}.")


def __generate_job_script(targets, case: input.MFCInputFile):
    env = {}
    if ARG('gpus') is not None:
        gpu_ids = ','.join([str(_) for _ in ARG('gpus')])
        env.update({
            'CUDA_VISIBLE_DEVICES': gpu_ids,
            'HIP_VISIBLE_DEVICES':  gpu_ids
        })

    # Compute GPU mode booleans for templates
    gpu_mode = ARG('gpu')

    # Validate gpu_mode is one of the expected values
    valid_gpu_modes = {e.value for e in gpuConfigOptions}
    if gpu_mode not in valid_gpu_modes:
        raise MFCException(
            f"Invalid GPU mode '{gpu_mode}'. Must be one of: {', '.join(sorted(valid_gpu_modes))}"
        )

    gpu_enabled = gpu_mode != gpuConfigOptions.NONE.value
    gpu_acc = gpu_mode == gpuConfigOptions.ACC.value
    gpu_mp = gpu_mode == gpuConfigOptions.MP.value

    content = __get_template().render(
        **{**ARGS(), 'targets': targets},
        ARG=ARG,
        env=env,
        case=case,
        MFC_ROOT_DIR=MFC_ROOT_DIR,
        SIMULATION=SIMULATION,
        qsystem=queues.get_system(),
        profiler=shlex.join(__profiler_prepend()),
        gpu_enabled=gpu_enabled,
        gpu_acc=gpu_acc,
        gpu_mp=gpu_mp
    )

    file_write(__job_script_filepath(), content)


def __generate_input_files(targets, case: input.MFCInputFile):
    for target in targets:
        cons.print(f"Generating input files for [magenta]{target.name}[/magenta]...")
        cons.indent()
        case.generate(target)
        cons.unindent()


def __execute_job_script(qsystem: queues.QueueSystem):
    # We CD to the case directory before executing the batch file so that
    # any files the queue system generates (like .err and .out) are created
    # in the correct directory.
    cmd = qsystem.gen_submit_cmd(__job_script_filepath())

    if system(cmd, cwd=os.path.dirname(ARG("input"))).returncode != 0:
        raise MFCException(f"Submitting batch file for {qsystem.name} failed. It can be found here: {__job_script_filepath()}. Please check the file for errors.")


def run(targets = None, case = None):
    targets = get_targets(list(REQUIRED_TARGETS) + (targets or ARG("targets")))
    case    = case or input.load(ARG("input"), ARG("--"))

    build(targets)

    cons.print("[bold]Run[/bold]")
    cons.indent()

    if ARG("clean"):
        cons.print("Cleaning up previous run...")
        cons.indent()
        case.clean(targets)
        cons.unindent()

    qsystem = queues.get_system()
    cons.print(f"Using queue system [magenta]{qsystem.name}[/magenta].")

    __generate_job_script(targets, case)
    __validate_job_options()
    __generate_input_files(targets, case)

    if not ARG("dry_run"):
        if ARG("output_summary") is not None:
            file_dump_yaml(ARG("output_summary"), {
                "invocation": sys.argv[1:],
                "lock":       dataclasses.asdict(CFG())
            })
        __execute_job_script(qsystem)
