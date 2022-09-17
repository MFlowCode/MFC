import typing
import dataclasses

from ..util import common

# Note: This file is now only used when running
#       in serial mode.

@dataclasses.dataclass
class MPIBinary:
    name: str
    bin:  str

    def is_present(self) -> bool:
        return common.does_command_exist(self.bin)

    def gen_params(self, args: list) -> str:
        raise common.MFCException(f"MPIBinary::gen_params <{self.name}> not implemented.")


class JSRUN(MPIBinary):
    def __init__(self):
        super().__init__("IBM's JSRUN", "jsrun")

    def gen_params(self, args: list) -> str:
        # ORNL Summit: https://docs.olcf.ornl.gov/systems/summit_user_guide.html?highlight=lsf#launching-a-job-with-jsrun
        # We create one resource-set per CPU(Core)/GPU pair.
        nrs=args["cpus_per_node"]*args["nodes"]
        cpus_per_rs=1
        gpus_per_rs=min(args["gpus_per_node"], 1)
        tasks_per_rs=1

        arguments=[
            f'--nrs {nrs}',
            f'--cpu_per_rs {cpus_per_rs}',
            f'--gpu_per_rs {gpus_per_rs}',
            f'--tasks_per_rs {tasks_per_rs}'
        ]

        if gpus_per_rs >= 1:
            arguments.append('--smpiargs="-gpu"')

        return ' '.join(arguments)


class SRUN(MPIBinary):
    def __init__(self):
        super().__init__("SLURM's SRUN", "srun")

    def gen_params(self, args: list) -> str:
        params = f' --ntasks-per-node {args["cpus_per_node"]}'

        if args["nodes"] != 1:
            params += f' -N {args["nodes"]}'

        # MFC binds its GPUs on its own, as long as they have been allocated
        # by the system's scheduler, or are present on your local machine,
        # if running in serial mode.
        #
        # if args["gpus_per_node"] != 0:
        #    params += f' -G {args["gpus_per_node"]}'


        if not common.isspace(args["account"]):
            params += f' -A "{args["account"]}"'

        if not common.isspace(args["partition"]):
            params += f' -p "{args["partition"]}"'

        return params


class MPIEXEC(MPIBinary):
    def __init__(self):
        super().__init__("MPIEXEC", "mpiexec")

    def gen_params(self, args: list) -> str:
        np = args["cpus_per_node"]*args["nodes"]

        return f"-np {np}"

class MPIRUN(MPIBinary):
    def __init__(self):
        super().__init__("MPIRUN", "mpirun")

    def gen_params(self, args: list) -> str:
        np = args["cpus_per_node"]*args["nodes"]

        return f"-np {np}"

# In descending order of priority (if no override present)
BINARIES: list = [ JSRUN(), SRUN(), MPIRUN(), MPIEXEC() ]

def get_binary(args: list, exclude: typing.List[str] = None) -> MPIBinary:
    if exclude is None:
        exclude = []

    binaries = [
        b for b in BINARIES if b.is_present() and b.bin not in exclude
    ]

    if len(binaries) == 0:
        raise common.MFCException("No MPI binary found.")

    # Handle user override
    if args["binary"] is not None:
        for binary in binaries:
            binary: MPIBinary

            if binary.bin == args["binary"]:
                return binary

        raise common.MFCException(f"MPI Binary <{args['binary']}> not found.")

    return binaries[0]
