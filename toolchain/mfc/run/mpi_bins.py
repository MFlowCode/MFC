import typing, dataclasses

from ..      import common
from ..state import ARG

# Note: This file is now only used when running
#       in serial mode.

@dataclasses.dataclass
class MPIBinary:
    name: str
    bin:  str

    def is_present(self) -> bool:
        return common.does_command_exist(self.bin)

    def gen_params(self) -> typing.List[str]:
        raise common.MFCException(f"MPIBinary::gen_params <{self.name}> not implemented.")


class JSRUN(MPIBinary):
    def __init__(self):
        super().__init__("IBM's JSRUN", "jsrun")

    def gen_params(self) -> typing.List[str]:
        # ORNL Summit: https://docs.olcf.ornl.gov/systems/summit_user_guide.html?highlight=lsf#launching-a-job-with-jsrun
        # We create one resource-set per CPU(Core)/GPU pair.
        nrs=ARG("tasks_per_node")*ARG("nodes")
        cores_per_rs=1
        gpus_per_rs=min(ARG("tasks_per_node"), 1)
        tasks_per_rs=1

        arguments=[
            '--nrs',          nrs,
            '--cpu_per_rs',   cores_per_rs,
            '--gpu_per_rs',   gpus_per_rs,
            '--tasks_per_rs', tasks_per_rs
        ]

        if gpus_per_rs >= 1:
            arguments.append('--smpiargs=-gpu')

        return arguments


class SRUN(MPIBinary):
    def __init__(self):
        super().__init__("SLURM's SRUN", "srun")

    def gen_params(self) -> typing.List[str]:
        params = ['--ntasks-per-node', ARG("tasks_per_node")]

        if ARG("nodes") != 1:
            params += ['-N', ARG("nodes")]

        # MFC binds its GPUs on its own, as long as they have been allocated
        # by the system's scheduler, or are present on your local machine,
        # if running in serial mode.

        if not common.isspace(ARG("account")):
            params += ['-A', ARG("account")]

        if not common.isspace(ARG("partition")):
            params += ['-p', ARG("partition")]

        return params


class MPIEXEC(MPIBinary):
    def __init__(self):
        super().__init__("MPIEXEC", "mpiexec")

    def gen_params(self) -> str:
        return ["-np", ARG("tasks_per_node")*ARG("nodes")]


class MPIRUN(MPIBinary):
    def __init__(self):
        super().__init__("MPIRUN", "mpirun")

    def gen_params(self) -> str:
        return ["-np", ARG("tasks_per_node")*ARG("nodes")]


# In descending order of priority (if no override present)
BINARIES: list = [ JSRUN(), SRUN(), MPIRUN(), MPIEXEC() ]

def get_binary(exclude: typing.List[str] = None) -> MPIBinary:
    if exclude is None:
        exclude = []

    binaries = [
        b for b in BINARIES if b.is_present() and b.bin not in exclude
    ]

    if len(binaries) == 0:
        raise common.MFCException("No MPI binary found.")

    # Handle user override
    if ARG("binary") is not None:
        for binary in binaries:
            binary: MPIBinary

            if binary.bin == ARG("binary"):
                return binary

        raise common.MFCException(f"MPI Binary <{ARG('binary')}> not found.")

    return binaries[0]
