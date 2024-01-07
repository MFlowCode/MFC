import os, typing, dataclasses

from mfc import common
from ..state import ARG

@dataclasses.dataclass
class QueueSystem:
    name:     str
    template: str

    def __init__(self, name: str, filename: str) -> None:
        self.name     = name
        self.template = common.file_read(os.sep.join(["toolchain", "templates", filename]))

    def is_active(self) -> bool:
        raise common.MFCException("QueueSystem::is_active: not implemented.")

    def gen_submit_cmd(self, filepath: str) -> typing.List[str]:
        raise common.MFCException("QueueSystem::gen_submit_cmd: not implemented.")


class PBSSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("PBS", "pbs.sh")

    def is_active(self) -> bool:
        return common.does_command_exist("qsub")

    def gen_submit_cmd(self, filepath: str) -> typing.List[str]:
        if ARG("wait"):
            raise common.MFCException("PBS Queue: Sorry, --wait is unimplemented.")

        return ["qsub", filepath]


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF", "lsf.sh")

    def is_active(self) -> bool:
        return common.does_command_exist("bsub") and common.does_command_exist("bqueues")

    def gen_submit_cmd(self, filepath: str) -> None:
        cmd = ["bsub"]

        if ARG("wait"):
            cmd += ["--wait"]

        return cmd + [filepath]


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM", "slurm.sh")

    def is_active(self) -> bool:
        return common.does_command_exist("sbatch")

    def gen_submit_cmd(self, filepath: str) -> None:
        cmd = ["sbatch"]

        if ARG("wait"):
            cmd += ["--wait"]

        return cmd + [filepath]


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]

def get_system() -> QueueSystem:
    for system in QUEUE_SYSTEMS:
        if system.is_active():
            return system

    raise common.MFCException("Failed to detect a queue system.")
