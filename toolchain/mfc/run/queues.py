import os, typing, dataclasses

from .. import common


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

    def gen_submit_cmd(self, filename: str) -> typing.List[str]:
        return ["qsub", filename]


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF", "lsf.sh")

    def is_active(self) -> bool:
        return common.does_command_exist("bsub") and common.does_command_exist("bqueues")

    def gen_submit_cmd(self, filename: str) -> None:
        return ["bsub", filename]


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM", "slurm.sh")

    def is_active(self) -> bool:
        return common.does_command_exist("sbatch")

    def gen_submit_cmd(self, filename: str) -> None:
        return ["sbatch", filename]


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]

def get_system() -> QueueSystem:
    for system in QUEUE_SYSTEMS:
        if system.is_active():
            return system

    raise common.MFCException(f"Failed to detect a queue system.")
