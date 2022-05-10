import os
import dataclasses

import common

@dataclasses.dataclass
class QueueSystem:
    name: str

    def is_active(self) -> bool:
        raise common.MFCException("QueueSystem::is_active: not implemented.")

    def gen_batch_header(self, args: dict, target_name: str) -> str:
        raise common.MFCException("QueueSystem::gen_batch_header: not implemented.")

    def gen_submit_cmd(self, filename: str) -> None:
        raise common.MFCException("QueueSystem::gen_submit_cmd: not implemented.")


class PBSSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("PBS")

    def is_active(self) -> bool:
        return 0 == os.system(f"qsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#PBS -A {args["account"]}
#PBS -l nodes={args["nodes"]}:ppn={args["cpus_per_node"]}
#PBS -l walltime={args["walltime"]}
#PBS -q {args["partition"]}
#PBS -N {job_name}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"qsub {filename}"


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")

    def is_active(self) -> bool:
        return 0 == os.system(f"bsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#BSUB -P {args["account"]}
#BSUB -W {args["walltime"][:-3]}
#BSUB -J {job_name}
#BSUB -nnodes {args["nodes"]}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"bsub {filename}"


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")

    def is_active(self) -> bool:
        return 0 == os.system(f"sbatch -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#SBATCH --job-name="{job_name}"
#SBATCH --time={args["walltime"]}
#SBATCH --nodes={args["nodes"]}
#SBATCH --ntasks-per-node={args["cpus_per_node"]}
#SBATCH --cpus-per-task={1}
#SBATCH --partition={args["partition"]}
#SBATCH --account={args["account"]}
#SBATCH --mail-user={args["email"]}
#SBATCH --mail-type=BEGIN, END, FAIL
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"sbatch {filename}"


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]

def get_system() -> QueueSystem:
    for system in QUEUE_SYSTEMS:
        if system.is_active():
            return system

    raise common.MFCException(f"Failed to detect a queue system.")
