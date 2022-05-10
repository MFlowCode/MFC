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
        header = f"""\
#PBS -N {job_name}
#PBS -l nodes={args["nodes"]}:ppn={args["cpus_per_node"]}
"""

        if not common.isspace(args["account"]):
            header += f'#PBS -A {args["account"]}\n'

        if not common.isspace(args["walltime"]):
            header += f'#PBS -l walltime={args["walltime"]}\n'

        if not common.isspace(args["walltime"]):
            header += f'#PBS -q {args["partition"]}\n'

        return header

    def gen_submit_cmd(self, filename: str) -> None:
        return f"qsub {filename}"


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")

    def is_active(self) -> bool:
        return 0 == os.system(f"bsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        header = f"""\
#BSUB -J {job_name}
#BSUB -nnodes {args["nodes"]}
"""

        if not common.isspace(args["account"]):
            header += f'#BSUB -P {args["account"]}\n'

        if not common.isspace(args["walltime"][:-3]):
            header += f'#BSUB -W {args["walltime"][:-3]}\n'

        return header

    def gen_submit_cmd(self, filename: str) -> None:
        return f"bsub {filename}"


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")

    def is_active(self) -> bool:
        return 0 == os.system(f"sbatch -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        header = f"""\
#SBATCH --job-name="{job_name}"
#SBATCH --nodes={args["nodes"]}
#SBATCH --ntasks-per-node={args["cpus_per_node"]}
#SBATCH --cpus-per-task={1}
"""

        if not common.isspace(args["walltime"]):
            header += f'#SBATCH --time={args["walltime"]}\n'

        if not common.isspace(args["partition"]):
            header += f'#SBATCH --partition={args["partition"]}\n'

        if not common.isspace(args["account"]):
            header += f'#SBATCH --account={args["account"]}\n'

        if not common.isspace(args["email"]):
            header += f"""\
#SBATCH --mail-user="{args["email"]}"
#SBATCH --mail-type="BEGIN, END, FAIL"
"""

        if args["gpus_per_node"] != 0:
            header += f'#SBATCH --gpus=v100-16:{args["gpus_per_node"]}\n'

        return header

    def gen_submit_cmd(self, filename: str) -> None:
        return f"sbatch {filename}"


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]

def get_system() -> QueueSystem:
    for system in QUEUE_SYSTEMS:
        if system.is_active():
            return system

    raise common.MFCException(f"Failed to detect a queue system.")
