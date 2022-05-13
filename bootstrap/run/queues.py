import os
import dataclasses

import common

def queue_helper(queue_pre: str, initial: list, dyn: list, args: dict) -> str:
    return "\n".join([ f'{queue_pre} {flag}' for flag in initial + [ flag for flag, predicate in dyn if predicate ] + args["flags"] ])


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
        return queue_helper(
            f"#PBS",
            [ f"-N {job_name}",
              f"-l nodes={args['nodes']}:ppn={args['cpus_per_node']}" ],
            [ (f'-A {args["account"]}',           not common.isspace(args["account"])),
              (f'-l walltime={args["walltime"]}', not common.isspace(args["walltime"])),
              (f'-q {args["partition"]}',         not common.isspace(args["partition"])),
              (f'-M {args["email"]}',             not common.isspace(args["email"])) ],
            args
        )

    def gen_submit_cmd(self, filename: str) -> None:
        return f"qsub {filename}"


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")

    def is_active(self) -> bool:
        return 0 == os.system(f"bsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return queue_helper(
            f"#BSUB",
            [ f"-J {job_name}",
              f'-nnodes {args["nodes"]}',
              f'-N' ],
            [ (f'-P {args["account"]}',       not common.isspace(args["account"])),
              (f'-W {args["walltime"][:-3]}', not common.isspace(args["walltime"][:-3])) ],
            args
        )

    def gen_submit_cmd(self, filename: str) -> None:
        return f"bsub {filename}"


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")

    def is_active(self) -> bool:
        return 0 == os.system(f"sbatch -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return queue_helper(
            f"#SBATCH",
            [ f'--job-name="{job_name}"',
              f'--nodes={args["nodes"]}',
              f'--ntasks-per-node={args["cpus_per_node"]}',
              f'--cpus-per-task={1}' ],
            [ (f'--time={args["walltime"]}',       not common.isspace(args["walltime"])), 
              (f'--partition={args["partition"]}', not common.isspace(args["partition"])),
              (f'--account={args["account"]}',     not common.isspace(args["account"])),
              (f'--mail-user="{args["email"]}"',   not common.isspace(args["email"])),
              (f'--mail-type="BEGIN, END, FAIL"',  not common.isspace(args["email"])) ],            
            args
        )

    def gen_submit_cmd(self, filename: str) -> None:
        return f"sbatch {filename}"


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]

def get_system() -> QueueSystem:
    for system in QUEUE_SYSTEMS:
        if system.is_active():
            return system

    raise common.MFCException(f"Failed to detect a queue system.")
