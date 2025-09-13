import re, os, csv, glob, statistics
from dataclasses import dataclass, fields

CDIR = os.path.abspath(os.path.join("examples", "scaling"))
LDIR = os.path.join(CDIR, "logs")


def get_num(s: str) -> float:
    try:
        return float(re.findall(r"[0-9]+\.[0-9]+(?:E[-+][0-9]+)?", s, re.MULTILINE)[0])
    except:
        return None


def get_nums(arr):
    return {get_num(_) for _ in arr if get_num(_)}


@dataclass(frozen=True, order=True)
class Configuration:
    nodes: int
    mem: int
    rdma_mpi: bool


@dataclass
class Result:
    ts_avg: float
    mpi_avg: float
    init_t: float
    sim_t: float


runs = {}

for logpath in glob.glob(os.path.join(LDIR, "run-*-sim*")):
    logdata = open(logpath, "r").read()

    tss = get_nums(re.findall(r"^ TS .+", logdata, re.MULTILINE))
    mpis = get_nums(re.findall(r"^ MPI .+", logdata, re.MULTILINE))
    try:
        perf = get_num(re.findall(r"^ Performance: .+", logdata, re.MULTILINE)[0])
    except:
        perf = "N/A"

    if len(tss) == 0:
        tss = [-1.0]
    if len(mpis) == 0:
        mpis = [-1.0]

    pathels = os.path.relpath(logpath, LDIR).split("-")

    runs[Configuration(nodes=int(pathels[1]), mem=int(pathels[2]), rdma_mpi=pathels[3] == "T")] = Result(
        ts_avg=statistics.mean(tss),
        mpi_avg=statistics.mean(mpis),
        init_t=get_num(re.findall(r"Init took .+", logdata, re.MULTILINE)[0]),
        sim_t=get_num(re.findall(r"sim_duration .+", logdata, re.MULTILINE)[0]),
    )

with open(os.path.join(CDIR, "export.csv"), "w") as f:
    writer = csv.writer(f, delimiter=",")
    writer.writerow([_.name for _ in fields(Configuration) + fields(Result)])

    for cfg in sorted(runs.keys()):
        writer.writerow([getattr(cfg, _.name) for _ in fields(Configuration)] + [getattr(runs[cfg], _.name) for _ in fields(Result)])

for rdma_mpi in (False, True):
    with open(os.path.join(CDIR, f"strong_scaling{'-rdma_mpi' if rdma_mpi else ''}.csv"), "w") as f:
        writer = csv.writer(f, delimiter=",")

        for nodes in sorted({_.nodes for _ in runs.keys() if _.rdma_mpi == rdma_mpi}):
            row = (nodes * 8,)
            for mem in sorted(
                {_.mem for _ in runs.keys() if _.nodes == nodes and _.rdma_mpi == rdma_mpi},
                reverse=True,
            ):
                ref = runs[
                    Configuration(
                        nodes=sorted({_.nodes for _ in runs.keys() if _.rdma_mpi == rdma_mpi})[0],
                        mem=mem,
                        rdma_mpi=rdma_mpi,
                    )
                ]
                run = runs[Configuration(nodes=nodes, mem=mem, rdma_mpi=rdma_mpi)]
                row = (*row, run.sim_t, ref.sim_t / nodes)

            writer.writerow(row)
