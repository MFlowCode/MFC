import re, os, glob, statistics

CDIR=os.path.abspath(os.path.join("examples", "scaling"))
LDIR=os.path.join(CDIR, "logs")

def get_num(s: str) -> float:
    try:
        return float(re.findall(r"[0-9]+\.[0-9]+(?:E[-+][0-9]+)?", s, re.MULTILINE)[0])
    except:
        return None

def get_nums(arr):
    return {get_num(_) for _ in arr if get_num(_)}

with open(os.path.join(CDIR, "export.csv"), "w") as f:
    f.write("Nodes,Mem (GB),cu_mpi,init (s),sim (s),ts avg,ts samples,mpi avg,mpi samples,perf (ns/gp/eq/rhs)\n")
    for logpath in sorted(glob.glob(os.path.join(LDIR, "run-*-sim*"))):
        logdata = open(logpath, "r").read()

        print(f" + {os.path.relpath(logpath, LDIR)}")
 
        tss  = get_nums(re.findall(r'^ TS .+', logdata, re.MULTILINE))
        mpis = get_nums(re.findall(r'^ MPI .+', logdata, re.MULTILINE))
        try:
            perf = get_num(re.findall(r"^ Performance: .+", logdata, re.MULTILINE)[0])
        except:
            perf = 'N/A'
        if len(tss)  == 0: tss  = [-1.0]
        if len(mpis) == 0: mpis = [-1.0]

        pathels = os.path.relpath(logpath, LDIR).split('-')
        nodes   = int(pathels[1])
        mem     = int(pathels[2])
        cu_mpi  = pathels[3] == 'T'
        ts_avg  = statistics.mean(tss)
        mpi_avg = statistics.mean(mpis)
        init_t  = get_num(re.findall(r"Init took .+", logdata, re.MULTILINE)[0])
        sim_t   = get_num(re.findall(r"sim_duration .+", logdata, re.MULTILINE)[0])
 
        f.write(f"{nodes},{mem},{cu_mpi},{init_t},{sim_t},{ts_avg},{len(tss)},{mpi_avg},{len(mpis)},{perf}\n")

