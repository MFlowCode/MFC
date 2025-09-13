import os, re
import pandas as pd
from io import StringIO


def parse_time_avg(path):
    last_val = None
    pattern = re.compile(r"Time Avg =\s*([0-9.E+-]+)")
    with open(path) as f:
        for line in f:
            match = pattern.search(line)
            if match:
                last_val = float(match.group(1))
    return last_val


def parse_grind_time(path):
    last_val = None
    pattern = re.compile(r"Performance: \s*([0-9.E+-]+)")
    with open(path) as f:
        for line in f:
            match = pattern.search(line)
            if match:
                last_val = float(match.group(1))
    return last_val


def parse_reference_file(filename):
    with open(filename) as f:
        content = f.read()

    records = []
    blocks = re.split(r"\n(?=Weak|Strong|Grind)", content.strip())

    for block in blocks:
        lines = block.strip().splitlines()
        header = lines[0].strip()
        body = "\n".join(lines[1:])

        df = pd.read_csv(StringIO(body), delim_whitespace=True)

        if header.startswith("Weak Scaling"):
            # Parse metadata from header
            mem_match = re.search(r"Memory: ~(\d+)GB", header)
            rdma_match = re.search(r"RDMA: (\w)", header)
            memory = int(mem_match.group(1)) if mem_match else None
            rdma = rdma_match.group(1) if rdma_match else None

            for _, row in df.iterrows():
                records.append({"scaling": "weak", "nodes": int(row["nodes"]), "memory": memory, "rdma": rdma, "phase": "sim", "time_avg": row["time_avg"], "efficiency": row["efficiency"]})

        elif header.startswith("Strong Scaling"):
            mem_match = re.search(r"Memory: ~(\d+)GB", header)
            rdma_match = re.search(r"RDMA: (\w)", header)
            memory = int(mem_match.group(1)) if mem_match else None
            rdma = rdma_match.group(1) if rdma_match else None

            for _, row in df.iterrows():
                records.append(
                    {
                        "scaling": "strong",
                        "nodes": int(row["nodes"]),
                        "memory": memory,
                        "rdma": rdma,
                        "phase": "sim",
                        "time_avg": row["time_avg"],
                        "speedup": row["speedup"],
                        "efficiency": row["efficiency"],
                    }
                )

        elif header.startswith("Grind Time"):
            for _, row in df.iterrows():
                records.append({"scaling": "grind", "memory": int(row["memory"]), "grind_time": row["grind_time"]})

    return pd.DataFrame(records)


# Get log files and filter for simulation logs
files = os.listdir("examples/scaling/logs/")
files = [f for f in files if "sim" in f]

records = []
for fname in files:
    # Remove extension
    parts = fname.replace(".out", "").split("-")
    scaling, nodes, memory, rdma, phase = parts
    records.append({"scaling": scaling, "nodes": int(nodes), "memory": int(memory), "rdma": rdma, "phase": phase, "file": fname})

df = pd.DataFrame(records)

ref_data = parse_reference_file("examples/scaling/reference.dat")

print()

weak_df = df[df["scaling"] == "weak"]
strong_df = df[df["scaling"] == "strong"]
grind_df = df[df["scaling"] == "grind"]

weak_ref_df = ref_data[ref_data["scaling"] == "weak"]
strong_ref_df = ref_data[ref_data["scaling"] == "strong"]
grind_ref_df = ref_data[ref_data["scaling"] == "grind"]

weak_scaling_mem = weak_df["memory"].unique()
weak_scaling_rdma = weak_df["rdma"].unique()

for mem in weak_scaling_mem:
    for rdma in weak_scaling_rdma:
        subset = weak_df[(weak_df["memory"] == mem) & (weak_df["rdma"] == rdma)]
        subset = subset.sort_values(by="nodes")
        ref = weak_ref_df[(weak_ref_df["memory"] == mem) & (weak_ref_df["rdma"] == rdma) & (weak_ref_df["nodes"].isin(subset["nodes"]))]
        ref = ref.sort_values(by="nodes")

        times = []
        for _, row in subset.iterrows():
            time_avg = parse_time_avg(os.path.join("examples/scaling/logs", row["file"]))
            times.append(time_avg)

        subset = subset.copy()
        ref = ref.copy()
        subset["time_avg"] = times
        base_time = subset.iloc[0]["time_avg"]

        subset["efficiency"] = base_time / subset["time_avg"]
        subset["rel_perf"] = subset["time_avg"] / ref["time_avg"].values
        print(f"Weak Scaling - Memory: ~{mem}GB, RDMA: {rdma}")
        print(subset[["nodes", "time_avg", "efficiency", "rel_perf"]].to_string(index=False))
        print()

strong_scaling_mem = strong_df["memory"].unique()
strong_scaling_rdma = strong_df["rdma"].unique()

for mem in strong_scaling_mem:
    for rdma in strong_scaling_rdma:
        subset = strong_df[(strong_df["memory"] == mem) & (strong_df["rdma"] == rdma)]
        subset = subset.sort_values(by="nodes")

        ref = strong_ref_df[(strong_ref_df["memory"] == mem) & (strong_ref_df["rdma"] == rdma) & (strong_ref_df["nodes"].isin(subset["nodes"]))]
        ref = ref.sort_values(by="nodes")

        times = []
        for _, row in subset.iterrows():
            time_avg = parse_time_avg(os.path.join("examples/scaling/logs", row["file"]))
            times.append(time_avg)

        subset = subset.copy()
        ref = ref.copy()
        subset["time_avg"] = times
        base_time = subset.iloc[0]["time_avg"]

        subset["speedup"] = base_time / subset["time_avg"]
        subset["efficiency"] = base_time / ((subset["nodes"] / subset.iloc[0]["nodes"]) * subset["time_avg"])
        subset["rel_perf"] = subset["time_avg"] / ref["time_avg"].values
        print(f"Strong Scaling - Memory: ~{mem}GB, RDMA: {rdma}")
        print(subset[["nodes", "time_avg", "speedup", "efficiency", "rel_perf"]].to_string(index=False))
        print()

if not grind_df.empty:
    grind_mem = grind_df["memory"].unique()
    subset = grind_df.sort_values(by="memory")
    ref = grind_ref_df[(grind_ref_df["memory"].isin(subset["memory"]))]
    ref = ref.sort_values(by="memory")

    times = []
    for _, row in subset.iterrows():
        grind_time = parse_grind_time(os.path.join("examples/scaling/logs", row["file"]))
        times.append(grind_time)

    subset = subset.copy()
    ref = ref.copy()

    subset["grind_time"] = times
    subset["rel_perf"] = subset["grind_time"] / ref["grind_time"].values
    print(f"Grind Time - Single Device")
    print(subset[["memory", "grind_time", "rel_perf"]].to_string(index=False))

print()
