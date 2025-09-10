import os, re
import pandas as pd

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

# Get log files and filter for simulation logs
files = os.listdir("logs/")
files = [f for f in files if "sim" in f]

records = []
for fname in files:
    # Remove extension
    parts = fname.replace(".out", "").split("-")
    scaling, nodes, memory, rdma, phase = parts
    records.append({
        "scaling": scaling,
        "nodes": int(nodes),
        "memory": int(memory),
        "rdma": rdma,
        "phase": phase,
        "file": fname
    })

df = pd.DataFrame(records)

print()

weak_df = df[df["scaling"] == "weak"]
strong_df = df[df["scaling"] == "strong"]
grind_df = df[df["scaling"] == "grind"]

weak_scaling_mem = weak_df["memory"].unique()
weak_scaling_rdma = weak_df["rdma"].unique()

for mem in weak_scaling_mem:
    for rmda in weak_scaling_rdma:
        subset = weak_df[(weak_df["memory"] == mem) & (weak_df["rdma"] == rmda)]
        subset = subset.sort_values(by="nodes")
        times = []
        for _, row in subset.iterrows():
            time_avg = parse_time_avg(os.path.join("logs", row["file"]))
            times.append(time_avg)
        subset = subset.copy()
        subset["time_avg"] = times
        base_time = subset.iloc[0]["time_avg"]
        subset["efficiency"] = base_time / subset["time_avg"]
        print(f"Weak Scaling - Memory: ~{mem}GB, RDMA: {rmda}")
        print(subset[["nodes", "time_avg", "efficiency"]].to_string(index=False))
        print()

strong_scaling_mem = strong_df["memory"].unique()
strong_scaling_rdma = strong_df["rdma"].unique()

for mem in strong_scaling_mem:
    for rdma in strong_scaling_rdma:
        subset = strong_df[(strong_df["memory"] == mem) & (strong_df["rdma"] == rdma)]
        subset = subset.sort_values(by="nodes")
        times = []
        for _, row in subset.iterrows():
            time_avg = parse_time_avg(os.path.join("logs", row["file"]))
            times.append(time_avg)
        subset = subset.copy()
        subset["time_avg"] = times
        base_time = subset.iloc[0]["time_avg"]
        subset["speedup"] = base_time / subset["time_avg"]
        subset["efficiency"] = base_time / (subset["nodes"] * subset["time_avg"])
        print(f"Strong Scaling - Memory: ~{mem}GB, RDMA: {rdma}")
        print(subset[["nodes", "time_avg", "speedup", "efficiency"]].to_string(index=False))
        print()

if not grind_df.empty:
    grind_mem = grind_df["memory"].unique()
    subset = grind_df.sort_values(by="memory")
    times = []
    for _, row in subset.iterrows():
        grind_time = parse_grind_time(os.path.join("logs", row["file"]))
        times.append(grind_time)
    subset = subset.copy()
    subset["grind_time"] = times
    print(f"Grind Time - Single Device")
    print(subset[["memory", "grind_time"]].to_string(index=False))

print()

