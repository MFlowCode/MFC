import csv
import numpy as np
import statistics
from tqdm import tqdm

import mfc.viz
from case import dt, sol_L as sol


case = mfc.viz.Case(".", dt)

for name in tqdm(sol.species_names, desc="Loading Variables"):
    case.load_variable(f"Y_{name}", f"prim.{5 + sol.species_index(name)}")
case.load_variable("rho", "prim.1")
case.load_variable("u", "prim.2")
case.load_variable("p", "prim.3")
case.load_variable("T", "prim.15")

steps = case.get_timesteps()

for step in [min(steps), max(steps)]:
    t = step * dt

    with open(f"mfc-{step}.csv", "w") as f:
        writer = csv.writer(f)
        keys = ['x'] + list(set(case.get_data()[0].keys()) - set(["x"]))
        writer.writerow(keys)
        for ix, x in enumerate(sorted(case.get_coords()[0])):
            row = [case.get_data()[step][key][ix] for key in keys]
            writer.writerow(row)
