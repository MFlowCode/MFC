import csv
import cantera as ct
from tqdm import tqdm

import mfc.viz
from case import dt, NS, Tend, SAVE_COUNT, sol


case = mfc.viz.Case(".", dt)

for name in tqdm(sol.species_names, desc="Loading Variables"):
    case.load_variable(f"Y_{name}", f"prim.{5 + sol.species_index(name)}")
case.load_variable("rho", "prim.1")

time_save = Tend / SAVE_COUNT

oh_idx = sol.species_index("OH")


def generate_ct_saves() -> tuple:
    reactor = ct.IdealGasReactor(sol)
    reactor_network = ct.ReactorNet([reactor])

    ct_time = 0.0
    ct_ts, ct_Ys, ct_rhos = [0.0], [reactor.thermo.Y], [reactor.thermo.density]

    while ct_time < Tend:
        reactor_network.advance(ct_time + time_save)
        ct_time += time_save
        ct_ts.append(ct_time)
        ct_Ys.append(reactor.thermo.Y)
        ct_rhos.append(reactor.thermo.density)

    return ct_ts, ct_Ys, ct_rhos


ct_ts, ct_Ys, ct_rhos = generate_ct_saves()

with open("mfc.csv", "w") as f:
    writer = csv.writer(f)
    keys = ["t"] + list(set(case.get_data()[0].keys()) - set(["x"]))
    writer.writerow(keys)
    for i, t_step in enumerate(sorted(case.get_timesteps())):
        t = t_step * dt
        row = [t] + [case.get_data()[t_step][key][0] for key in keys[1:]]
        writer.writerow(row)

with open("cantera.csv", "w") as f:
    writer = csv.writer(f)
    keys = ["t"] + [f"Y_{_}" for _ in list(sol.species_names)] + ["rho"]
    writer.writerow(keys)
    for step in range(len(ct_ts)):
        row = [ct_ts[step]] + [ct_Ys[step][i] for i in range(len(sol.species_names))] + [ct_rhos[step]]
        print([ct_ts[step]], row)
        writer.writerow(row)


def find_induction_time(ts: list, Ys: list, rhos: list) -> float:
    for t, y, rho in zip(ts, Ys, rhos):
        if (y * rho / sol.molecular_weights[oh_idx]) >= 1e-6:
            return t

    return None


skinner_induction_time = 0.052e-3
ct_induction_time = find_induction_time(ct_ts, [y[oh_idx] for y in ct_Ys], [rho for rho in ct_rhos])
mfc_induction_time = find_induction_time(
    sorted(case.get_timestamps()),
    [case.get_data()[step]["Y_OH"][0] for step in sorted(case.get_timesteps())],
    [case.get_data()[step]["rho"][0] for step in sorted(case.get_timesteps())],
)

print("Induction Times ([OH] >= 1e-6):")
print(f" + Skinner et al.: {skinner_induction_time} s")
print(f" + Cantera:        {ct_induction_time} s")
print(f" + (Che)MFC:       {mfc_induction_time} s")
