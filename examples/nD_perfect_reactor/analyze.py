import cantera as ct
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt

import mfc.viz
from case import dt, Tend, SAVE_COUNT, sol


case = mfc.viz.Case(".", dt)

sns.set_theme(style=mfc.viz.generate_cpg_style())

Y_MAJORS = set(["H", "O", "OH", "HO2"])
Y_MINORS = set(["H2O", "H2O2"])
Y_VARS = Y_MAJORS | Y_MINORS

for name in tqdm(Y_VARS, desc="Loading Variables"):
    case.load_variable(f"Y_{name}", f"prim.{5 + sol.species_index(name)}")
case.load_variable("rho", "prim.1")

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

mfc_plots = [[], []]
for y in Y_MAJORS:
    mfc_plots[0].append(case.plot_time(f"Y_{y}", ax=axes[0], label=f"${y}$"))
for y in Y_MINORS:
    mfc_plots[1].append(case.plot_time(f"Y_{y}", ax=axes[1], label=f"${y}$"))

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
for y in Y_VARS:
    sns.lineplot(
        x=ct_ts,
        y=[yt[sol.species_index(y)] for yt in ct_Ys],
        linestyle=":",
        ax=axes[0 if y in Y_MAJORS else 1],
        color="white",
        alpha=0.5,
        label=f"{y} (Cantera)",
    )


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
print(f" + Skinner et al.: {skinner_induction_time:.3e} s")
print(f" + Cantera:        {ct_induction_time:.3e} s")
print(f" + (Che)MFC:       {mfc_induction_time:.3e} s")

axes[0].add_artist(
    axes[0].legend(
        [
            axes[0].axvline(x=skinner_induction_time, color="r", linestyle="-"),
            axes[0].axvline(x=mfc_induction_time, color="b", linestyle="-."),
            axes[0].axvline(x=ct_induction_time, color="g", linestyle=":"),
        ],
        ["Skinner et al.", "(Che)MFC", "Cantera"],
        title="Induction Times",
        loc="lower right",
    )
)

for i in range(2):
    axes[i].legend(title="Species", ncol=2)
    axes[i].set_ylabel("$Y_k$")
    axes[i].set_xscale("log")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Time")

plt.tight_layout()
plt.savefig(f"plots.png", dpi=300)
plt.close()
