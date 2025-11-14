import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

N = [32, 64, 128, 256, 512, 1024]
Ord = [1, 3, 5]

errors = np.nan * np.zeros((len(N), len(Ord), 3))

TEND = 200000

for i in range(len(N)):
    for j in range(len(Ord)):

        sim_a1 = pd.read_csv(f"N{N[i]}_O{Ord[j]}/D/cons.5.00.{TEND}.dat", sep=r"\s+", header=None, names=["x", "y"])
        sim_a2 = pd.read_csv(f"N{N[i]}_O{Ord[j]}/D/cons.6.00.{TEND}.dat", sep=r"\s+", header=None, names=["x", "y"])

        exact_a1 = pd.read_csv(f"N{N[i]}_O{Ord[j]}/D/cons.5.00.000000.dat", sep=r"\s+", header=None, names=["x", "y"])
        exact_a2 = pd.read_csv(f"N{N[i]}_O{Ord[j]}/D/cons.6.00.000000.dat", sep=r"\s+", header=None, names=["x", "y"])

        ## 2 norm
        errors[i, j, 0] = np.linalg.norm(sim_a1.y - exact_a1.y) / np.sqrt(N[i])
        errors[i, j, 0] += np.linalg.norm(sim_a2.y - exact_a2.y) / np.sqrt(N[i])

        ## 1 norm
        errors[i, j, 1] = 1 / N[i] * np.sum(np.abs(sim_a1.y - exact_a1.y))
        errors[i, j, 1] += 1 / N[i] * np.sum(np.abs(sim_a2.y - exact_a2.y))

        ## Inf norm
        errors[i, j, 2] = np.max([np.nanmax(np.abs(sim_a1.y - exact_a1.y)), np.nanmax(np.abs(sim_a2.y - exact_a2.y))])

fig, ax = plt.subplots(1, 3, figsize=(12, 8), sharex=True)

colors = ["blue", "green", "red", "purple"]

ref = np.nan * np.zeros((len(N), len(Ord)))

for i in range(3):
    ax[i].plot(N, 30 / np.array(N) ** 1, label="Slope = -1", color=colors[0])
    ax[i].plot(N, 3000 / np.array(N) ** 3, label="Slope = -3", color=colors[1])
    ax[i].plot(N, 5000 / np.array(N) ** 5, label="Slope = -5", color=colors[2])

for j in range(len(Ord)):
    ax[0].plot(N, errors[:, j, 0], "o-", color=colors[j])
    ax[0].set_xscale("log", base=2)
    ax[0].set_yscale("log")
    ax[0].set_title("||error||_2")
    ax[0].legend()

    ax[1].plot(N, errors[:, j, 1], "o-", color=colors[j])
    ax[1].set_xscale("log", base=2)
    ax[1].set_yscale("log")
    ax[1].set_title("||error||_1")

    ax[2].plot(N, errors[:, j, 2], "o-", color=colors[j])
    ax[2].set_xscale("log", base=2)
    ax[2].set_yscale("log")
    ax[2].set_title("||error||_inf")

plt.tight_layout()
plt.show()

errors = np.column_stack((N, errors[:, :, 0]))
errors = np.column_stack((errors, 7 / np.array(N) ** 1))
errors = np.column_stack((errors, 700 / np.array(N) ** 3))
errors = np.column_stack((errors, 2000 / np.array(N) ** 5))

df = pd.DataFrame(errors, columns=["N", "1", "3", "5", "R1", "R3", "R5"], index=N)
df.to_csv("errors.csv", index=False)
