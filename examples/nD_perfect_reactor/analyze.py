#!/usr/bin/env python3
"""
Cantera vs. (Che)MFC validation for the perfectly stirred reactor example.

Reads MFC post-process Silo output, compares species mass fractions and
induction time against a Cantera 0-D ideal-gas reactor reference, and
reproduces the Skinner & Ringrose (1965) induction-time comparison.

Run from this directory after post_process has been executed:
    python analyze.py

All dependencies (cantera, matplotlib, tqdm, h5py) are installed automatically
by the MFC toolchain.

MFC writes species as Y_{name} (e.g. Y_OH) and density as alpha_rho1.
Run `./mfc.sh viz . --list-vars` to verify variable names in your output.
"""
from case import dt, Tend, SAVE_COUNT, sol
from mfc.viz.silo_reader import assemble_silo
from mfc.viz.reader import discover_timesteps
from tqdm import tqdm
import cantera as ct
import matplotlib.pyplot as plt
import sys

import matplotlib
matplotlib.use('Agg')


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CASE_DIR = '.'
Y_MAJORS = {'H', 'O', 'OH', 'HO2'}
Y_MINORS = {'H2O', 'H2O2'}
Y_VARS = Y_MAJORS | Y_MINORS
oh_idx = sol.species_index('OH')
skinner_induction_time = 0.052e-3   # Skinner & Ringrose (1965)


# ---------------------------------------------------------------------------
# Load MFC output
# ---------------------------------------------------------------------------
steps = discover_timesteps(CASE_DIR, 'silo')
if not steps:
    sys.exit('No silo timesteps found — did you run post_process?')

mfc_times = []
mfc_rhos = []
mfc_Ys = {y: [] for y in Y_VARS}

for step in tqdm(steps, desc='Loading MFC output'):
    assembled = assemble_silo(CASE_DIR, step)
    # Perfectly stirred reactor: spatially uniform — take the midpoint cell.
    mid = assembled.x_cc.size // 2
    mfc_times.append(step * dt)
    # alpha_rho1 = partial density of fluid 1; equals total density for single-fluid cases.
    mfc_rhos.append(float(assembled.variables['alpha_rho1'][mid]))
    for y in Y_VARS:
        mfc_Ys[y].append(float(assembled.variables[f'Y_{y}'][mid]))

# ---------------------------------------------------------------------------
# Cantera 0-D reference
# ---------------------------------------------------------------------------
time_save = Tend / SAVE_COUNT


def generate_ct_saves():
    reactor = ct.IdealGasReactor(sol, clone=True)
    net = ct.ReactorNet([reactor])
    phase = reactor.phase
    ct_time = 0.0
    ct_ts = [0.0]
    ct_Ys = [phase.Y.copy()]
    ct_rhos = [phase.density]
    while ct_time < Tend:
        net.advance(ct_time + time_save)
        ct_time += time_save
        ct_ts.append(ct_time)
        ct_Ys.append(phase.Y.copy())
        ct_rhos.append(phase.density)
    return ct_ts, ct_Ys, ct_rhos


ct_ts, ct_Ys, ct_rhos = generate_ct_saves()

# ---------------------------------------------------------------------------
# Induction time: first step where [OH] molar concentration >= 1e-6 mol/m^3
# ---------------------------------------------------------------------------


def find_induction_time(ts, Ys_OH, rhos):
    for t, y_oh, rho in zip(ts, Ys_OH, rhos):
        if y_oh * rho / sol.molecular_weights[oh_idx] >= 1e-6:
            return t
    return None


ct_induction = find_induction_time(ct_ts, [Y[oh_idx] for Y in ct_Ys], ct_rhos)
mfc_induction = find_induction_time(mfc_times, mfc_Ys['OH'], mfc_rhos)

print('Induction Times ([OH] >= 1e-6 mol/m^3):')
print(f'  Skinner et al.: {skinner_induction_time:.3e} s')
print(f'  Cantera:        {ct_induction:.3e} s'
      if ct_induction is not None else '  Cantera:        not reached')
print(f'  (Che)MFC:       {mfc_induction:.3e} s'
      if mfc_induction is not None else '  (Che)MFC:       not reached')

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
_color = {y: _colors[i % len(_colors)] for i, y in enumerate(sorted(Y_VARS))}

for ax, group in zip(axes, [sorted(Y_MAJORS), sorted(Y_MINORS)]):
    for y in group:
        ax.plot(mfc_times, mfc_Ys[y], color=_color[y], label=f'${y}$')
        ax.plot(ct_ts, [Y[sol.species_index(y)] for Y in ct_Ys],
                linestyle=':', color=_color[y], alpha=0.6, label=f'{y} (Cantera)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('$Y_k$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(title='Species', ncol=2)

# Mark induction times on both panels
induction_lines = [
    (skinner_induction_time, 'r', '-',  'Skinner et al.'),
    (mfc_induction,          'b', '-.', '(Che)MFC'),
    (ct_induction,           'g', ':',  'Cantera'),
]
for ax in axes:
    for t, c, ls, _ in induction_lines:
        if t is not None:
            ax.axvline(t, color=c, linestyle=ls)

axes[0].legend(
    handles=[plt.Line2D([0], [0], color=c, linestyle=ls)
             for t, c, ls, _lbl in induction_lines if t is not None],
    labels=[lbl for t, _c, _ls, lbl in induction_lines if t is not None],
    title='Induction Times',
    loc='lower right',
)

plt.tight_layout()
plt.savefig('plots.png', dpi=300)
plt.close()
print('Saved: plots.png')
