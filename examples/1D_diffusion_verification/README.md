# 1D Species Diffusion and Heat Conduction Verification

Mixture-averaged species diffusion (with the correction velocity), Fourier conduction and the transport time-step limit
(`chem_params%transport_dt`) on the `h2o2.yaml` mechanism, reactions and viscosity off. One `case.py`, four modes:

| `--mode` | Problem | Reference |
|---|---|---|
| `erf_species` | trace H2 in N2 (mass fraction step 1e-6), 300 K, 1 atm | exact erfc solution, D from Cantera |
| `erf_heat` | pure N2, a 0.01 K step in temperature | the same erfc with lambda/(rho c_p) |
| `mixing` | H2/N2 at 300 K meets O2/H2O/N2 at 1000 K, 3 mm, a wall at x = 0 | `ref_mixing_layer.py` (independent 1D solver) |
| `stab` | periodic, 64 cells, a 0.1 K spike in one cell; dt is the transport limit | growth of the grid-scale mode |

The full results, limits and references are on the @ref validation_species_transport page.

## Run

    ./mfc.sh run examples/1D_diffusion_verification/case.py -t pre_process simulation post_process -n 1 -- --mode erf_species --n 256
    python3 examples/1D_diffusion_verification/check_diffusion.py examples/1D_diffusion_verification --case-args '--mode erf_species --n 256'

Each run writes into the example directory, so check it before the next one. The erf runs end when the front width
sqrt(4 D t) is `--width` (default 0.05) of the domain, 6.4 cells at `--n 128`, so a sweep of `--n` is a convergence
study. The step starts on a cell face, so the initial cell averages are exact. `--tm 2` selects unity-Lewis transport.
Without arguments the case is the CI-sized run (128 cells, about 800 steps); the suite runs it as
`1D -> Chemistry -> Diffusion Verification -> Species Erf`, `-> Heat Erf` and `-> Unity-Lewis Transport-Limited dt`.

The mixing layer needs the reference run once (`python3 ref_mixing_layer.py 1500 ref1500.npz`, a few minutes on one core) and
`python3 check_mixing_layer.py <run directory> ref1500.npz --case-args '--mode mixing --n 300'`.

## Results

GNU 13, one rank, `master` at 0cdfdf02a. For the erf cases L1 is sum |e| dx / (A s), s = sqrt(4 D t_end), against the cell-averaged
erfc solution; the front is 1.6 to 25.6 cells wide.

| N | species L1 | order | heat L1 | order |
|---:|---:|---:|---:|---:|
| 32 | 7.07e-2 | | 7.07e-2 | |
| 64 | 1.58e-2 | 2.16 | 1.58e-2 | 2.16 |
| 128 | 3.85e-3 | 2.04 | 3.84e-3 | 2.04 |
| 256 | 9.57e-4 | 2.01 | 9.51e-4 | 2.01 |
| 512 | 2.39e-4 | 2.00 | 2.33e-4 | 2.03 |

Species diffusion and conduction use the same two-point operator and agree to three digits at equal cells per front width.

Mixing layer against the 3000-cell reference, L1 at 20 us (mass coordinate, `--bc_out -3`):

| N | H2 | O2 | H2O | N2 | T |
|---:|---:|---:|---:|---:|---:|
| 150 | 8.67e-4 | 6.71e-4 | 8.04e-4 | 7.15e-4 | 3.33e-4 |
| 300 | 2.84e-4 | 2.14e-4 | 2.57e-4 | 2.28e-4 | 8.05e-5 |
| 600 | 9.68e-5 | 6.82e-5 | 8.14e-5 | 7.20e-5 | 1.95e-5 |
| 1200 | 3.85e-5 | 2.35e-5 | 2.79e-5 | 2.44e-5 | 1.22e-5 |

The species errors fall at 1.3 to 1.7 and T at 2.0 up to 600 cells, then stalls (see the page).

Stability of `transport_dt` (growth of the grid-scale temperature mode after 400 steps; 1 is neutral):

| mixture | model | cfl_target 0.5 | 0.6 | 0.65 | 0.7 |
|---|---|---:|---:|---:|---:|
| N2, 300 K | 1 mixture-averaged | 3.5e-5 | 1.1e-4 | 7.4e3 | 1.7e4 |
| N2, 300 K | 2 unity Lewis | 3.5e-5 | 1.1e-4 | 9.9e3 | 2.8e4 |
| H2/O2/H2O/N2, 1050 K | 1 | 1.5e-4 | 8.9e-5 | 1.4e-4 | 1.1e-4 |
| H2/O2/H2O/N2, 1050 K | 2 | 7.8e-5 | 3.5e-5 | 1.1e4 | 7.1e4 |

The limit is the SSP-RK3 real-axis bound (diffusion number 0.628) on max(D_k, lambda/(rho c_v)). In the mixture the step is set by H2, twice the conduction limit, and the spike in T does not excite the H2 mode.
