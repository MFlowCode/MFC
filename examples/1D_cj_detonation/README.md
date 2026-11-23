# 1D Chapman-Jouguet Detonation Against the Equilibrium CJ Solution

A detonation in 2H2 + O2 + 7Ar at 100 kPa and 298 K travels down a tube from a closed end (x = 0, reflecting wall) at the speed an equilibrium Chapman-Jouguet (CJ) analysis predicts, and is followed by a Taylor-Zeldovich expansion. `case.py` runs it with the 5-equation model (one fluid), mapped WENO5 with `mp_weno`, HLLC, SSP-RK3 and the split alpha-QSS reaction integrator with adaptive substeps (`h2o2.yaml`, no diffusion or viscosity), at a fixed `dt = 0.3 dx / (u_vN + c_vN)` from the frozen von Neumann (vN) state.

    python3 examples/1D_cj_detonation/cj_ref.py          # rewrites cj_ref.json and taylor_eq.txt (Cantera only)
    ./mfc.sh run examples/1D_cj_detonation/case.py -n 1 -t pre_process simulation post_process -- --ic tz --ppi 5 --L 0.5
    python3 examples/1D_cj_detonation/check_cj.py examples/1D_cj_detonation --case-args "--ic tz --ppi 5 --L 0.5"

`cj_ref.py` computes the equilibrium CJ state (Rankine-Hugoniot with the sonic condition against the equilibrium sound speed), the frozen vN state, the ZND induction length Delta_i (distance to the maximum temperature gradient behind the frozen shock at 1.001 D_CJ) and the Taylor-Zeldovich fan, all from the mechanism file MFC uses. The committed `cj_ref.json` and `taylor_eq.txt` are its output for 2H2 + O2 + 7Ar (`--p1 4e4` gives the 40 kPa mixture; the checker and `case.py` read the files, so rerun it only for that).

`case.py` takes `--ic` (`tz`: the Taylor-Zeldovich solution at 40 Delta_i as ten constant states plus a vN slab `--vnw` induction lengths thick; `driver`: a slug of equilibrium products at `--pdrive` times the explosion pressure), `--ppi` (cells per Delta_i), `--L` (tube length [m]), `--steps` (stop after this many steps), `--nsave` and the integrator options `--substeps` and `--substeps-max`. `check_cj.py` reads the Silo output and prints the front speed, the peak pressure and (for `tz`) the plateau pressure against the reference. It needs the Silo output (`post_process`).

Measured on `master` at 0cdfdf02a (the CJ solution: D_CJ = 1692.225 m/s, p_CJ = 1684.9 kPa, T_CJ = 3075 K, Delta_i = 0.0763 mm, plateau 622.9 kPa and 2570 K):

| Tube | cells per Delta_i | cells | D / D_CJ - 1 | D_end / D_CJ - 1 | local sd | p_max/p_vN | plateau p | earlier (D, D_end, plateau) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 50 cm | 2.5 | 16380 | -0.24 % | -0.19 % | 0.06 % | 0.943 | +0.12 % | -0.24 %, -0.17 %, +0.15 % |
| 50 cm | 5 | 32760 | -0.25 % | -0.18 % | 0.06 % | 0.985 | +0.10 % | -0.25 %, -0.16 %, +0.12 % |
| 10 cm | 2.5 | 3276 | -1.51 % | -1.14 % | 0.33 % | 0.920 | -0.20 % | -1.53 %, -1.14 %, -0.15 % |
| 10 cm | 5 | 6552 | -1.52 % | -1.14 % | 0.30 % | 0.965 | -0.31 % | -1.51 %, -1.15 %, -0.25 % |

The grids agree to 0.02 % in the front speed at the same time and the shortfall shrinks with distance (the 10 cm tubes end at 52 us, the 50 cm tubes at 264 us), so it is the approach to CJ and not discretization error. The plateau pressure is within 0.3 % and the front peak pressure reaches 0.92 to 0.99 of the vN value at 2.5 to 5 cells per Delta_i. NVHPC 25.11 on an A100 and on the CPU agree to the digits shown (10 cm, 2.5 cells per Delta_i); all runs are double precision.

The suite registers `1D -> Chemistry -> CJ Detonation -> Front Formation` (`--ic driver`, 200 cells, 300 steps). Results and limits: @ref validation_detonation.
