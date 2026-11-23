@page validation_detonation Detonation

# Detonation validation

**Status:** 1D Chapman-Jouguet front speed and post-detonation state: validated for 2H2 + O2 + 7Ar. 2D cell width: partially validated (the cell width is not unique and depends on the perturbation seed). H2/air: open discrepancy (resolved waves decouple). **Code:** `master` at 0cdfdf02a. **Examples:** `examples/1D_cj_detonation`. See also @ref validation.

A reacting-Euler detonation couples shock capturing with the stiff split reaction integrator. The 1D check compares the front speed and the expansion behind it with the equilibrium Chapman-Jouguet (CJ) solution and the Taylor-Zeldovich (TZ) similarity solution for a closed end \cite Taylor1950, computed from the mechanism file MFC uses (`h2o2.yaml`), so it tests the solver and not the mechanism. Tables marked "earlier" were measured on a pre-merge build with the settings listed and were not repeated.

## One-dimensional CJ wave {#sec-validation-detonation-cj}

| Item | Value |
| :--- | :--- |
| Geometry | 1D tube, reflecting wall at x = 0, length 10 cm or 50 cm, extrapolation at the far end (the front never reaches it) |
| Gas | 2H2 + O2 + 7Ar (X: H2 2, O2 1, Ar 7), 298 K, 100 kPa, at rest; Euler (no diffusion or viscosity) |
| Reference | equilibrium CJ state and frozen vN state (`cj_ref.py`): \f$D_{CJ}\f$ = 1692.225 m/s, \f$p_{CJ}\f$ = 1684.9 kPa, \f$T_{CJ}\f$ = 3075 K, \f$u_{CJ}\f$ = 735.1 m/s, \f$p_{vN}\f$ = 2874 kPa, \f$T_{vN}\f$ = 2053 K; ZND induction length \f$\Delta_i\f$ = 0.0763 mm; TZ plateau 622.9 kPa and 2570 K |
| Start | the TZ solution at \f$t_0\f$ = 40 \f$\Delta_i/D_{CJ}\f$ (the u = 0 plateau, the fan as 7 constant steps, then a slab of fresh gas at the vN state 5 \f$\Delta_i\f$ thick) |
| Numerics | 5-equation model, mapped WENO5 with `mp_weno`, HLLC, SSP-RK3, split alpha-QSS reactions with adaptive substeps (2 to 16), fixed \f$\Delta t = 0.3\,\Delta x/(u_{vN} + c_{vN})\f$ |
| Measured | `check_cj.py`: shock position (rightmost crossing of \f$p_1 + (p_{CJ} - p_1)/2\f$), least-squares slope over the second half of the run (D), mean local speed over the last 5 saves (D_end), mean peak pressure within 20 \f$\Delta_i\f$ of the shock over \f$p_{vN}\f$, mean plateau pressure |

Measured with `--ic tz` on `master` at 0cdfdf02a (NVHPC 25.11, one NVIDIA A100, double precision). D is the least-squares slope of the shock position over the second half of the run; "earlier" is the pre-merge study (same case on a different build, tube lengths and cells per \f$\Delta_i\f$ as here).

| Tube | cells per \f$\Delta_i\f$ | cells | D / \f$D_{CJ}\f$ - 1 | D_end / \f$D_{CJ}\f$ - 1 | local sd | \f$p_{max}/p_{vN}\f$ | plateau \f$p\f$ | earlier (D, D_end, plateau) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 50 cm | 2.5 | 16380 | -0.24 % | -0.19 % | 0.06 % | 0.943 | +0.12 % | -0.24 %, -0.17 %, +0.15 % |
| 50 cm | 5 | 32760 | -0.25 % | -0.18 % | 0.06 % | 0.985 | +0.10 % | -0.25 %, -0.16 %, +0.12 % |
| 10 cm | 2.5 | 3276 | -1.51 % | -1.14 % | 0.33 % | 0.920 | -0.20 % | -1.53 %, -1.14 %, -0.15 % |
| 10 cm | 5 | 6552 | -1.52 % | -1.14 % | 0.30 % | 0.965 | -0.31 % | -1.51 %, -1.15 %, -0.25 % |

The grids agree to 0.02 % in the front speed at the same time and the shortfall shrinks with distance (the 10 cm tubes end at 52 us, the 50 cm tubes at 264 us), so it is the approach to CJ and not discretization error. The plateau pressure is within 0.3 % and the front peak pressure reaches 0.92 to 0.99 of the vN value at 2.5 to 5 cells per \f$\Delta_i\f$. NVHPC 25.11 on an A100 and on the CPU agree to the digits shown (10 cm, 2.5 cells per \f$\Delta_i\f$); all runs are double precision.

## H2/air at 1 atm {#sec-validation-detonation-h2air}

Measured earlier on a pre-merge build with `h2o2.yaml` (stoichiometric H2/air, 300 K, 100 kPa, \f$D_{CJ}\f$ = 1976.04 m/s, \f$\Delta_i\f$ = 0.2305 mm, a 2 \f$\Delta_i\f$ vN slab on the TZ start), not repeated here. At 2.5 and 5 cells per \f$\Delta_i\f$ the front speed is within 0.6 % of \f$D_{CJ}\f$ (plateau pressure within 0.5 %) but the wave is held together by numerical radical diffusion: the first post-shock cell already reacts. At 10, 20 and 40 cells per \f$\Delta_i\f$ the wave decouples into an inert shock within 3 to 14 us (open discrepancy). The cause is open. The constant-volume ignition delay behind a frozen shock lengthens sharply with a small speed deficit in this mechanism (2.6 times at 0.95 \f$D_{CJ}\f$ and 18.7 times at 0.9 \f$D_{CJ}\f$, against 1.5 and 2.5 times for 2H2 + O2 + 7Ar), so a resolved wave with a slight deficit decouples; whether the 1D H2/air wave is stable at CJ in this model is not settled here. A more dissipative reconstruction (first order or WENO3) keeps the 10 cells per \f$\Delta_i\f$ wave coupled to 37 us, and the default WENO5-M with `mp_weno` decouples at 13 to 14 us.

## Two-dimensional cell width {#sec-validation-detonation-cells}

Measured earlier on a pre-merge build, not repeated (7.5 cells per \f$\Delta_i\f$ in a 34.65 mm periodic channel to 15 cm took about 3 hours on 16 GPU devices; `examples/2D_detonation_cell` runs the same mixture at 6.67 kPa from a driver and is not compared with the measured cell width). Mixture: 2H2 + O2 + 7Ar at 40 kPa (\f$D_{CJ}\f$ = 1667.30 m/s, \f$\Delta_i\f$ = 0.2103 mm; `cj_ref.py --p1 4e4`), the TZ start with a transverse velocity perturbation of 60 m/s on the vN slab, a numerical soot foil from the running maximum of pressure, and the cell width taken from the dominant transverse mode of the foil. The measured cell width at 40 kPa is 6.3 mm (the detonation database of \cite KaneshigeShepherd1997).

| Grid | seed (modes in the channel) | cells after 15 cm | cell width |
| :--- | :--- | :--- | ---: |
| 7.5 cells per \f$\Delta_i\f$ | 7 | 7 | 4.95 mm (0.79 of measured) |
| 7.5 cells per \f$\Delta_i\f$ | 3 | 3 | 11.55 mm (1.83 of measured) |
| 5 cells per \f$\Delta_i\f$ | 3 and 9 | 7 to 9, weakening after 14 cm | 3.9 to 5.0 mm |

The 2D front runs 4 to 8 % below \f$D_{CJ}\f$ at 4 to 15 cm and 2 % below at 20 to 32 cm, the same deficit as the 1D runs of this mixture. Status: partially validated. MFC forms regular, strong, self-sustained cellular detonations, but the cell width in this periodic channel is not unique (each seed keeps its own mode over 15 cm, bracketing the measured value), so it is not a single predicted value. A unique width would need a much longer and wider channel or a random multi-mode seed.


## Known limitations {#sec-validation-detonation-limits}

- The reference uses the same mechanism file as MFC, so the check is of the solver. The CJ solution of `cj_ref.py` agrees with an independent equilibrium code (the Shock and Detonation Toolbox) to 1.2e-6 relative in speed at 40 kPa (1667.298 and 1667.300 m/s) and reproduces its published H2/air value with the Mevel mechanism to 5e-6 (earlier; not repeated here).
- The front approaches \f$D_{CJ}\f$ slowly: the recombination tail behind the front is long (the H2O mass fraction reaches 99 % of its equilibrium value 7 mm, about 92 \f$\Delta_i\f$, behind the shock). A 10 cm tube therefore ends 1.1 to 1.2 % below \f$D_{CJ}\f$; this is the approach to CJ, not discretization error, since the grids agree.
- The start needs a vN slab thick enough to ignite before the rarefaction from its rear reaches the shock; with 1 \f$\Delta_i\f$ the 10 cells per \f$\Delta_i\f$ run fails, and an overdriven slab (5 \f$\Delta_i\f$ in H2/air) overshoots \f$D_{CJ}\f$ by 27 %. This is a property of the start.
- A resolved ZND profile cannot be prescribed through a case file (an initial condition is at most 10 constant patches, `num_patches_max`); the TZ start uses 10.
- The cell-width runs of the earlier study are not repeated: they need tens of millions of cell-steps (multi-node, hours).
