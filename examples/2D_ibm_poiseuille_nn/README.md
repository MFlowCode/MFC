# 2D IBM-Walled Power-Law Poiseuille Channel

Validates the **immersed boundary (IBM) + non-Newtonian viscosity interaction**:
the channel's no-slip walls are two rectangular IB slabs instead of domain
boundary conditions, so the flow exercises the Herschel-Bulkley viscous fluxes at
the IB walls, which `s_accumulate_ib_face_fluxes` in `m_ibm.fpp` also sums into
the IB force. Companion to
`examples/2D_poiseuille_thickening_nn`, which validates the same fluid against the
same analytic profile with BC walls.

## Geometry and parameters

Domain `x in [0, 0.2]` (periodic), `y in [0, 0.3]`, grid `m = 24`, `n = 95`
(`dy = 0.003125`). Two rectangle IB slabs (`patch_ib%geometry = 3`, no-slip) bound
the flow gap `y in [0.05, 0.25]` (half-height `H = 0.1`, centerline `y_c = 0.15`,
64 cells across the gap). Each slab extends beyond the domain in `x` and mostly
outside the domain in `y`, so the only IB surface seen by the flow is its flat gap
face; the slab centroids sit just *inside* the domain (a centroid exactly on the
boundary is owned by no rank and its `ib_state` force record is never written).
The domain BCs behind the slabs are no-slip walls (`bc_y = -16`), which keep the
body-forced dead fluid inside the slabs benign. `patch_ib%mass = 0` so the
reported IB force is the fluid force alone (no
`bf_x*mass` bookkeeping term); `ib_state_wrt = T` writes it at every save.

Fluid and forcing match the BC-walled template (single Papanastasiou-regularized
Herschel-Bulkley fluid, no yield stress):

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(1)%K`   | `5.0e-2` | consistency index |
| `fluid_pp(1)%nn`  | `1.5`    | flow index `> 1` -> shear-thickening |
| `fluid_pp(1)%tau0`| `0.0`    | no yield stress (pure power law) |
| `fluid_pp(1)%hb_m`| `1000.0` | Papanastasiou regularization parameter |
| `fluid_pp(1)%mu_min`/`mu_max` | `1e-6` / `0.035` | viscosity clamp (`~1.5x mu_wall = 0.0232`, clamp inactive) |

Driven by `g_x = 5e-2`, `rho = 1`, `pres = 10` (Mach ~3e-3); `cfl_adap_dt` with
`cfl_target = 0.3` to `t_stop = 0.9` (~2 wall-viscous diffusion times).

## Analytic solution

In the gap the steady fully-developed power-law profile is

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - y_c|^((n+1)/n) )

and the steady x-force per wall per unit depth is `tau_w * L_x` with
`tau_w = rho*g*H = 5e-3`.

## How to run

    ./mfc.sh run examples/2D_ibm_poiseuille_nn/case.py -n 2
    ./build/venv/bin/python3 examples/2D_ibm_poiseuille_nn/compare_analytic.py

(~2.5 min on 2 CPU ranks.) For the n = 1 equivalence check, run the `newtonian`
and `nn1` modes into two scratch copies and compare:

    for MODE in newtonian nn1; do
      mkdir -p build/ibm_nn_equiv/$MODE
      cp examples/2D_ibm_poiseuille_nn/case.py build/ibm_nn_equiv/$MODE/
      IBM_NN_MODE=$MODE ./mfc.sh run build/ibm_nn_equiv/$MODE/case.py -n 2
    done
    ./build/venv/bin/python3 examples/2D_ibm_poiseuille_nn/check_equivalence.py \
        build/ibm_nn_equiv/newtonian build/ibm_nn_equiv/nn1

## Validation results

**A — n = 1 Newtonian equivalence (IBM mu_eff degeneracy).** A power-law fluid
with `nn = 1, tau0 = 0, K = 0.02` is analytically the same fluid as a Newtonian
one with `mu = 0.02`. Running both modes with the same fixed `dt = 6e-5` to
`t = 0.3` gives **bitwise identical** velocity fields (max abs and rel L2
difference `0.0`) *and* bitwise identical IB wall forces — the
non-Newtonian IBM path reduces exactly to the Newtonian one.

**B — analytic power-law profile (n = 1.5).** Relative L2 error of the steady
x-averaged gap profile vs. the analytic solution: **5.8%** with the nominal
`H = 0.1` (slab faces), **2.6%** with `H = 0.1016` fitted from the `u -> 0`
crossings (IBM walls are sharp only to ~half a cell; the fitted walls sit
`~dy/2 = 0.0016` outside the faces). Steady-state drift between the last two
saves 2.0e-3; profile bluntness (mean/peak) 0.640 vs. the `n = 1.5` theory
0.625 (parabola 0.667), confirming the pointed shear-thickening profile.
The BC-walled template achieves 1.46% on the same fluid; the extra error is the
diffuse-wall representation, not the viscosity model.

**C — IB wall force.** The IB force, the momentum flux the solver applies across
the faces between slab and fluid cells, reaches **9.97e-4** per wall at
`t = 0.9` (both walls identical by symmetry) vs. the analytic
`tau_w*L_x = 1.0e-3`, a ratio of **0.997**. The earlier volume-integration
estimator, which differenced ghost-cell states inside the body, gave 0.971.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
