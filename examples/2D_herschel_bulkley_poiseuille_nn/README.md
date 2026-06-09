# 2D General Herschel-Bulkley Poiseuille Channel

Validates the **combined** non-Newtonian terms of MFC's Herschel-Bulkley viscosity
against a closed-form analytic Poiseuille profile: a shear-thinning power law
(`nn = 0.5 < 1`) **and** a yield stress (`tau0 > 0`) acting together. The companion
examples isolate each effect — `2D_poiseuille_nn` / `2D_poiseuille_thickening_nn`
(power-law only, `tau0 = 0`) and `2D_bingham_poiseuille_nn` (yield only, `nn = 1`).
The signature of a correct general Herschel-Bulkley model is a rigid **plug** near the
centerline (where `|tau| < tau0`) joined to a shear-thinning sheared profile at the walls.

## Regime and parameters

Single Papanastasiou-regularized Herschel-Bulkley fluid with both a sub-unity flow
index and a finite yield stress:

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(1)%K`   | `1.5e-2` | consistency index |
| `fluid_pp(1)%nn`  | `0.5`    | flow index `< 1` -> shear-thinning |
| `fluid_pp(1)%tau0`| `3.5e-3` | yield stress -> plug half-width `y0 = 0.35 H` |
| `fluid_pp(1)%hb_m`| `1.0e4`  | sharp Papanastasiou yield regularization |
| `fluid_pp(1)%mu_min`/`mu_max` | `1e-6` / `0.3` | viscosity clamp (rigid plug) |

Driven by `g_x = 0.1`, `rho = 1`, `pres = 10`, giving `tau_w = rho*g*H = 1e-2`,
`u_plug ~ 4e-3` and Mach ~1e-3. Channel `L_y = 0.2`, `H = 0.1`, no-slip walls,
periodic in `x`. Grid `m = 24` (x), `n = 63` (y).

The plug viscosity diverges as the shear rate `-> 0`, so the clamp `mu_max` sets the
plug rigidity. `mu_max = 0.3` (~6x the wall effective viscosity) keeps a clear plug
while keeping the explicit viscous timestep `dt ~ dy^2 rho/mu_max` tractable: the
Bingham-scale `mu_max = 1.0` drove `dt` to ~3e-6 (ETA ~20 min) *and* relaxed to steady
state too slowly to finish in minutes.

## Governing physics and analytic solution

The shear stress is `tau = rho*g*(H - y)`; the fluid only flows where `|tau| > tau0`.
With `tau = tau0 + K*|du/dy|^n` and `tau_w = rho*g*H > tau0`:

    plug half-width   : y0   = tau0/(rho*g)
    sheared region    : u(y) = (n/((n+1)*rho*g)) * K^(-1/n) *
                               [ (tau_w - tau0)^((n+1)/n)
                                 - (rho*g*(H-y) - tau0)^((n+1)/n) ]   (|y-H| >= y0)
    plug              : u_plug = (n/((n+1)*rho*g)) * K^(-1/n) *
                                 (tau_w - tau0)^((n+1)/n)             (|y-H| <  y0)

(upper half mirrors about `y = H`). Requires `tau_w > tau0` for any flow.

## How to run

    ./mfc.sh run examples/2D_herschel_bulkley_poiseuille_nn/case.py -n 2
    python examples/2D_herschel_bulkley_poiseuille_nn/compare_analytic.py

The run reaches `t_stop = 0.4` in ~5 min on 2 CPU ranks with `dt ~ 1e-5`.

## Validation result

Relative L2 error vs. the analytic Herschel-Bulkley profile: **3.8%** (2-rank run,
steady-state drift between the last two saves ~1.1%). The sheared-region momentum
balance `K*|du/dy|^n + tau0 = rho*g*(H-y)` holds to ~1% near the walls. A flat plug
forms at the centerline: because the Papanastasiou plug is regularized (not perfectly
rigid), the strict `>=99% u_max` band understates it, but the **`>=95% u_max` plug
half-width = 0.0359 = 1.03 y0**, matching the analytic `y0 = tau0/(rho*g) = 0.035 = 0.35 H`.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
