# 2D Bingham (Yield-Stress) Poiseuille Channel

Validates the **yield-stress term** of MFC's Herschel-Bulkley non-Newtonian
viscosity against a closed-form analytic Poiseuille profile. Demonstrates the
Bingham regime (`nn = 1`, `tau0 > 0`): a rigid **plug** of uniform velocity forms
near the centerline, where the shear stress falls below the yield stress.

## Regime and parameters

Single Papanastasiou-regularized Herschel-Bulkley fluid with unit flow index, so
`K = mu` is a plain Newtonian consistency and the only non-Newtonian effect is the
yield stress:

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(1)%K`   | `5.0e-2` | `n = 1` -> plain dynamic viscosity `mu` |
| `fluid_pp(1)%nn`  | `1.0`    | flow index = 1 (Bingham, no power-law) |
| `fluid_pp(1)%tau0`| `4.0e-3` | yield stress -> plug half-width `y0 = 0.4 H` |
| `fluid_pp(1)%hb_m`| `1.0e4`  | sharp Papanastasiou yield regularization |
| `fluid_pp(1)%mu_min`/`mu_max` | `1e-6` / `1.0` | viscosity clamp (rigid plug) |

Driven by `g_x = 0.1`, `rho = 1`, `pres = 10`, giving `tau_w = rho*g*H = 1e-2`,
`u_plug ~ 3.6e-3` and Mach ~1e-3. Channel `L_y = 0.2`, `H = 0.1`, no-slip walls,
periodic in `x`.

## Governing physics and analytic solution

The shear stress is `tau = rho*g*(H - y)`; the fluid only flows where `|tau| > tau0`.
With `n = 1`, `K = mu`, `tau_w = rho*g*H > tau0`:

    plug half-width   : y0   = tau0/(rho*g)
    sheared region    : u(y) = (1/(2*mu*rho*g)) *
                               [ (tau_w - tau0)^2 - (rho*g*(H-y) - tau0)^2 ]   (|y-H| >= y0)
    plug              : u_plug = (1/(2*mu*rho*g)) * (tau_w - tau0)^2            (|y-H| <  y0)

The signature of a correct yield term is the flat plug of uniform velocity within
`|y - H| < y0`. Requires `tau_w > tau0` for any flow.

## How to run

    ./mfc.sh run examples/2D_bingham_poiseuille_nn/case.py -n 2
    python examples/2D_bingham_poiseuille_nn/compare_analytic.py

## Validation result

Relative L2 error vs. the analytic Bingham profile: **2.5%** (2-rank run, steady
state confirmed). A plug forms at the centerline with half-width matching the
analytic `y0 = tau0/(rho*g) = 0.4 H`.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
