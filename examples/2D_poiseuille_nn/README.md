# 2D Power-Law (Shear-Thinning) Poiseuille Channel

Validates the **power-law term** of MFC's Herschel-Bulkley non-Newtonian viscosity
against a closed-form analytic Poiseuille profile. Demonstrates the shear-thinning
regime (`nn < 1`): the velocity profile is blunter (flatter-topped) than a parabola.

## Regime and parameters

Single Papanastasiou-regularized Herschel-Bulkley fluid with **no yield stress**
(`tau0 = 0`), so the effective viscosity is the pure power law `mu = K*gamma_dot^(n-1)`,
clamped to `[mu_min, mu_max]`:

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(1)%K`   | `2.0e-2` | consistency index |
| `fluid_pp(1)%nn`  | `0.7`    | flow index `< 1` -> shear-thinning |
| `fluid_pp(1)%tau0`| `0.0`    | no yield stress (pure power law) |
| `fluid_pp(1)%hb_m`| `1000.0` | Papanastasiou regularization parameter |
| `fluid_pp(1)%mu_min`/`mu_max` | `1e-6` / `10.0` | viscosity clamp |

Driven by a constant body acceleration `g_x = 8e-2`, `rho = 1`, `pres = 10`
(sound speed ~3.74), giving `u_max ~ 0.011` and Mach ~3e-3 (effectively incompressible).
Channel `L_y = 0.2`, half-height `H = L_y/2 = 0.1`, no-slip walls at `y = 0, L_y`,
periodic in `x`.

## Governing physics and analytic solution

Fully-developed steady channel flow balances the body force against the shear
stress, `tau = rho*g*(H - y)`. With `tau = K*|du/dy|^n` (power law) the closed form is

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - H|^((n+1)/n) )

(`n < 1` blunt/flat-topped; `n = 1` parabola; `n > 1` pointed). For `n < 1` the
viscosity diverges at the shear-free centerline, so any regularized solver caps it
there; the near-wall momentum balance `K*|du/dy|^n = rho*g*(H-y)` is the cleanest
pointwise correctness test and holds regardless of the cap.

## How to run

    ./mfc.sh run examples/2D_poiseuille_nn/case.py -n 2
    python examples/2D_poiseuille_nn/compare_analytic.py

## Validation result

Relative L2 error vs. the analytic power-law profile: **0.68%** (2-rank run, steady
state confirmed). The local momentum balance holds to ~0.1% across the channel, and
the profile bluntness (mean/peak ~0.706) matches the `n = 0.7` theory.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
