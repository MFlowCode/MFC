# 2D Power-Law (Shear-Thickening) Poiseuille Channel

Validates the **power-law term** of MFC's Herschel-Bulkley non-Newtonian viscosity
against a closed-form analytic Poiseuille profile, in the shear-**thickening** regime
(`nn > 1`): the velocity profile is more **pointed** (sharper-topped) than a parabola.
Companion to `examples/2D_poiseuille_nn` (shear-thinning, `nn = 0.7`).

## Regime and parameters

Single Papanastasiou-regularized Herschel-Bulkley fluid with **no yield stress**
(`tau0 = 0`), so the effective viscosity is the pure power law `mu = K*gamma_dot^(n-1)`,
clamped to `[mu_min, mu_max]`:

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(1)%K`   | `5.0e-2` | consistency index |
| `fluid_pp(1)%nn`  | `1.5`    | flow index `> 1` -> shear-thickening |
| `fluid_pp(1)%tau0`| `0.0`    | no yield stress (pure power law) |
| `fluid_pp(1)%hb_m`| `1000.0` | Papanastasiou regularization parameter |
| `fluid_pp(1)%mu_min`/`mu_max` | `1e-6` / `0.035` | viscosity clamp |

Driven by a constant body acceleration `g_x = 5e-2`, `rho = 1`, `pres = 10`
(sound speed ~3.74), giving `u_max ~ 0.013` and Mach ~3e-3 (effectively incompressible).
Channel `L_y = 0.2`, half-height `H = L_y/2 = 0.1`, no-slip walls at `y = 0, L_y`,
periodic in `x`. Grid `m = 24` (x), `n = 63` (y).

For `n > 1` the maximum physical viscosity is at the **wall** (highest shear rate),
`mu_wall = K^(1/n) * (rho*g*H)^((n-1)/n) = 0.0232`. `mu_max = 0.035 ~ 1.5*mu_wall` sits
just above that maximum, so the clamp **never activates** (the analytic profile stays
exact everywhere) while keeping the explicit viscous timestep large — the timestep
scales as `1/mu_max`, so set `mu_max` just above the physical maximum viscosity.

## Governing physics and analytic solution

Fully-developed steady channel flow balances the body force against the shear
stress, `tau = rho*g*(H - y)`. With `tau = K*|du/dy|^n` (power law) the closed form is

    u(y) = (n/(n+1)) * (rho*g/K)^(1/n) * ( H^((n+1)/n) - |y - H|^((n+1)/n) )

(`n < 1` blunt/flat-topped; `n = 1` parabola; `n > 1` pointed). For `n > 1` the
effective viscosity `mu = K*gamma_dot^(n-1) -> 0` at the shear-free centerline (rather
than diverging as for `n < 1`), so no regularization cap is needed and the analytic
profile is an exact reference everywhere.

## How to run

    ./mfc.sh run examples/2D_poiseuille_thickening_nn/case.py -n 2
    python examples/2D_poiseuille_thickening_nn/compare_analytic.py

The run reaches `t_stop = 0.9` (~2 wall-viscous diffusion times) in ~1 min on 2 CPU
ranks with `dt = 8.4e-5`.

## Validation result

Relative L2 error vs. the analytic power-law profile: **1.46%** (2-rank run,
steady-state drift between the last two saves 3.9e-3). The local momentum balance
`K*|du/dy|^n = rho*g*(H-y)` holds to ~1.4% across the channel, and the profile
bluntness (mean/peak = **0.626**) matches the `n = 1.5` theory `(n+1)/(2n+1) = 0.625`,
confirming the pointed shear-thickening profile. `u_max = 1.288e-2` matches the
analytic `1.291e-2`, confirming the clamp stays inactive.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
