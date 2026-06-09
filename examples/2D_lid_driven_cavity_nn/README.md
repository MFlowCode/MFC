# 2D Shear-Thinning Lid-Driven Cavity

Qualitative demonstration of MFC's Herschel-Bulkley non-Newtonian viscosity in a
recirculating flow. A shear-thinning fluid fills a unit square cavity driven by a
moving top lid. Unlike the Poiseuille examples, **this case has no closed-form
analytic solution** — it is a qualitative demonstration of the expected
shear-thinning trend (a primary vortex center shifted toward the moving lid, with
stronger near-wall velocity gradients relative to the Newtonian case).

## Regime and parameters

Two identical Papanastasiou-regularized Herschel-Bulkley fluids (a two-fluid setup
sharing one rheology), pure power law (`tau0 = 0`):

| Parameter | Value | Role |
|-----------|-------|------|
| `fluid_pp(i)%K`   | `1.0e-2` | consistency index; `Re_eff = 1/K = 100` at unit shear rate |
| `fluid_pp(i)%nn`  | `0.5`    | flow index `< 1` -> shear-thinning |
| `fluid_pp(i)%tau0`| `0.0`    | no yield stress (pure power law) |
| `fluid_pp(i)%hb_m`| `1000.0` | Papanastasiou regularization parameter |
| `fluid_pp(i)%mu_min`/`mu_max` | `1e-6` / `1.0` | viscosity clamp |

Unit square `[0,1]^2`, `m = n = 99` (coarse smoke-run grid), all walls no-slip with
the top lid moving at `bc_y%ve1 = 0.5`. Effective Reynolds number `Re_eff = 1/K = 100`
at unit shear rate; the conventional lid-based Reynolds number
`rho*U*L/mu(1) = 0.5/1e-2 = 50` with `U = 0.5` and `mu(1) = K`.

The auto-registered CI test of this example is truncated to 50 time steps and serves
as smoke coverage only, not a physics anchor.

## Governing physics

Incompressible recirculating cavity flow with the shear-dependent power-law
viscosity `mu = K*gamma_dot^(n-1)`. With `n < 1` the fluid thins under the strong
shear beneath the lid and in the corner boundary layers, while the slow cavity core
stays comparatively viscous.

**What to look for** (qualitative, no analytic match): a primary recirculating
vortex whose center, relative to a Newtonian cavity at the same Reynolds number, is
shifted toward the moving lid, with stronger near-wall velocity gradients — the
expected shear-thinning trend. Do not expect a quantitative
error; the committed grid (`m = n = 99`) is intentionally coarse. For quantitative
comparison use `m = n = 499` or finer with a longer `t_step_stop`.

## How to run

    ./mfc.sh run examples/2D_lid_driven_cavity_nn/case.py -n 2

Post-process and inspect the velocity / vorticity (`omega_wrt(3)`) fields; there is
no `compare_analytic.py` for this qualitative case.

## References

- Papanastasiou, T. C. (1987). Flows of materials with yield. *J. Rheol.* 31, 385.
