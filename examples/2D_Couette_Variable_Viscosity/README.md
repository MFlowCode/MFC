# 2D Couette flow with temperature-dependent viscosity

A unit-level validation of MFC's Arrhenius viscosity model
`fluid_pp%visc_model = 1`, `mu(T) = exp(C + D/T)`, against an **exact** steady
solution. This is the analytical gold standard for the variable-viscosity
feature: it pins down the spatial order of accuracy on a problem with a known
closed-form answer, so a sign error or discretization bug in the variable-`mu`
viscous flux shows up directly.

## The physics, and why it tests mu(T)

A single fluid fills the gap `0 <= y <= H` between two no-slip walls:

- bottom wall (`y = 0`): fixed, held at temperature `T0 = 300 K` (cold),
- top wall (`y = H`): sliding at speed `U`, held at temperature `T1 = 400 K` (hot),
- periodic in the streamwise (`x`) direction.

Conduction sustains a temperature gradient across the gap, and because
`mu(T) = exp(C + D/T)` with `D > 0`, the viscosity **falls as temperature rises**
(the liquid/Arrhenius law). The hot fluid near the moving wall is therefore
*thinner* than the cold fluid near the fixed wall.

In steady state the shear stress `mu(T) du/dy` is uniform across the gap (no
pressure gradient, no body force). Uniform stress with a non-uniform `mu` forces
`du/dy` to vary inversely with `mu`: the velocity climbs slowly through the cold,
viscous lower half and steeply through the hot, thin upper half. **The velocity
profile is curved.** A constant-viscosity solver gives a straight line
`u = U y/H` with `u(H/2) = U/2`; here the exact mid-gap velocity is
`u(H/2) ≈ 0.352 U`. That ~30% deficit *is* the `mu(T)` signal, and it must match
the exact profile to within the discretization error.

## The exact solution (`reference.py`)

The incompressible steady state solves the coupled two-point BVP

```
momentum:  d/dy [ mu(T) du/dy ] = 0
energy:    d/dy [ k  dT/dy ] + mu(T) (du/dy)^2 = 0
```

with `u(0)=0, u(H)=U, T(0)=T0, T(H)=T1`. `reference.py` solves it with
`scipy.solve_bvp` to ~1e-10. The viscous-heating term `mu(T)(du/dy)^2` couples
the energy equation back to the momentum equation, so the reference also exercises
the energy side of `mu(T)`. At this case's Brinkman number (`Br ≈ 0.014`) the
heating bump is small (~0.17 K above the linear profile), and the closed-form
`Br -> 0` limit (linear `T`, `u = U * integral_0^y dy'/mu / integral_0^H dy'/mu`)
agrees with the coupled solver to ~0.02 % of `U` -- an independent check on the
solver.

## Compressible-solver setup

MFC is a compressible code; temperature is not stored but follows from the
stiffened-gas EOS `T = (p + p_inf)/((gam-1) rho cv)`. Two consequences shape the
deck (`case.py`, constants in `couette_config.py`):

- **Soft sound speed.** A deliberately small `cv` puts the sound speed near
  `c ≈ 10` so the flow sits at `Ma ≈ 0.1`. Low Mach keeps the incompressible
  reference faithful and the acoustic time step affordable.
- **Temperature via a density profile.** The target linear `T(y)` is imposed by
  initializing a `y`-varying density `rho(y) = (p0+p_inf)/((gam-1) cv T(y))` at
  uniform pressure -- a valid mechanical equilibrium (no gravity).
- **Modes.** The deck uses `model_eqns = 3` (6-equation) and `riemann_solver = 2`
  (HLLC); `visc_model = 1` also runs with `model_eqns = 2` and HLL. MFC evaluates
  `mu(T)` at the mixture temperature of the reconstructed left and right states at
  each cell face, so the face viscosity follows the reconstructed temperature.
  Fourier conduction (`fluid_pp(1)%k_therm > 0`) with isothermal no-slip walls
  (`bc_y = -16`, `isothermal_in/out`) holds the walls at `T0` and `T1`.
- **Compact viscous gradient (`weno_Re_flux = F`, `weno_avg = F`).** `weno_Re_flux`
  is designed for *discontinuous* viscosity at material interfaces. For a single
  fluid with smooth `mu(T)` it leaves the `2*dy` (odd-even) mode of the
  wall-normal velocity undamped -- and because this flow has no advection
  (`v ~ 0`), the HLLC upwinding has nothing to damp either -- so `u(y)` grows a
  grid-scale checkerboard (~5 % of `U` at `Ny = 33`). The compact central
  gradient damps it; the resolved profile then matches the exact solution to the
  discretization error. (The two flags must be set together: `weno_avg = F` with
  `weno_Re_flux = T` disables the viscous coupling entirely.)

The flow starts from rest (`u = 0`) and develops toward the curved steady profile,
so the match is *not* baked into the initial condition.

Dimensionless parameters: `Re ≈ 20`, `Ma ≈ 0.1`, `Pr = 1`, `Br ≈ 0.014`,
viscosity contrast `mu(T0)/mu(T1) ≈ 3.5`.

## Running it

```bash
# 1. Run the grid-refinement sweep. The first grid is a normal building run
#    (the analytic density IC compiles into pre_process, so the build is
#    case-specific); the remaining grids then reuse that binary (--no-build)
#    and run concurrently, each in its own runs/n<N>/ on one rank.
#    Default grids: Ny = 33, 65, 97.
cd examples/2D_Couette_Variable_Viscosity
./run_suite.sh                 # or: ./run_suite.sh 32 64 96 128

# 2. Compare to the exact solution; writes summary.json + figures/
python3 validate.py
```

A single grid can be run directly, choosing the wall-normal resolution with an
environment variable:

```bash
COUETTE_N=64 ./mfc.sh run examples/2D_Couette_Variable_Viscosity/case.py -n 1
```

## Files

| file | role |
|------|------|
| `couette_config.py` | single source of truth for all physical constants |
| `reference.py`      | exact coupled-BVP solution (`scipy.solve_bvp`) |
| `case.py`           | MFC input deck, parametrized by `COUETTE_N` |
| `run_suite.sh`      | runs the grid sweep into `runs/n<N>/` |
| `validate.py`       | compares MFC to the exact solution; figures + `summary.json` |

## Results

The MFC steady velocity matches the exact `mu(T)` Couette profile and converges at
second order; temperature is exact to ~1e-6 and the runs are fully steady
(`unsteadiness ~ 1e-9`). The viscosity contrast bends `u(H/2)` to **0.352** versus
the constant-viscosity **0.500** -- the signal a constant-`mu` code cannot produce.

| Ny | L2(u)/U  | Linf(u)/U | L2(T)/dT |
|----|----------|-----------|----------|
| 33 | 1.33e-4  | 1.83e-4   | 1.58e-6  |
| 65 | 3.41e-5  | 4.85e-5   | 3.93e-7  |
| 97 | 1.53e-5  | 2.20e-5   | 1.75e-7  |

Observed spatial order (L2 velocity) = **2.00**. `validate.py` writes to `figures/`: profile
overlay (`couette_profiles.png`), L2 convergence (`couette_convergence.png`), and
the spatial relative error (`couette_error.png`) -- the latter shows the
pointwise `(u - u_exact)/U` shrinking with refinement and collapsing under an
`Ny^2` rescaling, i.e. the second-order accuracy holds everywhere in the gap, not
just in the L2 norm, and the full data to `summary.json`.

## Scope and the companion benchmark

This case validates the **exponential/Arrhenius (liquid) law** that MFC actually
implements -- the right gold standard for temperature-dependent droplet
viscosity. It is a unit test: `mu(T)` is isolated on an exact answer, free of the
EOS/buoyancy coupling that entangles system-level benchmarks.

The natural *system-level* companion is the low-Mach differentially-heated cavity
(Le Quéré / Vierendeels), where boundary-layer asymmetry is driven by the
viscosity contrast and validated against tabulated wall Nusselt numbers. That
benchmark uses **Sutherland's law (a gas law)**, which MFC does not yet
implement, so it is a separate follow-up (add the Sutherland model, then the
cavity case).

### Cost note

The time step is explicit-stability limited: acoustic (`~dy`) on coarse grids and
**viscous (`~dy^2`)** once `Ny` exceeds ~60. The default sweep stays in the cheap
range; `Ny = 129` is viscous-CFL limited (~221k steps, hours on a single rank) and
is optional.
