# 2D JWL Blast With Moving Particle Bed

This case is an idealized 2D blast interaction problem for the current MFC JWL implementation. A compact TNT-like JWL products region is placed near the lower boundary, the lower boundary acts as a ghost-cell extrapolation, and the blast loads a semicircular bed of small moving immersed-boundary particles.

The goal is not to reproduce a specific experiment. The goal is to stress the pieces that matter for this branch: mixed JWL/air thermodynamics, strong shock propagation, moving IBM bodies, particle-particle collision settings, and IBM state output.

## Geometry

The domain is larger than the earlier smoke test so that particles have space to move:

```text
domain     2.0 x 1.2
grid       200 x 120
dx = dy    about 1.0e-2
```

That grid is intentionally moderate so the case is easier to run while debugging the moving-IBM behavior. For a sharper production-style check with roughly 10 cells across each particle diameter, use about `m = 1000` and `n = 600`.

The blast is initialized near the bottom wall:

```text
blast center    (1.0, 0.06)
blast radius    0.025
blast pressure  2.0e10 Pa
```

The particles are 2D circular immersed boundaries. In a 2D simulation, this is best interpreted as a cross-section through long cylinders, not as true 3D spherical particles.

```text
number of particles  40
particle radius      0.01
particle diameter D  0.02
particle mass        1.5e-3
```

The particle bed is arranged as an upper half-cylinder around the blast. The center-to-center spacing is about `1.2D`, so the clear surface gap is about `0.2D`. That spacing is tight enough to look like a bed, but it avoids initial overlap.

## Boundaries

The bottom boundary is a ghost-cell extrapolation:

```python
"bc_y%beg": -3
```

The left, right, and top boundaries are ghost-cell extrapolation/open boundaries:

```python
"bc_x%beg": -3
"bc_x%end": -3
"bc_y%end": -3
```

Important: `-3` is a fluid boundary condition. It does not delete IBM particles and it does not stop the run when a particle leaves the domain. If particle escape should terminate or remove particles, that needs separate IBM logic.

## Physics Model

The case uses two fluids:

```text
fluid 1  JWL explosive products
fluid 2  ideal-gas air
```

The JWL constants are TNT-like:

```text
A      3.712e11 Pa
B      3.231e9 Pa
R1     4.15
R2     0.95
omega  0.30
rho0   1630 kg/m3
E0     1.0089e10 J/kg
```

The small `eps = 1.0e-8` volume fraction is intentional. It prevents exactly pure-fluid cells in the mixture model, which helps avoid singular mass-fraction and EOS states at sharp JWL/air interfaces.

The exact equation is described in the root `README-JWL-EOS.md`. This case uses the current piecewise implementation: ideal-gas air fallback for very small JWL mass fraction, effective JWL constants in mixed cells, and the standard JWL pressure form in JWL-rich cells.

## IBM Settings

Each particle starts at rest:

```text
velocity        0
angular velocity 0
moving_ibm      2
slip            T
```

Collision settings are enabled:

```python
"collision_model": 1
"ib_coefficient_of_friction": 0.03
"collision_time": 5.0e-7
"coefficient_of_restitution": 0.8
```

These settings make the case useful for testing particle-bed response, but they should not be treated as calibrated material data.

## Run

From the MFC repository root:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

If the build is stale:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py
```

The case writes primitive variables, pressure, density, sound speed, energy, and IBM state files.

## What To Look For

A healthy run should show:

- a strong pressure wave expanding upward from the bottom blast,
- wave reflection from the bottom slip wall,
- shock interaction with the semicircular particle bed,
- nonzero particle forces and velocities in `restart_data/ib_state_*.dat`,
- no negative pressure, NaNs, or immediate IBM instability.

This is a strong regression test for “does the coupled JWL + moving IBM path work?” It is still not a full validation-quality blast-particle calculation.
