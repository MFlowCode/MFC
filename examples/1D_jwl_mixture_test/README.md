# 1D JWL Mixture Smoke Test

This is the smallest practical test for the mixed JWL/air EOS path. It keeps the geometry simple so that failures are easy to diagnose: a JWL-rich high-pressure driver expands into mostly ideal-gas air.

The case is intentionally modest. It is not a calibrated detonation problem and it is not meant to match explosive test data. It is a quick check that the new EOS path can initialize pressure, convert between primitive and conservative variables, compute sound speed, march a few steps, and post-process output without producing bad thermodynamic states.

## What It Exercises

The case uses:

```text
model_eqns  2
num_fluids  2
fluid 1     JWL explosive products
fluid 2     ideal-gas air
```

The JWL fluid is selected with:

```python
"fluid_pp(1)%eos": 2
```

Air remains the standard ideal/stiffened-gas path:

```python
"fluid_pp(2)%eos": 1
```

The mixed cells use a tiny volume-fraction floor, `eps = 1.0e-8`, to avoid exactly zero volume fraction for either material. That is a numerical safety device for the mixture model, not a physical interface thickness.

The exact EOS used by this case is documented in the root `README-JWL-EOS.md`. In short, cells with `Y_JWL <= 1.0e-2` use the ideal-gas air fallback, partially mixed cells use effective JWL constants, and JWL-rich cells use the standard JWL pressure form.

## Run

From the MFC repository root:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
```

If the build is stale:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py
```

## Expected Result

A successful run completes:

```text
syscheck
pre_process
simulation
post_process
```

and exits with code `0`.

Because this is a smoke test, the important result is not a particular shock location. The important result is that pressure, density, energy, and sound speed stay finite and the run completes.

## Parameters

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

The mixed JWL/air helper also uses reference air values:

```text
air_e0      2.5575e5 J/kg
air_rho0    1.225 kg/m3
air_gamma   0.4
```

In this branch, `jwl_air_gamma` is the `gamma - 1` coefficient used by the current mixed JWL/air helper, matching the convention used in the implementation.
