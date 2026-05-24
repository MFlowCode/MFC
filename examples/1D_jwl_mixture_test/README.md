# 1D JWL Mixture Test

This is a small sanity test for the mixed-material JWL EOS path in MFC.

The case uses two fluids:

- `fluid_pp(1)%eos = 2`: JWL explosive products
- `fluid_pp(2)%eos = 1`: ideal-gas air

The domain starts as mostly air with a JWL-rich, high-pressure driver on the left. The purpose is not to model a calibrated detonation, but to verify that pressure, energy conversion, sound speed, and post-processing all work for a simple mixed JWL/air setup.

## Run

From the MFC repository root:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
```

If the current build is stale, omit `--no-build`:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py
```

## Expected Result

A successful run completes:

- `syscheck`
- `pre_process`
- `simulation`
- `post_process`

and exits with code `0`.

The case is intentionally tiny: `79x0x0` cells and `20` simulation steps.

## Key Parameters

- JWL constants are TNT-like:
  - `A = 3.712e11`
  - `B = 3.231e9`
  - `R1 = 4.15`
  - `R2 = 0.95`
  - `omega = 0.30`
  - `rho0 = 1630.0`
- Mixed JWL support is exercised through:
  - `fluid_pp(1)%jwl_E0`
  - `fluid_pp(1)%jwl_air_e0`
  - `fluid_pp(1)%jwl_air_rho0`
  - `fluid_pp(1)%jwl_air_gamma`

## Notes

This is a numerical smoke test. It checks that the implementation runs and produces output without NaNs or EOS conversion failures. It is not a validation case against experimental explosive data.
