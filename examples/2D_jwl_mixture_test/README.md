# 2D JWL Mixture Test

This is a small 2D sanity test for the mixed-material JWL EOS path in MFC.

The case uses two fluids:

- `fluid_pp(1)%eos = 2`: JWL explosive products
- `fluid_pp(2)%eos = 1`: ideal-gas air

The domain starts as mostly air with a circular JWL-rich, high-pressure driver placed very close to a small circular immersed boundary. The IB body uses the two-way moving-IBM path with zero initial translational and angular velocity, and a moderate low mass so the blast pressure force can accelerate it without immediately destabilizing the IB pressure correction. This exercises the same mixed JWL pressure/energy/sound-speed path as the 1D case, plus transverse fluxes, 2D IBM marker generation, moving IBM state handling, and IB-state post-processing.

## Run

From the MFC repository root:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

If the current build is stale, omit `--no-build`:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py
```

## Expected Result

A successful run completes:

- `syscheck`
- `pre_process`
- `simulation`
- `post_process`

and exits with code `0`.

The case is intentionally small: `63x31x0` cells with adaptive-CFL time stepping from a fresh start, `n_start = 0`. The initial timestep is `2.5e-8` s, `cfl_target = 0.5`, and the final physical time is `9.15e-5` s, which matches the previous `3660` fixed-step window. Output is requested every `2.5e-7` s. The IB starts from rest with zero translational and angular velocity, while the nearby JWL driver is initialized at `5.0e9` Pa.

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
- Moving IBM support is exercised through:
  - `ib = T`
  - `num_ibs = 1`
  - `patch_ib(1)%geometry = 2`
  - `patch_ib(1)%moving_ibm = 2`
  - `patch_ib(1)%slip = T`
  - `patch_ib(1)%vel(1) = 0.0`
  - `patch_ib(1)%angular_vel(3) = 0.0`
  - `patch_ib(1)%mass = 0.25`
  - `ib_state_wrt = T`

## Notes

This is a numerical smoke test. It is meant to answer, “does the mixed JWL implementation run in 2D with a nearby light immersed body using the moving-IBM setup and IB-state output?” It is not a calibrated detonation, blast, or fluid-structure validation case.

On macOS, an occasional dynamic-loader failure can occur before `pre_process` starts. If that happens, rerun the same command; the case itself should run cleanly.
