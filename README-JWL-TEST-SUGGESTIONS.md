# JWL Testing And Implementation Roadmap

This file is the working roadmap for bringing the MFC JWL implementation closer to ROCFLU-level robustness. The main README documents what exists today. This file explains what should be tested next, why those tests matter, and what implementation work is still missing.

## Current Baseline

Two examples are available now:

```text
examples/1D_jwl_mixture_test/case.py
examples/2D_jwl_mixture_test/case.py
```

The 1D case is a smoke test for the EOS. It should be used whenever the JWL thermodynamic formulas are touched.

The 2D case is an integrated blast/IBM test. It uses a compact bottom-wall JWL blast and a 40-particle moving IBM bed. It is the better test for checking that the implementation survives a strong shock, mixed-material conversion, moving bodies, collisions, and post-processing.

Recommended quick checks:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

For the 2D case, inspect `restart_data/ib_state_*.dat` if particle motion is the question. The fluid boundary condition `-3` will not stop or remove particles that leave the domain.

## What We Have Implemented

The current implementation includes:

- JWL EOS selection with `fluid_pp(i)%eos = 2`,
- JWL material parameters in `physical_parameters`,
- pressure from density and energy,
- energy from pressure and density,
- temperature estimate,
- JWL and mixed JWL/air sound-speed support,
- primitive/conservative/flux conversion integration,
- MPI broadcast of JWL parameters,
- pre-process, simulation, and post-process defaults,
- Python parameter definitions and descriptions,
- validation that rejects unsupported model combinations,
- 1D mixed JWL/air smoke test,
- 2D JWL blast with moving IBM particles.

The moving IBM force path was adjusted to follow the working `m_ibm-4.fpp` pressure-gradient force accumulation style for moving particles.

## What Is Still Missing

The main gaps relative to ROCFLU are:

```text
1. Standalone Yfix-style mass-fraction correction
2. Full bad-state guard and recovery workflow
3. Dedicated JWL energy floor
4. Production NaN recovery paths
5. Multiple JWL materials per case
6. PICL-style particle-in-cell coupling
7. Validation against known JWL blast or shock benchmarks
```

These are not small polish items. They are the difference between a useful research implementation and a production-hardened explosive-flow model.

## Priority 1: Yfix-Style Mass-Fraction Correction

ROCFLU has a `RFLU_JWL_Yfix()` style correction that keeps mass fractions away from singular values. MFC currently uses small `eps` values in the cases, but that is not the same as a dedicated correction pass.

What to implement:

- a helper that computes the JWL mass fraction from `alpha_rho`,
- lower and upper clipping for the JWL mass fraction,
- consistent correction of affected mixture quantities,
- diagnostics when clipping happens too often.

Suggested tests:

```text
examples/1D_jwl_yfix_tiny_fraction/case.py
examples/1D_jwl_yfix_near_one/case.py
```

Pass criteria:

- no NaNs,
- no negative pressure,
- no negative sound-speed squared,
- finite pressure in cells with tiny JWL fraction,
- mostly-air cells behave like air.

## Priority 2: Low-JWL-Fraction Air Fallback

Cells with almost no explosive products should not behave like stiff explosive products. ROCFLU has more mature fallback behavior around this.

What to implement:

- a clear threshold for “effectively no JWL,”
- air-dominated pressure and energy behavior below that threshold,
- warnings or counters when fallback is used.

Suggested test:

```text
examples/1D_jwl_low_y_air_fallback/case.py
```

Compare it against an equivalent ideal-gas air case. The pressure and sound speed should be close when `Y_jwl` is near zero.

## Priority 3: Energy Floor

The current formulas can reject some invalid states, but there is not yet a dedicated ROCFLU-equivalent JWL energy floor.

What to implement:

- compute a minimum safe internal energy for JWL-rich states,
- apply the floor before pressure and sound-speed calls,
- track how often the floor is active,
- avoid hiding real numerical instability by silently flooring everything.

Suggested tests:

```text
examples/1D_jwl_energy_floor/case.py
examples/1D_jwl_energy_floor_shock/case.py
```

Pass criteria:

- pressure remains finite,
- sound speed remains real,
- the shock still propagates,
- the floor activates only in intended cells.

## Priority 4: Bad-State And NaN Recovery

ROCFLU's production path has extra recovery behavior around bad states. MFC should eventually have similar guardrails.

What to implement:

- checks around pressure recovery,
- checks around energy inversion,
- checks around sound-speed squared,
- local fallback when a state is recoverable,
- clear aborts when a state is not recoverable.

Suggested tests:

```text
examples/1D_jwl_bad_state_recovery/case.py
examples/2D_jwl_bad_interface_recovery/case.py
```

These cases should intentionally create difficult mixed cells near a JWL/air interface. The pass condition is not “nothing ever goes wrong.” The pass condition is that the code either recovers cleanly or fails with a useful diagnostic instead of producing silent NaNs.

## Priority 5: ROCFLU-Style Wrapper Workflow

ROCFLU exposes higher-level mixture routines such as `ComputePressureMixt()` and `ComputeEnergyMixt()`. MFC currently has helper routines wired into existing conversion paths, but not an exact wrapper layer.

What to implement:

- a pressure wrapper that owns all JWL mixture pressure guard logic,
- an energy wrapper that owns pressure-to-energy guard logic,
- one place for fallback thresholds,
- one place for diagnostics,
- focused unit-style tests around these wrappers.

This would make the implementation easier to reason about and easier to compare against ROCFLU.

## Priority 6: Multiple JWL Materials

Right now, the code supports one JWL fluid per case and rejects more than one. That is intentional.

To support multiple JWL materials, MFC would need:

- material-index-aware JWL helper calls,
- unambiguous mixture rules when more than one explosive product exists,
- separate reference energies and densities per JWL material,
- validation cases with two distinct JWL products.

Do not remove the current rejection until this is designed and tested.

## Priority 7: Particle Coupling Beyond IBM

ROCFLU's PICL particle-in-cell coupling has not been ported. The current examples use MFC's moving IBM particles. That is useful, but it is not the same model.

Possible future work:

- decide whether JWL should couple to MFC Lagrangian/PIC-style particles,
- define particle force and heat-transfer models,
- add particle escape/removal behavior for open boundaries,
- add restart and post-processing tests for particle state.

## Validation Ideas

Once the code is stable, add comparisons against known solutions or reference calculations:

- 1D JWL expansion into air,
- shock tube with a JWL driver,
- pressure decay from a compact JWL products region,
- 2D blast reflection from a slip wall,
- blast loading on a fixed cylinder before enabling moving IBM,
- moving particle-bed response after the fixed-body force path is trusted.

The current examples are smoke and integration tests. They are not validation data.

## Practical Testing Notes

Use these checks often:

```bash
./mfc.sh build -j 2
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

When debugging the 2D case:

- look at pressure output to confirm the shock forms,
- look at `restart_data/ib_state_*.dat` for force and velocity,
- remember that open fluid boundaries do not remove particles,
- keep the domain large enough that particle escape does not confuse the result,
- use smaller `cfl_target` for very strong blast pressures.

The right development rhythm is: first make the 1D EOS path boring and reliable, then make the 2D fixed-body blast reliable, and only then use the full moving particle-bed case as the stress test.
