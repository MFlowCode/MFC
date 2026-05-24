# JWL Test Suggestions

This file lists suggested tests for improving the MFC JWL implementation toward ROCFLU-style robustness. It is separate from `README-JWL-EOS.md`, which documents the current implementation.

## Test Priorities

Recommended order:

```text
1. Existing 1D and 2D smoke tests
2. Yfix and low-JWL-fraction fallback tests
3. Energy floor tests
4. NaN and bad-state recovery tests
5. ComputePressureMixt / ComputeEnergyMixt wrapper tests
6. Multiple-JWL-fluid validation tests
7. Unsupported-physics rejection tests
8. PICL coupling tests, only if/when particle coupling is added
```

## Current Baseline Tests

### 1D Mixed JWL/Air Smoke Test

Path:

```text
examples/1D_jwl_mixture_test/case.py
```

Purpose:

```text
Verify that mixed JWL/air pressure, energy, and sound-speed paths run in 1D.
```

Expected result:

```text
pre_process, simulation, and post_process complete with exit code 0.
```

Run:

```bash
./mfc.sh run examples/1D_jwl_mixture_test/case.py --no-build
```

### 2D JWL Blast With Light IBM Body

Path:

```text
examples/2D_jwl_mixture_test/case.py
```

Purpose:

```text
Verify mixed JWL/air in 2D with a high-pressure blast, moving-IBM setup, IB-state output, and a light body initialized from rest.
```

Current important values:

```text
driver pressure   5.0e9 Pa
IB mass           0.01
IB initial vel    0
steps             1000
save period       250
```

Expected result:

```text
simulation reaches t_step = 1000
post_process writes output at 0, 250, 500, 750, 1000
IB files are written in restart_data and p_all/p0
```

Run:

```bash
./mfc.sh run examples/2D_jwl_mixture_test/case.py --no-build
```

## Yfix Tests

These tests should be added after implementing a standalone JWL mass-fraction correction helper.

### Tiny JWL Fraction

Suggested case:

```text
examples/1D_jwl_yfix_tiny_fraction/case.py
```

Initial condition:

```text
mostly air
JWL alpha and alpha_rho near machine floor
small high-pressure region
```

Purpose:

```text
Confirm that Y_jwl is clipped to a safe lower bound and the air fallback path is used when Y_jwl is nearly zero.
```

Pass criteria:

```text
no NaN pressure
no negative energy
simulation completes
pressure remains finite and positive
```

### Nearly Pure JWL Fraction

Suggested case:

```text
examples/1D_jwl_yfix_near_one/case.py
```

Initial condition:

```text
JWL alpha = 1 - eps
air alpha = eps
```

Purpose:

```text
Confirm that Y_jwl is clipped below 1 and mixed formulas do not divide by zero or produce NaNs.
```

Pass criteria:

```text
simulation completes
pressure and sound speed remain finite
```

## Low-JWL-Fraction Air Fallback Tests

These tests should prove that cells with almost no explosive products behave like ideal-gas air.

Suggested case:

```text
examples/1D_jwl_low_y_air_fallback/case.py
```

Initial condition:

```text
uniform mostly-air state
Y_jwl <= 1.0e-8
moderate pressure perturbation
```

Expected behavior:

```text
pressure from JWL helper approximately matches air pressure
energy from pressure approximately matches ideal-gas air energy
sound speed remains finite and close to ideal-gas value
```

Useful checks:

```text
compare pressure output against an equivalent stiffened-gas/ideal-gas case
verify relative pressure difference stays small
```

## Energy Floor Tests

These tests should be added after implementing a dedicated JWL energy floor.

### Low Internal Energy Recovery

Suggested case:

```text
examples/1D_jwl_energy_floor/case.py
```

Initial condition:

```text
JWL-rich region with pressure chosen close to the minimum safe pressure
mostly-air region with low internal energy
```

Purpose:

```text
Confirm that energy is raised to a physically safe floor before pressure and sound-speed calculations.
```

Pass criteria:

```text
no negative pressure
no negative sound-speed squared
no NaNs
simulation completes
```

### Shock With Energy Floor

Suggested case:

```text
examples/1D_jwl_energy_floor_shock/case.py
```

Initial condition:

```text
strong pressure jump
thin mixed cells between JWL and air
```

Purpose:

```text
Confirm the floor does not break shock propagation and only activates in problematic cells.
```

Pass criteria:

```text
shock remains visible
pressure remains finite
run completes without fallback dominating the whole domain
```

## NaN And Bad-State Recovery Tests

These tests intentionally create difficult states to check recovery logic.

### Bad Mixture Cell

Suggested case:

```text
examples/1D_jwl_bad_state_recovery/case.py
```

Initial condition:

```text
one or more cells near the JWL/air interface have extreme mass-fraction imbalance
large pressure jump
```

Purpose:

```text
Confirm recovery replaces bad pressure, energy, temperature, or sound-speed values with safe finite values.
```

Pass criteria:

```text
simulation does not abort
no NaN appears in output
log or counter reports fallback activation if such diagnostics are added
```

### Negative Pressure Guard

Suggested unit-style check:

```text
directly call JWL pressure helper with an energy low enough to produce negative pressure
```

Purpose:

```text
Confirm pressure is floored or recovered consistently.
```

Pass criteria:

```text
returned pressure is finite and positive
```

## ComputePressureMixt / ComputeEnergyMixt Wrapper Tests

After adding wrapper helpers, test them directly and through full MFC cases.

### Round-Trip Pressure/Energy

Suggested unit-style check:

```text
given rho, pressure, and Y
compute e = EnergyMixt(rho, pressure, Y)
then compute p2 = PressureMixt(rho, e, Y)
```

Pass criteria:

```text
p2 approximately equals input pressure
temperature is finite
sound speed is finite and positive
```

Test several values:

```text
Y = 1.0e-10
Y = 1.0e-6
Y = 0.5
Y = 1.0 - 1.0e-10
```

### Full Conversion Round Trip

Suggested case:

```text
examples/1D_jwl_conversion_round_trip/case.py
```

Purpose:

```text
Exercise primitive-to-conservative and conservative-to-primitive conversion repeatedly during a short run.
```

Pass criteria:

```text
pressure remains finite and positive
mass fractions remain bounded
energy remains finite
```

## Multiple JWL Fluid Tests

These tests are for a future implementation that removes the single `jwl_idx` assumption.

### Rejection Test

Current expected behavior:

```text
case with two fluids where both have eos = 2 should abort cleanly
```

Suggested case:

```text
examples/negative_tests/2D_two_jwl_fluids_should_fail/case.py
```

Pass criteria:

```text
clear error message:
JWL EOS currently supports one JWL fluid per case.
```

### Identical-Parameter Multi-JWL Test

Future expected behavior:

```text
two JWL fluids with identical JWL constants may be allowed
total JWL mass fraction is sum of both JWL mass fractions
```

Pass criteria:

```text
result matches equivalent one-JWL-fluid case within tolerance
```

## Unsupported Physics Rejection Tests

These should remain negative tests until real coupling is implemented.

Suggested cases:

```text
JWL + MHD should fail
JWL + chemistry should fail
JWL + Eulerian bubbles should fail
JWL + model_eqns = 4 should fail
```

Pass criteria:

```text
case aborts before simulation with a clear error message
no segmentation fault
no silent fallback to wrong physics
```

Expected messages should match the guards in:

```text
src/common/m_variables_conversion.fpp
```

## PICL Coupling Tests

Only add these if MFC gets a PICL or particle-in-cell coupling path for JWL.

Suggested future tests:

```text
single JWL particle in uniform air
JWL particle cloud depositing energy into air
particle-to-grid mass conservation test
particle-to-grid energy conservation test
grid-to-particle pressure update test
```

Pass criteria:

```text
total mass conserved within tolerance
total energy conserved or changes only by documented source terms
particle pressure remains finite
grid pressure remains finite
no NaNs in particle or grid state
```

## Diagnostics To Add

These diagnostics would make the tests easier to judge:

```text
number of Yfix activations
number of low-Y air fallbacks
number of energy-floor activations
number of bad-state recoveries
minimum pressure per time step
minimum internal energy per time step
minimum sound-speed squared before clipping
maximum IBM force
IBM center-of-mass displacement
IBM velocity history
```

For the 2D IBM blast case, especially useful outputs are:

```text
ib_state_*.dat
ib_data.dat
pressure field
density field
JWL volume fraction or mass fraction
```

## Suggested Acceptance Levels

Use three levels of confidence:

```text
Smoke
  Case runs to completion and writes outputs.

Robustness
  Case survives extreme but expected mixed cells without NaNs, negative pressure, or negative sound-speed squared.

Regression
  Key outputs are compared against a stored reference or analytic/near-analytic expectation.
```

The current 1D and 2D cases are smoke tests. The next goal should be robustness tests for Yfix, low-JWL-fraction fallback, energy floors, and bad-state recovery.

## Implementation Roadmap

This section describes how to implement the missing ROCFLU-style robustness features in MFC. The goal is to add the stability layer in small, testable pieces instead of making one large thermodynamics change.

## ROCFLU-To-MFC Porting Map

The ROCFLU JWL implementation is centered around a dedicated thermodynamics module:

```text
ROCFLU:
modflu/RFLU_ModJWL.F90
```

The matching MFC implementation lives mainly in:

```text
MFC:
src/common/m_variables_conversion.fpp
```

The MFC file is the right place because pressure, energy, sound speed, and primitive/conservative conversion already pass through this module. The remaining work should continue there unless the feature requires a different physics subsystem.

### Routine Mapping

```text
ROCFLU routine                    MFC equivalent/status
-----------------------------------------------------------------------------------------------
RFLU_JWL_P_ER                     s_jwl_pressure_er
                                  Implemented for mixed JWL/air, but needs stronger Yfix,
                                  energy-floor, and bad-state recovery wrapping.

RFLU_JWL_E_PR                     s_jwl_energy_pr
                                  Implemented for pressure-to-energy conversion, but needs
                                  centralized recovery and floor enforcement.

RFLU_JWL_T_PR                     s_jwl_temperature_pr
                                  Implemented as a helper. Needs to be called consistently from
                                  a higher-level ComputeEnergyMixt/ComputePressureMixt wrapper.

RFLU_JWL_C_ER                     s_jwl_sound_speed_squared
                                  s_jwl_mixture_sound_speed_squared
                                  Implemented. Needs robust clipping/recovery for negative c^2.

RFLU_JWL_ComputePressureMixt      Not yet implemented as one ROCFLU-style wrapper.
                                  Current MFC calls lower-level helpers directly from conversion
                                  paths. Add s_jwl_compute_pressure_mixt.

RFLU_JWL_ComputeEnergyMixt        Not yet implemented as one ROCFLU-style wrapper.
                                  Current MFC calls s_jwl_energy_pr directly. Add
                                  s_jwl_compute_energy_mixt.

RFLU_JWL_Yfix                     Not yet implemented as a standalone helper.
                                  Add s_jwl_yfix and call it everywhere Y_jwl is used.

ROCFLU NaN/recovery guards        Partially covered by parameter validation only.
                                  Add explicit finite-value checks and fallback/recovery helpers.

ROCFLU energy floor behavior      Partially covered by formulas and positive parameter checks.
                                  Add s_jwl_energy_floor and use it before/after JWL EOS calls.

ROCFLU low-JWL-fraction fallback  Not yet complete.
                                  Add air fallback when Y_jwl is near zero.

ROCFLU PICL coupling              Not ported.
                                  Requires MFC particle/PICL infrastructure design first.
```

### What Has Already Been Ported

MFC already has the core thermodynamic and case-parameter plumbing:

```text
1. EOS selector
   fluid_pp(i)%eos = 2 selects JWL.

2. JWL constants
   fluid_pp(i)%jwl_A
   fluid_pp(i)%jwl_B
   fluid_pp(i)%jwl_R1
   fluid_pp(i)%jwl_R2
   fluid_pp(i)%jwl_omega
   fluid_pp(i)%jwl_rho0

3. Mixed JWL/air constants
   fluid_pp(i)%jwl_E0
   fluid_pp(i)%jwl_air_e0
   fluid_pp(i)%jwl_air_rho0
   fluid_pp(i)%jwl_air_gamma

4. Parameter defaults
   pre_process, simulation, and post_process global-parameter modules.

5. MPI broadcast
   pre_process, simulation, and post_process MPI proxy modules.

6. GPU-visible arrays
   jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s,
   jwl_E0s, jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas.

7. Core helpers
   s_jwl_pcold
   s_jwl_dpcold_drho
   s_jwl_sound_speed_squared
   s_jwl_mixture_sound_speed_squared
   s_jwl_pressure_er
   s_jwl_energy_pr
   s_jwl_temperature_pr

8. MFC conversion integration
   s_compute_pressure
   s_convert_conservative_to_primitive_variables
   s_convert_primitive_to_conservative_variables
   s_convert_primitive_to_flux_variables
   s_compute_speed_of_sound

9. Toolchain/schema support
   toolchain/mfc/params/definitions.py
   toolchain/mfc/params/descriptions.py

10. Smoke tests
   examples/1D_jwl_mixture_test
   examples/2D_jwl_mixture_test
```

### What Is Still Left

The remaining work is mostly robustness and full ROCFLU workflow parity:

```text
1. Standalone Yfix helper
   Current status: missing.
   Needed because every mixed-cell JWL formula should receive a bounded mass fraction.

2. Low-JWL-fraction air fallback
   Current status: incomplete.
   Needed because mostly-air cells should not be forced through explosive-product math.

3. Dedicated energy floor
   Current status: incomplete.
   Needed because pressure recovery from energy can become nonphysical in difficult mixed cells.

4. NaN and bad-state recovery
   Current status: incomplete.
   Needed because a single bad cell can otherwise poison pressure, energy, or sound speed.

5. ComputePressureMixt wrapper
   Current status: missing.
   Needed to match ROCFLU's centralized pressure workflow.

6. ComputeEnergyMixt wrapper
   Current status: missing.
   Needed to match ROCFLU's centralized pressure-to-energy workflow.

7. Diagnostics counters
   Current status: missing.
   Needed to know when Yfix, fallback, flooring, or recovery actually activate.

8. Multiple JWL fluids
   Current status: rejected.
   Needed only if cases require more than one JWL material.

9. MHD, chemistry, Eulerian bubbles, model_eqns = 4
   Current status: rejected.
   Needed only after separate thermodynamic coupling work.

10. PICL coupling
   Current status: not ported.
   Needed only if MFC particle-in-cell JWL coupling is required.
```

## Detailed Porting Approach

The safest way to port the remaining ROCFLU behavior is to make MFC's JWL path look like ROCFLU's layered design:

```text
raw formulas
  -> Yfix
  -> low-Y fallback
  -> energy floor
  -> pressure/energy computation
  -> temperature and sound speed
  -> bad-state recovery
  -> diagnostics
```

### Target MFC Call Structure

Current MFC flow is roughly:

```text
conversion routine
  -> compute Y_jwl locally
  -> call s_jwl_pressure_er or s_jwl_energy_pr directly
```

Target flow should be:

```text
conversion routine
  -> compute raw Y_jwl locally
  -> call s_jwl_compute_pressure_mixt or s_jwl_compute_energy_mixt
       -> Yfix
       -> low-Y air fallback
       -> energy floor
       -> JWL pressure or energy
       -> temperature
       -> sound speed
       -> recovery
```

This keeps the safety behavior identical whether the caller is primitive initialization, conservative recovery, flux conversion, or sound-speed calculation.

### Porting Rule

Do not spread ROCFLU-style recovery logic across many call sites. Put the protection inside the JWL wrapper helpers, then make all MFC conversion paths call those helpers.

Good:

```text
s_compute_pressure
  -> s_jwl_compute_pressure_mixt
```

Riskier:

```text
s_compute_pressure
  -> local Yfix
  -> local energy floor
  -> local pressure recovery
  -> s_jwl_pressure_er
```

The wrapper approach is easier to audit and much closer to ROCFLU's `ComputePressureMixt()` and `ComputeEnergyMixt()` style.

### MFC Data Translation

ROCFLU and MFC store the same physics in different shapes.

In MFC:

```text
rho = mixture density
Y_jwl = JWL partial density / mixture density
e = specific internal energy
rho*e = total energy per volume - kinetic energy per volume
pressure is primitive variable eqn_idx%E in many MFC primitive arrays
```

When entering JWL helpers:

```fortran
rho_safe = max(rho, sgm_eps)
Y_raw = alpha_rho_jwl/max(rho_safe, sgm_eps)
```

Then the wrapper should call:

```fortran
call s_jwl_yfix(Y_raw, Y_safe)
```

When converting primitive pressure to conservative energy:

```fortran
call s_jwl_compute_energy_mixt(rho, pres, Y_raw, ..., e, T, c2)
energy = rho*e + kinetic_energy
```

When converting conservative energy to primitive pressure:

```fortran
e = (total_energy - kinetic_energy)/rho
call s_jwl_compute_pressure_mixt(rho, e, Y_raw, ..., pres, T, c2)
```

### Parameter Translation

ROCFLU-style JWL parameter names map to MFC like this:

```text
ROCFLU prepRealVal14  -> fluid_pp(i)%jwl_A
ROCFLU prepRealVal15  -> fluid_pp(i)%jwl_B
ROCFLU prepRealVal17  -> fluid_pp(i)%jwl_omega
ROCFLU prepRealVal18  -> fluid_pp(i)%jwl_R1
ROCFLU prepRealVal19  -> fluid_pp(i)%jwl_R2
ROCFLU prepRealVal24  -> fluid_pp(i)%jwl_rho0
ROCFLU ETNT/E0        -> fluid_pp(i)%jwl_E0
```

MFC also stores these air-reference parameters for the mixed JWL/air path:

```text
fluid_pp(i)%jwl_air_e0
fluid_pp(i)%jwl_air_rho0
fluid_pp(i)%jwl_air_gamma
```

### Done Criteria For ROCFLU Parity

The MFC port should be considered ROCFLU-style robust when:

```text
1. All JWL mass fractions pass through s_jwl_yfix.
2. Mostly-air cells use the air fallback.
3. Energy is floored before pressure evaluation and after pressure-to-energy conversion.
4. Pressure, energy, temperature, and sound speed are recovered if NaN or nonphysical.
5. Pressure and energy wrappers are the only public JWL entry points used by conversion routines.
6. Diagnostics can report how often each guard activates.
7. Existing 1D/2D smoke tests still pass.
8. New robustness tests pass for tiny-Y, near-one-Y, low-energy, and bad-state cases.
```

Recommended implementation order:

```text
1. Add reusable JWL safety constants and helper routines
2. Add Yfix
3. Add low-JWL-fraction air fallback
4. Add energy floor
5. Add bad-state and NaN recovery
6. Add ComputePressureMixt / ComputeEnergyMixt wrappers
7. Add diagnostics counters
8. Add tests for each feature
9. Consider multiple JWL fluids
10. Consider PICL and unsupported-physics coupling
```

## Step 1: JWL Safety Constants

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add constants near the existing module variables:

```fortran
real(wp), parameter :: jwl_y_min = 1.e-8_wp
real(wp), parameter :: jwl_y_max = 1._wp - 1.e-8_wp
real(wp), parameter :: jwl_p_min = 1.e-12_wp
real(wp), parameter :: jwl_e_min = 1.e-12_wp
real(wp), parameter :: jwl_c2_min = 1.e-12_wp
```

Implementation notes:

```text
jwl_y_min protects air-dominated cells.
jwl_y_max protects pure-JWL cells from divide-by-zero style singularities.
jwl_p_min prevents negative or zero pressure after recovery.
jwl_e_min prevents nonphysical internal energy.
jwl_c2_min prevents negative sound-speed squared from reaching sqrt().
```

Validation:

```text
Build MFC.
Run existing 1D and 2D JWL cases.
No behavior should change yet if constants are unused.
```

## Step 2: Yfix Helper

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add a small helper:

```fortran
subroutine s_jwl_yfix(Y_in, Y_out)
    real(wp), intent(in)  :: Y_in
    real(wp), intent(out) :: Y_out

    if (Y_in /= Y_in) then
        Y_out = jwl_y_min
    else
        Y_out = min(jwl_y_max, max(jwl_y_min, Y_in))
    end if
end subroutine s_jwl_yfix
```

If this helper is called inside GPU-decorated routines, add the same GPU routine macro style used by the existing JWL helpers:

```fortran
$:GPU_ROUTINE(function_name='s_jwl_yfix',parallelism='[seq]', cray_noinline=True)
```

Call sites:

```text
s_compute_pressure
s_compute_speed_of_sound
s_convert_conservative_to_primitive_variables
s_convert_primitive_to_conservative_variables
s_convert_primitive_to_flux_variables
s_jwl_pressure_er
s_jwl_energy_pr
s_jwl_temperature_pr
s_jwl_mixture_sound_speed_squared
```

Pattern:

```fortran
call s_jwl_yfix(Y_jwl, Y_jwl)
```

Validation tests:

```text
examples/1D_jwl_yfix_tiny_fraction/case.py
examples/1D_jwl_yfix_near_one/case.py
```

Pass criteria:

```text
Y never leaves [jwl_y_min, jwl_y_max] inside JWL helpers.
Pressure, energy, temperature, and sound speed remain finite.
Existing smoke tests still pass.
```

## Step 3: Low-JWL-Fraction Air Fallback

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add helper routines:

```fortran
subroutine s_jwl_air_pressure_er(rho, e, air_gamma, pres)
    real(wp), intent(in)  :: rho, e, air_gamma
    real(wp), intent(out) :: pres

    pres = max(jwl_p_min, air_gamma*rho*max(e, jwl_e_min))
end subroutine s_jwl_air_pressure_er

subroutine s_jwl_air_energy_pr(rho, pres, air_gamma, e)
    real(wp), intent(in)  :: rho, pres, air_gamma
    real(wp), intent(out) :: e

    e = max(jwl_e_min, max(pres, jwl_p_min)/(air_gamma*max(rho, sgm_eps)))
end subroutine s_jwl_air_energy_pr
```

Then use them inside:

```text
s_jwl_pressure_er
s_jwl_energy_pr
s_jwl_temperature_pr
s_jwl_mixture_sound_speed_squared
```

Pattern:

```fortran
call s_jwl_yfix(Y, Y_safe)

if (Y_safe <= jwl_y_min) then
    call s_jwl_air_pressure_er(rho_safe, e_safe, air_gamma, pres)
    return
end if
```

For energy:

```fortran
if (Y_safe <= jwl_y_min) then
    call s_jwl_air_energy_pr(rho_safe, pres, air_gamma, e)
    return
end if
```

Sound-speed fallback:

```fortran
if (Y_safe <= jwl_y_min) then
    c2 = max(jwl_c2_min, (1._wp + air_gamma)*max(pres, jwl_p_min)/max(rho_safe, sgm_eps))
    return
end if
```

Validation tests:

```text
examples/1D_jwl_low_y_air_fallback/case.py
```

Pass criteria:

```text
mostly-air JWL-mixed cells behave close to ideal gas
pressure remains positive
energy remains positive
sound speed remains positive
```

## Step 4: Energy Floor

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add helper:

```fortran
subroutine s_jwl_energy_floor(e_in, Y, E0, air_e0, e_out)
    real(wp), intent(in)  :: e_in, Y, E0, air_e0
    real(wp), intent(out) :: e_out
    real(wp) :: e_ref

    e_ref = (1._wp - Y)*air_e0 + Y*max(E0, air_e0)
    e_out = max(e_in, jwl_e_min*max(e_ref, 1._wp))
end subroutine s_jwl_energy_floor
```

Call it before pressure evaluation:

```fortran
call s_jwl_energy_floor(e, Y_safe, E0, air_e0, e_safe)
call s_jwl_pressure_er(..., e_safe, ...)
```

Call it after pressure-to-energy conversion:

```fortran
call s_jwl_energy_pr(..., e)
call s_jwl_energy_floor(e, Y_safe, E0, air_e0, e)
```

Implementation note:

```text
The first floor can be conservative. Start with a tiny floor to avoid changing normal blast behavior.
Increase it only if the robustness tests still produce bad states.
```

Validation tests:

```text
examples/1D_jwl_energy_floor/case.py
examples/1D_jwl_energy_floor_shock/case.py
```

Pass criteria:

```text
no negative pressure
no negative internal energy
no NaNs
shock remains visible in the pressure field
```

## Step 5: Bad-State And NaN Recovery

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add logical helpers:

```fortran
logical function f_jwl_is_bad_scalar(a)
    real(wp), intent(in) :: a

    f_jwl_is_bad_scalar = (a /= a)
end function f_jwl_is_bad_scalar

logical function f_jwl_bad_state(rho, e, pres, Y)
    real(wp), intent(in) :: rho, e, pres, Y

    f_jwl_bad_state = rho <= 0._wp .or. e <= 0._wp .or. pres <= 0._wp &
        .or. rho /= rho .or. e /= e .or. pres /= pres .or. Y /= Y
end function f_jwl_bad_state
```

Add recovery helper:

```fortran
subroutine s_jwl_recover_pressure(rho, e, air_gamma, pres)
    real(wp), intent(in)    :: rho, e, air_gamma
    real(wp), intent(inout) :: pres

    if (pres /= pres .or. pres <= jwl_p_min) then
        pres = max(jwl_p_min, air_gamma*max(rho, sgm_eps)*max(e, jwl_e_min))
    end if
end subroutine s_jwl_recover_pressure
```

Use after pressure calculation:

```fortran
call s_jwl_pressure_er(...)
call s_jwl_recover_pressure(rho, e_safe, air_gamma, pres)
```

Use similar recovery after energy and sound-speed calculations:

```fortran
if (e /= e .or. e <= jwl_e_min) e = jwl_e_min
if (c2 /= c2 .or. c2 <= jwl_c2_min) c2 = jwl_c2_min
```

Validation tests:

```text
examples/1D_jwl_bad_state_recovery/case.py
unit-style negative pressure guard check
```

Pass criteria:

```text
bad cells recover to finite positive pressure and energy
simulation completes
fallback activation can be counted if diagnostics are added
```

## Step 6: ComputePressureMixt / ComputeEnergyMixt Wrappers

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Add wrapper helpers that mirror ROCFLU's workflow more closely:

```fortran
subroutine s_jwl_compute_pressure_mixt(rho, e, Y, A, B, R1, R2, omega, rho0, E0, air_e0, air_rho0, air_gamma, cv, pres, T, c2)
```

Workflow:

```text
1. clamp rho to a safe value
2. call Yfix
3. apply energy floor
4. if Y is tiny, use air fallback
5. compute pressure
6. recover bad pressure
7. compute temperature
8. compute sound-speed squared
9. recover bad temperature or sound speed
```

Energy wrapper:

```fortran
subroutine s_jwl_compute_energy_mixt(rho, pres, Y, A, B, R1, R2, omega, rho0, E0, air_e0, air_rho0, air_gamma, cv, e, T, c2)
```

Workflow:

```text
1. clamp rho and pressure
2. call Yfix
3. if Y is tiny, use air fallback
4. compute JWL energy from pressure
5. apply energy floor
6. compute temperature
7. compute sound-speed squared
8. recover bad values
```

Then replace direct lower-level calls in conversion routines with the wrappers:

```text
s_compute_pressure
s_convert_primitive_to_conservative_variables
s_convert_primitive_to_flux_variables
s_compute_speed_of_sound
```

This makes all conversion paths use the same safeguards.

Validation tests:

```text
round-trip pressure/energy test
existing 1D and 2D smoke tests
Yfix and energy-floor tests
```

Pass criteria:

```text
PressureMixt(EnergyMixt(p)) approximately returns p for representative Y values.
All full cases run cleanly.
```

## Step 7: Diagnostics Counters

Primary file:

```text
src/common/m_variables_conversion.fpp
```

Optional counters:

```fortran
integer :: jwl_yfix_count
integer :: jwl_air_fallback_count
integer :: jwl_energy_floor_count
integer :: jwl_recovery_count
```

Add increments inside the helpers:

```fortran
if (Y_out /= Y_in) jwl_yfix_count = jwl_yfix_count + 1
```

Implementation caution:

```text
For GPU and MPI runs, simple module counters may not be enough.
Start with CPU-only diagnostics or debug-only counters.
Then add MPI reductions if the counters become part of normal output.
```

Possible output locations:

```text
run_time_info
time_data.dat
dedicated jwl_diagnostics.dat
```

Validation:

```text
Yfix test should produce nonzero Yfix count.
Low-Y fallback test should produce nonzero air fallback count.
Normal smoke tests should produce small or zero recovery count.
```

## Step 8: Multiple JWL Fluids

Current implementation uses one global JWL fluid index:

```text
jwl_idx
```

Safe first extension:

```text
allow multiple JWL fluids only if their JWL constants are identical
```

Implementation approach:

```text
1. Keep the current rejection as the default.
2. Add a parameter such as allow_multiple_identical_jwl if desired.
3. During initialization, find all eos = 2 fluids.
4. Verify their A, B, R1, R2, omega, rho0, E0, air_e0, air_rho0, and air_gamma match within tolerance.
5. Replace Y_jwl = q(jwl_idx)/rho with Y_jwl_total = sum(q(i)/rho for all JWL fluids).
6. Use the shared JWL constants.
```

Do not mass-average unlike JWL constants unless the mixture closure is clearly defined.

Validation tests:

```text
negative two-JWL-fluid test still fails by default
identical two-JWL-fluid case matches equivalent single-JWL case if extension is enabled
```

## Step 9: Unsupported Physics Coupling

Keep these as explicit rejections until each coupling is designed.

### MHD

Implementation needs:

```text
subtract magnetic energy before passing internal energy to JWL
include magnetic pressure consistently in total pressure paths
verify MHD fluxes use JWL thermodynamic pressure correctly
```

Suggested first test:

```text
uniform JWL material with weak magnetic field
no velocity
pressure remains constant
```

### Chemistry

Implementation needs:

```text
define whether JWL E0 replaces or adds to chemical heat release
avoid double-counting explosive energy
couple species reaction source terms to JWL energy consistently
```

Suggested first test:

```text
chemistry disabled but species present
then simple one-step energy release with known final pressure
```

### Eulerian Bubbles

Implementation needs:

```text
new pressure-equilibrium closure between JWL mixture and bubble phase
new sound-speed model
careful volume-fraction bounds
```

Suggested first test:

```text
static uniform JWL/air/bubble mixture remains stable
```

### model_eqns = 4

Implementation needs:

```text
wire JWL into the four-equation model pressure closure
audit volume-fraction and partial-density assumptions
adapt primitive/conservative conversion for that model
```

Suggested first test:

```text
1D two-material interface with no velocity should remain stable
```

## Step 10: PICL Coupling

Only start this if MFC has or receives a particle-in-cell path that needs JWL.

Implementation approach:

```text
1. Identify particle state variables: mass, density, velocity, energy, pressure, material id, JWL constants.
2. Add JWL EOS evaluation for particle pressure and temperature.
3. Add grid-to-particle interpolation for density, pressure, velocity, and energy as needed.
4. Add particle-to-grid deposition for mass, momentum, and energy.
5. Ensure particle and grid JWL parameters use the same units and reference values.
6. Add conservation diagnostics for mass, momentum, and energy.
```

Validation tests:

```text
single stationary JWL particle in uniform air
JWL particle cloud in air
particle energy deposition into one grid cell
particle energy deposition into multiple grid cells
```

Pass criteria:

```text
mass conserved within tolerance
energy conserved or changes only by documented source terms
particle pressure finite
grid pressure finite
no NaNs in particle or grid state
```

## Suggested Development Branches

Keep the work split into small branches:

```text
jwl-yfix
jwl-air-fallback
jwl-energy-floor
jwl-bad-state-recovery
jwl-mixture-wrappers
jwl-diagnostics
jwl-multi-fluid
```

Each branch should add or update at least one focused test case before moving to the next feature.

## Feature-By-Feature Implementation Plan

This section restates each missing ROCFLU capability as an implementation task.

### 1. RFLU_JWL_Yfix Equivalent

Status:

```text
Not implemented as a standalone helper.
```

Goal:

```text
Every JWL helper receives a finite mass fraction clipped to a safe interval.
```

Files to edit:

```text
src/common/m_variables_conversion.fpp
```

Implementation:

```text
1. Add s_jwl_yfix(Y_in, Y_out).
2. Add GPU routine decoration if it is called from GPU routines.
3. Replace direct use of raw Y_jwl with Y_safe.
4. Call s_jwl_yfix inside the new pressure/energy wrappers.
5. Keep local call-site Y computations simple; let the wrapper handle correction.
```

Important call sites:

```text
conservative-to-primitive pressure recovery
primitive-to-conservative energy initialization
primitive-to-flux energy conversion
sound-speed calculation
temperature calculation
```

Tests:

```text
tiny JWL fraction
nearly pure JWL fraction
mixed-cell shock
```

### 2. ROCFLU-Style ComputePressureMixt

Status:

```text
Lower-level pressure helper exists, but no centralized mixture-pressure workflow exists yet.
```

Goal:

```text
All pressure recovery from JWL states should go through one routine.
```

Suggested helper:

```fortran
subroutine s_jwl_compute_pressure_mixt(rho, e, Y, A, B, R1, R2, omega, rho0, E0, air_e0, air_rho0, air_gamma, cv, pres, T, c2)
```

Internal order:

```text
1. rho_safe = max(rho, sgm_eps)
2. call s_jwl_yfix(Y, Y_safe)
3. call s_jwl_energy_floor(e, Y_safe, E0, air_e0, e_safe)
4. if Y_safe is tiny, compute air pressure and skip JWL pressure
5. otherwise call s_jwl_pressure_er
6. recover pressure if bad
7. compute temperature
8. compute c2
9. recover T and c2 if bad
```

MFC call sites to replace:

```text
s_compute_pressure
s_convert_conservative_to_primitive_variables
```

Test:

```text
Use a strong mixed JWL/air shock and confirm no direct lower-level pressure calls bypass the wrapper.
```

### 3. ROCFLU-Style ComputeEnergyMixt

Status:

```text
s_jwl_energy_pr exists, but no centralized mixture-energy workflow exists yet.
```

Goal:

```text
All pressure-to-energy initialization and flux-state energy conversion should go through one routine.
```

Suggested helper:

```fortran
subroutine s_jwl_compute_energy_mixt(rho, pres, Y, A, B, R1, R2, omega, rho0, E0, air_e0, air_rho0, air_gamma, cv, e, T, c2)
```

Internal order:

```text
1. rho_safe = max(rho, sgm_eps)
2. pres_safe = max(pres, jwl_p_min)
3. call s_jwl_yfix(Y, Y_safe)
4. if Y_safe is tiny, compute air energy
5. otherwise call s_jwl_energy_pr
6. apply energy floor
7. compute temperature
8. compute c2
9. recover e, T, and c2 if bad
```

MFC call sites to replace:

```text
s_convert_primitive_to_conservative_variables
s_convert_primitive_to_flux_variables
```

Test:

```text
Pressure-energy-pressure round trip for Y = tiny, moderate, and nearly one.
```

### 4. NaN Recovery And Bad-State Handling

Status:

```text
Only basic parameter validation exists.
```

Goal:

```text
No NaN or negative pressure/energy/sound-speed value should leave a JWL wrapper.
```

Implementation:

```text
1. Add scalar finite checks using a /= a for NaN.
2. Add f_jwl_bad_state for rho, e, pressure, and Y.
3. Add pressure recovery using air fallback pressure.
4. Add energy recovery using air fallback energy or energy floor.
5. Add c2 recovery using jwl_c2_min.
6. Add optional counters to measure each recovery.
```

Where recovery should happen:

```text
inside s_jwl_compute_pressure_mixt
inside s_jwl_compute_energy_mixt
inside s_jwl_mixture_sound_speed_squared if it remains callable directly
```

Tests:

```text
bad-state recovery case
negative pressure guard
strong mixed-cell shock
```

### 5. Energy Floor Enforcement

Status:

```text
No dedicated ROCFLU-style energy floor helper yet.
```

Goal:

```text
Prevent JWL pressure from being evaluated with nonphysical internal energy.
```

Implementation:

```text
1. Add s_jwl_energy_floor.
2. Use a conservative floor first, based on air_e0 and E0.
3. Apply before pressure-from-energy.
4. Apply after energy-from-pressure.
5. Add a diagnostic counter for floor activations.
```

Suggested first floor:

```text
e_floor = jwl_e_min * max((1 - Y)*air_e0 + Y*max(E0, air_e0), 1)
e = max(e, e_floor)
```

Tests:

```text
low-energy JWL state
strong shock with thin mixed interface
```

### 6. Low-JWL-Fraction Ideal-Gas Fallback

Status:

```text
Not complete.
```

Goal:

```text
Mostly-air cells should behave like air, not explosive products.
```

Implementation:

```text
1. Choose jwl_y_min.
2. If Y_safe <= jwl_y_min, use air pressure, energy, temperature, and sound speed.
3. Count fallback activations.
4. Keep fallback inside wrappers so all call sites behave the same.
```

Air formulas:

```text
p = air_gamma * rho * e
e = p / (air_gamma * rho)
c^2 = (1 + air_gamma) * p / rho
T = e / cv, if cv is positive
```

Tests:

```text
mostly-air JWL case
comparison with equivalent ideal-gas air case
```

### 7. PICL Particle-In-Cell Coupling

Status:

```text
Not ported.
```

Reason:

```text
PICL is not just an EOS formula. It needs particle state, grid-particle interpolation, particle-grid deposition, and conservation diagnostics.
```

Implementation approach:

```text
1. First identify whether MFC already has the required particle/PICL infrastructure.
2. Add JWL material id and JWL constants to particle state.
3. Add particle pressure/temperature update using JWL wrappers.
4. Add grid-to-particle interpolation for required fields.
5. Add particle-to-grid deposition for mass, momentum, and energy.
6. Ensure deposited energy is compatible with MFC mixture energy.
7. Add conservation diagnostics.
```

Do this only after:

```text
Yfix, fallback, energy floor, and recovery wrappers are complete.
```

### 8. More Than One JWL Fluid

Status:

```text
Currently rejected.
```

Reason:

```text
The current MFC implementation stores one active JWL index in jwl_idx.
```

Safe first implementation:

```text
Allow multiple JWL fluids only if all JWL constants are identical.
```

Implementation:

```text
1. Keep the current rejection as default behavior.
2. Add optional support for identical-parameter JWL fluids.
3. During initialization, collect all eos = 2 fluid indices.
4. Verify their JWL constants match within tolerance.
5. Replace Y_jwl with Y_jwl_total.
6. Use the shared JWL constants.
```

Do not do this first:

```text
Do not mass-average different A, B, R1, R2, omega, rho0, or E0 values unless a valid closure is defined.
```

Tests:

```text
two JWL fluids should fail by default
two identical JWL fluids should match one equivalent JWL fluid if optional support is enabled
```

### 9. MHD Support

Status:

```text
Currently rejected.
```

Reason:

```text
JWL pressure must use internal energy after subtracting kinetic and magnetic energy. MHD fluxes also need thermodynamic pressure consistently separated from magnetic pressure.
```

Implementation:

```text
1. Audit MFC MHD energy definition.
2. Compute internal energy as total - kinetic - magnetic.
3. Pass only thermodynamic internal energy to JWL.
4. Ensure MHD pressure terms use JWL thermodynamic pressure correctly.
5. Add uniform-state and shock tests with weak magnetic fields.
```

### 10. Chemistry Support

Status:

```text
Currently rejected.
```

Reason:

```text
JWL E0 and chemical heat release can double-count explosive energy if coupled naively.
```

Implementation:

```text
1. Define the energy model: JWL E0, chemistry heat release, or both with a clear reference.
2. Audit species energy source terms.
3. Make pressure recovery use the chemically updated internal energy.
4. Add one-step reaction tests with known pressure response.
```

### 11. Eulerian Bubble Support

Status:

```text
Currently rejected.
```

Reason:

```text
Bubble models use a different mixture pressure and sound-speed closure.
```

Implementation:

```text
1. Derive pressure equilibrium between JWL mixture and bubble phase.
2. Define the correct sound speed for the coupled system.
3. Add volume-fraction bounds and recovery.
4. Start with static uniform tests before shocks.
```

### 12. model_eqns = 4 Support

Status:

```text
Currently rejected.
```

Reason:

```text
The four-equation model has different volume-fraction and pressure-equilibrium assumptions than the current JWL mixed-material path.
```

Implementation:

```text
1. Audit model_eqns = 4 primitive and conservative variable definitions.
2. Decide where JWL pressure equilibrium belongs.
3. Wire JWL wrappers into that model's conversion routines.
4. Add static-interface tests before shock tests.
```
