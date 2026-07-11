# 3D JWL++ Reactive Detonation Bursting into Air

A full 3D **self-propagating detonation**: a spherical TNT charge (ρ₀ = 1630
kg/m³, r = 5 cm) at rest in free air is ignited by a small central hot spot. A
JWL++ reactive burn (Souers 2000),

```
dλ/dt = jwl_G · p^jwl_b_exp · (1 − λ)
```

sweeps outward through the charge at the TNT detonation (CJ) speed, converting
unreacted explosive to hot JWL products, which then vent and mix into the
surrounding air. It is the 3D analogue of `examples/2D_jwl_reactive_air_mixing`,
using the same TNT JWL products EOS and the `Y·(1−λ)` reaction gate that keeps
the air inert.

## Run it (≈ real-world detonation in a few minutes on a laptop)

The case is sized for a **multicore CPU with MPI** — it targets < 10 min on an
8-core machine:

```bash
./mfc.sh build --mpi -j 8                              # one-time (MPI build)
./mfc.sh run examples/3D_jwl_reactive_air_burst/case.py --mpi -n 8
./examples/3D_jwl_reactive_air_burst/make_movie.py examples/3D_jwl_reactive_air_burst
```

`make_movie.py` extracts the **z = 0 mid-plane** from the per-rank 3D output and
renders a two-panel contour movie, `detonation_burst.mp4`:

- **left** — pressure [GPa] (log, `inferno`): the detonation front and the blast
  vented into air;
- **right** — TNT/products mass fraction (`magma`) with the λ = ½ reaction
  contour: the burn front consuming the charge and the products/air mixing.

Requires `ffmpeg` for the mp4 (falls back to an animated `.gif`).

## Configuration

| | |
|---|---|
| mesh | 120³ uniform, cube [−0.1, 0.1] m (dx = 1.7 mm) |
| solver | 5-eq (`model_eqns=2`), HLLC (`riemann_solver=2`, required for reactive burn), WENO3, RK3 |
| explosive | TNT JWL products (`fluid_pp(1)%eos=2`), full-density charge |
| ambient | air, ideal gas (`fluid_pp(2)%eos=1`) |
| burn | `jwl_reactive=T`, `jwl_G=5e-14`, `jwl_b_exp=2` |
| time | dt = 3×10⁻⁸ s, 400 steps (12 µs), 20 saved frames |

The charge is initialized as **unreacted explosive at ambient pressure**; the
central hot spot (25 GPa, r = 6 mm) provides the ignition kernel. The detonation
is self-sustaining thereafter — no imposed front. All six domain boundaries are
non-reflecting so the vented blast leaves cleanly.

**Fast-rate note.** On a laptop-affordable mesh (dx = 1.7 mm) the *resolved* JWL++
reaction zone would be sub-grid, so a slow rate decouples the burn from the shock
and the detonation fizzles. The rate coefficient `jwl_G = 5e-14` is chosen so the
reaction completes within roughly a cell — the standard fast-reaction limit that
gives a **grid-independent Chapman–Jouguet detonation** rather than a resolved ZND
structure. That is the right trade for a 3D demo; resolved ZND lives in the 1D case.

## What you see (verified in this run)

The detonation front locks onto the **TNT CJ speed (~6.7 km/s)** and sweeps the
5 cm charge in ~7 µs, holding **13–15 GPa** at the front (below the ~20 GPa planar
CJ value — the expected reduction for a *diverging spherical* front), then breaks
out and vents the hot products into the air. λ rises from 0 (unreacted) to ~1
(burnt) behind the front, and the products finger into the air along a
Richtmyer–Meshkov-unstable contact — the same physics validated quantitatively in
1D (`examples/1D_jwl_znd_detonation`, analytic CJ/ZND) and 2D
(`examples/2D_jwl_reactive_air_mixing`).

> This example is a **qualitative, real-world-appearance** demonstration of the
> reactive JWL capability in 3D, not a calibrated benchmark. For quantitative
> validation see the 1D ZND, 2D Sedov–Taylor, and underwater-explosion (Cole)
> cases.
