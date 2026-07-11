# 3D JWL Spherical TNT Free-Air Benchmark

This benchmark is a 3D JWL/TNT free-air blast candidate for MFC. A spherical
region of TNT detonation products expands into ambient air, and native point
probes record pressure-time histories for arrival time, peak incident
overpressure, and positive-phase impulse extraction.

Reference values are intentionally not fabricated. The MFC probe extraction is
implemented and the case runs, but the published Giam et al. reference values
below still need to be filled in before comparison against that specific paper
should be described as completed validation.

A separate, already-working far-field check exists alongside this candidate
table: `validate_kingery.py` compares the same gauge output against the
**Kinney & Graham (1985)** closed-form free-air overpressure correlation
(despite the script's filename, it does not implement Kingery-Bulmash/CONWEP —
see the script's own docstring). At last run this passed at 20% tolerance in
the far field (Z >= 0.30): 15.32% and 8.74% error at the two farthest gauges;
near-field gauges (Z < 0.30) are excluded, as expected for the coarse-mesh
burst-IC setup this benchmark currently uses (see "Physics Scope" below).

## Citation

Giam, Toh, and Tan, "Numerical Review of Jones-Wilkins-Lee Parameters for
Trinitrotoluene Explosive in Free-Air Blast," Journal of Applied Mechanics,
2020. DOI: `10.1115/1.4046243`.

The case is motivated by that free-air TNT/JWL review, but this repository does
not currently include an accessible table of the paper's exact benchmark
values. Until those values, Kingery-Bulmash/CONWEP values, or another trusted
numeric free-air TNT reference are added, this is a benchmark candidate rather
than completed validation.

## Physics Scope

This benchmark initializes TNT detonation products directly and validates 3D JWL product-air blast expansion. It does not model detonation initiation, reaction-zone structure, afterburn, or structural coupling.

The benchmark exercises:

- 3D JWL product-air expansion.
- Radial free-air blast propagation in a Cartesian grid.
- Native MFC pressure probes and post-processing.
- A comparison framework for published or standard free-air blast quantities.

The default (`--mode burst`) case does not exercise detonation initiation,
reaction-zone physics, or afterburn -- the products region starts already
detonated. An optional `--mode prog_burn` (see "Initiation Mode" below) adds a
kinematic detonation front plus afterburn mixing on top of the same geometry,
still against the same free-air Kinney-Graham comparison. Neither mode
validates structural coupling, confined blast, tunnel blast, or ground
reflection (see "Roadmap: Ground/Surface Burst" below).

## Why This Is Genuinely 3D

The default case resolves an octant of a spherical products region in a 3D
Cartesian domain. The shock expands in x, y, and z, with active transverse
reconstruction and Riemann fluxes. This is not a 1D shock tube extruded through
the 3D solver path.

Symmetry planes at `x = 0`, `y = 0`, and `z = 0` recover the full spherical
solution while reducing the cell count by a factor of about eight. The probes
lie near the positive x radial direction, at the nearest cell centers to the
target radii. Keeping probes on cell centers avoids boundary and
MPI-decomposition interpolation artifacts on the coarse default grid.

## Geometry

| Quantity | Value |
|---|---:|
| Charge center | `(0, 0, 0) m` in the octant representation |
| Charge radius | `0.05 m` |
| TNT/product reference density | `1630 kg/m^3` |
| Full-sphere TNT mass | `0.853466 kg` |
| Domain | `x, y, z in [0, 0.5] m` |
| Final time | `2.0e-4 s` |

The simulated octant contains one eighth of the full sphere. The symmetry
boundaries make the resolved field equivalent to the full charge, so scaled
distance uses the full TNT mass:

```text
W = (4/3) pi (0.05 m)^3 (1630 kg/m^3) = 0.853466 kg
Z = r / W^(1/3)
```

## Grid

Default local grid:

| Quantity | Value |
|---|---:|
| `m = n = p` | `63` |
| Cells in octant | `64^3` |
| Cell size | `0.0078125 m` |
| Reconstruction | mapped WENO3 |
| Riemann solver | HLLC |
| Time stepper | RK3 |
| CFL target | `0.3` |

Lower-resolution smoke run:

```bash
./mfc.sh run benchmarks/3D_jwl_spherical_tnt_free_air_validation/case.py -n 4 -- --grid 31
python3 benchmarks/3D_jwl_spherical_tnt_free_air_validation/gauges.py --grid 31
```

The `64^3` octant default gives a visibly cleaner spherical shock than the
older `32^3` smoke grid, but it takes several minutes and writes larger output
files.

## Initiation Mode (`--mode`)

```bash
./mfc.sh run benchmarks/3D_jwl_spherical_tnt_free_air_validation/case.py -n 4 -- --mode prog_burn
python3 benchmarks/3D_jwl_spherical_tnt_free_air_validation/gauges.py
python3 benchmarks/3D_jwl_spherical_tnt_free_air_validation/validate_kingery.py
```

`--mode burst` (default) is the original baseline above: the charge starts
already at the JWL reference (`V=1`) detonation-products pressure. `--mode
prog_burn` starts the charge unreacted at ambient pressure and ignites it with
a kinematic `prog_burn` front at the TNT CJ velocity (6930 m/s) from the
symmetry corner, plus `jwl_afterburn` mixing-rate combustion as products
contact the ambient air -- the same source combination validated end-to-end in
`examples/2D_jwl_prog_burn_afterburn_mixing/case.py`. Both modes are checked
against the same `validate_kingery.py` far-field comparison; running both and
comparing tests whether a real detonation front (rather than an instantaneous
burst) changes the near-field accuracy, which the burst-IC baseline's excluded
near-field error band has been suspected to reflect.

This case uses `cfl_adap_dt`, so `pb_D_cj*dt <= pb_width`'s startup check
(only enforced for fixed-dt runs) does not apply here; `pb_width` is sized to
3 grid cells, matching the ratio already validated in the 2D case, to keep the
front from outrunning its deposition band under the adaptive `dt` in practice.

## Boundary Conditions

| Boundary | MFC value | Meaning |
|---|---:|---|
| `bc_x%beg`, `bc_y%beg`, `bc_z%beg` | `-2` | symmetry planes |
| `bc_x%end`, `bc_y%end`, `bc_z%end` | `-3` | non-reflecting/open boundaries |

## EOS And Initial Conditions

Ambient air fills the full octant before the products sphere is overlaid:

| Parameter | Value |
|---|---:|
| Air EOS | ideal gas/stiffened gas, `eos = 1` |
| Air pressure | `101325 Pa` |
| Air density | `1.225 kg/m^3` |
| Air physical gamma | `1.4` |
| MFC air `fluid_pp(2)%gamma` | `2.5` |
| Air sound speed | `340.3 m/s` |
| Air `cv` | `717.5 J/(kg K)` |

TNT products use the JWL constants already present in the MFC JWL examples and
benchmarks. They are documented here as the repo-local parameter source, not as
a claim that these are the exact Giam et al. tabulated values:

| Parameter | Value |
|---|---:|
| Products EOS | JWL, `eos = 2` |
| Products density/reference density `rho0` | `1630 kg/m^3` |
| Products volume-fraction floor | `1.0e-6` |
| JWL `A` | `3.712e11 Pa` |
| JWL `B` | `3.231e9 Pa` |
| JWL `R1` | `4.15` |
| JWL `R2` | `0.95` |
| JWL `omega` | `0.30` |
| Initial specific energy equivalent | `E0/rho0 = 6.1908e6 J/kg` |
| Initial internal energy density `E0` | `1.0089e10 J/m^3` |
| JWL cold pressure at `rho = rho0` | `6.2837e9 Pa` |
| Initial products pressure | `9.3104e9 Pa` |
| Products `cv` | `613.5 J/(kg K)` |
| `jwl_air_e0` | `2.5575e5 J/kg` |
| `jwl_air_rho0` | `1.225 kg/m^3` |
| air Grüneisen (derived from `fluid_pp(2)%gamma`) | `1/2.5 = 0.4` |

The initialized products pressure is computed in `case.py` from:

```text
p = A (1 - omega/(R1 V)) exp(-R1 V)
  + B (1 - omega/(R2 V)) exp(-R2 V)
  + omega E0,
V = rho0 / rho = 1.
```

## Probe Locations

The default `--grid 63` run uses the nearest cell centers to target x positions
`0.15, 0.25, 0.35, 0.45 m`:

| Gauge | x (m) | y (m) | z (m) | r (m) | Z (m/kg^(1/3)) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.152344 | 0.003906 | 0.003906 | 0.152444 | 0.160712 |
| 2 | 0.253906 | 0.003906 | 0.003906 | 0.253966 | 0.267741 |
| 3 | 0.347656 | 0.003906 | 0.003906 | 0.347700 | 0.366558 |
| 4 | 0.449219 | 0.003906 | 0.003906 | 0.449253 | 0.473618 |

MFC's native probe interpolation can write a local initial pressure that differs
from the analytic ambient state on this coarse octant grid. The reducer therefore
uses each probe's first sample as its local baseline and reports pressure rise
relative to that baseline.

Arrival time is the first time the pressure rise exceeds 5 percent of ambient:

```text
p(t) - p(t=0) > 0.05 p_ambient
```

Peak incident overpressure is estimated from the same baseline-corrected trace:

```text
dp_peak = max(p(t) - p(t=0))
```

Positive-phase impulse is integrated over the baseline-corrected history:

```text
I = integral max(p(t) - p(t=0), 0) dt
```

## Run Commands

From the repository root:

```bash
./mfc.sh run benchmarks/3D_jwl_spherical_tnt_free_air_validation/case.py -n 4
python3 benchmarks/3D_jwl_spherical_tnt_free_air_validation/gauges.py
./mfc.sh precheck
```

The probe reducer reads:

```text
benchmarks/3D_jwl_spherical_tnt_free_air_validation/D/probe<i>_prim.dat
```

and writes:

```text
benchmarks/3D_jwl_spherical_tnt_free_air_validation/gauge_results.csv
```

## Comparison Table

| Gauge | r (m) | Z (m/kg^(1/3)) | MFC arrival (s) | Ref arrival (s) | MFC peak dp (Pa) | Ref peak dp (Pa) | MFC impulse (Pa s) | Ref impulse (Pa s) | Error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.152444 | 0.160712 | 8.000000e-06 | pending reference | 9.189414e+06 | pending reference | 1.142617e+03 | pending reference | pending |
| 2 | 0.253966 | 0.267741 | 3.700000e-05 | pending reference | 7.494431e+06 | pending reference | 3.127588e+02 | pending reference | pending |
| 3 | 0.347700 | 0.366558 | 7.200000e-05 | pending reference | 5.750937e+06 | pending reference | 2.453974e+02 | pending reference | pending |
| 4 | 0.449253 | 0.473618 | 1.150000e-04 | pending reference | 4.023872e+06 | pending reference | 3.029903e+02 | pending reference | pending |

Reference values may be filled from Giam et al. if numeric values are available,
from Kingery-Bulmash/CONWEP free-air TNT data if used for comparison, or from
another clearly cited spherical free-air TNT reference. If values are digitized
from a figure, label them as digitized estimates.

## Pass/Fail Expectations

This candidate benchmark passes its local sanity checks when:

- The case completes without NaNs.
- Density, pressure, and internal energy remain positive.
- Shock arrival time increases monotonically with radius.
- Peak incident overpressure decreases monotonically with radius.
- The pressure histories show a radially reasonable outward blast.

Quantitative error should be reported only after trusted reference arrival,
peak incident overpressure, and impulse values are added.

## Limitations And Remaining Work

- This is not completed validation until the reference columns are populated
  (the Giam et al. paper table above -- distinct from the already-working
  `validate_kingery.py` far-field check, see "Physics Scope").
- A `32^3` octant run is available as a faster smoke test, but the default is now `64^3`.
- `--mode burst` (default) does not model a detonation wave or reaction zone;
  `--mode prog_burn` adds a kinematic front and afterburn (see "Initiation
  Mode") but is still free-air, unconfined, and has no ground reflection.
- A higher-resolution run should be used before drawing quantitative
  conclusions.

## Roadmap: Ground/Surface Burst

Free-air (this benchmark) is a convenient validation geometry, but most
real-world full-scale TNT test data is closer to a near-surface or
hemispherical burst than an idealized free-air sphere. A follow-up benchmark
sketch, not yet built:

- **Geometry**: a hemispherical charge at the domain corner, `bc_x%beg = -2`
  as the physical ground-reflection plane, 2D `cyl_coord` axisymmetric with
  `bc_y%beg = -2` as the (unrelated, purely numerical) axisymmetric-axis
  condition `cyl_coord` requires -- the same BC code serving two different
  physical roles at the two edges that meet at the charge's corner.
- **Feasibility gate first**: no existing MFC example combines `cyl_coord`
  with `eos_jwl`. A small smoke case validating that combination in isolation
  should come before investing in a full validated benchmark, the same way
  `prog_burn` and `jwl_reactive` were each proven out small before this 3D
  free-air case existed.
- **Reference correlation**: comparing a surface burst against free-air
  Kinney-Graham would be methodologically wrong -- surface/hemispherical
  bursts reflect off the ground and roughly double the effective charge mass
  at the shock. This needs either a distinct surface-burst correlation or the
  standard ground-reflection charge-mass correction, which is itself an
  additional empirical approximation on top of the JWL EOS approximation --
  separate validation-methodology work, not just a geometry change.
