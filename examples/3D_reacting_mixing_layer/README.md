# 3D Temporal Reacting Mixing Layer (H2/N2 - air, Mc = 1.5)

A temporally-evolving supersonic reacting shear layer between a hot air stream and an
N2-diluted hydrogen stream. The base state comes from a 1-D flamelet solve (Cantera +
Pyrometheus + JAX) extruded into 3D by `hcid=371`. This is the supersonic counterpart to
`examples/2D_reacting_mixing_layer`, which runs the same flamelet machinery at `Mc = 0.3`.

## Configuration

| Parameter | Value |
|---|---|
| Oxidizer stream | air, `X_O2 = 0.21`, 500 K |
| Fuel stream | `X_H2 = 0.5`, balance N2, 300 K |
| Pressure | 101325 Pa |
| Vorticity thickness `delta_omega` | 1.0e-3 m |
| Convective Mach number `Mc` | 1.5 |
| Domain | `[15, 20, 10] delta_omega` in `(x, y, z)` |
| Grid | 14 pts/`delta_omega` in x and z, 28 in y, so 210 x 560 x 140 |
| Time step | 1e-9 s, RK3 |
| Numerics | WENO5 (mapped, monotonicity-preserving), HLLC |
| Boundaries | periodic in x and z, ghost-cell extrapolation in y |
| Chemistry | San Diego mechanism, 9 species, unity-Lewis transport |

x is streamwise, y is cross-stream (the flamelet profile axis), z is spanwise. The
resolution density follows the temporal mixing-layer DNS of Wang et al. (*Combustion and
Flame*, 2024), with the box reduced to a quarter of theirs in each direction.

The fuel stream is diluted because pure H2's sound speed is over 4x the oxidizer's, so
reaching `Mc = 1.5` would demand a velocity split of roughly 2650 m/s. `X_H2 = 0.5` brings
that to about 1394 m/s and moves the stoichiometric mixture fraction off the domain edge.

## Initial condition

`case.py` calls `flamelet_ic.py` to solve a 1-D flamelet on the cross-stream grid and write
it as `prim.<n>.00.000000.dat` under `IC/`. `hcid=370` would extrude those files uniformly
across z and leave the flow z-invariant, so this case uses `hcid=371`: it scales the file's
cross-stream velocity by `1 + 0.5*cos(k_z z)` and sets the spanwise component from the
result, giving the IC 3D content at step 0. `k_z = 2*pi/L_z` is taken over the global
domain, so the IC does not depend on the MPI decomposition.

The in-plane `(x, y)` perturbation is baked into the files by `perturb_xy`, seeded by
`perturb_seed` in `case.py`. `IC/` is gitignored and regenerated on a fresh checkout, so
the fixed seed is what keeps the case reproducible. Regeneration is skipped when `IC/`
already matches the grid and the physical parameters, tracked in `.cache_key.json`.

The file spacing must match the run grid. A mismatch aborts in `pre_process`, so delete
`IC/` and let it regenerate after changing the grid or `--scale`.

## Running

    ./mfc.sh run examples/3D_reacting_mixing_layer/case.py -n 8

`--scale` shrinks the grid for cheap runs; `--scale 0.05` gives 32^3, which is what the
`3D -> Chemistry -> Reacting Mixing Layer` regression test uses. `--hot` runs the full
flamelet Newton/BDF solve instead of the default cold mollified profile.

The mechanism ships alongside the case as `sandiego.yaml` (UC San Diego Combustion
Research Group, <https://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html>).
