# Boundary-condition patch geometry has to match the dimensionality

`s_apply_boundary_patches` (`src/pre_process/m_boundary_conditions.fpp`) dispatches by dimensionality:

```fortran
if (p > 0) then                                  ! 3D
    if      (patch_bc(i)%geometry == 2) call s_circle_bc(i, bc_type)
    else if (patch_bc(i)%geometry == 3) call s_rectangle_bc(i, bc_type)
else if (n > 0) then                             ! 2D
    if      (patch_bc(i)%geometry == 1) call s_line_segment_bc(i, bc_type)
```

There is no `else`. A geometry belonging to the other dimensionality falls straight through: **the patch is
never applied, nothing is printed, and the face silently keeps whatever `bc_[xyz]` gave it.**

That is quiet in the worst way. A nozzle cut into a no-slip wall with a 3D geometry in a 2D case simply stays
a solid wall — the run completes, writes output, and the jet has a velocity of exactly zero for all time with
no indication anything was ignored.

## Running it

```
./mfc.sh validate examples/2D_bc_patch_geometry/case.py                  # geometry 1, valid in 2D
GEOMETRY=3 ./mfc.sh validate examples/2D_bc_patch_geometry/case.py       # a 3D geometry in a 2D case
```

| | before | after |
| --- | --- | --- |
| `GEOMETRY=1` | passes | passes |
| `GEOMETRY=3` | **passes, then does nothing at run time** | `patch_bc(1)%geometry must be 1 (line segment) in 2D; geometry 3 is never applied` |

The 3D direction is symmetric: geometry 1 in a 3D case is equally ignored, and is now equally refused.

## A second trap in the same corner

The case carries a second initial-condition patch a few cells thick at the inlet, and it is load-bearing:
the Dirichlet buffer is filled by pre_process from the **initial condition at the boundary face**. A domain
initialised at rest stores rest in that buffer, and the "inflow" then delivers zero for all time — again with
no warning. This is not what the validator change addresses; it is noted here because the two failures look
identical from the outside, and knowing that saves working out which one is in play.

## Scope

Validator only; no source change and no golden files. All 182 example cases still validate.
