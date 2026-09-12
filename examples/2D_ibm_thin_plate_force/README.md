# Immersed-boundary force on a thin plate (issue #1849)

A 2D flat plate pitching about its leading edge, 0 → 45° on an Eldredge smoothed ramp (a = 21), K = π/8
(case C1), Re_c = 300, Ma 0.2, plate thickness 2.5 % of chord. There is a published measurement to compare
against:

> Jantzen, Taira, Granlund & Ol, *Phys. Fluids* **26**, 053606 (2014), Fig. 10, 2D panel, curve C1.

`C_L = F_y / (½ ρ U² c) = 2 F_y` here, with ρ = U = c = 1 and MFC's 2D force being per unit depth.

This case exists to **measure** the defect in #1849, not to fix it. Whoever does fix it needs a number to move.

## Running the sweep

`NCELL` sets how many cells lie across the plate thickness. The physical problem is identical at every level —
same chord, same thickness, same domain, same times — so the sweep isolates the resolution requirement from
any change of geometry, which a thickness sweep would confound.

```
NCELL=2  ./mfc.sh run examples/2D_ibm_thin_plate_force/case.py     # dx = 0.0125 c,   0.22 M cells
NCELL=4  ./mfc.sh run examples/2D_ibm_thin_plate_force/case.py     # dx = 0.00625 c,  0.90 M cells
NCELL=8  ./mfc.sh run examples/2D_ibm_thin_plate_force/case.py     # dx = 0.003125 c, 3.58 M cells
NCELL=16 ./mfc.sh run examples/2D_ibm_thin_plate_force/case.py     # dx = 0.0015625 c, 14.3 M cells
```

`SUMMARY=1 python3 case.py` prints the grid and step count without running anything.

## What the sweep shows

| cells across the thickness | peak C_L | C_L at t = 4 | rms difference from the reference |
| --- | --- | --- | --- |
| 2 | 6.413 | 1.209 | 20.9 % |
| 4 | 4.682 | 1.271 | 33.5 % |
| 8 | 4.458 | 1.111 | 27.9 % |
| 16 | 4.833 | 1.005 | 26.9 % |
| reference | 6.998 | 1.521 | — |

Two separate readings:

**The peak converges from four cells up.** 4.682, 4.458 and 4.833 across a fourfold refinement — a spread of
8 % with no trend. Two cells is genuinely under-resolved and 35 % out; four is enough. So a thin body does not
need ten or more cells before the immersed boundary resolves it, which is worth knowing on its own for cost
estimates.

**The disagreement with the reference does not shrink.** The last column sits at 27–34 % at every resolution
and shows no sign of falling as the grid refines. **That is the measurement that matters for #1849**: it rules
out under-resolution as the explanation and points at the force computation itself.

## Why a zero-reference check is not available here

The obvious cheap test — a symmetric body whose true force is exactly zero — does not apply to a pitching
plate, whose lift is large and unknown. For that style of check see
`examples/2D_ibm_force_decomposition` (a cylinder, whose lift must be zero) and
`examples/3D_ibm_neighborhood_radius`. This case trades the exact reference for a published one.

## Related

- #1849 — the force integral over body cells on a thin plate, which this measures
- #1859 — an out-of-bounds coefficient read in the same force path, fixed
- #1863 — a spurious transverse force on multi-rank runs, open
