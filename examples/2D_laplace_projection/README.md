# Static Drop (2D, all-Mach pressure projection, surface tension)

A quarter water drop (R = 0.15, sigma = 8) at rest in air, as in
`examples/2D_laplace_pressure_jump` but with the five-equation model. The exact state is
rest with a pressure jump sigma/R = 53.33; any velocity is a spurious current.

```shell
./mfc.sh run examples/2D_laplace_projection/case.py                          # well-balanced
./mfc.sh run examples/2D_laplace_projection/case.py -- --model conservative
python3 examples/2D_laplace_projection/analyze.py examples/2D_laplace_projection
```

`surface_tension_model` selects how surface tension enters the projection:

- `conservative` (1): the divergence of MFC's capillary stress tensor, added to the cell momentum as
  in the explicit solver. It reaches the cells while the pressure that should cancel it acts on faces,
  so the two do not balance.
- `well_balanced` (2): a Brackbill CSF force, sigma*kappa_f*(alpha_nb - alpha)/dx, on faces next to
  the pressure gradient, with the same face density as the pressure operator. A constant curvature is
  then balanced exactly by a pressure jump, and the spurious currents come only from curvature error.
  The indicator is the volume fraction (the force is unchanged under alpha -> 1 - alpha), which stays
  sharp; the color function the stress-tensor model uses is upwinded and smears.

Maximum spurious speed (the pressure jump is within ~1% of sigma/R = 53.33 in every case):

|                                   | N = 50, t <= 2 | N = 100, t <= 1 |
|-----------------------------------|----------------|-----------------|
| projection, `well_balanced`       | 0.08-0.10 m/s  | 0.20-0.34 m/s   |
| projection, `conservative`        | 0.78-1.1 m/s   | 1.0-2.2 m/s     |
| explicit HLLC (`--explicit`)      | 0.02 m/s       |                 |

The explicit run needs 100,000 steps to t = 0.5 at N = 50; the projection needs 500. The remaining
well-balanced currents are the curvature error of kappa = -div(grad alpha/|grad alpha|) on a 2-3 cell
interface, which grows as the interface thins in cells.
