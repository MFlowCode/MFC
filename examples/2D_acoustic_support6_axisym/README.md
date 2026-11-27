# 2D axisymmetric focused transducer

A spherical-cap transducer (`acoustic(1)%support = 6`, `cyl_coord = T`) in water: radius of curvature F = 0.54 mm,
aperture 0.8 mm, apex on the axis at x = 0.06 mm, one Gaussian pulse ($\sigma_t$ = 20 ns, delay 100 ns, `mag` = 1 Pa).
The converging wave focuses on the axis at x = 0.6 mm, so the case tests the axisymmetric source, the geometric source
terms and the axis together.

On the axis the exact linear field of a cap with uniform normal velocity (O'Neil 1949) is

$$ p(z, t) = \frac{F}{F - z}\left[p_0(t - z/c) - p_0(t - R_e/c)\right], \quad R_e = \sqrt{a^2 + (z - h)^2}, \quad h = F - \sqrt{F^2 - a^2}, $$

which at the focus is $(h/c)\,dp_0/dt$: two peaks of $\pm\,\mathrm{mag}\,h/(c\sigma_t\sqrt{e})$ = ±3.60 Pa here.
`check_focus.py` reads the pressure of a probe in the cell just past the focus (written every step with `--probe`) and
prints the peak values and times against this solution and the L2 error of the signal.

## Run

    ./mfc.sh run examples/2D_acoustic_support6_axisym/case.py -t pre_process simulation -n 1 -- --refine 4 --probe
    python3 examples/2D_acoustic_support6_axisym/check_focus.py examples/2D_acoustic_support6_axisym --case-args "--refine 4 --probe"

`--refine k` multiplies the 100 x 50 cells by k and divides dt by k. The peaks converge to the O'Neil value (−0.4 % at
k = 8); see @ref validation_coordinates.
