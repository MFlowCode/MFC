# 1D Gaussian acoustic pulse from a planar source

A planar acoustic source (`acoustic(1)%support = 1`) at x = 0.1 mm in water (stiffened gas, c = 1492.8 m/s) emits one
Gaussian pulse (`pulse = 2`, $\sigma_t$ = 20 ns, delay 100 ns, `mag` = 1 Pa) in +x. The source adds
$S = \mathrm{mag}\,g(t) f(x)$ to the momentum, $S/c$ to the mass and $S c/(\gamma - 1)$ to the energy, so in
linear acoustics only the right-going invariant $p + \rho c u$ is driven and, away from the support,

$$ p'(x, t) = \mathrm{mag} \sum_j f_j \Delta x_j\, g\big(t - t_d - (x - x_j)/c\big), $$

with $f_j$ the discrete Gaussian weights of the support (standard deviation 1.25 $\Delta x_j$).

`check_acoustic.py` compares the pressure of a save with this expression (by default the last save with the whole signal
inside the domain) and prints the relative L2 and max errors outside the support, and the amplitude A and time shift
$\tau$ of a least-squares fit $p'_\mathrm{MFC} = A\,p'_\mathrm{lin}(t + \tau)$. The same script reads the
sibling examples `1D_acoustic_gauss_sigmadist`, `1D_acoustic_sine_frequency`, `1D_acoustic_sine_wavelength` and
`1D_acoustic_dipole` (`dipole = T`: two pulses of the same sign, one each way).

## Run

    ./mfc.sh run examples/1D_acoustic_gauss_sigmatime/case.py -t pre_process simulation post_process -n 1
    python3 examples/1D_acoustic_gauss_sigmatime/check_acoustic.py examples/1D_acoustic_gauss_sigmatime

`case.py` takes `-N` (cells; dt and the steps scale with it, c dt/dx = 0.6) and `--stretch a_x` (a stretched grid
clustered on 0.3-0.6 mm). For a sweep, run copies of the directory and pass the case arguments of each run:

    python3 check_acoustic.py run200 run400 run800 --case-args "; -N 400; -N 800"

## Result

The amplitude converges to the linear value and the shape after the shift converges at about third order, but the pulse lags
by $\Delta t/2$ at every resolution: the source is evaluated at the step-start time in all Runge-Kutta stages
([#2462](https://github.com/MFlowCode/MFC/issues/2462)). The full tables are on @ref validation_acoustics.
