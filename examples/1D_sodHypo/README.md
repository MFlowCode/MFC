# 1D hypoelastic water shock tube

Water (stiffened gas, $\gamma$ = 4.4, $\pi_\infty$ = 6e8 Pa, 1000 kg/m³) with a shear modulus G = 10 GPa
(`hypoelasticity = T`), at rest, with a pressure jump at x = 0.5 m: 1e8 Pa on the left (`--pL`), 1e5 Pa on the right;
WENO5, HLL, SSP-RK3, 64 µs. The jump in the normal stress $\sigma = -p + \tau_{xx}$ splits into two longitudinal waves at
$c_L = \sqrt{(K + 4G/3)/\rho}$, $K = \gamma(p + \pi_\infty)$, about 4000 m/s here against a sound speed of 1625 m/s.

`check_elastic.py` compares the last save with the linear elastodynamic Riemann solution: the plateau velocity
$u^* = (\sigma_R - \sigma_L)/(Z_l + Z_r)$, $Z = \rho c_L$, the plateau pressure and $\tau_{xx}$ on each side of the contact
($dp = -(K/M)\,d\sigma$, $d\tau_{xx} = (4G/3)/M\,d\sigma$, $M = K + 4G/3$), and the front speeds (where u crosses $u^*/2$).
The same script reads `2D_hypo_shear_contact`, whose $\tau_{xy}$ jump launches shear waves at $c_S = \sqrt{G/\rho}$.
The prediction is first order in the jump, so a finite jump adds $O(\Delta\sigma/M)$ nonlinear errors: the default 1e8 Pa
jump is 0.6 % of M, `--pL 1e6` is 6e-5.

## Run

    ./mfc.sh run examples/1D_sodHypo/case.py -t pre_process simulation post_process -n 1 -- --pL 1e6
    python3 examples/1D_sodHypo/check_elastic.py examples/1D_sodHypo --case-args "--pL 1e6"

`case.py` also takes `-N` (cells, default 202; dt scales with 202/N). Results: @ref validation_elastic_waves.
