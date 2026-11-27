# Sedov-Taylor blast with AMR

A point blast at the center of \f$[-0.6, 0.6]^3\f$: energy \f$E = 1\f$ in a sphere of radius `--r0` = 0.04, gas at rest with
\f$\rho_0 = 1\f$, \f$p_0 = 10^{-3}\f$, \f$\gamma = 1.4\f$, extrapolation boundaries (the shock stays inside). The whole sphere is
simulated: an AMR block must lie `buff_size` cells inside the domain, so an octant with symmetry walls through the blast
cannot refine the shock where it meets them. The self-similar solution has the shock at \f$R = \xi_0 (E t^2/\rho_0)^{1/5}\f$ with
\f$\xi_0 = 1.0328\f$, and a post-shock density of \f$(\gamma + 1)/(\gamma - 1) = 6\,\rho_0\f$.

- `--mode amr`: one 2:1 level regridded onto the shock every 4 steps (absolute density-jump tags), lock-step or `--subcycle`.
- `--mode uniform`: `--N` cells per direction, or 2N with `--fine`.

MUSCL-2 with the minmod limiter, HLLC, SSP-RK3, CFL-adaptive \f$\Delta t\f$ (0.3), saves every `--t-save` = 0.02 to `--t-end` =
0.16 (\f$R = 0.50\f$). WENO5 with `mp_weno = T` went NaN next to the blast center at t = 0.004 on a uniform 256^3 octant (no AMR).

`check_sedov.py` reads every save (`parallel_io = T`; level 0 and the fine blocks), takes the blast energy from the data
(the total energy above the ambient, which is conserved), averages the density over spherical shells one
finest cell thick, and takes the shock radius at the outermost half-jump between \f$\rho_0\f$ and the peak. It integrates the
similarity equations for \f$\xi_0\f$ and prints the radius error and the peak density ratio. Save times come from
`run_time.inf` (`run_time_info = T`).

```
./mfc.sh run examples/3D_amr_sedov/case.py -t pre_process simulation -n 1 -- --N 96 --mode amr --subcycle
python3 examples/3D_amr_sedov/check_sedov.py examples/3D_amr_sedov
```

The checkers import `mfc.viz`, so run them in MFC's Python environment (`source build/venv/bin/activate` after any `./mfc.sh` command).

Results and limits: the AMR section of the validation page.
