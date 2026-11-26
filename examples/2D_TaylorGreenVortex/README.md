# 2D Taylor-Green Vortex (viscous decay)

$u = U_0\sin(x/L_0)\cos(y/L_0)$, $v = -U_0\cos(x/L_0)\sin(y/L_0)$ in a periodic box $[-\pi, \pi]^2$ (`patch_icpp%geometry = 20`), at
$U_0 L_0/\nu = 1.2\times10^{-5}$ and Mach $3\times10^{-4}$. The exact solution decays as $e^{-2\nu t/L_0^2}$, so the kinetic energy
decays at the rate $4\nu/L_0^2$. `-N` sets the cells per direction (default 200).

`check_tgv.py` fits the decay rate of the kinetic energy over all saves and prints its error, and the error of the final kinetic
energy, against the exact solution:

    ./mfc.sh run examples/2D_TaylorGreenVortex/case.py -t pre_process simulation -n 1 -- -N 64
    python3 examples/2D_TaylorGreenVortex/check_tgv.py <run N=32> <run N=64> [--tol 1e-4]

The rate error is 7e-7 at 32 cells and 1e-10 at 128. Results and limits: @ref validation_vortex_order.
