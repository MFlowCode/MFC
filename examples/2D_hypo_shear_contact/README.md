# 2D hypoelastic shear contact

A density contact (1 | 2 at x = 0.5) in a hypoelastic ideal gas ($\gamma$ = 1.4, G = 1, p = 1), periodic in y, with a
uniform tangential velocity $v_0 = 1$ and a jump of the shear stress $\tau_{xy}$ = 0 | A. The jump launches shear waves at
$c_S = \sqrt{G/\rho}$ = 1 and 0.707. First-order MUSCL, HLLD, forward Euler.

The case serves two checks:

- `Convergence -> HypoShearContact -> HLLD / HLLC` in `./mfc.sh test` sweeps A over one step and checks the order of the
  pressure response (see the docstring of `case.py`);
- `examples/1D_sodHypo/check_elastic.py` compares a longer run with the linear shear Riemann solution
  $v^* - v_0 = \Delta\tau_{xy}/(Z_l + Z_r)$, $\tau_{xy}^* = Z_l (v^* - v_0)$, $Z = \sqrt{\rho G}$, and the front speeds.

## Run

    ./mfc.sh run examples/2D_hypo_shear_contact/case.py -t pre_process simulation post_process -n 1 -- -N 128 --steps 106
    python3 examples/1D_sodHypo/check_elastic.py examples/2D_hypo_shear_contact --case-args "-N 128 --steps 106"

`--steps` = ceil(0.8265 N) runs to t = 0.25. Results: @ref validation_elastic_waves.
