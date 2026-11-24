# Premixed H2/Air Flame and a Vortex Ring (3D)

The 3D companion of `2D_flame_vortex_pair`: the same flame, mechanism, boundaries and fixed-dt numerics in a 3.2 mm cube, hit by a vortex ring (radius 0.4 mm, Gaussian core 0.15 mm, circulation 0.022 m^2/s, axis along x, centre at x = 0.9 mm). The ring's velocity is the exact axisymmetric Stokes stream function of its core. `pre_process` writes a uniform state and `ic3d.py` overwrites the restart file with the flame and the ring, so the run takes three commands and the case is not in the test suite. The restart file needs `parallel_io` (an MPI build).

    ./mfc.sh run examples/3D_flame_vortex_ring/case.py -n 1 -t pre_process -- --N 64
    python3 examples/3D_flame_vortex_ring/ic3d.py examples/3D_flame_vortex_ring
    ./mfc.sh run examples/3D_flame_vortex_ring/case.py -n 1 -t simulation post_process -- --N 64
    python3 examples/3D_flame_vortex_ring/analyze_fv3d.py RUN64 RUN128                      # copies of this directory, one per grid

`analyze_fv3d.py` prints the H2O mass, the flame area (T = 1200 K isosurface) and the heat release at 40, 70 and 100 us, with the observed order and Richardson value for three grids (it needs scikit-image). `ic3d.py --gamma 0` with `--nyz 32` gives the flat-flame control.

Run on an A100 at 32^3 to 20 us on `master` at 0cdfdf02a to check the three-command procedure (H2O mass 4.73474e-10 kg, flame area 1.02521e-5 m^2, heat release 59.59 W at 20 us). Convergence of the pre-merge study (64^3, 128^3, 256^3; not re-measured here) is in the page below: H2O mass orders 1.8 to 2.6, flame area order 2.0 at 70 us.

Results and limits: @ref validation_flame_vortex.
