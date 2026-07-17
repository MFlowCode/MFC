# 3D spherical PETN ZND detonation, octant symmetry

A grid-convergence study of a genuinely 3D spherical detonation, using octant
symmetry (reflective `bc_*%beg = -2` at `x = y = z = 0`) so the full sphere is
recovered at 1/8 the mesh cost. MFC has no native spherical (1/r^2) source
term, so this uses the code's real 3D machinery rather than an analytic
1D radial approximation -- the same technique already validated in
`benchmarks/3D_jwl_spherical_tnt_free_air_validation` and
`examples/2D_sedov_taylor_blast`.

The explosive is PETN with Kuhl's own JWL fit (Kuhl, Bell, Beckner,
Khasainov, "Simulation of Turbulent Combustion Fields of Shock-Dispersed
Aluminum Using the AMR Code," UCRL-PROC-225822, 2006, Appendix; `rho0 = 1.0
g/cc`, Kuhl's pressed-charge reference density, not full-density PETN --
see `toolchain/mfc/jwl_products.py`). Unlike the TNT burst benchmark, the
detonation here is a genuine self-propagating reactive burn (`jwl_reactive`,
Souers 2000 JWL++) with the Garno et al. (2020) Eq. 17 reactant energy
offset, so it carries real ZND structure rather than an instantaneous or
kinematic release.

## The CJ state and the overdriven 1D calibration

The self-sustained Chapman-Jouguet state of this material, from the
Hugoniot/Rayleigh tangency (rate-independent, EOS only; the same
construction reproduces literature CJ for full-density PETN, TNT and Comp B
to a few percent, and is corroborated by empirical low-density PETN at
1 g/cc, D ~ 5.5 km/s):

    D_CJ ~ 5384 m/s
    P_CJ ~ 7.35 GPa
    P_vN ~ 13.6 GPa (single-fluid offset model; uncalibrated vs experiment)
    jwl_delta_e = -7.98e6 J/kg  (Garno Eq. 17 offset, closed-form -- not fit)

Kuhl reports only T_CJ ~ 4600 K for this PETN; P_CJ and D_CJ are not in his
papers, so the pair above is computed, not quoted.

The shipped 1D calibration run is **overdriven**, not at CJ: it propagates at
D ~ 6245 m/s with a ~16.8 GPa state behind a ~19.8 GPa spike. That pair is
exactly the overdriven lambda = 1 Hugoniot/Rayleigh strong branch at
D = 6245 m/s (16.7 GPa) -- the front is supported by the oversized 30 GPa
`rxn_val = 1` booster and the closed rear boundary acting as a piston, so a
planar wave with no rear rarefaction never relaxes to CJ. A `jwl_G`
invariance check does **not** prove this front is at CJ: a piston-supported
overdriven front is equally G-independent (its speed is set by the piston,
not the reaction rate). Re-calibrating to the true CJ (softer booster,
rarefiable rear boundary) is follow-up work.

A pressure-only hot spot fizzles for this material (P4/P6: sub-critical
initiation) rather than transitioning to a self-sustained detonation, so the
booster patch is seeded already-reacted (`rxn_val = 1`) for robust direct
initiation, in both the 1D calibration case and here.

## Running the resolution study

    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 39   # dx = 1.00 mm (40^3 octant cells)
    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 79   # dx = 0.50 mm (80^3 octant cells)
    ./mfc.sh run examples/3D_jwl_znd_spherical_octant/case.py -n 8 -- --grid 159  # dx = 0.25 mm (160^3 octant cells)

Each resolution should be run from its own subdirectory copy of `case.py`
(`grid_39/`, `grid_79/`, `grid_159/`) so the probe/field output does not
collide -- MFC writes output relative to the case file's own directory.
Five radial probes (`probe_wrt`) sit at r = 0.3, 0.6, 1.0, 1.5, 2.0 times the
charge radius (4.5, 9.0, 15.0, 22.5, 30.0 mm), spanning the charge interior,
edge, and breakout region, so the front velocity, the lambda 0->1 rise behind
the front, and breakout can each be read off the probe histories.

    ./analyze_convergence.py grid_39 grid_79 grid_159

fits a front-arrival-time speed through the probes and reports the observed
order of grid convergence (P8: 3 grids x factor 2, `p =
log2((f_c-f_m)/(f_m-f_f))`) against the measured `D_CJ`.

## Status: propagation confirmed, front-speed diagnostic is NOT yet trustworthy

Both `grid_39` (dx = 1.00 mm, ~30 s) and `grid_79` (dx = 0.50 mm, ~150 s) run
clean (no NaN/abort) and show a compression wave sweeping outward from the
booster and breaking out into the surrounding air -- i.e. this is not simply
a decaying blast, the qualitative first check in the validation checklist
passes.

The naive `analyze_convergence.py` front-arrival fit (pressure crossing
2x ambient) does NOT yet give a trustworthy speed: it reports ~5.0 km/s at
dx = 1.00 mm and ~2.3 km/s at dx = 0.50 mm -- moving the WRONG way under
refinement (a genuine converging quantity should approach a fixed value
monotonically, not roughly halve). Two compounding diagnostic weaknesses
are the likely cause, not necessarily a physics failure:

1. **MFC's ASCII probe output prints only 6 decimal digits of time**
   (microsecond resolution) -- coarse relative to the ~1-2 us transit time
   probed here, so the "arrival time" at each probe is quantized to values
   like exactly 0, 1e-6, 2e-6 s regardless of the true sub-microsecond
   crossing time.
2. **A low arrival threshold (2x ambient) mixes the true shock front with a
   smeared numerical precursor** whose width is itself resolution-dependent
   (more numerical dissipation at coarse dx spreads the precursor further
   ahead, which can make the coarse run look FASTER by this metric even if
   the true front is not).

Both a 2D calibration check (short cylindrical PETN runs, same jwl_G, at
this case's resolution) and this octant case show pressure relaxing from the
overdriven booster state through several GPa before charge-edge breakout,
consistent with normal ZND relaxation toward CJ plus a real breakout
pressure drop -- not obviously a quench -- but this has not been confirmed
against a $\lambda$-based (reaction-completion) front tracker, which would be
immune to both weaknesses above. **Do not yet cite the front-speed numbers
from this case as a validated grid-convergence result.**

Before trusting a GCI number here: (1) track the front via the cell where
$\lambda$ crosses 0.5 (rides the HLLC contact flux, not smeared like a raw
pressure threshold) instead of a pressure threshold, (2) increase probe
output time precision or reduce `t_save`/probe cadence so arrival times are
resolved well below the ~1 us transit-time increments seen here, and (3) add
the dx = 0.25 mm (`grid_159`) run once (1)-(2) are in place -- a 2-point
"trend" is not sufficient to claim convergence order regardless.
