import itertools
import math
import os
import typing

from mfc import common

from ..state import ARG
from .case import CaseGeneratorStack, Nt, TestCaseBuilder, define_case_d, define_case_f, define_convergence_case
from .convergence import ConvergenceSpec, run_amp_sweep, run_dt_sweep, run_h_sweep, run_sod_l1

# Convergence test specs.
# One TestCase per (problem, scheme) pair. Trace prefix "Convergence ->" is
# the filter handle (`./mfc.sh test --only Convergence`); convergence cases
# are skipped by default.

# Advection convergence cases. Cell-shift mode: T = K*h per resolution, compare
# q(T) to np.roll(q(0), +K) per dim. Cost is O(1) in N (Nt = K*c/CFL independent
# of resolution) — wins ~10-100x vs period mode (full advection period).
# WENO7/TENO7 stay in period mode: at typical N their cell-shift signal sinks
# below machine precision (h^8 < 1e-15 at N=64) before any rate develops.
#
# expected_order is always the scheme's spatial order p. The runner subtracts
# 1 from the displayed rate in cell-shift mode (where raw rate = p+1) so the
# reported "spatial order" matches expected_order in both modes.
_CONS_VARS_1D = [("density", 1), ("x-momentum", 2), ("energy", 3)]
_CONS_VARS_2D = [("density", 1), ("energy", 4)]
_CONS_VARS_3D = [("density", 1), ("energy", 5)]

# (label, extra_args, expected_order, tol, resolutions)
# WENO3-JS at smooth extrema empirically achieves ~1.5 in MFC (Henrick mapping
# enabled). MUSCL2 unlimited central → effective spatial order 2.
_CONVERGENCE_1D_SCHEMES = [
    ("WENO5", ["--order", "5", "--cfl", "0.02"], 5, 0.3, [32, 64, 128]),
    ("WENO3", ["--order", "3", "--cfl", "0.02"], 1.5, 0.3, [64, 128, 256]),
    ("WENO1", ["--order", "1", "--cfl", "0.02"], 1, 0.2, [64, 128, 256]),
    ("MUSCL2", ["--muscl", "--muscl-lim", "0", "--cfl", "0.02"], 2, 0.3, [64, 128, 256]),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6", "--cfl", "0.02"], 5, 0.3, [32, 64, 128]),
]
# WENO7/TENO7 in 1D: period mode (full period T=1.0, see 1D case.py).
_CONVERGENCE_1D_PERIOD_SCHEMES = [
    ("WENO7", ["--order", "7", "--cfl", "0.005"], 7, 0.5, [64, 128]),
    ("TENO7", ["--order", "7", "--teno", "--teno-ct", "1e-9", "--cfl", "0.005"], 7, 0.5, [64, 128]),
]

_CONVERGENCE_2D_SCHEMES = [
    ("WENO5", ["--order", "5", "--cfl", "0.02"], 5, 0.3, [32, 64, 96]),
    ("WENO3", ["--order", "3", "--cfl", "0.02"], 1.5, 0.3, [32, 64, 128]),
    ("WENO1", ["--order", "1", "--cfl", "0.02"], 1, 0.2, [32, 64, 128]),
    ("MUSCL2", ["--muscl", "--muscl-lim", "0", "--cfl", "0.02"], 2, 0.3, [32, 64, 128]),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6", "--cfl", "0.02"], 5, 0.3, [32, 64, 96]),
]
_CONVERGENCE_2D_PERIOD_SCHEMES = [
    ("WENO7", ["--order", "7", "--cfl", "0.005"], 7, 0.5, [80, 96]),
    ("TENO7", ["--order", "7", "--teno", "--teno-ct", "1e-9", "--cfl", "0.005"], 7, 0.5, [80, 96]),
]

# 3D diagonal advection: only cell-shift mode (period T=1/3 with N^3 cells
# would dominate CI even at N=64). WENO7/TENO7 skipped — at N=64 with K=1
# the spatial error signal is below machine precision.
_CONVERGENCE_3D_SCHEMES = [
    ("WENO5", ["--order", "5", "--cfl", "0.02"], 5, 0.3, [32, 64]),
    ("WENO3", ["--order", "3", "--cfl", "0.02"], 1.5, 0.3, [32, 64]),
    ("WENO1", ["--order", "1", "--cfl", "0.02"], 1, 0.2, [32, 64]),
    ("MUSCL2", ["--muscl", "--muscl-lim", "0", "--cfl", "0.02"], 2, 0.3, [32, 64]),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6", "--cfl", "0.02"], 5, 0.3, [32, 64]),
]

# Sod L1 self-convergence: any conservative monotone scheme converges at L1
# rate ~1 (Godunov). SUPERBEE is over-compressive; min_N=128 skips its
# pre-asymptotic point.
_RES_SOD_DEFAULT = [128, 256, 512, 1024]
_CONVERGENCE_SOD_SCHEMES = [
    ("WENO1", ["--order", "1"], 1, 0.5, None),
    ("WENO3", ["--order", "3"], 1, 0.3, None),
    ("WENO5", ["--order", "5"], 1, 0.3, None),
    ("WENO7", ["--order", "7"], 1, 0.3, None),
    ("MUSCL-minmod", ["--muscl", "--muscl-lim", "1"], 1, 0.3, None),
    ("MUSCL-MC", ["--muscl", "--muscl-lim", "2"], 1, 0.3, None),
    ("MUSCL-VanLeer", ["--muscl", "--muscl-lim", "4"], 1, 0.3, None),
    ("MUSCL-SUPERBEE", ["--muscl", "--muscl-lim", "5"], 1, 0.5, 128),
    ("TENO5", ["--order", "5", "--teno", "--teno-ct", "1e-6"], 1, 0.3, None),
    ("TENO7", ["--order", "7", "--teno", "--teno-ct", "1e-9"], 1, 0.3, None),
]

# Temporal order: fixed N=512 / WENO5; vary CFL.
_CONVERGENCE_TEMPORAL_SCHEMES = [
    ("RK1", ["--order", "5", "--time-stepper", "1"], 1, 0.1, [0.10, 0.05]),
    ("RK2", ["--order", "5", "--time-stepper", "2"], 2, 0.2, [0.50, 0.25]),
    ("RK3", ["--order", "5", "--time-stepper", "3"], 3, 0.3, [0.50, 0.25]),
]


def add_convergence_cases(cases):
    num_ranks = 4

    def _h_sweep(case_path, ndim, cons_vars, extra_args, expected, tol, resolutions, cell_shift):
        return ConvergenceSpec(
            runner=run_h_sweep,
            case_path=case_path,
            extra_args=extra_args,
            expected_order=expected,
            tol=tol,
            cons_vars=cons_vars,
            resolutions=resolutions,
            ndim=ndim,
            cell_shift=cell_shift,
            num_ranks=num_ranks,  # ignored by run_h_sweep when cell_shift > 0
        )

    advection_groups = [
        (_CONVERGENCE_1D_SCHEMES, "1D", "examples/1D_euler_convergence/case.py", 1, _CONS_VARS_1D, 1, 1),
        (_CONVERGENCE_1D_PERIOD_SCHEMES, "1D", "examples/1D_euler_convergence/case.py", 1, _CONS_VARS_1D, 0, num_ranks),
        (_CONVERGENCE_2D_SCHEMES, "2D", "examples/2D_advection_convergence/case.py", 2, _CONS_VARS_2D, 1, 1),
        (_CONVERGENCE_2D_PERIOD_SCHEMES, "2D", "examples/2D_advection_convergence/case.py", 2, _CONS_VARS_2D, 0, num_ranks),
        (_CONVERGENCE_3D_SCHEMES, "3D", "examples/3D_advection_convergence/case.py", 3, _CONS_VARS_3D, 1, 1),
    ]
    for schemes, dim_label, case_path, ndim, cons_vars, cell_shift, ppn in advection_groups:
        for label, extra_args, expected, tol, resolutions in schemes:
            cases.append(
                define_convergence_case(
                    f"Convergence -> {dim_label} -> {label}",
                    spec=_h_sweep(case_path, ndim, cons_vars, extra_args, expected, tol, resolutions, cell_shift),
                    ppn=ppn,
                )
            )

    for label, extra_args, expected, tol, min_N in _CONVERGENCE_SOD_SCHEMES:
        resolutions = [N for N in _RES_SOD_DEFAULT if min_N is None or N >= min_N]
        cases.append(
            define_convergence_case(
                f"Convergence -> Sod -> {label}",
                spec=ConvergenceSpec(
                    runner=run_sod_l1,
                    case_path="examples/1D_sod_convergence/case.py",
                    extra_args=extra_args,
                    expected_order=expected,
                    tol=tol,
                    resolutions=resolutions,
                    num_ranks=num_ranks,
                ),
                ppn=num_ranks,
            )
        )

    for label, extra_args, expected, tol, cfls in _CONVERGENCE_TEMPORAL_SCHEMES:
        cases.append(
            define_convergence_case(
                f"Convergence -> Temporal -> {label}",
                spec=ConvergenceSpec(
                    runner=run_dt_sweep,
                    case_path="examples/1D_euler_convergence/case.py",
                    extra_args=extra_args,
                    expected_order=expected,
                    tol=tol,
                    cons_vars=_CONS_VARS_1D,
                    cfls=cfls,
                    fixed_N=512,
                    num_ranks=num_ranks,
                ),
                ppn=num_ranks,
            )
        )

    # Hypoelastic shear-contact amplitude order (see examples/2D_hypo_shear_contact):
    # HLLD's paired tangential momentum/energy star fluxes give a quadratic pressure
    # response to a tau_xy jump; HLLC's mismatched weights give linear. The
    # HLLC leg is the control that the case still discriminates (with v0 = 0 both
    # solvers are quadratic, so a drifted case would silently pass HLLD alone).
    for label, solver, expected, tol in [("HLLD", "hlld", 2.0, 0.1), ("HLLC", "hllc", 1.0, 0.2)]:
        cases.append(
            define_convergence_case(
                f"Convergence -> HypoShearContact -> {label}",
                spec=ConvergenceSpec(
                    runner=run_amp_sweep,
                    case_path="examples/2D_hypo_shear_contact/case.py",
                    extra_args=["--solver", solver],
                    expected_order=expected,
                    tol=tol,
                    amps=[1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3, 1.6e-2],
                ),
                ppn=1,
            )
        )


def make_3d_box_patches(
    z_centroids=(0.05, 0.45, 0.9),
    z_lengths=(0.1, 0.7, 0.2),
    geometry=9,
) -> dict:
    """3-patch 3D box IC: uniform xy plane (centroid=0.5, length=1), z spacing given."""
    d = {}
    for pid in range(1, 4):
        d[f"patch_icpp({pid})%geometry"] = geometry
        for vel in (1, 2, 3):
            d[f"patch_icpp({pid})%vel({vel})"] = 0.0
        d[f"patch_icpp({pid})%x_centroid"] = 0.5
        d[f"patch_icpp({pid})%length_x"] = 1
        d[f"patch_icpp({pid})%y_centroid"] = 0.5
        d[f"patch_icpp({pid})%length_y"] = 1
        d[f"patch_icpp({pid})%z_centroid"] = z_centroids[pid - 1]
        d[f"patch_icpp({pid})%length_z"] = z_lengths[pid - 1]
    return d


def get_bc_mods(bc: int, dimInfo):
    params = {}
    for dimCmp in dimInfo[0]:
        params.update({f"bc_{dimCmp}%beg": bc, f"bc_{dimCmp}%end": bc})

    return params


def get_dimensions():
    r = []

    for dimInfo in [(["x"], {"m": 299, "n": 0, "p": 0}, {"geometry": 1}), (["x", "y"], {"m": 49, "n": 39, "p": 0}, {"geometry": 3}), (["x", "y", "z"], {"m": 24, "n": 24, "p": 24}, {"geometry": 9})]:
        dimParams = {**dimInfo[1]}

        for dimCmp in dimInfo[0]:
            dimParams.update({f"{dimCmp}_domain%beg": 0.0e00, f"{dimCmp}_domain%end": 1.0e00})

        dimParams.update(get_bc_mods(-3, dimInfo))

        for patchID in range(1, 3 + 1):
            dimParams[f"patch_icpp({patchID})%geometry"] = dimInfo[2].get("geometry")

            if "z" in dimInfo[0]:
                dimParams.update(
                    {
                        f"patch_icpp({1})%z_centroid": 0.05,
                        f"patch_icpp({1})%length_z": 0.1,
                        f"patch_icpp({2})%z_centroid": 0.45,
                        f"patch_icpp({2})%length_z": 0.7,
                        f"patch_icpp({3})%z_centroid": 0.9,
                        f"patch_icpp({3})%length_z": 0.2,
                        f"patch_icpp({patchID})%y_centroid": 0.5,
                        f"patch_icpp({patchID})%length_y": 1,
                        f"patch_icpp({patchID})%x_centroid": 0.5,
                        f"patch_icpp({patchID})%length_x": 1,
                    }
                )

            elif "y" in dimInfo[0]:
                dimParams.update(
                    {
                        f"patch_icpp({1})%y_centroid": 0.05,
                        f"patch_icpp({1})%length_y": 0.1,
                        f"patch_icpp({2})%y_centroid": 0.45,
                        f"patch_icpp({2})%length_y": 0.7,
                        f"patch_icpp({3})%y_centroid": 0.9,
                        f"patch_icpp({3})%length_y": 0.2,
                        f"patch_icpp({patchID})%x_centroid": 0.5,
                        f"patch_icpp({patchID})%length_x": 1,
                    }
                )
            else:
                dimParams.update(
                    {
                        f"patch_icpp({1})%x_centroid": 0.05,
                        f"patch_icpp({1})%length_x": 0.1,
                        f"patch_icpp({2})%x_centroid": 0.45,
                        f"patch_icpp({2})%length_x": 0.7,
                        f"patch_icpp({3})%x_centroid": 0.9,
                        f"patch_icpp({3})%length_x": 0.2,
                    }
                )

            if "x" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(1)"] = 0.0

            if "y" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(2)"] = 0.0

            if "z" in dimInfo[0]:
                dimParams[f"patch_icpp({patchID})%vel(3)"] = 0.0

        r.append((dimInfo, dimParams))

    return r


# Always-run "canary" smoke set: one cheap, feature-dominant regression case per major
# physics module. Tagged canary=True in list_cases() so coverage-based selection
# (toolchain/mfc/test/coverage.py:select_tests) can never skip them on any lane -- a silent
# regression that disables a feature then trips on every PR. Validated in list_cases(), so a
# renamed/removed trace fails loudly instead of silently un-tagging the canary.
_CANARY_TRACES = frozenset(
    {
        "1D -> 1 Fluid(s) -> Viscous",  # m_viscous (Newtonian, Re=1e-4)
        "1D -> 1 Fluid(s) -> Non-Newtonian -> nn=0.5",  # m_hb_function (Herschel-Bulkley)
        "2D -> 2 Fluid(s) -> capillary=T -> model_eqns=3",  # m_surface_tension
        "1D -> Bubbles -> QBMM",  # m_qbmm / m_bubbles_EE
        "2D -> Lagrange Bubbles -> One-way Coupling",  # m_bubbles_EL
        "1D -> MHD -> HLLD",  # m_mhd / m_riemann_solver_hlld
        "1D -> Hypoelasticity -> 1 Fluid(s)",  # m_hypoelastic
        "1D -> Chemistry -> Perfect Reactor",  # chemistry
        "2D -> 1 Fluid(s) -> IBM -> Circle -> slip",  # m_ibm
        "1D -> Phase Change model 5 -> 2 Fluid(s) -> model equation -> 3",  # m_pressure_relaxation (6-eq)
        "1D -> Acoustic Source -> Sine -> Frequency",  # m_acoustic_src
        "1D -> Bodyforces",  # m_body_forces
    }
)


def list_cases() -> typing.List[TestCaseBuilder]:
    stack, cases = CaseGeneratorStack(), []

    def alter_bcs(dimInfo):
        for bc in [-1, -2, -4, -5, -6, -7, -8, -9, -10, -11, -12, -3, -15, -16, -17]:
            cases.append(define_case_d(stack, f"bc={bc}", get_bc_mods(bc, dimInfo)))

    def alter_grcbc(dimInfo):
        if len(dimInfo[0]) == 1:
            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(1)": 1.0,
                    "patch_icpp(2)%vel(1)": 1.0,
                    "patch_icpp(3)%vel(1)": 1.0,
                    "bc_x%beg": -7,
                    "bc_x%end": -8,
                    "bc_x%grcbc_in": "T",
                    "bc_x%grcbc_out": "T",
                    "bc_x%grcbc_vel_out": "T",
                    "bc_x%vel_in(1)": 1.0,
                    "bc_x%vel_in(2)": 0.0,
                    "bc_x%vel_in(3)": 0.0,
                    "bc_x%vel_out(1)": 1.0,
                    "bc_x%vel_out(2)": 0.0,
                    "bc_x%vel_out(3)": 0.0,
                    "bc_x%pres_in": 1.0,
                    "bc_x%pres_out": 1.0,
                    "bc_x%alpha_in(1)": 1.0,
                    "bc_x%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc x"], {}))
            stack.pop()
        elif len(dimInfo[0]) == 2:
            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(1)": 1.0,
                    "patch_icpp(2)%vel(1)": 1.0,
                    "patch_icpp(3)%vel(1)": 1.0,
                    "bc_x%beg": -7,
                    "bc_x%end": -8,
                    "bc_x%grcbc_in": "T",
                    "bc_x%grcbc_out": "T",
                    "bc_x%grcbc_vel_out": "T",
                    "bc_x%vel_in(1)": 1.0,
                    "bc_x%vel_in(2)": 0.0,
                    "bc_x%vel_in(3)": 0.0,
                    "bc_x%vel_out(1)": 1.0,
                    "bc_x%vel_out(2)": 0.0,
                    "bc_x%vel_out(3)": 0.0,
                    "bc_x%pres_in": 1.0,
                    "bc_x%pres_out": 1.0,
                    "bc_x%alpha_in(1)": 1.0,
                    "bc_x%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc x"], {}))
            stack.pop()

            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(2)": 1.0,
                    "patch_icpp(2)%vel(2)": 1.0,
                    "patch_icpp(3)%vel(2)": 1.0,
                    "bc_y%beg": -7,
                    "bc_y%end": -8,
                    "bc_y%grcbc_in": "T",
                    "bc_y%grcbc_out": "T",
                    "bc_y%grcbc_vel_out": "T",
                    "bc_y%vel_in(1)": 0.0,
                    "bc_y%vel_in(2)": 1.0,
                    "bc_y%vel_in(3)": 0.0,
                    "bc_y%vel_out(1)": 0.0,
                    "bc_y%vel_out(2)": 1.0,
                    "bc_y%vel_out(3)": 0.0,
                    "bc_y%pres_in": 1.0,
                    "bc_y%pres_out": 1.0,
                    "bc_y%alpha_in(1)": 1.0,
                    "bc_y%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc y"], {}))
            stack.pop()
        elif len(dimInfo[0]) == 3:
            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(1)": 1.0,
                    "patch_icpp(2)%vel(1)": 1.0,
                    "patch_icpp(3)%vel(1)": 1.0,
                    "bc_x%beg": -7,
                    "bc_x%end": -8,
                    "bc_x%grcbc_in": "T",
                    "bc_x%grcbc_out": "T",
                    "bc_x%grcbc_vel_out": "T",
                    "bc_x%vel_in(1)": 1.0,
                    "bc_x%vel_in(2)": 0.0,
                    "bc_x%vel_in(3)": 0.0,
                    "bc_x%vel_out(1)": 1.0,
                    "bc_x%vel_out(2)": 0.0,
                    "bc_x%vel_out(3)": 0.0,
                    "bc_x%pres_in": 1.0,
                    "bc_x%pres_out": 1.0,
                    "bc_x%alpha_in(1)": 1.0,
                    "bc_x%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc x"], {}))
            stack.pop()

            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(2)": 1.0,
                    "patch_icpp(2)%vel(2)": 1.0,
                    "patch_icpp(3)%vel(2)": 1.0,
                    "bc_y%beg": -7,
                    "bc_y%end": -8,
                    "bc_y%grcbc_in": "T",
                    "bc_y%grcbc_out": "T",
                    "bc_y%grcbc_vel_out": "T",
                    "bc_y%vel_in(1)": 0.0,
                    "bc_y%vel_in(2)": 1.0,
                    "bc_y%vel_in(3)": 0.0,
                    "bc_y%vel_out(1)": 0.0,
                    "bc_y%vel_out(2)": 1.0,
                    "bc_y%vel_out(3)": 0.0,
                    "bc_y%pres_in": 1.0,
                    "bc_y%pres_out": 1.0,
                    "bc_y%alpha_in(1)": 1.0,
                    "bc_y%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc y"], {}))
            stack.pop()

            stack.push(
                "",
                {
                    "patch_icpp(1)%vel(3)": 1.0,
                    "patch_icpp(2)%vel(3)": 1.0,
                    "patch_icpp(3)%vel(3)": 1.0,
                    "bc_z%beg": -7,
                    "bc_z%end": -8,
                    "bc_z%grcbc_in": "T",
                    "bc_z%grcbc_out": "T",
                    "bc_z%grcbc_vel_out": "T",
                    "bc_z%vel_in(1)": 0.0,
                    "bc_z%vel_in(2)": 0.0,
                    "bc_z%vel_in(3)": 1.0,
                    "bc_z%vel_out(1)": 0.0,
                    "bc_z%vel_out(2)": 0.0,
                    "bc_z%vel_out(3)": 1.0,
                    "bc_z%pres_in": 1.0,
                    "bc_z%pres_out": 1.0,
                    "bc_z%alpha_in(1)": 1.0,
                    "bc_z%alpha_rho_in(1)": 1.0,
                },
            )
            cases.append(define_case_d(stack, ["grcbc z"], {}))
            stack.pop()

    def alter_capillary():
        stack.push("", {"patch_icpp(1)%cf_val": 1, "patch_icpp(2)%cf_val": 0, "patch_icpp(3)%cf_val": 1, "sigma": 1, "model_eqns": 3, "surface_tension": "T"})
        cases.append(define_case_d(stack, ["capillary=T", "model_eqns=3"], {}))
        stack.pop()

    def alter_weno(dimInfo):
        for weno_order in [3, 5, 7]:
            stack.push(f"weno_order={weno_order}", {"weno_order": weno_order})
            for mapped_weno, wenoz, teno, mp_weno in itertools.product("FT", repeat=4):
                if sum(var == "T" for var in [mapped_weno, wenoz, teno, mp_weno]) > 1:
                    continue
                if mp_weno == "T" and weno_order != 5:
                    continue
                if teno == "T" and weno_order == 3:
                    continue

                trace = [f"{var}={val}" for var, val in zip(["mapped_weno", "wenoz", "teno", "mp_weno"], [mapped_weno, wenoz, teno, mp_weno]) if val == "T"]
                data = {var: "T" for var, val in zip(["mapped_weno", "wenoz", "teno", "mp_weno"], [mapped_weno, wenoz, teno, mp_weno]) if val == "T"}

                if "teno" in data:
                    data["teno_CT"] = 1e-6
                if "wenoz" in data and weno_order == 7:
                    data["wenoz_q"] = 3.0

                if weno_order == 7:
                    data = {**data, "weno_eps": 1e-6}  # increase damping for stability

                    if "z" in dimInfo[0]:
                        data = {**data, "m": 35, "n": 35, "p": 35}

                cases.append(define_case_d(stack, trace, data))

            stack.pop()

    def alter_igr():
        stack.push("IGR", {"igr": "T", "alf_factor": 10, "num_igr_iters": 10, "elliptic_smoothing": "T", "elliptic_smoothing_iters": 10, "num_igr_warm_start_iters": 10})

        for order in [3, 5]:
            stack.push(f"igr_order={order}", {"igr_order": order})

            cases.append(define_case_d(stack, "Jacobi", {"igr_iter_solver": 1}))
            if order == 5:
                cases.append(define_case_d(stack, "Gauss Seidel", {"igr_iter_solver": 2}))

            stack.pop()

        stack.pop()

    def alter_muscl():
        for muscl_order in [1, 2]:
            stack.push(f"muscl_order={muscl_order}", {"muscl_order": muscl_order, "recon_type": 2, "weno_order": 0, "weno_eps": None, "wenoz_q": None, "teno_CT": None})

            if muscl_order == 2:
                for muscl_lim in [2, 3, 4, 5]:
                    cases.append(define_case_d(stack, f"muscl_lim={muscl_lim}", {"muscl_lim": muscl_lim}))
                stack.push("muscl_eps=0", {"muscl_eps": 0})
                for muscl_lim in [1, 2, 3, 4, 5]:
                    cases.append(define_case_d(stack, f"muscl_lim={muscl_lim}", {"muscl_lim": muscl_lim}))
                stack.pop()
            stack.pop()

    def alter_eos():
        # BASE_CFG already sets fluid_pp(1)%pi_inf = 0, so the fluid is an ideal gas either way and
        # naming it one must not change a single digit. That is the point: it exercises the selector
        # end to end without moving the physics.
        cases.append(define_case_d(stack, "eos=ideal_gas", {"fluid_pp(1)%eos": "ideal_gas"}))

    def alter_riemann_solvers(num_fluids):
        for riemann_solver in [1, 5, 2]:
            stack.push(f"riemann_solver={riemann_solver}", {"riemann_solver": riemann_solver})

            cases.append(define_case_d(stack, "mixture_err", {"mixture_err": "T"}))

            if riemann_solver in (1, 2):
                cases.append(define_case_d(stack, "avg_state=1", {"avg_state": 1}))
                cases.append(define_case_d(stack, "wave_speeds=2", {"wave_speeds": 2}))

                # The averaged state is only read by the pressure-based wave speeds, so neither
                # case above reaches the Roe average: one computes it and discards it, the other
                # takes the arithmetic branch. Combining them is the only coverage it gets.
                if num_fluids == 1:
                    stack.push("avg_state=1", {"avg_state": 1})
                    cases.append(define_case_d(stack, "wave_speeds=2", {"wave_speeds": 2}))
                    stack.pop()

                if riemann_solver == 2:
                    cases.append(define_case_d(stack, "model_eqns=3", {"model_eqns": 3}))

                if num_fluids == 2:
                    if riemann_solver == 2:
                        cases.append(define_case_d(stack, "alt_soundspeed", {"alt_soundspeed": "T"}))

                    cases.append(define_case_d(stack, "mpp_lim", {"mpp_lim": "T"}))

            stack.pop()

    def add_hll_u_interface_cases(trace_prefix: str, u_interface_mods: typing.Optional[dict] = None):
        cases.append(
            define_case_d(
                stack,
                f"{trace_prefix} -> u-interface",
                {"riemann_solver": 1, "hll_u_interface": "T", **(u_interface_mods or {})},
            )
        )
        cases.append(define_case_d(stack, f"{trace_prefix} -> u-interface -> alt_soundspeed", {"riemann_solver": 1, "hll_u_interface": "T", "alt_soundspeed": "T", **(u_interface_mods or {})}))

    def alter_low_Mach_correction():
        stack.push("", {"fluid_pp(1)%gamma": 0.16, "fluid_pp(1)%eos": "stiffened_gas", "fluid_pp(1)%pi_inf": 3515.0, "dt": 1e-7})

        stack.push("riemann_solver=1", {"riemann_solver": 1})
        cases.append(define_case_d(stack, "low_Mach=1", {"low_Mach": 1}))
        stack.pop()
        stack.push("riemann_solver=2", {"riemann_solver": 2})
        cases.append(define_case_d(stack, "low_Mach=1", {"low_Mach": 1}))
        cases.append(define_case_d(stack, "low_Mach=2", {"low_Mach": 2}))
        stack.pop()

        stack.pop()

    def alter_int_comp(dimInfo):
        eps = 1e-6
        sharp_ic = {
            "patch_icpp(1)%alpha_rho(1)": 1.0 - eps,
            "patch_icpp(1)%alpha(1)": 1.0 - eps,
            "patch_icpp(1)%alpha_rho(2)": eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(2)%alpha_rho(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha_rho(2)": eps,
            "patch_icpp(2)%alpha(2)": eps,
            "patch_icpp(3)%alpha_rho(1)": eps,
            "patch_icpp(3)%alpha(1)": eps,
            "patch_icpp(3)%alpha_rho(2)": 1.0 - eps,
            "patch_icpp(3)%alpha(2)": 1.0 - eps,
        }

        stack.push("", sharp_ic)

        stack.push("weno_order=5", {"weno_order": 5})
        cases.append(define_case_d(stack, "int_comp=1", {"int_comp": 1}))
        if "y" in dimInfo[0]:  # Only test MTHINC in 2D and 3D
            cases.append(define_case_d(stack, "int_comp=2", {"int_comp": 2}))
            stack.push(
                "surface_tension=T",
                {
                    "surface_tension": "T",
                    "sigma": 1,
                    "patch_icpp(1)%cf_val": 1,
                    "patch_icpp(2)%cf_val": 0,
                    "patch_icpp(3)%cf_val": 1,
                },
            )
            cases.append(define_case_d(stack, "int_comp=1", {"int_comp": 1}))
            stack.pop()
        stack.pop()

        stack.push("muscl_order=2", {"muscl_order": 2, "recon_type": 2, "weno_order": 0, "weno_eps": None, "wenoz_q": None, "teno_CT": None})
        stack.push("int_comp=1", {"int_comp": 1})
        cases.append(define_case_d(stack, "muscl_lim=1", {"muscl_lim": 1}))
        stack.pop()
        if "y" in dimInfo[0]:  # Only test MTHINC in 2D and 3D
            stack.push("int_comp=2", {"int_comp": 2})
            cases.append(define_case_d(stack, "muscl_lim=1", {"muscl_lim": 1}))
            stack.pop()
        stack.pop()

        stack.pop()  # sharp IC

    def alter_num_fluids(dimInfo):
        for num_fluids in [1, 2]:
            stack.push(f"{num_fluids} Fluid(s)", {"num_fluids": num_fluids})

            if num_fluids == 2:
                stack.push(
                    "",
                    {
                        "fluid_pp(2)%gamma": 2.5,
                        "fluid_pp(2)%eos": "ideal_gas",
                        "fluid_pp(2)%pi_inf": 0.0,
                        "patch_icpp(1)%alpha_rho(1)": 0.81,
                        "patch_icpp(1)%alpha(1)": 0.9,
                        "patch_icpp(1)%alpha_rho(2)": 0.19,
                        "patch_icpp(1)%alpha(2)": 0.1,
                        "patch_icpp(2)%alpha_rho(1)": 0.25,
                        "patch_icpp(2)%alpha(1)": 0.5,
                        "patch_icpp(2)%alpha_rho(2)": 0.25,
                        "patch_icpp(2)%alpha(2)": 0.5,
                        "patch_icpp(3)%alpha_rho(1)": 0.08,
                        "patch_icpp(3)%alpha(1)": 0.2,
                        "patch_icpp(3)%alpha_rho(2)": 0.0225,
                        "patch_icpp(3)%alpha(2)": 0.8,
                    },
                )

                if len(dimInfo[0]) > 1:
                    alter_capillary()

            alter_riemann_solvers(num_fluids)
            if num_fluids == 2 and len(dimInfo[0]) > 1:
                # The existing 2D fluid-only HLL Method-2 row also guards the characteristic-CBC shared-velocity representation.
                cbc_mods = {"bc_y%end": -6} if len(dimInfo[0]) == 2 else None
                add_hll_u_interface_cases("riemann_solver=1", cbc_mods)
            alter_low_Mach_correction()
            if num_fluids == 1:
                alter_eos()
            alter_ib(dimInfo)
            if len(dimInfo[0]) > 1:
                alter_igr()

            if num_fluids == 2:
                alter_int_comp(dimInfo)

            if num_fluids == 1:
                stack.push("Viscous", {"fluid_pp(1)%Re(1)": 0.0001, "dt": 1e-11, "patch_icpp(1)%vel(1)": 1.0, "viscous": "T"})

                alter_ib(dimInfo, six_eqn_model=True, viscous=True)

                if len(dimInfo[0]) > 1:
                    alter_igr()

                cases.append(define_case_d(stack, "", {"weno_Re_flux": "F"}))
                cases.append(define_case_d(stack, "weno_Re_flux", {"weno_Re_flux": "T"}))
                cases.append(define_case_d(stack, "riemann_solver=5", {"riemann_solver": 5}))

                for weno_Re_flux in ["T"]:
                    stack.push("weno_Re_flux" if weno_Re_flux == "T" else "", {"weno_Re_flux": "T"})
                    cases.append(define_case_d(stack, "weno_avg", {"weno_avg": "T"}))
                    stack.pop()

                stack.pop()

                if len(dimInfo[0]) <= 2:
                    stack.push(
                        "Non-Newtonian",
                        {
                            "dt": 1e-11,
                            "patch_icpp(1)%vel(1)": 1.0,
                            "viscous": "T",
                            "riemann_solver": 2,
                            "model_eqns": 2,
                            "fluid_pp(1)%Re(1)": 1.0e4,
                            "fluid_pp(1)%non_newtonian": "T",
                            "fluid_pp(1)%tau0": 0.0,
                            "fluid_pp(1)%K": 1e-4,
                            "fluid_pp(1)%mu_max": 0.1,
                            "fluid_pp(1)%mu_min": 1e-6,
                            "fluid_pp(1)%hb_m": 1000.0,
                        },
                    )
                    cases.append(define_case_d(stack, "nn=0.5", {"fluid_pp(1)%nn": 0.5}))
                    cases.append(define_case_d(stack, "nn=1.5", {"fluid_pp(1)%nn": 1.5}))
                    cases.append(define_case_d(stack, "tau0=0.001", {"fluid_pp(1)%nn": 0.5, "fluid_pp(1)%tau0": 1.0e-3, "fluid_pp(1)%hb_m": 1.0e3}))
                    if len(dimInfo[0]) == 2:
                        # IBM + non-Newtonian: ib_state_wrt also exercises the
                        # per-stencil-sample HB viscosity in the IB force integration
                        cases.append(
                            define_case_d(
                                stack,
                                "IBM -> nn=0.5",
                                {
                                    "fluid_pp(1)%nn": 0.5,
                                    "ib": "T",
                                    "num_ibs": 1,
                                    "fd_order": 2,
                                    "ib_state_wrt": "T",
                                    "patch_ib(1)%geometry": 3,
                                    "patch_ib(1)%x_centroid": 0.5,
                                    "patch_ib(1)%y_centroid": 0.5,
                                    "patch_ib(1)%length_x": 0.05,
                                    "patch_ib(1)%length_y": 0.05,
                                    "patch_ib(1)%slip": "F",
                                },
                            )
                        )
                    stack.pop()

            if num_fluids == 2:
                stack.push(
                    "Viscous",
                    {"fluid_pp(1)%Re(1)": 0.001, "fluid_pp(1)%Re(2)": 0.001, "fluid_pp(2)%Re(1)": 0.001, "fluid_pp(2)%Re(2)": 0.001, "dt": 1e-11, "patch_icpp(1)%vel(1)": 1.0, "viscous": "T"},
                )

                alter_ib(dimInfo, six_eqn_model=True, viscous=True)

                if len(dimInfo[0]) > 1:
                    alter_igr()

                cases.append(define_case_d(stack, "", {"weno_Re_flux": "F"}))
                cases.append(define_case_d(stack, "weno_Re_flux", {"weno_Re_flux": "T"}))
                cases.append(define_case_d(stack, "riemann_solver=5", {"riemann_solver": 5}))
                for weno_Re_flux in ["T"]:
                    stack.push("weno_Re_flux" if weno_Re_flux == "T" else "", {"weno_Re_flux": "T"})
                    cases.append(define_case_d(stack, "weno_avg", {"weno_avg": "T"}))
                    stack.pop()

                stack.pop()

                if len(dimInfo[0]) == 2:
                    # Mixed non-Newtonian (fluid 1) / Newtonian (fluid 2) case
                    cases.append(
                        define_case_d(
                            stack,
                            "Non-Newtonian",
                            {
                                "dt": 1e-11,
                                "patch_icpp(1)%vel(1)": 1.0,
                                "viscous": "T",
                                "riemann_solver": 2,
                                "model_eqns": 2,
                                "fluid_pp(1)%Re(1)": 1.0e4,
                                "fluid_pp(1)%non_newtonian": "T",
                                "fluid_pp(1)%tau0": 0.0,
                                "fluid_pp(1)%K": 1e-4,
                                "fluid_pp(1)%nn": 0.5,
                                "fluid_pp(1)%mu_max": 0.1,
                                "fluid_pp(1)%mu_min": 1e-6,
                                "fluid_pp(1)%hb_m": 1000.0,
                                "fluid_pp(2)%Re(1)": 1.0e4,
                            },
                        )
                    )

                stack.pop()

            stack.pop()

    def alter_2d():
        stack.push(
            "Axisymmetric",
            {
                "num_fluids": 2,
                "bc_y%beg": -2,
                "cyl_coord": "T",
                "fluid_pp(2)%gamma": 2.5,
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                "patch_icpp(1)%alpha_rho(1)": 0.81,
                "patch_icpp(1)%alpha(1)": 0.9,
                "patch_icpp(1)%alpha_rho(2)": 0.19,
                "patch_icpp(1)%alpha(2)": 0.1,
                "patch_icpp(2)%alpha_rho(1)": 0.25,
                "patch_icpp(2)%alpha(1)": 0.5,
                "patch_icpp(2)%alpha_rho(2)": 0.25,
                "patch_icpp(2)%alpha(2)": 0.5,
                "patch_icpp(3)%alpha_rho(1)": 0.08,
                "patch_icpp(3)%alpha(1)": 0.2,
                "patch_icpp(3)%alpha_rho(2)": 0.0225,
                "patch_icpp(3)%alpha(2)": 0.8,
                "patch_icpp(1)%vel(1)": 0.0,
            },
        )

        cases.append(define_case_d(stack, "model_eqns=2", {"model_eqns": 2}))
        cases.append(define_case_d(stack, "model_eqns=3", {"model_eqns": 3}))
        cases.append(define_case_d(stack, "HLL", {"riemann_solver": 1}))
        add_hll_u_interface_cases("HLL")

        stack.push("Viscous", {"fluid_pp(1)%Re(1)": 0.0001, "fluid_pp(1)%Re(2)": 0.0001, "fluid_pp(2)%Re(1)": 0.0001, "fluid_pp(2)%Re(2)": 0.0001, "dt": 1e-11, "viscous": "T"})

        cases.append(define_case_d(stack, "", {"weno_Re_flux": "F"}))
        cases.append(define_case_d(stack, "weno_Re_flux", {"weno_Re_flux": "T"}))
        for weno_Re_flux in ["T"]:
            stack.push("weno_Re_flux" if weno_Re_flux == "T" else "", {"weno_Re_flux": "T"})
            cases.append(define_case_d(stack, "weno_avg", {"weno_avg": "T"}))
            stack.pop()

        stack.pop()
        stack.pop()

    def alter_3d():
        stack.push(
            "Cylindrical",
            {
                "bc_y%beg": -14,
                "bc_z%beg": -1,
                "bc_z%end": -1,
                "cyl_coord": "T",
                "x_domain%beg": 0.0e00,
                "x_domain%end": 5.0e00,
                "y_domain%beg": 0.0e00,
                "y_domain%end": 1.0e00,
                "z_domain%beg": 0.0e00,
                "z_domain%end": 2.0 * 3.141592653589793e00,
                "m": 29,
                "n": 29,
                "p": 29,
                "patch_icpp(1)%geometry": 10,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.0e00,
                "patch_icpp(1)%z_centroid": 0.0e00,
                "patch_icpp(1)%radius": 1.0,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": -1e6,
                "patch_icpp(1)%length_z": -1e6,
                "patch_icpp(2)%geometry": 10,
                "patch_icpp(2)%x_centroid": 2.5,
                "patch_icpp(2)%y_centroid": 0.0e00,
                "patch_icpp(2)%z_centroid": 0.0e00,
                "patch_icpp(2)%radius": 1.0,
                "patch_icpp(2)%length_x": 3.0,
                "patch_icpp(2)%length_y": -1e6,
                "patch_icpp(2)%length_z": -1e6,
                "patch_icpp(3)%geometry": 10,
                "patch_icpp(3)%x_centroid": 4.5,
                "patch_icpp(3)%y_centroid": 0.0e00,
                "patch_icpp(3)%z_centroid": 0.0e00,
                "patch_icpp(3)%radius": 1.0,
                "patch_icpp(3)%length_x": 1.0,
                "patch_icpp(3)%length_y": -1e6,
                "patch_icpp(3)%length_z": -1e6,
                "patch_icpp(1)%vel(1)": 0.0,
                "num_fluids": 2,
                "fluid_pp(2)%gamma": 2.5,
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                "patch_icpp(1)%alpha_rho(1)": 0.81,
                "patch_icpp(1)%alpha(1)": 0.9,
                "patch_icpp(1)%alpha_rho(2)": 0.19,
                "patch_icpp(1)%alpha(2)": 0.1,
                "patch_icpp(2)%alpha_rho(1)": 0.25,
                "patch_icpp(2)%alpha(1)": 0.5,
                "patch_icpp(2)%alpha_rho(2)": 0.25,
                "patch_icpp(2)%alpha(2)": 0.5,
                "patch_icpp(3)%alpha_rho(1)": 0.08,
                "patch_icpp(3)%alpha(1)": 0.2,
                "patch_icpp(3)%alpha_rho(2)": 0.0225,
                "patch_icpp(3)%alpha(2)": 0.8,
            },
        )

        cases.append(define_case_d(stack, "model_eqns=2", {"model_eqns": 2}))

        stack.push("cfl_adap_dt=T", {"cfl_adap_dt": "T", "cfl_target": 0.08, "t_save": 0.1, "n_start": 0, "t_stop": 0.1})
        cases.append(define_case_d(stack, "", {}))

        stack.pop()

        stack.push("Viscous", {"fluid_pp(1)%Re(1)": 0.0001, "fluid_pp(1)%Re(2)": 0.0001, "fluid_pp(2)%Re(1)": 0.0001, "fluid_pp(2)%Re(2)": 0.0001, "dt": 1e-10, "viscous": "T"})

        cases.append(define_case_d(stack, "", {"weno_Re_flux": "F"}))
        cases.append(define_case_d(stack, "weno_Re_flux", {"weno_Re_flux": "T"}))
        for weno_Re_flux in ["T"]:
            stack.push("weno_Re_flux" if weno_Re_flux == "T" else "", {"weno_Re_flux": "T"})
            cases.append(define_case_d(stack, "weno_avg", {"weno_avg": "T"}))
            stack.pop()

        stack.pop()
        stack.pop()

    def alter_ppn(dimInfo):
        if len(dimInfo[0]) == 3:
            cases.append(define_case_d(stack, "2 MPI Ranks", {"m": 29, "n": 29, "p": 49}, ppn=2))
            if ARG("rdma_mpi"):
                cases.append(define_case_d(stack, "2 MPI Ranks -> RDMA MPI", {"m": 29, "n": 29, "p": 49, "rdma_mpi": "T"}, ppn=2))
            cases.append(
                define_case_d(
                    stack,
                    "2 MPI Ranks -> IBM Sphere",
                    {
                        "m": 29,
                        "n": 29,
                        "p": 49,
                        "ib": "T",
                        "num_ibs": 1,
                        "fd_order": 2,
                        "patch_ib(1)%geometry": 8,
                        "patch_ib(1)%x_centroid": 0.5,
                        "patch_ib(1)%y_centroid": 0.5,
                        "patch_ib(1)%z_centroid": 0.5,
                        "patch_ib(1)%radius": 0.1,
                        "patch_icpp(1)%vel(1)": 0.001,
                        "patch_icpp(2)%vel(1)": 0.001,
                        "patch_icpp(3)%vel(1)": 0.001,
                        "patch_ib(1)%slip": "F",
                    },
                    ppn=2,
                )
            )
        else:
            cases.append(define_case_d(stack, "2 MPI Ranks", {}, ppn=2))
            if ARG("rdma_mpi"):
                cases.append(define_case_d(stack, "2 MPI Ranks -> RDMA MPI", {"rdma_mpi": "T"}, ppn=2))

    def alter_ib(dimInfo, six_eqn_model=False, viscous=False):
        for slip in [True, False]:
            stack.push(
                "IBM",
                {
                    "ib": "T",
                    "num_ibs": 1,
                    "fd_order": 2,
                    "patch_ib(1)%x_centroid": 0.5,
                    "patch_ib(1)%y_centroid": 0.5,
                    "patch_ib(1)%radius": 0.1,
                    "patch_icpp(1)%vel(1)": 0.001,
                    "patch_icpp(2)%vel(1)": 0.001,
                    "patch_icpp(3)%vel(1)": 0.001,
                    "patch_ib(1)%slip": "T" if slip else "F",
                },
            )

            suffix = " -> slip" if slip else ""

            if len(dimInfo[0]) == 3:
                cases.append(
                    define_case_d(
                        stack,
                        f"Sphere{suffix}",
                        {
                            "patch_ib(1)%z_centroid": 0.5,
                            "patch_ib(1)%geometry": 8,
                        },
                    )
                )

                cases.append(
                    define_case_d(
                        stack,
                        f"Cuboid{suffix}",
                        {
                            "patch_ib(1)%z_centroid": 0.5,
                            "patch_ib(1)%length_x": 0.1,
                            "patch_ib(1)%length_y": 0.1,
                            "patch_ib(1)%length_z": 0.1,
                            "patch_ib(1)%geometry": 9,
                        },
                    )
                )

                cases.append(
                    define_case_d(
                        stack,
                        f"Cylinder{suffix}",
                        {
                            "patch_ib(1)%z_centroid": 0.5,
                            "patch_ib(1)%length_x": 0.1,
                            "patch_ib(1)%geometry": 10,
                        },
                    )
                )

            elif len(dimInfo[0]) == 2:
                cases.append(
                    define_case_d(
                        stack,
                        f"Rectangle{suffix}",
                        {
                            "patch_ib(1)%length_x": 0.05,
                            "patch_ib(1)%length_y": 0.05,
                            "patch_ib(1)%geometry": 3,
                        },
                    )
                )
                cases.append(define_case_d(stack, f"Circle{suffix}", {"patch_ib(1)%geometry": 2, "n": 49}))
                if six_eqn_model:
                    cases.append(
                        define_case_d(
                            stack,
                            f"model_eqns=3{suffix}",
                            {
                                "patch_ib(1)%geometry": 2,
                                "model_eqns": 3,
                                "n": 49,  # there is a machine-level precision sensitivity to circles with n=39
                            },
                        )
                    )

            stack.pop()

        if len(dimInfo[0]) == 2 and not viscous:
            cases.append(
                define_case_d(
                    stack,
                    "IBM -> Periodic Circle",
                    {
                        "ib": "T",
                        "num_ibs": 1,
                        "fd_order": 2,
                        "bc_x%beg": -1,
                        "bc_x%end": -1,
                        "bc_y%beg": -1,
                        "bc_y%end": -1,
                        "patch_ib(1)%geometry": 2,
                        "patch_ib(1)%x_centroid": 0.0,
                        "patch_ib(1)%y_centroid": 0.0,
                        "patch_ib(1)%radius": 0.1,
                        "patch_icpp(1)%vel(1)": 0.001,
                        "patch_icpp(2)%vel(1)": 0.001,
                        "patch_icpp(3)%vel(1)": 0.001,
                        "patch_ib(1)%slip": "F",
                        "n": 49,
                    },
                )
            )
            cases.append(
                define_case_d(
                    stack,
                    "IBM -> Particle Cloud -> Hemisphere Shell",
                    {
                        "ib": "T",
                        "num_ibs": 0,
                        "num_particle_clouds": 1,
                        "fd_order": 2,
                        "ib_state_wrt": "T",
                        "n": 49,
                        "particle_cloud(1)%cloud_geometry": 2,
                        "particle_cloud(1)%packing_method": 1,
                        "particle_cloud(1)%x_centroid": 0.5,
                        "particle_cloud(1)%y_centroid": 0.1,
                        "particle_cloud(1)%num_particles": 4,
                        "particle_cloud(1)%radius": 0.02,
                        "particle_cloud(1)%mass": 1.0,
                        "particle_cloud(1)%min_spacing": 0.005,
                        "particle_cloud(1)%shell_inner_radius": 0.1,
                        "particle_cloud(1)%shell_outer_radius": 0.3,
                        "particle_cloud(1)%moving_ibm": 0,
                        "particle_cloud(1)%seed": 12345,
                        "patch_icpp(1)%vel(1)": 0.001,
                        "patch_icpp(2)%vel(1)": 0.001,
                        "patch_icpp(3)%vel(1)": 0.001,
                    },
                )
            )
            cases.append(
                define_case_d(
                    stack,
                    "IBM -> Particle Cloud -> Box",
                    {
                        "ib": "T",
                        "num_ibs": 0,
                        "num_particle_clouds": 1,
                        "fd_order": 2,
                        "ib_state_wrt": "T",
                        "n": 49,
                        "particle_cloud(1)%cloud_geometry": 1,
                        "particle_cloud(1)%packing_method": 1,
                        "particle_cloud(1)%x_centroid": 0.5,
                        "particle_cloud(1)%y_centroid": 0.5,
                        "particle_cloud(1)%length_x": 0.6,
                        "particle_cloud(1)%length_y": 0.6,
                        "particle_cloud(1)%num_particles": 50,
                        "particle_cloud(1)%radius": 0.02,
                        "particle_cloud(1)%mass": 1.0,
                        "particle_cloud(1)%min_spacing": 0.005,
                        "particle_cloud(1)%moving_ibm": 0,
                        "particle_cloud(1)%seed": 12345,
                        "patch_icpp(1)%vel(1)": 0.001,
                        "patch_icpp(2)%vel(1)": 0.001,
                        "patch_icpp(3)%vel(1)": 0.001,
                    },
                )
            )

        if len(dimInfo[0]) == 3 and not viscous:
            cases.append(
                define_case_d(
                    stack,
                    "IBM -> Particle Cloud -> Hemisphere Shell",
                    {
                        "ib": "T",
                        "num_ibs": 0,
                        "num_particle_clouds": 1,
                        "fd_order": 2,
                        "ib_state_wrt": "T",
                        "m": 29,
                        "n": 29,
                        "p": 29,
                        "particle_cloud(1)%cloud_geometry": 2,
                        "particle_cloud(1)%packing_method": 1,
                        "particle_cloud(1)%x_centroid": 0.5,
                        "particle_cloud(1)%y_centroid": 0.5,
                        "particle_cloud(1)%z_centroid": 0.1,
                        "particle_cloud(1)%num_particles": 4,
                        "particle_cloud(1)%radius": 0.06,
                        "particle_cloud(1)%mass": 1.0,
                        "particle_cloud(1)%min_spacing": 0.005,
                        "particle_cloud(1)%shell_inner_radius": 0.1,
                        "particle_cloud(1)%shell_outer_radius": 0.35,
                        "particle_cloud(1)%moving_ibm": 0,
                        "particle_cloud(1)%seed": 12345,
                        "patch_icpp(1)%vel(1)": 0.001,
                        "patch_icpp(2)%vel(1)": 0.001,
                        "patch_icpp(3)%vel(1)": 0.001,
                    },
                )
            )

    def ibm_stl():
        common_mods = {
            "t_step_stop": Nt,
            "t_step_save": Nt,
            "fd_order": 2,
            "num_stl_models": 1,
            "patch_ib(1)%model_id": 1,
            "stl_models(1)%model_scale(1)": 5.0,
            "stl_models(1)%model_scale(2)": 5.0,
            "stl_models(1)%model_scale(3)": 5.0,
            "stl_models(1)%model_threshold": 0.5,
        }

        for ndim in range(2, 4):
            cases.append(define_case_f(f"{ndim}D -> IBM -> STL", f"examples/{ndim}D_ibm_stl_test/case.py", ["--ndim", str(ndim)], mods=common_mods))

        # ICPP STL: the same flat-array winding-number model path as IBM, exercised as a constant-IC patch (geometry 21)
        cases.append(define_case_f("3D -> ICPP -> STL", "examples/3D_icpp_stl_cube/case.py", [], mods={"t_step_stop": Nt, "t_step_save": Nt}))
        cases.append(define_case_f("2D -> ICPP -> STL", "examples/2D_icpp_stl_circle/case.py", [], mods={"t_step_stop": Nt, "t_step_save": Nt}))

    ibm_stl()

    def alter_acoustic_src(dimInfo):
        stack.push("Acoustic Source", {"acoustic_source": "T", "acoustic(1)%support": 1, "dt": 1e-3, "t_step_stop": 50, "t_step_save": 50})

        transducer_params = {"acoustic(1)%loc(1)": 0.2, "acoustic(1)%foc_length": 0.4, "acoustic(1)%aperture": 0.6}

        if len(dimInfo[0]) == 1:
            for pulse_type in ["Sine", "Square"]:
                stack.push(pulse_type, {"acoustic(1)%pulse": 1 if pulse_type == "Sine" else 3})
                cases.append(define_case_d(stack, "Frequency", {"acoustic(1)%frequency": 50}))
                cases.append(define_case_d(stack, "Wavelength", {"acoustic(1)%wavelength": 0.02}))
                cases.append(define_case_d(stack, "Delay", {"acoustic(1)%delay": 0.02, "acoustic(1)%wavelength": 0.02}))
                cases.append(define_case_d(stack, "Number of Pulses", {"acoustic(1)%npulse": 2, "acoustic(1)%wavelength": 0.01}))
                stack.pop()

            stack.push("Gaussian", {"acoustic(1)%pulse": 2, "acoustic(1)%delay": 0.02})
            cases.append(define_case_d(stack, "Sigma Time", {"acoustic(1)%gauss_sigma_time": 0.01}))
            cases.append(define_case_d(stack, "Sigma Dist", {"acoustic(1)%gauss_sigma_dist": 0.01}))
            cases.append(define_case_d(stack, "Dipole", {"acoustic(1)%gauss_sigma_dist": 0.01, "acoustic(1)%dipole": "T"}))
            stack.pop()

        elif len(dimInfo[0]) == 2:
            stack.push("", {"acoustic(1)%loc(2)": 0.5, "acoustic(1)%wavelength": 0.02})

            stack.push("Planar", {})
            stack.push("support=2", {"acoustic(1)%support": 2})
            cases.append(define_case_d(stack, "", {}))
            cases.append(define_case_d(stack, "Dipole", {"acoustic(1)%dipole": "T"}))
            stack.pop()
            stack.pop()

            stack.push("Transducer", transducer_params)
            for support in [5, 6]:
                stack.push(f"support={support}", {"acoustic(1)%support": support, "cyl_coord": "T" if support == 6 else "F", "bc_y%beg": -2 if support == 6 else -3})
                cases.append(define_case_d(stack, "Sine", {}))
                cases.append(define_case_d(stack, "Gaussian", {"acoustic(1)%pulse": 2, "acoustic(1)%delay": 0.02, "acoustic(1)%gauss_sigma_dist": 0.01}))
                cases.append(define_case_d(stack, "Delay", {"acoustic(1)%delay": 0.02}))
                stack.pop()
            stack.pop()

            stack.push("Transducer Array", {**transducer_params, "acoustic(1)%num_elements": 4, "acoustic(1)%element_spacing_angle": 0.05, "acoustic(1)%element_on": 0})
            stack.push("support=9", {"acoustic(1)%support": 9})
            cases.append(define_case_d(stack, "All Elements", {}))
            cases.append(define_case_d(stack, "One element", {"acoustic(1)%element_on": 1}))
            stack.pop()
            cases.append(define_case_d(stack, "support=10", {"acoustic(1)%support": 10, "cyl_coord": "T", "bc_y%beg": -2}))
            stack.pop()

            stack.pop()

        elif len(dimInfo[0]) == 3:
            stack.push("", {"acoustic(1)%loc(2)": 0.5, "acoustic(1)%loc(3)": 0.5, "acoustic(1)%wavelength": 0.02})

            stack.push("Planar", {})
            stack.push("support=3", {"acoustic(1)%support": 3, "acoustic(1)%height": 0.25})
            cases.append(define_case_d(stack, "", {}))
            cases.append(define_case_d(stack, "Dipole", {"acoustic(1)%dipole": "T"}))
            stack.pop()
            stack.pop()

            stack.push("Transducer", transducer_params)
            cases.append(define_case_d(stack, "support=7", {"acoustic(1)%support": 7}))
            stack.pop()

            stack.push("Transducer Array", {**transducer_params, "acoustic(1)%num_elements": 6, "acoustic(1)%element_polygon_ratio": 0.7})
            stack.push("support=11", {"acoustic(1)%support": 11})
            cases.append(define_case_d(stack, "All Elements", {}))
            cases.append(define_case_d(stack, "One element", {"acoustic(1)%element_on": 1}))
            stack.pop()
            stack.pop()

            stack.pop()

        stack.pop()

    def alter_bubbles(dimInfo):
        if len(dimInfo[0]) > 0:
            stack.push("Bubbles", {"bubbles_euler": "T"})

            stack.push(
                "",
                {
                    "nb": 3,
                    "fluid_pp(1)%gamma": 0.16,
                    "fluid_pp(1)%eos": "stiffened_gas",
                    "fluid_pp(1)%pi_inf": 3515.0,
                    "bub_pp%R0ref": 1.0,
                    "bub_pp%p0ref": 1.0,
                    "bub_pp%rho0ref": 1.0,
                    "bub_pp%T0ref": 1.0,
                    "bub_pp%ss": 0.07179866765358993,
                    "bub_pp%pv": 0.02308216136195411,
                    "bub_pp%vd": 0.2404125083932959,
                    "bub_pp%mu_l": 0.009954269975623244,
                    "bub_pp%mu_v": 8.758168074360729e-05,
                    "bub_pp%mu_g": 0.00017881922111898042,
                    "bub_pp%gam_v": 1.33,
                    "bub_pp%gam_g": 1.4,
                    "bub_pp%M_v": 18.02,
                    "bub_pp%M_g": 28.97,
                    "bub_pp%k_v": 0.5583395141263873,
                    "bub_pp%k_g": 0.7346421281308791,
                    "bub_pp%R_v": 1334.8378710170155,
                    "bub_pp%R_g": 830.2995663005393,
                    "patch_icpp(1)%alpha_rho(1)": 0.96,
                    "patch_icpp(1)%alpha(1)": 4e-02,
                    "patch_icpp(2)%alpha_rho(1)": 0.96,
                    "patch_icpp(2)%alpha(1)": 4e-02,
                    "patch_icpp(3)%alpha_rho(1)": 0.96,
                    "patch_icpp(3)%alpha(1)": 4e-02,
                    "patch_icpp(1)%pres": 1.0,
                    "patch_icpp(2)%pres": 1.0,
                    "patch_icpp(3)%pres": 1.0,
                    "acoustic(1)%support": 1,
                    "acoustic(1)%wavelength": 0.25,
                },
            )

            stack.push("", {"acoustic_source": "T"})

            if len(dimInfo[0]) >= 2:
                stack.push("", {"acoustic(1)%loc(2)": 0.5, "acoustic(1)%support": 2})

            if len(dimInfo[0]) >= 3:
                stack.push("", {"acoustic(1)%support": 3, "acoustic(1)%height": 1e10})

            for polytropic in ["T", "F"]:
                stack.push("Polytropic" if polytropic == "T" else "", {"polytropic": polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {"bubble_model": bubble_model})

                    if not (polytropic == "F" and bubble_model == 3):
                        cases.append(define_case_d(stack, "", {}))

                    stack.pop()

                stack.pop()

            stack.push("", {"polytropic": "T", "bubble_model": 2})
            cases.append(define_case_d(stack, "nb=1", {"nb": 1}))

            stack.push("adv_n=T", {"adv_n": "T"})
            cases.append(define_case_d(stack, "", {}))
            cases.append(define_case_d(stack, "adap_dt=T", {"adap_dt": "T"}))
            stack.pop()

            stack.push("", {"fluid_pp(1)%eos": "stiffened_gas", "fluid_pp(1)%pi_inf": 351.5})
            cases.append(define_case_d(stack, "artificial_Ma", {"pi_fac": 0.1}))

            stack.pop()

            cases.append(define_case_d(stack, "low_Mach=1", {"low_Mach": 1}))
            cases.append(define_case_d(stack, "low_Mach=2", {"low_Mach": 2}))

            stack.push("QBMM", {"qbmm": "T"})
            cases.append(define_case_d(stack, "", {}))

            stack.push("Non-polytropic", {"polytropic": "F"})
            cases.append(define_case_d(stack, "", {}))

            stack.pop()

            stack.push("bubble_model=3", {"bubble_model": 3, "polytropic": "T"})
            cases.append(define_case_d(stack, "", {}))

            stack.push("Non-polytropic", {"polytropic": "F"})
            cases.append(define_case_d(stack, "", {}))

            for _ in range(7):
                stack.pop()

            if len(dimInfo[0]) >= 2:
                stack.pop()

            if len(dimInfo[0]) >= 3:
                stack.pop()

    def alter_hypoelasticity(dimInfo):
        # Hypoelasticity checks
        for num_fluids in [1, 2]:
            stack.push(
                f"Hypoelasticity -> {num_fluids} Fluid(s)",
                {
                    "hypoelasticity": "T",
                    "num_fluids": num_fluids,
                    "riemann_solver": 1,
                    "fd_order": 4,
                    "fluid_pp(1)%gamma": 0.3,
                    "fluid_pp(1)%eos": "stiffened_gas",
                    "fluid_pp(1)%pi_inf": 7.8e05,
                    "patch_icpp(1)%pres": 1.0e06,
                    "patch_icpp(1)%alpha_rho(1)": 1000.0e00,
                    "patch_icpp(2)%pres": 1.0e05,
                    "patch_icpp(2)%alpha_rho(1)": 1000.0e00,
                    "patch_icpp(3)%pres": 5.0e05,
                    "patch_icpp(3)%alpha_rho(1)": 1000.0e00,
                    "patch_icpp(1)%tau_e(1)": 0.0e-00,
                    "patch_icpp(2)%tau_e(1)": 0.0e-00,
                    "patch_icpp(3)%tau_e(1)": 0.0e-00,
                    "fluid_pp(1)%G": 1.0e05,
                },
            )

            if num_fluids == 2:
                stack.push(
                    "",
                    {
                        "fluid_pp(2)%gamma": 0.3,
                        "fluid_pp(2)%eos": "stiffened_gas",
                        "fluid_pp(2)%pi_inf": 7.8e05,
                        "patch_icpp(1)%alpha_rho(1)": 900.0e00,
                        "patch_icpp(1)%alpha(1)": 0.9,
                        "patch_icpp(1)%alpha_rho(2)": 100,
                        "patch_icpp(1)%alpha(2)": 0.1,
                        "patch_icpp(2)%alpha_rho(1)": 100,
                        "patch_icpp(2)%alpha(1)": 0.1,
                        "patch_icpp(2)%alpha_rho(2)": 900,
                        "patch_icpp(2)%alpha(2)": 0.9,
                        "patch_icpp(3)%alpha_rho(1)": 900,
                        "patch_icpp(3)%alpha(1)": 0.9,
                        "patch_icpp(3)%alpha_rho(2)": 100,
                        "patch_icpp(3)%alpha(2)": 0.1,
                        "fluid_pp(2)%G": 5.0e04,
                    },
                )

            if len(dimInfo[0]) >= 2:
                stack.push(
                    "",
                    {
                        "patch_icpp(1)%tau_e(2)": 0.0e00,
                        "patch_icpp(1)%tau_e(3)": 0.0e00,
                        "patch_icpp(2)%tau_e(2)": 0.0e00,
                        "patch_icpp(2)%tau_e(3)": 0.0e00,
                        "patch_icpp(3)%tau_e(2)": 0.0e00,
                        "patch_icpp(3)%tau_e(3)": 0.0e00,
                    },
                )

            if len(dimInfo[0]) == 3:
                stack.push(
                    "",
                    {
                        "patch_icpp(1)%tau_e(4)": 0.0e00,
                        "patch_icpp(1)%tau_e(5)": 0.0e00,
                        "patch_icpp(1)%tau_e(6)": 0.0e00,
                        "patch_icpp(2)%tau_e(4)": 0.0e00,
                        "patch_icpp(2)%tau_e(5)": 0.0e00,
                        "patch_icpp(2)%tau_e(6)": 0.0e00,
                        "patch_icpp(3)%tau_e(4)": 0.0e00,
                        "patch_icpp(3)%tau_e(5)": 0.0e00,
                        "patch_icpp(3)%tau_e(6)": 0.0e00,
                    },
                )

            cases.append(define_case_d(stack, "", {}))

            if len(dimInfo[0]) == 2 and num_fluids == 1:
                cases.append(
                    define_case_d(
                        stack,
                        "probe -> nonuniform stress",
                        {
                            "probe_wrt": "T",
                            "num_probes": 1,
                            "probe(1)%x": 0.5,
                            "probe(1)%y": 0.5,
                            "patch_icpp(1)%tau_e(1)": 1.0e04,
                            "patch_icpp(1)%tau_e(2)": 5.0e03,
                            "patch_icpp(1)%tau_e(3)": 2.0e03,
                            "patch_icpp(2)%tau_e(1)": -1.0e04,
                            "patch_icpp(2)%tau_e(2)": 3.0e03,
                            "patch_icpp(2)%tau_e(3)": -2.0e03,
                        },
                    )
                )

            reflective_params = {"bc_x%beg": -2, "bc_x%end": -2, "bc_y%beg": -2, "bc_y%end": -2}
            if len(dimInfo[0]) == 3:
                reflective_params.update({"bc_z%beg": -2, "bc_z%end": -2})

            if num_fluids == 1:
                cases.append(define_case_d(stack, "cont_damage", {"cont_damage": "T", "tau_star": 0.0, "cont_damage_s": 2.0, "alpha_bar": 1e-4}))
                if len(dimInfo[0]) == 2:
                    cases.append(
                        define_case_d(
                            stack,
                            "cont_damage -> HLLC -> nonuniform stress",
                            {
                                "riemann_solver": 2,
                                "cont_damage": "T",
                                "tau_star": 0.0,
                                "cont_damage_s": 2.0,
                                "alpha_bar": 1e-4,
                                "patch_icpp(1)%tau_e(1)": 100.0,
                                "patch_icpp(1)%tau_e(2)": 25.0,
                                "patch_icpp(1)%tau_e(3)": -100.0,
                                "patch_icpp(2)%tau_e(1)": 200.0,
                                "patch_icpp(2)%tau_e(2)": 50.0,
                                "patch_icpp(2)%tau_e(3)": -200.0,
                                "patch_icpp(3)%tau_e(1)": 300.0,
                                "patch_icpp(3)%tau_e(2)": 75.0,
                                "patch_icpp(3)%tau_e(3)": -300.0,
                            },
                        )
                    )
                if len(dimInfo[0]) >= 2:
                    cases.append(define_case_d(stack, "bc=-2", reflective_params))
                if len(dimInfo[0]) == 2:
                    cases.append(define_case_d(stack, "Axisymmetric", {**reflective_params, "cyl_coord": "T"}))

            stack.pop()

            if num_fluids == 2:
                stack.pop()

            if len(dimInfo[0]) == 2:
                stack.pop()

            if len(dimInfo[0]) == 3:
                for _ in range(2):
                    stack.pop()

    def alter_body_forces(dimInfo):
        ndims = len(dimInfo[0])

        stack.push("Bodyforces", {"bf_x": "T", "k_x": 1, "w_x": 1, "p_x": 1, "g_x": 10})

        if ndims >= 2:
            stack.push("", {"bf_y": "T", "k_y": 1, "w_y": 1, "p_y": 1, "g_y": 10})

        if ndims == 3:
            stack.push("", {"bf_z": "T", "k_z": 1, "w_z": 1, "p_z": 1, "g_z": 10})

        cases.append(define_case_d(stack, "", {}))

        stack.push("cfl_adap_dt=T", {"cfl_adap_dt": "T", "cfl_target": 0.08, "t_save": 0.025, "n_start": 0, "t_stop": 0.025})
        cases.append(define_case_d(stack, "", {}))

        stack.pop()

        stack.pop()

        if ndims >= 2:
            stack.pop()

        if ndims == 3:
            stack.pop()

        # Spatially supported body force (Wei & Freund, JFM 2005) in isolation: no
        # bf_x/y/z, no chemistry, so this golden pins the forcing kernel (momentum
        # sources + u*f energy work term) at the default 1e-12 tolerance. The only
        # other coverage (Spatial Reacting Mixing Layer) entangles it with stiff
        # chemistry at 1e-6. 2D-only by construction (m_checker rejects 1D/3D).
        if ndims == 2:
            stack.push(
                "SpatialBodyforces",
                {
                    "bf_spatial_support": "T",
                    "spatial_bf%amp": 1.0,
                    "spatial_bf%x_centroid": 0.5,
                    "spatial_bf%y_centroid": 0.5,
                    "spatial_bf%conv_vel": 1.0,
                    "spatial_bf%sigma": 100.0,
                    **{f"spatial_bf%freq({i})": 2.0 * i for i in range(1, 9)},
                    **{f"spatial_bf%phase({i})": 0.3 * i for i in range(1, 9)},
                },
            )
            cases.append(define_case_d(stack, "", {}))
            stack.pop()

    def alter_synthetic_turbulence(dimInfo):
        # 3-D solenoidal synthetic-turbulence forcing (m_body_forces): a quiescent,
        # triply-periodic single-fluid box driven by a deterministic (compiler-independent)
        # random-Fourier-mode volume force over two energy shells. Forced turbulence is
        # chaotic: tiny cross-compiler libm differences (cos/exp) amplify to O(field) within
        # ~20 steps (the reason the airfoil example is skipped). Run only 3 steps, where the
        # solution is still ~dt*F -- the deterministic forcing itself -- so the double-precision
        # golden holds to <1e-12 across CI compilers while fully exercising mode generation +
        # assembly + application. Single precision needs a looser tolerance (see
        # compute_tolerance): roundoff-driven divergence reaches ~1e-3 by t_step 3.
        if len(dimInfo[0]) == 3:
            cases.append(
                define_case_d(
                    stack,
                    "synthetic_turbulence",
                    {
                        "m": 24,
                        "n": 24,
                        "p": 24,
                        "dt": 1e-2,
                        "t_step_stop": 3,
                        "t_step_save": 3,
                        "num_fluids": 1,
                        "x_domain%beg": 0.0,
                        "x_domain%end": 1.0,
                        "y_domain%beg": 0.0,
                        "y_domain%end": 1.0,
                        "z_domain%beg": 0.0,
                        "z_domain%end": 1.0,
                        "bc_x%beg": -1,
                        "bc_x%end": -1,
                        "bc_y%beg": -1,
                        "bc_y%end": -1,
                        "bc_z%beg": -1,
                        "bc_z%end": -1,
                        # Keep the base's three box-tiling patches (all active, vel=0) but set
                        # them to a single uniform quiescent state, so the IC has no
                        # discontinuity -- the forcing alone drives the field.
                        "patch_icpp(1)%pres": 1.0,
                        "patch_icpp(1)%alpha_rho(1)": 1.0,
                        "patch_icpp(1)%alpha(1)": 1.0,
                        "patch_icpp(2)%pres": 1.0,
                        "patch_icpp(2)%alpha_rho(1)": 1.0,
                        "patch_icpp(2)%alpha(1)": 1.0,
                        "patch_icpp(3)%pres": 1.0,
                        "patch_icpp(3)%alpha_rho(1)": 1.0,
                        "patch_icpp(3)%alpha(1)": 1.0,
                        "synthetic_turbulence": "T",
                        "synth_seed": 1,
                        "synth_n_shells": 2,
                        "synth_U_inf": 1.0,
                        "num_turbulent_sources": 1,
                        "synth_k_shell(1)": 2 * math.pi / 0.5,
                        "synth_k_shell(2)": 2 * math.pi / 0.25,
                        "synth_amp_shell(1)": 0.02,
                        "synth_amp_shell(2)": 0.01,
                        "synth_n_waves_per_shell(1)": 2,
                        "synth_n_waves_per_shell(2)": 3,
                        "turb_pos(1,1)": 0.5,
                        "turb_pos(1,2)": 0.5,
                        "turb_pos(1,3)": 0.5,
                        "synth_L(1,1)": 4.0,
                        "synth_L(1,2)": 4.0,
                        "synth_L(1,3)": 4.0,
                    },
                )
            )

    def alter_mixlayer_perturb(dimInfo):
        if len(dimInfo[0]) == 3:
            cases.append(
                define_case_d(
                    stack,
                    "mixlayer_perturb",
                    {
                        "m": 24,
                        "n": 64,
                        "p": 24,
                        "dt": 1e-2,
                        "num_patches": 1,
                        "num_fluids": 1,
                        "x_domain%beg": 0.0,
                        "x_domain%end": 20.0,
                        "bc_x%beg": -1,
                        "bc_x%end": -1,
                        "y_domain%beg": -10.0,
                        "y_domain%end": 10.0,
                        "bc_y%beg": -6,
                        "bc_y%end": -6,
                        "z_domain%beg": 0.0,
                        "z_domain%end": 20.0,
                        "bc_z%beg": -1,
                        "bc_z%end": -1,
                        "mixlayer_vel_profile": "T",
                        "mixlayer_perturb": "T",
                        "weno_Re_flux": "F",
                        "weno_avg": "T",
                        "wenoz": "T",
                        "fluid_pp(1)%gamma": 2.5,
                        "fluid_pp(1)%eos": "ideal_gas",
                        "fluid_pp(1)%pi_inf": 0.0,
                        "fluid_pp(1)%Re(1)": 1.6881644098979287,
                        "viscous": "T",
                        "patch_icpp(1)%geometry": 9,
                        "patch_icpp(1)%x_centroid": 10.0,
                        "patch_icpp(1)%length_x": 20.0,
                        "patch_icpp(1)%y_centroid": 0.0,
                        "patch_icpp(1)%length_y": 20.0,
                        "patch_icpp(1)%z_centroid": 10.0,
                        "patch_icpp(1)%length_z": 20.0,
                        "patch_icpp(1)%vel(1)": 1.0,
                        "patch_icpp(1)%vel(2)": 0.0,
                        "patch_icpp(1)%vel(3)": 0.0,
                        "patch_icpp(1)%pres": 17.8571428571,
                        "patch_icpp(1)%alpha_rho(1)": 1.0,
                        "patch_icpp(1)%alpha(1)": 1.0,
                        "patch_icpp(1)%r0": -1e6,
                        "patch_icpp(1)%v0": -1e6,
                        "patch_icpp(2)%geometry": -100,
                        "patch_icpp(2)%x_centroid": -1e6,
                        "patch_icpp(2)%length_x": -1e6,
                        "patch_icpp(2)%y_centroid": -1e6,
                        "patch_icpp(2)%length_y": -1e6,
                        "patch_icpp(2)%z_centroid": -1e6,
                        "patch_icpp(2)%length_z": -1e6,
                        "patch_icpp(2)%vel(1)": -1e6,
                        "patch_icpp(2)%vel(2)": -1e6,
                        "patch_icpp(2)%vel(3)": -1e6,
                        "patch_icpp(2)%r0": -1e6,
                        "patch_icpp(2)%v0": -1e6,
                        "patch_icpp(3)%geometry": -100,
                        "patch_icpp(3)%x_centroid": -1e6,
                        "patch_icpp(3)%length_x": -1e6,
                        "patch_icpp(3)%y_centroid": -1e6,
                        "patch_icpp(3)%length_y": -1e6,
                        "patch_icpp(3)%z_centroid": -1e6,
                        "patch_icpp(3)%length_z": -1e6,
                        "patch_icpp(3)%vel(1)": -1e6,
                        "patch_icpp(3)%vel(2)": -1e6,
                        "patch_icpp(3)%vel(3)": -1e6,
                        "patch_icpp(3)%r0": -1e6,
                        "patch_icpp(3)%v0": -1e6,
                    },
                )
            )

    def alter_phasechange(dimInfo):
        ndims = len(dimInfo[0])

        # Phase Change checks
        for relax_model in [5] + ([6] if ndims <= 2 else []):
            for num_fluids in ([2] if ndims == 1 or relax_model == 5 else []) + [3]:
                for model_eqns in [3, 2]:
                    stack.push(
                        f"Phase Change model {relax_model} -> {num_fluids} Fluid(s) -> model equation -> {model_eqns}",
                        {
                            "relax": "T",
                            "relax_model": relax_model,
                            "model_eqns": model_eqns,
                            "palpha_eps": 1e-02,
                            "ptgalpha_eps": 1e-02,
                            "num_fluids": num_fluids,
                            "riemann_solver": 2,
                            "fluid_pp(1)%gamma": 0.7409,
                            "fluid_pp(1)%eos": "stiffened_gas",
                            "fluid_pp(1)%pi_inf": 1.7409e09,
                            "fluid_pp(1)%cv": 1816,
                            "fluid_pp(1)%qv": -1167000,
                            "fluid_pp(1)%qvp": 0.0,
                            "fluid_pp(2)%gamma": 2.3266,
                            "fluid_pp(2)%eos": "ideal_gas",
                            "fluid_pp(2)%pi_inf": 0.0e00,
                            "fluid_pp(2)%cv": 1040,
                            "fluid_pp(2)%qv": 2030000,
                            "fluid_pp(2)%qvp": -23400,
                            "patch_icpp(1)%pres": 4.3755e05,
                            "patch_icpp(1)%alpha(1)": 8.7149e-06,
                            "patch_icpp(1)%alpha_rho(1)": 9.6457e02 * 8.7149e-06,
                            "patch_icpp(1)%alpha(2)": 1 - 8.7149e-06,
                            "patch_icpp(1)%alpha_rho(2)": 2.3132 * (1 - 8.7149e-06),
                            "patch_icpp(2)%pres": 9.6602e04,
                            "patch_icpp(2)%alpha(1)": 3.6749e-05,
                            "patch_icpp(2)%alpha_rho(1)": 1.0957e03 * 3.6749e-05,
                            "patch_icpp(2)%alpha(2)": 1 - 3.6749e-05,
                            "patch_icpp(2)%alpha_rho(2)": 0.5803 * (1 - 3.6749e-05),
                            "patch_icpp(3)%pres": 9.6602e04,
                            "patch_icpp(3)%alpha(1)": 3.6749e-05,
                            "patch_icpp(3)%alpha_rho(1)": 1.0957e03 * 3.6749e-05,
                            "patch_icpp(3)%alpha(2)": 1 - 3.6749e-05,
                            "patch_icpp(3)%alpha_rho(2)": 0.5803 * (1 - 3.6749e-05),
                        },
                    )

                    if num_fluids == 3:
                        stack.push(
                            "",
                            {
                                "fluid_pp(3)%gamma": 2.4870,
                                "fluid_pp(3)%eos": "ideal_gas",
                                "fluid_pp(3)%pi_inf": 0.0e00,
                                "fluid_pp(3)%cv": 717.5,
                                "fluid_pp(3)%qv": 0.0e00,
                                "fluid_pp(3)%qvp": 0.0,
                                "patch_icpp(1)%alpha(2)": 2.5893e-02,
                                "patch_icpp(1)%alpha_rho(2)": 2.3132 * 2.5893e-02,
                                "patch_icpp(2)%alpha(2)": 2.8728e-02,
                                "patch_icpp(2)%alpha_rho(2)": 0.5803 * 2.8728e-02,
                                "patch_icpp(3)%alpha(2)": 2.8728e-02,
                                "patch_icpp(3)%alpha_rho(2)": 0.5803 * 2.8728e-02,
                                "patch_icpp(1)%alpha(3)": 1 - 8.7149e-06 - 2.5893e-02,
                                "patch_icpp(1)%alpha_rho(3)": 3.5840 * (1 - 8.7149e-06 - 2.5893e-02),
                                "patch_icpp(2)%alpha(3)": 1 - 3.6749e-05 - 2.8728e-02,
                                "patch_icpp(2)%alpha_rho(3)": 0.8991 * (1 - 3.6749e-05 - 2.8728e-02),
                                "patch_icpp(3)%alpha(3)": 1 - 3.6749e-05 - 2.8728e-02,
                                "patch_icpp(3)%alpha_rho(3)": 0.8991 * (1 - 3.6749e-05 - 2.8728e-02),
                            },
                        )

                    if ndims == 1:
                        stack.push("", {"patch_icpp(1)%vel(1)": 606.15, "patch_icpp(2)%vel(1)": 10.0, "patch_icpp(3)%vel(1)": 10.0})
                    elif ndims == 2:
                        stack.push(
                            "",
                            {
                                "patch_icpp(1)%vel(1)": 0.0,
                                "patch_icpp(2)%vel(1)": 0.0,
                                "patch_icpp(3)%vel(1)": 0.0,
                                "patch_icpp(1)%vel(2)": 606.15,
                                "patch_icpp(2)%vel(2)": 10.0,
                                "patch_icpp(3)%vel(2)": 10.0,
                            },
                        )
                    elif ndims == 3:
                        stack.push(
                            "",
                            {
                                "patch_icpp(1)%vel(1)": 0.0,
                                "patch_icpp(2)%vel(1)": 0.0,
                                "patch_icpp(3)%vel(1)": 0.0,
                                "patch_icpp(1)%vel(2)": 0.0,
                                "patch_icpp(2)%vel(2)": 0.0,
                                "patch_icpp(3)%vel(2)": 0.0,
                                "patch_icpp(1)%vel(3)": 606.15,
                                "patch_icpp(2)%vel(3)": 10.0,
                                "patch_icpp(3)%vel(3)": 10.0,
                            },
                        )

                    cases.append(define_case_d(stack, "", {}))

                    stack.pop()
                    stack.pop()

                    if num_fluids == 3:
                        stack.pop()

    def alter_viscosity(dimInfo):
        # Viscosity & bubbles checks
        if len(dimInfo[0]) > 0:
            stack.push("Viscosity -> Bubbles", {"fluid_pp(1)%Re(1)": 50, "bubbles_euler": "T", "viscous": "T"})

            stack.push(
                "",
                {
                    "nb": 1,
                    "fluid_pp(1)%gamma": 0.16,
                    "fluid_pp(1)%eos": "stiffened_gas",
                    "fluid_pp(1)%pi_inf": 3515.0,
                    "bub_pp%R0ref": 1.0,
                    "bub_pp%p0ref": 1.0,
                    "bub_pp%rho0ref": 1.0,
                    "bub_pp%T0ref": 1.0,
                    "bub_pp%ss": 0.07179866765358993,
                    "bub_pp%pv": 0.02308216136195411,
                    "bub_pp%vd": 0.2404125083932959,
                    "bub_pp%mu_l": 0.009954269975623244,
                    "bub_pp%mu_v": 8.758168074360729e-05,
                    "bub_pp%mu_g": 0.00017881922111898042,
                    "bub_pp%gam_v": 1.33,
                    "bub_pp%gam_g": 1.4,
                    "bub_pp%M_v": 18.02,
                    "bub_pp%M_g": 28.97,
                    "bub_pp%k_v": 0.5583395141263873,
                    "bub_pp%k_g": 0.7346421281308791,
                    "bub_pp%R_v": 1334.8378710170155,
                    "bub_pp%R_g": 830.2995663005393,
                    "patch_icpp(1)%alpha_rho(1)": 0.96,
                    "patch_icpp(1)%alpha(1)": 4e-02,
                    "patch_icpp(2)%alpha_rho(1)": 0.96,
                    "patch_icpp(2)%alpha(1)": 4e-02,
                    "patch_icpp(3)%alpha_rho(1)": 0.96,
                    "patch_icpp(3)%alpha(1)": 4e-02,
                    "patch_icpp(1)%pres": 1.0,
                    "patch_icpp(2)%pres": 1.0,
                    "patch_icpp(3)%pres": 1.0,
                },
            )

            for polytropic in ["T", "F"]:
                stack.push("Polytropic" if polytropic == "T" else "", {"polytropic": polytropic})

                for bubble_model in [3, 2]:
                    stack.push(f"bubble_model={bubble_model}", {"bubble_model": bubble_model})

                    if not (polytropic == "F" and bubble_model == 3):
                        cases.append(define_case_d(stack, "", {}))

                    stack.pop()

                stack.pop()

            stack.push("", {"polytropic": "T", "bubble_model": 2})
            cases.append(define_case_d(stack, "nb=1", {"nb": 1}))

            stack.push("QBMM", {"qbmm": "T"})
            cases.append(define_case_d(stack, "", {}))

            stack.push("bubble_model=3", {"bubble_model": 3})
            cases.append(define_case_d(stack, "", {}))

            stack.push("cfl_adap_dt=T", {"cfl_adap_dt": "T", "cfl_target": 0.8, "t_save": 0.01, "n_start": 0, "t_stop": 0.01, "m": 24})
            cases.append(define_case_d(stack, "", {}))

            stack.pop()

            stack.push("cfl_const_dt=T", {"cfl_const_dt": "T", "cfl_target": 0.8, "t_save": 0.01, "n_start": 0, "t_stop": 0.01, "m": 24})
            cases.append(define_case_d(stack, "", {}))

            for _ in range(6):
                stack.pop()

    def alter_lag_bubbles(dimInfo):
        # Lagrangian bubbles
        if len(dimInfo[0]) > 1:
            for adap_dt in ["F", "T"]:
                for couplingMethod in [1, 2]:
                    stack.push(
                        "Lagrange Bubbles",
                        {
                            "bubbles_lagrange": "T",
                            "dt": 1e-06,
                            "lag_params%pressure_corrector": "T",
                            "bubble_model": 2,
                            "num_fluids": 2,
                            "lag_params%heatTransfer_model": "T",
                            "lag_params%massTransfer_model": "T",
                            "fluid_pp(1)%gamma": 0.16,
                            "fluid_pp(1)%eos": "stiffened_gas",
                            "fluid_pp(1)%pi_inf": 3515.0,
                            "fluid_pp(2)%gamma": 2.5,
                            "fluid_pp(2)%eos": "ideal_gas",
                            "fluid_pp(2)%pi_inf": 0.0,
                            "patch_icpp(1)%alpha_rho(1)": 0.96,
                            "patch_icpp(1)%alpha(1)": 4e-02,
                            "patch_icpp(1)%alpha_rho(2)": 0.0,
                            "patch_icpp(1)%alpha(2)": 0.0,
                            "patch_icpp(2)%alpha_rho(1)": 0.96,
                            "patch_icpp(2)%alpha(1)": 4e-02,
                            "patch_icpp(2)%alpha_rho(2)": 0.0,
                            "patch_icpp(2)%alpha(2)": 0.0,
                            "patch_icpp(3)%alpha_rho(1)": 0.96,
                            "patch_icpp(3)%alpha(1)": 4e-02,
                            "patch_icpp(3)%alpha_rho(2)": 0.0,
                            "patch_icpp(3)%alpha(2)": 0.0,
                            "patch_icpp(1)%pres": 1.0,
                            "patch_icpp(2)%pres": 1.0,
                            "patch_icpp(3)%pres": 1.0,
                            "acoustic_source": "T",
                            "acoustic(1)%loc(2)": 0.5,
                            "acoustic(1)%wavelength": 0.25,
                            "acoustic(1)%mag": 2e04,
                            "t_step_start": 0,
                            "t_step_stop": 50,
                            "t_step_save": 50,
                            "lag_txt_wrt": "T",
                            "lag_header": "T",
                            "lag_db_wrt": "T",
                            "lag_id_wrt": "T",
                            "lag_pos_wrt": "T",
                            "lag_pos_prev_wrt": "T",
                            "lag_vel_wrt": "T",
                            "lag_rad_wrt": "T",
                            "lag_rvel_wrt": "T",
                            "lag_r0_wrt": "T",
                            "lag_rmax_wrt": "T",
                            "lag_rmin_wrt": "T",
                            "lag_dphidt_wrt": "T",
                            "lag_pres_wrt": "T",
                            "lag_mv_wrt": "T",
                            "lag_mg_wrt": "T",
                            "lag_betaT_wrt": "T",
                            "lag_betaC_wrt": "T",
                            "lag_params%write_bubbles": "T",
                            "lag_params%write_bubbles_stats": "T",
                            "lag_params%write_void_evol": "T",
                            "lag_params%valmaxvoid": 0.99,
                            "lag_params%nBubs_glb": 1,
                            "polytropic": "F",
                            "bub_pp%R0ref": 1.0,
                            "bub_pp%p0ref": 1.0,
                            "bub_pp%rho0ref": 1.0,
                            "bub_pp%T0ref": 1.0,
                            "bub_pp%ss": 7.131653759435349e-07,
                            "bub_pp%pv": 0.02292716400352907,
                            "bub_pp%vd": 2.4752475247524753e-06,
                            "bub_pp%mu_l": 9.920792079207921e-08,
                            "bub_pp%gam_v": 1.33,
                            "bub_pp%gam_g": 1.4,
                            "bub_pp%M_v": 18.02,
                            "bub_pp%M_g": 28.97,
                            "bub_pp%k_v": 5.618695895665441e-06,
                            "bub_pp%k_g": 7.392868685947116e-06,
                            "bub_pp%R_v": 1347.810235139403,
                            "bub_pp%R_g": 838.3686723235085,
                            "bub_pp%cp_g": 2921.2822272326243,
                            "bub_pp%cp_v": 6134.692677188511,
                        },
                    )

                    if len(dimInfo[0]) == 2:
                        stack.push("", {"acoustic(1)%support": 2, "lag_params%charwidth": 2, "lag_params%charNz": 25})
                    else:
                        stack.push("", {"acoustic(1)%support": 3, "acoustic(1)%height": 1e10})

                    if couplingMethod == 1:
                        stack.push("One-way Coupling", {"lag_params%solver_approach": 1})
                    else:
                        stack.push("Two-way Coupling", {"lag_params%solver_approach": 2})

                    if adap_dt == "F":
                        stack.push("", {})
                    else:
                        stack.push("adap_dt=T", {"adap_dt": "T"})

                    cases.append(define_case_d(stack, "", {}))

                    # Regression guard for issue #1706: the Lagrange bubble initial gas
                    # pressure open-coded the stiffened-gas inversion and dropped the qv
                    # term, so seeding qv /= 0 is the only way to catch it. Every other
                    # lag-bubble case leaves qv at its 0 default, which is why the bug
                    # survived. Restricted to one configuration to add a single golden.
                    if len(dimInfo[0]) == 2 and couplingMethod == 1 and adap_dt == "F":
                        stack.push("qv_nonzero", {"fluid_pp(1)%qv": 0.01})
                        cases.append(define_case_d(stack, "", {}))
                        stack.pop()

                    if len(dimInfo[0]) == 3 and couplingMethod == 2:
                        stack.push("Tracer Bubbles", {"lag_params%vel_model": 1, "fd_order": 2})
                        cases.append(define_case_d(stack, "", {}))
                        stack.pop()

                        stack.push(
                            "Inertial Bubbles",
                            {"lag_params%vel_model": 2, "viscous": "T", "fluid_pp(1)%Re(1)": 100.0, "fluid_pp(2)%Re(1)": 100.0},
                        )
                        if adap_dt == "F":
                            inertial_matrix = [(d, f) for d in [0, 1, 2] for f in [1, 2, 4]]
                        else:
                            inertial_matrix = [(0, 1), (1, 2), (2, 4)]
                        for dragModel, fdOrder in inertial_matrix:
                            stack.push(f"drag_model={dragModel}", {"lag_params%drag_model": dragModel})
                            stack.push(f"fd_order={fdOrder}", {"fd_order": fdOrder})
                            cases.append(define_case_d(stack, "", {}))
                            stack.pop()
                            stack.pop()
                        stack.pop()

                    if len(dimInfo[0]) == 2 and couplingMethod == 1:
                        stack.push("Tracer Bubbles", {"lag_params%vel_model": 1, "fd_order": 2})
                        cases.append(define_case_d(stack, "", {}))
                        stack.pop()

                        stack.push(
                            "Inertial Bubbles",
                            {"lag_params%vel_model": 2, "lag_params%drag_model": 2, "fd_order": 2, "viscous": "T", "fluid_pp(1)%Re(1)": 100.0, "fluid_pp(2)%Re(1)": 100.0},
                        )
                        cases.append(define_case_d(stack, "", {}))
                        stack.pop()

                    stack.pop()

                    stack.pop()

                    stack.pop()

                    stack.pop()

    def alter_elliptic_smoothing():
        # Elliptic Smoothing

        stack.push("Smoothing", {"elliptic_smoothing": "T", "elliptic_smoothing_iters": 10})

        cases.append(define_case_d(stack, "", {}))

        stack.pop()

    def alter_bc_patches(dimInfo):
        # BC_Patches

        stack.push("BC Patches", {"num_bc_patches": 1})

        if len(dimInfo[0]) > 2:
            for direc in [1, 2, 3]:
                stack.push(
                    "Circle",
                    {
                        "patch_bc(1)%geometry": 2,
                        "patch_bc(1)%dir": direc,
                        "patch_bc(1)%type": -17,
                        "patch_bc(1)%loc": -1,
                    },
                )

                if direc == 1:
                    stack.push(
                        "X",
                        {
                            "patch_bc(1)%centroid(2)": 0,
                            "patch_bc(1)%centroid(3)": 0,
                            "patch_bc(1)%radius": 0.000125,
                        },
                    )
                elif direc == 2:
                    stack.push(
                        "Y",
                        {
                            "patch_bc(1)%centroid(1)": 0,
                            "patch_bc(1)%centroid(3)": 0,
                            "patch_bc(1)%radius": 0.000125,
                        },
                    )
                else:
                    stack.push(
                        "Z",
                        {
                            "patch_bc(1)%centroid(1)": 0,
                            "patch_bc(1)%centroid(2)": 0,
                            "patch_bc(1)%radius": 0.000125,
                        },
                    )

                cases.append(define_case_d(stack, "", {}))

                stack.pop()

                stack.pop()

        elif len(dimInfo[0]) > 1:
            for direc in [1, 2]:
                stack.push("Line Segment", {"patch_bc(1)%geometry": 1, "patch_bc(1)%dir": direc, "patch_bc(1)%type": -17, "patch_bc(1)%loc": -1})

                if direc == 1:
                    stack.push("X", {"patch_bc(1)%centroid(2)": 0.0, "patch_bc(1)%length(2)": 0.0025})
                else:
                    stack.push("Y", {"patch_bc(1)%centroid(1)": 0.0, "patch_bc(1)%length(1)": 0.0025})

                cases.append(define_case_d(stack, "", {}))

                stack.pop()

                stack.pop()

        stack.pop()

    def mhd_cases():
        params = {
            "1D": {"m": 200, "dt": 0.001, "t_step_stop": 200, "t_step_save": 200},
            "2D": {"m": 50, "n": 50, "dt": 0.002, "t_step_stop": 500, "t_step_save": 500},
            "3D": {"m": 25, "n": 25, "p": 25, "dt": 0.005, "t_step_stop": 200, "t_step_save": 200},
        }

        case_specs = [
            ("1D -> MHD -> HLL", "examples/1D_brio_wu/case.py", params["1D"]),
            ("1D -> MHD -> HLLD", "examples/1D_brio_wu_hlld/case.py", params["1D"]),
            ("1D -> RMHD", "examples/1D_brio_wu_rmhd/case.py", params["1D"]),
            ("2D -> MHD -> HLL", "examples/2D_orszag_tang/case.py", params["2D"]),
            ("2D -> MHD -> HLLD", "examples/2D_orszag_tang/case.py", {**params["2D"], "riemann_solver": 4}),
            ("2D -> MHD -> hyper_cleaning", "examples/2D_orszag_tang_hyper_cleaning/case.py", params["2D"]),
            ("2D -> RMHD", "examples/2D_shock_cloud_rmhd/case.py", params["2D"]),
            ("3D -> MHD", "examples/3D_brio_wu/case.py", params["3D"]),
            ("3D -> RMHD", "examples/3D_brio_wu/case.py", {**params["3D"], "relativity": "T"}),
        ]

        for name, path, param in case_specs:
            cases.append(define_case_f(name, path, mods=param))

    def hypo_example_cases():
        # Inline Riemann problem configs for hypoelastic solver regression testing.
        # Two-material pressure discontinuity: liquid (fluid 1) vs solid (fluid 2).
        # Enriched ICs: non-zero transverse velocity and initial stress to exercise
        # all solver code paths including shear and geometry source terms.
        _eps = 1e-8
        _fl_g = 1.0e00 / (4.4e00 - 1.0e00)
        _fl_p = 4.4e00 * 5.57e08 / (4.4e00 - 1.0e00)
        _fluids = {
            "fluid_pp(1)%gamma": _fl_g,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%pi_inf": _fl_p,
            "fluid_pp(1)%G": 0.0,
            "fluid_pp(2)%gamma": _fl_g,
            "fluid_pp(2)%eos": "stiffened_gas",
            "fluid_pp(2)%pi_inf": _fl_p,
            "fluid_pp(2)%G": 1e7,
        }
        _common = {
            "run_time_info": "T",
            "t_step_start": 0,
            "t_step_stop": 10,
            "t_step_save": 10,
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "F",
            "time_stepper": 1,
            "weno_order": 1,
            "weno_eps": 1.0e-20,
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 4,
            "wave_speeds": 1,
            "avg_state": 2,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "rho_wrt": "T",
            "parallel_io": "T",
            "hypoelasticity": "T",
            "fd_order": 4,
            **_fluids,
        }
        _patch1_2d = {
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 10.0,
            "patch_icpp(1)%pres": 1e6,
            "patch_icpp(1)%tau_e(1)": 1e4,
            "patch_icpp(1)%alpha_rho(1)": 1000 * (1.0 - _eps),
            "patch_icpp(1)%alpha(1)": 1.0 - _eps,
            "patch_icpp(1)%alpha_rho(2)": 1000 * _eps,
            "patch_icpp(1)%alpha(2)": _eps,
        }
        _patch2_2d = {
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": -10.0,
            "patch_icpp(2)%pres": 1e5,
            "patch_icpp(2)%tau_e(1)": -1e4,
            "patch_icpp(2)%alpha_rho(1)": 1000 * _eps,
            "patch_icpp(2)%alpha(1)": _eps,
            "patch_icpp(2)%alpha_rho(2)": 1000 * (1.0 - _eps),
            "patch_icpp(2)%alpha(2)": 1.0 - _eps,
        }
        base_configs = {
            "2D -> Hypoelasticity": {
                **_common,
                "m": 24,
                "n": 24,
                "p": 0,
                "dt": 6.0e-6,
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                **_patch1_2d,
                **_patch2_2d,
            },
            "2D -> Axisymmetric -> Hypoelasticity": {
                **_common,
                "m": 24,
                "n": 24,
                "p": 0,
                "dt": 4.0e-6,
                "cyl_coord": "T",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -2,
                "bc_x%end": -2,
                "bc_y%beg": -2,
                "bc_y%end": -2,
                **_patch1_2d,
                **_patch2_2d,
                # Axis regularity: vel(2) is the radial velocity here, and a nonzero
                # value at r = 0 is not a smooth axisymmetric state (v_r must be O(r)).
                # Move the two-stream dynamics to the axial component instead.
                "patch_icpp(1)%vel(1)": 10.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(2)%vel(1)": -10.0,
                "patch_icpp(2)%vel(2)": 0.0,
            },
            "3D -> Hypoelasticity": {
                **_common,
                "m": 24,
                "n": 24,
                "p": 24,
                "dt": 6.0e-6,
                "cyl_coord": "F",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "z_domain%beg": 0.0,
                "z_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "bc_z%beg": -3,
                "bc_z%end": -3,
                "patch_icpp(1)%geometry": 9,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%z_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%length_z": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 10.0,
                "patch_icpp(1)%vel(3)": 5.0,
                "patch_icpp(1)%pres": 1e6,
                "patch_icpp(1)%tau_e(1)": 1e4,
                "patch_icpp(1)%alpha_rho(1)": 1000 * (1.0 - _eps),
                "patch_icpp(1)%alpha(1)": 1.0 - _eps,
                "patch_icpp(1)%alpha_rho(2)": 1000 * _eps,
                "patch_icpp(1)%alpha(2)": _eps,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%geometry": 9,
                "patch_icpp(2)%x_centroid": 0.75,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(2)%z_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.5,
                "patch_icpp(2)%length_y": 1.0,
                "patch_icpp(2)%length_z": 1.0,
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(2)%vel(2)": -10.0,
                "patch_icpp(2)%vel(3)": -5.0,
                "patch_icpp(2)%pres": 1e5,
                "patch_icpp(2)%tau_e(1)": -1e4,
                "patch_icpp(2)%alpha_rho(1)": 1000 * _eps,
                "patch_icpp(2)%alpha(1)": _eps,
                "patch_icpp(2)%alpha_rho(2)": 1000 * (1.0 - _eps),
                "patch_icpp(2)%alpha(2)": 1.0 - _eps,
            },
        }

        solver_specs = [
            {"trace": "HLLD", "mods": {"riemann_solver": 4}},
            {"trace": "HLLD -> ADC", "mods": {"riemann_solver": 4, "riemann_hypo_ADC": "T", "ADC_kappa": 1.0}},
            {"trace": "HLLC", "mods": {"riemann_solver": 2}},
            {"trace": "HLLC -> ADC", "mods": {"riemann_solver": 2, "riemann_hypo_ADC": "T", "ADC_kappa": 1.0}},
            {"trace": "HLL -> Interface RHS", "mods": {"riemann_solver": 1, "hypo_hll_interface_rhs": "T"}},
            {"trace": "HLL -> u-interface -> Interface RHS", "mods": {"riemann_solver": 1, "hypo_hll_interface_rhs": "T", "hll_u_interface": "T"}},
        ]

        def apply_solver(case: dict, base_cfg: dict, solver_mods: dict, alt_soundspeed: str):
            case.update(base_cfg)
            for key in ["riemann_hypo_ADC", "ADC_kappa", "hypo_hll_interface_rhs", "hll_u_interface"]:
                case.pop(key, None)
            case["alt_soundspeed"] = alt_soundspeed
            case.update(solver_mods)

        for base_trace, base_cfg in base_configs.items():
            for solver_spec in solver_specs:
                solver_trace = solver_spec["trace"]
                solver_mods = solver_spec["mods"]
                for alt_soundspeed in ["F", "T"]:
                    # 2D axisymmetric HLL Method 1 + alt_soundspeed remains intentionally unsupported.
                    if base_trace == "2D -> Axisymmetric -> Hypoelasticity" and solver_trace == "HLL -> Interface RHS" and alt_soundspeed == "T":
                        continue

                    # The shear-stress rows of this case stay near zero, so the absolute
                    # tolerance is the binding comparison there and compiler/backend roundoff
                    # can exceed the suite default. override_tol is case-wide.
                    is_axisym_hlld_no_alt_soundspeed = base_trace == "2D -> Axisymmetric -> Hypoelasticity" and solver_trace in {"HLLD", "HLLD -> ADC"} and alt_soundspeed == "F"
                    tol = 1e-5 if is_axisym_hlld_no_alt_soundspeed else None

                    trace = f"{base_trace} -> {solver_trace} -> alt_soundspeed={alt_soundspeed}"
                    cases.append(
                        define_case_f(
                            trace,
                            "",
                            mods={},
                            override_tol=tol,
                            functor=lambda case, bc=base_cfg, sm=solver_mods, alt_ss=alt_soundspeed: apply_solver(case, bc, sm, alt_ss),
                        )
                    )

        multifluid_cfg = {
            **base_configs["2D -> Hypoelasticity"],
            "num_fluids": 3,
            "alt_soundspeed": "F",
            "fluid_pp(3)%gamma": _fl_g,
            "fluid_pp(3)%eos": "stiffened_gas",
            "fluid_pp(3)%pi_inf": _fl_p,
            "fluid_pp(3)%G": 5e6,
            "patch_icpp(1)%alpha_rho(1)": 600.0,
            "patch_icpp(1)%alpha(1)": 0.6,
            "patch_icpp(1)%alpha_rho(2)": 300.0,
            "patch_icpp(1)%alpha(2)": 0.3,
            "patch_icpp(1)%alpha_rho(3)": 100.0,
            "patch_icpp(1)%alpha(3)": 0.1,
            "patch_icpp(2)%alpha_rho(1)": 100.0,
            "patch_icpp(2)%alpha(1)": 0.1,
            "patch_icpp(2)%alpha_rho(2)": 300.0,
            "patch_icpp(2)%alpha(2)": 0.3,
            "patch_icpp(2)%alpha_rho(3)": 600.0,
            "patch_icpp(2)%alpha(3)": 0.6,
        }
        multifluid_solvers = [
            ("HLL -> u-interface -> Interface RHS", {"riemann_solver": 1, "hypo_hll_interface_rhs": "T", "hll_u_interface": "T"}),
            ("HLLC", {"riemann_solver": 2}),
        ]
        for solver_trace, solver_mods in multifluid_solvers:
            cases.append(
                define_case_f(
                    f"2D -> Hypoelasticity -> 3 Fluid(s) -> {solver_trace} -> alt_soundspeed=F",
                    "",
                    mods={**multifluid_cfg, **solver_mods},
                )
            )

        edge_common = {
            "run_time_info": "F",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "n": 7,
            "p": 0,
            "t_step_start": 0,
            "t_step_stop": 1,
            "t_step_save": 1,
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 1,
            "recon_type": 1,
            "weno_order": 1,
            "weno_eps": 1.0e-16,
            "riemann_solver": 4,
            "wave_speeds": 1,
            "avg_state": 2,
            "riemann_hypo_ADC": "F",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "format": 1,
            "precision": 2,
            "parallel_io": "F",
            "hypoelasticity": "T",
            "fd_order": 4,
        }

        # A solid-anchored, light-gas face puts the raw right shear wave outside
        # the outer fan. Tangential slip makes HLL fallback differ from clipping.
        gamma = 1.4
        p0 = 1.0
        G_solid = 1.0
        rho_solid = 1.0
        rho_gas = 1.0e-2
        eps_alpha = 1.0e-12
        nx = 32
        dx = 1.0 / nx
        c_outer = math.sqrt(gamma * p0 / rho_gas)
        fan_order_cfg = {
            **edge_common,
            "m": nx - 1,
            "dt": 0.2 * dx / c_outer,
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%G": G_solid,
            "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%G": 0.0,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%tau_e(1)": 0.0,
            "patch_icpp(1)%tau_e(2)": 0.0,
            "patch_icpp(1)%tau_e(3)": 0.0,
            "patch_icpp(1)%alpha_rho(1)": rho_solid * (1.0 - eps_alpha),
            "patch_icpp(1)%alpha(1)": 1.0 - eps_alpha,
            "patch_icpp(1)%alpha_rho(2)": rho_solid * eps_alpha,
            "patch_icpp(1)%alpha(2)": eps_alpha,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 1.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%tau_e(1)": 0.0,
            "patch_icpp(2)%tau_e(2)": 0.0,
            "patch_icpp(2)%tau_e(3)": 0.0,
            "patch_icpp(2)%alpha_rho(1)": rho_gas * eps_alpha,
            "patch_icpp(2)%alpha(1)": eps_alpha,
            "patch_icpp(2)%alpha_rho(2)": rho_gas * (1.0 - eps_alpha),
            "patch_icpp(2)%alpha(2)": 1.0 - eps_alpha,
        }
        cases.append(
            define_case_f(
                "2D -> Hypoelasticity -> HLLD -> Fan-order fallback",
                "",
                mods=fan_order_cfg,
            )
        )

        # Exercise the two changed G_eff = G_hat + tau_nn_hat outcomes: resolved
        # positive HLLD and negative-state HLL fallback. The exact-zero HLLC
        # limit is unchanged by this patch.

        def make_shear_guard_cfg(regime):
            gamma = 1.4
            if regime == "positive":
                G = 1.0e-8
                rho0 = 1.0e-5
                p0 = 1.0e-5
                tau_xx = 0.0
            else:
                G = 1.0
                rho0 = 1.0
                p0 = 1.0
                tau_xx = -1.1 * G
            tau_yy = -tau_xx

            nx = 16
            dx = 1.0 / nx
            c_outer = math.sqrt((gamma * p0 + 4.0 * G / 3.0 + tau_xx) / rho0)
            return {
                **edge_common,
                "m": nx - 1,
                "dt": 0.2 * dx / c_outer,
                "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
                "fluid_pp(1)%eos": "ideal_gas",
                "fluid_pp(1)%pi_inf": 0.0,
                "fluid_pp(1)%G": G,
                "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%G": G,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 1.0,
                "patch_icpp(1)%pres": p0,
                "patch_icpp(1)%tau_e(1)": tau_xx,
                "patch_icpp(1)%tau_e(2)": 0.0,
                "patch_icpp(1)%tau_e(3)": tau_yy,
                "patch_icpp(1)%alpha_rho(1)": 0.5 * rho0,
                "patch_icpp(1)%alpha(1)": 0.5,
                "patch_icpp(1)%alpha_rho(2)": 0.5 * rho0,
                "patch_icpp(1)%alpha(2)": 0.5,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%geometry": 3,
                "patch_icpp(2)%x_centroid": 0.75,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.5,
                "patch_icpp(2)%length_y": 1.0,
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(2)%vel(2)": 2.0,
                "patch_icpp(2)%pres": p0,
                "patch_icpp(2)%tau_e(1)": tau_xx,
                "patch_icpp(2)%tau_e(2)": 0.0,
                "patch_icpp(2)%tau_e(3)": tau_yy,
                "patch_icpp(2)%alpha_rho(1)": 0.5 * rho0,
                "patch_icpp(2)%alpha(1)": 0.5,
                "patch_icpp(2)%alpha_rho(2)": 0.5 * rho0,
                "patch_icpp(2)%alpha(2)": 0.5,
            }

        for label, regime, override_tol in [
            # The positive-case signal is O(1e-8), below the generic hypoelastic
            # absolute tolerance; 1e-9 keeps the regression discriminating while
            # remaining loose relative to single-precision roundoff at O(1e-5).
            ("Resolved positive", "positive", 1e-9),
            ("Negative fallback", "negative", None),
        ]:
            cases.append(
                define_case_f(
                    f"2D -> Hypoelasticity -> HLLD -> Shear guard -> {label}",
                    "",
                    mods=make_shear_guard_cfg(regime),
                    override_tol=override_tol,
                )
            )

        # Axisymmetric fallback-face trace consistency: tau_rr = -1.1*G makes the
        # radial-sweep fan invalid everywhere (G_eff < 0), so the radial geometric
        # source and NC velocity traces must come from the HLL/one-sided face state
        # that matches the fallback flux -- not from the rejected fan's contact
        # speed and star states. The stress IC is deviatoric with tau_rr = tau_qq
        # (axis-regular); the radial band contrast (rho, p, axial velocity) makes
        # the two trace families measurably different at the band interface. The
        # ADC twin pins that ADC on/off is identical on fallback-selected faces.
        def make_axisym_fallback_cfg(adc):
            G = 1.0
            cfg = {
                **edge_common,
                "m": 15,
                "cyl_coord": "T",
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -2,
                "bc_y%end": -3,
                "dt": 2.0e-3,
                "fluid_pp(1)%gamma": 2.5,
                "fluid_pp(1)%eos": "ideal_gas",
                "fluid_pp(1)%pi_inf": 0.0,
                "fluid_pp(1)%G": G,
                "fluid_pp(2)%gamma": 2.5,
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%G": G,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%tau_e(1)": 2.2 * G,
                "patch_icpp(1)%tau_e(2)": 0.0,
                "patch_icpp(1)%tau_e(3)": -1.1 * G,
                "patch_icpp(1)%tau_e(4)": -1.1 * G,
                "patch_icpp(1)%alpha_rho(1)": 0.5,
                "patch_icpp(1)%alpha(1)": 0.5,
                "patch_icpp(1)%alpha_rho(2)": 0.5,
                "patch_icpp(1)%alpha(2)": 0.5,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%geometry": 3,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%y_centroid": 0.75,
                "patch_icpp(2)%length_x": 1.0,
                "patch_icpp(2)%length_y": 0.5,
                "patch_icpp(2)%vel(1)": 0.5,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(2)%pres": 1.5,
                "patch_icpp(2)%tau_e(1)": 2.2 * G,
                "patch_icpp(2)%tau_e(2)": 0.0,
                "patch_icpp(2)%tau_e(3)": -1.1 * G,
                "patch_icpp(2)%tau_e(4)": -1.1 * G,
                "patch_icpp(2)%alpha_rho(1)": 0.25,
                "patch_icpp(2)%alpha(1)": 0.5,
                "patch_icpp(2)%alpha_rho(2)": 0.25,
                "patch_icpp(2)%alpha(2)": 0.5,
            }
            if adc:
                cfg.update({"riemann_hypo_ADC": "T", "ADC_kappa": 1.0})
            return cfg

        for label, adc in [("", False), (" -> ADC", True)]:
            cases.append(
                define_case_f(
                    f"2D -> Axisymmetric -> Hypoelasticity -> HLLD -> Fallback traces{label}",
                    "",
                    mods=make_axisym_fallback_cfg(adc),
                )
            )

        # Uniform-state preservation across a characteristic boundary: the CBC flux
        # rebuild must reproduce HLLD's folded volume-fraction representation
        # (adv_src_mode_none; -/+ K*u_n under alt_soundspeed). The WENO5 path is the
        # one that calls the primitive-to-flux converter at the boundary. Pre-fix,
        # this alt_soundspeed case drifted alpha by O(1e-3) within five steps; the
        # golden pins exact uniformity. A uniform 50/50 mixture makes K large enough
        # for a robust regression. In simulations, HLLD characteristic boundaries
        # are limited to fluid regions; material interfaces must remain away.
        def make_cbc_uniform_cfg(alt_soundspeed):
            nx, ny = 64, 32
            alpha_w = 0.5
            return {
                **edge_common,
                "m": nx - 1,
                "n": ny - 1,
                "t_step_stop": 5,
                "dt": 1.0e-6,
                "time_stepper": 3,
                "weno_order": 5,
                "mapped_weno": "T",
                "alt_soundspeed": alt_soundspeed,
                "num_patches": 1,
                "bc_x%beg": -6,
                "bc_x%end": -6,
                "fluid_pp(1)%gamma": 1.0 / (4.4 - 1.0),
                "fluid_pp(1)%eos": "stiffened_gas",
                "fluid_pp(1)%pi_inf": 4.4 * 6.0e8 / (4.4 - 1.0),
                "fluid_pp(1)%G": 0.0,
                "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%G": 0.0,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 10.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0e5,
                "patch_icpp(1)%tau_e(1)": 0.0,
                "patch_icpp(1)%tau_e(2)": 0.0,
                "patch_icpp(1)%tau_e(3)": 0.0,
                "patch_icpp(1)%alpha_rho(1)": alpha_w * 1000.0,
                "patch_icpp(1)%alpha(1)": alpha_w,
                "patch_icpp(1)%alpha_rho(2)": (1.0 - alpha_w) * 1.2,
                "patch_icpp(1)%alpha(2)": 1.0 - alpha_w,
            }

        cases.append(
            define_case_f(
                "2D -> Hypoelasticity -> HLLD -> CBC uniform preservation -> alt_soundspeed",
                "",
                mods=make_cbc_uniform_cfg("T"),
            )
        )

        # Pins the recommended WENO5-M + RK3 configuration on interior two-patch
        # dynamics (the CBC cases share the discretization but test the boundary).
        weno5_cfg = {
            **make_cbc_uniform_cfg("T"),
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "t_step_stop": 10,
            "num_patches": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%y_centroid": 0.5,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%length_y": 1.0,
            "patch_icpp(2)%vel(1)": 10.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": 2.0e5,
            "patch_icpp(2)%tau_e(1)": 1.0e3,
            "patch_icpp(2)%tau_e(2)": 5.0e2,
            "patch_icpp(2)%tau_e(3)": -1.0e3,
            "patch_icpp(2)%alpha_rho(1)": 0.9 * 1000.0,
            "patch_icpp(2)%alpha(1)": 0.9,
            "patch_icpp(2)%alpha_rho(2)": 0.1 * 1.2,
            "patch_icpp(2)%alpha(2)": 0.1,
            "fluid_pp(1)%G": 1.0e6,
            "fluid_pp(2)%G": 0.0,
        }
        cases.append(
            define_case_f(
                "2D -> Hypoelasticity -> HLLD -> WENO5 RK3",
                "",
                mods=weno5_cfg,
            )
        )

        # Rigid axial translation in axisymmetric coordinates: every stress rate must
        # vanish identically (the pre-repair operator produced d(tau_xr)/dt = G*U0/r).
        # The golden pins the invariant as a permanent semantic regression.
        rigid_cfg = {
            **edge_common,
            "m": 23,
            "n": 23,
            "t_step_stop": 10,
            "dt": 1.0e-6,
            "cyl_coord": "T",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -2,
            "num_patches": 1,
            "fluid_pp(1)%gamma": 1.0 / (4.4 - 1.0),
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%pi_inf": 4.4 * 6.0e8 / (4.4 - 1.0),
            "fluid_pp(1)%G": 1.0e6,
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%G": 1.0e6,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 10.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0e5,
            "patch_icpp(1)%tau_e(1)": 0.0,
            "patch_icpp(1)%tau_e(2)": 0.0,
            "patch_icpp(1)%tau_e(3)": 0.0,
            "patch_icpp(1)%alpha_rho(1)": 0.5 * 1000.0,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": 0.5 * 1.2,
            "patch_icpp(1)%alpha(2)": 0.5,
        }
        cases.append(
            define_case_f(
                "2D -> Axisymmetric -> Hypoelasticity -> HLLD -> Rigid axial translation",
                "",
                mods=rigid_cfg,
            )
        )

    def foreach_dimension():
        for dimInfo, dimParams in get_dimensions():
            stack.push(f"{len(dimInfo[0])}D", dimParams)
            alter_bcs(dimInfo)
            alter_grcbc(dimInfo)
            alter_weno(dimInfo)
            alter_muscl()
            alter_num_fluids(dimInfo)
            if len(dimInfo[0]) == 2:
                alter_2d()
            if len(dimInfo[0]) == 3:
                alter_3d()
            alter_lag_bubbles(dimInfo)
            alter_ppn(dimInfo)
            stack.push("", {"dt": [1e-07, 1e-06, 1e-06][len(dimInfo[0]) - 1]})
            alter_acoustic_src(dimInfo)
            alter_bubbles(dimInfo)
            alter_hypoelasticity(dimInfo)
            alter_phasechange(dimInfo)
            alter_viscosity(dimInfo)
            alter_elliptic_smoothing()
            alter_body_forces(dimInfo)
            alter_synthetic_turbulence(dimInfo)
            alter_mixlayer_perturb(dimInfo)
            alter_bc_patches(dimInfo)
            stack.pop()
            stack.pop()

    def foreach_example():
        for path in os.listdir(common.MFC_EXAMPLE_DIRPATH):
            if path == "scaling":
                continue

            # # List of all example cases that will be skipped during testing
            casesToSkip = [
                "2D_ibm_cfl_dt",
                "1D_sodHypo",
                "2D_viscous",
                "2D_laplace_pressure_jump",
                "2D_bubbly_steady_shock",
                "2D_advection",
                "2D_hardcoded_ic",
                # File-based IC (hcid=273/274) sized to the full grid; the Example
                # suite's m/n cap breaks it. Covered by the Chemistry golden tests.
                "2D_reacting_mixing_layer",
                "2D_spatial_reacting_mixing_layer",
                "2D_ibm_multiphase",
                "2D_acoustic_broadband",
                "1D_inert_shocktube",
                "1D_reactive_shocktube",
                "1D_reactive_shocktube_adaptive",
                "2D_detonation_cell",
                "2D_ibm_steady_shock",
                "3D_performance_test",
                "3D_ibm_stl_ellipsoid",
                "3D_sphbubcollapse",
                "2D_ibm_stl_wedge",
                "3D_ibm_stl_pyramid",
                "3D_ibm_bowshock",
                "3D_turb_mixing",
                "2D_mixing_artificial_Ma",
                "2D_lagrange_bubblescreen",
                "3D_lagrange_bubblescreen",
                "2D_triple_point",
                "1D_shuosher_analytical",
                "1D_titarevtorro_analytical",
                "2D_acoustic_pulse_analytical",
                "2D_isentropicvortex_analytical",
                "1D_euler_convergence",
                "1D_advection_convergence",
                "1D_sod_convergence",
                "2D_advection_convergence",
                "3D_advection_convergence",
                "2D_hypo_shear_contact",  # exercised by the convergence suite
                "2D_zero_circ_vortex_analytical",
                "3D_TaylorGreenVortex_analytical",
                "3D_IGR_TaylorGreenVortex_nvidia",
                "2D_backward_facing_step",
                "2D_forward_facing_step",
                "1D_convergence",
                "3D_IGR_33jet",
                "1D_multispecies_diffusion",
                "2D_ibm_stl_MFCCharacter",
                "1D_qbmm",  # formatted I/O field overflow on gfortran 12
                "2D_moving_lag_bubs",  # adap_dt hangs on reduced grid
                "3D_moving_lag_particles",  # adap_dt hangs on reduced grid
                "2D_premixed_landau_insta",
                "1D_flamelet",
                "2D_premixed_flame_vortex",
                "2D_Thermal_Flatplate",  # formatted I/O field overflow on gfortran 12
                "2D_hypo_hlld",  # acoustic demo case, not a regression test
                "3D_hypo_hlld",  # acoustic demo case, not a regression test
                "2D_axisym_hypo_hlld",  # acoustic demo case, not a regression test
                "2D_lagrange_rising_bubble",
                "2D_lagrange_in_crossflow",
                # Non-Newtonian validation cases whose cfl_adap_dt run is viscous-CFL limited
                # by a large mu_max: even on the downsized grid the step count to reach t_stop
                # is too large for the CI smoke suite. The faster NN examples remain tested.
                "2D_poiseuille_nn",
                "2D_bingham_poiseuille_nn",
                # The CI grid cap (~25x25) thins this case's immersed-boundary wall slabs
                # to ~2 cells, an under-resolved IB whose body-forced dead-fluid dynamics
                # is platform-marginal (CPU goldens fail on most GPU lanes). The fast
                # "Non-Newtonian -> IBM" suite case covers IBM+NN portably at 1e-12.
                "2D_ibm_poiseuille_nn",
                # Synthetic turbulence now uses a deterministic (compiler-independent) PRNG,
                # but the 50-step forced run with a moving airfoil IB is FP-sensitive enough
                # that Intel's aggressive FP model (FMA/fast trig) diverges from the golden on
                # one variable past the 1e-3 Example tolerance, while all other lanes pass.
                # The forcing physics is correct; the golden is just cross-compiler-marginal.
                "2D_synthetic_turbulence",
                # Two immersed boundaries colliding across a periodic boundary via a stiff
                # soft-sphere spring. The spring acts as a strong amplifier: it turns the
                # ~1e-13 CPU/GPU floating-point difference in the hydrodynamic force (the
                # order-dependent atomic surface-pressure integral) into an exponentially
                # growing trajectory divergence, so the sharp-interface field fails the
                # golden tolerance on GPU lanes even though both runs are individually
                # reproducible. Not a correctness bug -- the case is genuinely chaotic at
                # this stiffness, so it is not a portable regression target.
                "3D_mibm_periodic_collision",
                # The violently stiff 3D bubble collapse amplifies compiler/arch floating-point
                # differences past the 1e-3 Example tolerance under the always-pTg phase-change
                # solver (a Newton equilibrium solve per cavitating cell, replacing the old
                # shortcut). CI's own value is compiler-version-dependent -- nvhpc 24.5 and 26.1
                # disagree by ~1.1e-3 (> tol) -- so no single golden passes every lane regardless
                # of where it is generated. The 2D bubble and all 18 phase-change unit tests remain
                # portable and CPU/GPU machine-zero; only this stiff 3D collapse is non-portable.
                "3D_phasechange_bubble",
                # Finite-rate propellant flame: the flame-front position after 50 steps is set by
                # accumulated stiff-kinetics roundoff, so it drifts past the 1e-3 Example tolerance
                # across compilers (nvhpc 25.11 golden vs 24.3/GNU/Intel/CCE/AMD disagree by ~1-5e-3).
                # No single golden is portable -- same class as the other flame examples above. The
                # reactive_burn detonation golden tests remain portable (machine-zero across lanes).
                "1D_propellant_flame",
                # Same finite-rate-flame non-portability as 1D_propellant_flame above: the diffusion
                # flame over the fuel slab is set by accumulated stiff-kinetics roundoff, so by step 50
                # the transverse momentum drifts past the 1e-3 Example tolerance across compilers
                # (nvhpc passes; Intel and CCE disagree by ~2e-3 absolute). No single golden is portable.
                "2D_hybrid_slab",
            ]
            if path in casesToSkip:
                continue
            name = f"{path.split('_')[0]} -> Example -> {'_'.join(path.split('_')[1:])}"
            case_path = os.path.join(common.MFC_EXAMPLE_DIRPATH, path, "case.py")
            if not os.path.isfile(case_path):
                continue

            def modify_example_case(case: dict):
                case["parallel_io"] = "F"
                if "t_step_stop" in case and case["t_step_stop"] >= 50:
                    case["t_step_start"] = 0
                    case["t_step_stop"] = 50
                    case["t_step_save"] = 50

                if case.get("recon_type") == 2:
                    for k in ("weno_order", "weno_eps", "wenoz_q", "teno_CT"):
                        case[k] = None

                caseSize = case["m"] * max(case["n"], 1) * max(case["p"], 1)
                if caseSize > 25 * 25:
                    if case["n"] == 0 and case["p"] == 0:
                        case["m"] = 25 * 25
                    elif case["p"] == 0:
                        case["m"] = 25
                        case["n"] = 25
                    elif caseSize > 25 * 25 * 25:
                        case["m"] = 25
                        case["n"] = 25
                        case["p"] = 25

            cases.append(define_case_f(name, case_path, [], {}, functor=modify_example_case))

    def chemistry_cases():
        common_mods = {"t_step_stop": Nt, "t_step_save": Nt}
        for ndim in range(1, 4):
            cases.append(define_case_f(f"{ndim}D -> Chemistry -> Perfect Reactor", "examples/nD_perfect_reactor/case.py", ["--ndim", str(ndim)], mods=common_mods))

        # Operator-split, sub-stepped reaction integration (chem_params%reaction_substeps > 0): the
        # autoigniting reactor exercises the split path that keeps stiff mechanisms stable.
        cases.append(
            define_case_f(
                "1D -> Chemistry -> Perfect Reactor -> Sub-stepped Reactions",
                "examples/nD_perfect_reactor/case.py",
                ["--ndim", "1"],
                mods={**common_mods, "chem_params%reaction_substeps": 10},
            )
        )

        for riemann_solver, gamma_method in itertools.product([1, 2], [1, 2]):
            cases.append(
                define_case_f(
                    f"1D -> Chemistry -> Inert Shocktube -> Riemann Solver {riemann_solver} -> Gamma Method {gamma_method}",
                    "examples/1D_inert_shocktube/case.py",
                    mods={**common_mods, "riemann_solver": riemann_solver, "chem_params%gamma_method": gamma_method, "weno_order": 3, "mapped_weno": "F", "mp_weno": "F"},
                    override_tol=10 ** (-10),
                )
            )

        # The reacting Roe sound speed - the chemistry average state, c_sum_Yi_Phi, and the
        # c = sqrt(c_c - (gamma - 1)*(vel_sum - H)) branch of s_compute_speed_of_sound_avg - is
        # reached only with avg_state = 1 AND wave_speeds = 2. Every other chemistry case sets
        # wave_speeds = 1, so none of it had coverage. Both solvers: HLLC used to pass a literal 0
        # here and take the frozen branch instead, which is why the gap went unseen (#1774).
        cases.append(
            define_case_f(
                "1D -> Chemistry -> Inert Shocktube -> Reacting Roe Average",
                "examples/1D_inert_shocktube/case.py",
                mods={**common_mods, "riemann_solver": 1, "avg_state": 1, "wave_speeds": 2, "weno_order": 3, "mapped_weno": "F", "mp_weno": "F"},
                override_tol=10 ** (-10),
            )
        )
        cases.append(
            define_case_f(
                "1D -> Chemistry -> Inert Shocktube -> Reacting Roe Average -> HLLC",
                "examples/1D_inert_shocktube/case.py",
                mods={**common_mods, "riemann_solver": 2, "avg_state": 1, "wave_speeds": 2, "weno_order": 3, "mapped_weno": "F", "mp_weno": "F"},
                override_tol=10 ** (-10),
            )
        )

        # 1D -> Chemistry -> Flamelet: temporarily removed from the suite. The stiff flamelet
        # integration is the most FP-sensitive chemistry case; on the Frontier CCE OpenMP-offload
        # backend it diverges from the single-reference golden by ~1e-9 (rel) -- compiler
        # FMA-contraction noise, not a physics difference (every other backend matches to <1e-10,
        # the pinned override_tol). Surfaced by the Riemann device-helper refactor (#1572). Re-enable
        # once goldens are regenerated per-backend or the tolerance model gains backend awareness.
        # cases.append(define_case_f("1D -> Chemistry -> Flamelet", "examples/1D_flamelet/case.py", mods={"t_step_stop": 1, "t_step_save": 1}, override_tol=10 ** (-10)))

        # Both reacting mixing-layer goldens drift across compiler/arch: on the temporal
        # case, a GCC-15 arm64 build (the macOS CI lane's architecture) diverges from the
        # GNU-13-Linux golden by up to 6.65e-8 (rel) on a momentum component after 50 stiff
        # chemistry steps, so the 1e-12 default is not portable. 1e-6 clears that with margin
        # (and headroom for the other compiler lanes) while still catching real regressions,
        # which shift these fields by far more. The spatial case adds bf_spatial_support's
        # Fourier forcing but stays in the same band on this short run.
        cases.append(
            define_case_f(
                "2D -> Chemistry -> Reacting Mixing Layer",
                "examples/2D_reacting_mixing_layer/case.py",
                ["--scale", "0.1"],  # cold (non-reacting) profile by default; see case.py
                mods=common_mods,
                override_tol=10 ** (-6),
            )
        )

        cases.append(
            define_case_f(
                "2D -> Chemistry -> Spatial Reacting Mixing Layer",
                "examples/2D_spatial_reacting_mixing_layer/case.py",
                ["--scale", "0.1"],  # cold (non-reacting) profile by default; see case.py
                mods=common_mods,
                override_tol=10 ** (-6),
            )
        )

        stack.push(
            "1D -> Chemistry -> Dual Isothermal Wall Gradient",
            {
                "m": 49,
                "n": 0,  # 1D case
                "p": 0,
                "dt": 8.0e-08,
                "num_patches": 1,
                "num_fluids": 1,
                "x_domain%beg": 0.0,
                "x_domain%end": 0.05,
                "bc_x%beg": -16,  # Left Isothermal Wall
                "bc_x%end": -16,  # Right Isothermal Wall
                "bc_x%isothermal_in": "T",
                "bc_x%Twall_in": 600.0,
                "bc_x%isothermal_out": "T",
                "bc_x%Twall_out": 900.0,
                "weno_order": 5,
                "weno_eps": 1e-16,
                "mapped_weno": "T",
                "mp_weno": "T",
                "riemann_solver": 2,
                "wave_speeds": 1,
                "avg_state": 2,
                "time_stepper": 3,
                "chemistry": "T",
                "chem_params%diffusion": "T",
                "chem_params%reactions": "F",
                "chem_wrt_T": "T",
                "cantera_file": "h2o2.yaml",
                "viscous": "T",
                "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
                "fluid_pp(1)%eos": "ideal_gas",
                "fluid_pp(1)%pi_inf": 0.0,
                "fluid_pp(1)%Re(1)": 100000,
                "patch_icpp(1)%geometry": 1,
                "patch_icpp(1)%hcid": 191,
                "patch_icpp(1)%x_centroid": 0.025,
                "patch_icpp(1)%length_x": 0.05,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%pres": 101325.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(1)%Y(1)": 1.0,
                "t_step_start": 0,
                "t_step_stop": 1000,
                "t_step_save": 1000,
            },
        )

        cases.append(define_case_d(stack, "", {}, override_tol=10 ** (-10)))

        stack.pop()

        stack.push(
            "2D -> Chemistry -> Isothermal Wall",
            {
                "m": 49,
                "n": 49,
                "dt": 4.0e-08,
                "num_patches": 1,
                "num_fluids": 1,
                "x_domain%beg": 0.0,
                "x_domain%end": 0.05,
                "y_domain%beg": 0.0,
                "y_domain%end": 0.05,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -16,
                "bc_y%end": -3,
                "bc_y%isothermal_in": "T",
                "bc_y%Twall_in": 600.0,
                "weno_order": 5,
                "weno_eps": 1e-16,
                "mapped_weno": "T",
                "mp_weno": "T",
                "riemann_solver": 2,
                "wave_speeds": 1,
                "avg_state": 2,
                "time_stepper": 3,
                "chemistry": "T",
                "chem_params%diffusion": "T",
                "chem_params%reactions": "F",
                "chem_wrt_T": "T",
                "cantera_file": "h2o2.yaml",
                "viscous": "T",
                "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
                "fluid_pp(1)%eos": "ideal_gas",
                "fluid_pp(1)%pi_inf": 0.0,
                "fluid_pp(1)%Re(1)": 100000,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%hcid": 291,
                "patch_icpp(1)%x_centroid": 0.025,
                "patch_icpp(1)%y_centroid": 0.025,
                "patch_icpp(1)%length_x": 0.05,
                "patch_icpp(1)%length_y": 0.05,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 101325.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(1)%Y(1)": 1.0,
                "t_step_start": 0,
                "t_step_stop": 50,
                "t_step_save": 50,
            },
        )
        cases.append(define_case_d(stack, "", {}, override_tol=10 ** (-10)))
        stack.pop()

        stack.push(
            "1D -> Chemistry -> MultiComponent Diffusion",
            {
                "m": 200,
                "dt": 0.1e-06,
                "num_patches": 1,
                "num_fluids": 1,
                "x_domain%beg": 0.0,
                "x_domain%end": 0.05,
                "bc_x%beg": -1,
                "bc_x%end": -1,
                "weno_order": 5,
                "weno_eps": 1e-16,
                "weno_avg": "F",
                "mapped_weno": "T",
                "mp_weno": "T",
                "weno_Re_flux": "F",
                "riemann_solver": 2,
                "wave_speeds": 1,
                "avg_state": 1,
                "chemistry": "T",
                "chem_params%diffusion": "T",
                "chem_params%reactions": "F",
                "chem_wrt_T": "T",
                "patch_icpp(1)%geometry": 1,
                "patch_icpp(1)%hcid": 182,
                "patch_icpp(1)%x_centroid": 0.05 / 2,
                "patch_icpp(1)%length_x": 0.05,
                "patch_icpp(1)%vel(1)": "0",
                "patch_icpp(1)%pres": 1.01325e5,
                "patch_icpp(1)%alpha(1)": 1,
                "fluid_pp(1)%gamma": 1.0e00 / (1.9326e00 - 1.0e00),
                "fluid_pp(1)%eos": "ideal_gas",
                "fluid_pp(1)%pi_inf": 0,
                "cantera_file": "h2o2.yaml",
                "t_step_start": 0,
                "t_step_stop": 50,
                "t_step_save": 50,
            },
        )
        cases.append(define_case_d(stack, "", {}, override_tol=10 ** (-10)))

        stack.pop()

    foreach_dimension()

    mhd_cases()
    hypo_example_cases()

    foreach_example()

    chemistry_cases()

    def reactive_burn_cases():
        # Condensed-phase programmed burn (m_reactive_burn): a uniform two-fluid stiffened-gas
        # reactant above the ignition pressure, so the burn source converts reactant (fluid 1)
        # to product (fluid 2) and releases energy through the qv difference. Uniform IC isolates
        # the source term (no shock/gradient) for a stable, reproducible golden; it also exercises
        # the reactive_burn precondition checks (num_fluids=2, shared EOS, qv(1) > qv(2)).
        stack.push(
            "1D -> Reactive Burn -> Condensed Programmed Detonation",
            {
                "m": 99,
                "n": 0,
                "p": 0,
                "dt": 1.0e-9,
                "num_patches": 1,
                "num_fluids": 2,
                "model_eqns": 2,
                "x_domain%beg": 0.0,
                "x_domain%end": 0.005,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "weno_order": 5,
                "weno_eps": 1e-16,
                "mapped_weno": "F",
                "mp_weno": "F",
                "riemann_solver": 2,
                "wave_speeds": 1,
                "avg_state": 2,
                "time_stepper": 3,
                "reactive_burn": "T",
                "rburn%k": 1.0e7,
                "rburn%pign": 1.0e9,
                "rburn%pref": 2.0e9,
                "rburn%n": 1.0,
                "fluid_pp(1)%gamma": 1.0e00 / (3.0e00 - 1.0e00),
                "fluid_pp(1)%eos": "stiffened_gas",
                "fluid_pp(1)%pi_inf": 9.0e8,
                "fluid_pp(1)%qv": 4.0e6,
                "fluid_pp(2)%gamma": 1.0e00 / (3.0e00 - 1.0e00),
                "fluid_pp(2)%eos": "stiffened_gas",
                "fluid_pp(2)%pi_inf": 9.0e8,
                "fluid_pp(2)%qv": 0.0,
                "patch_icpp(1)%geometry": 1,
                "patch_icpp(1)%x_centroid": 0.0025,
                "patch_icpp(1)%length_x": 0.005,
                # Uniform advection velocity: the field stays spatially uniform (sub-cell drift over
                # 40 steps), so the source term is still isolated, but the momentum output is now a
                # well-resolved O(1e5) quantity instead of near-zero roundoff -- otherwise the golden
                # compares ~0 momentum at 1e-12 tolerance, which is not portable across recompiles/backends.
                "patch_icpp(1)%vel(1)": 100.0,
                "patch_icpp(1)%pres": 2.0e9,
                "patch_icpp(1)%alpha_rho(1)": 1900.0 * 0.95,
                "patch_icpp(1)%alpha_rho(2)": 1900.0 * 0.05,
                "patch_icpp(1)%alpha(1)": 0.95,
                "patch_icpp(1)%alpha(2)": 0.05,
                "t_step_start": 0,
                "t_step_stop": 40,
                "t_step_save": 40,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        # Same burn on the 6-equation model (model_eqns=3): the reactant->product qv release
        # must manifest through the qv-consistent phasic-pressure relaxation. Guards that the
        # 5-eq source term is correct on the 6-eq model and that the qv threading holds up.
        stack.push("6eq", {"model_eqns": 3})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        # Arrhenius temperature-driven rate: rburn%ta > 0 multiplies the rate by exp(-rburn%ta/T_r),
        # with T_r the reactant phasic temperature (needs cv > 0). rburn%ta/cv are chosen so the
        # factor is O(0.4) at the IC temperature -- exercises the branch instead of leaving it ~1.
        stack.push("Arrhenius", {"rburn%ta": 500.0, "fluid_pp(1)%cv": 1500.0, "fluid_pp(2)%cv": 1500.0})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        # Same burn on 2 MPI ranks: the rburn parameters must be broadcast to non-root ranks, or
        # rank 1's half of the domain burns with the sentinel default and diverges. The single-rank
        # goldens cannot catch a broken rburn broadcast; this one does.
        stack.push("2 ranks", {})
        cases.append(define_case_d(stack, "", {}, ppn=2))
        stack.pop()
        stack.pop()

    reactive_burn_cases()

    def ibm_burn_rate_cases():
        """Vieille's-law pressure-coupled IB burn rate (patch_ib%burn_rate_exp/pref).

        The ibm_burning_grain/ibm_flameholder Example goldens cover surface injection
        (v_blow, inj_species) but both leave burn_rate_exp at 0, so the pressure-coupled
        scaling v_blow*(p/pref)^n has no golden. This registers the same example with the
        coupling on, at reduced resolution to stay cheap.
        """
        cases.append(
            define_case_f(
                "2D -> IBM -> Vieille Burn Rate",
                "examples/2D_ibm_burning_grain/case.py",
                ["--burn_exp", "0.5"],
                mods={"m": 25, "n": 25, "t_step_stop": 50, "t_step_save": 50, "parallel_io": "F"},
                # Same 1e-3 as the ibm_burning_grain Example golden this mirrors: a 50-step IBM +
                # finite-rate-chemistry run accumulates cross-compiler roundoff well past the 1e-10
                # the ib branch of compute_tolerance would otherwise pick. Turning the coupling on
                # moves the solution by ~5e-3 relative, so the burn-rate path stays covered at 1e-3.
                override_tol=1e-3,
            )
        )

    ibm_burn_rate_cases()

    def direction_symmetry_tests():
        """3D tests with shock propagating in x and y directions.

        Default 3D tests have the shock along z. These test x and y
        code paths to catch direction-specific bugs in reconstruction,
        Riemann solvers, and gradient calculations.
        """
        for direction in ["x", "y"]:
            others = [d for d in ["x", "y", "z"] if d != direction]
            mods = {
                "m": 24,
                "n": 24,
                "p": 24,
                "x_domain%beg": 0.0e00,
                "x_domain%end": 1.0e00,
                "y_domain%beg": 0.0e00,
                "y_domain%end": 1.0e00,
                "z_domain%beg": 0.0e00,
                "z_domain%end": 1.0e00,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "bc_z%beg": -3,
                "bc_z%end": -3,
            }

            centroids = [0.05, 0.45, 0.9]
            lengths = [0.1, 0.7, 0.2]

            for patchID in range(1, 4):
                mods[f"patch_icpp({patchID})%geometry"] = 9
                mods[f"patch_icpp({patchID})%vel(1)"] = 0.0
                mods[f"patch_icpp({patchID})%vel(2)"] = 0.0
                mods[f"patch_icpp({patchID})%vel(3)"] = 0.0
                mods[f"patch_icpp({patchID})%{direction}_centroid"] = centroids[patchID - 1]
                mods[f"patch_icpp({patchID})%length_{direction}"] = lengths[patchID - 1]
                for od in others:
                    mods[f"patch_icpp({patchID})%{od}_centroid"] = 0.5
                    mods[f"patch_icpp({patchID})%length_{od}"] = 1.0

            stack.push(f"3D Direction Symmetry -> Shock in {direction.upper()}", mods)
            cases.append(define_case_d(stack, "", {}))
            stack.pop()

    direction_symmetry_tests()

    def mpi_consistency_tests():
        """ppn=2 tests for physics sensitive to MPI decomposition.

        Exercises bubble dynamics, viscous flows, and hypoelasticity
        with 2 MPI ranks to catch broadcast/reduction bugs.
        """
        base_3d = {
            "m": 29,
            "n": 29,
            "p": 49,
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 1.0e00,
            "z_domain%beg": 0.0e00,
            "z_domain%end": 1.0e00,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
        }

        base_3d.update(make_3d_box_patches())

        # Bubbles with 2 MPI ranks
        stack.push(
            "MPI Consistency -> 3D -> Bubbles",
            {
                **base_3d,
                "dt": 1e-06,
                "bubbles_euler": "T",
                "nb": 1,
                "polytropic": "T",
                "bubble_model": 2,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%eos": "stiffened_gas",
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%T0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%vd": 0.2404125083932959,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%mu_v": 8.758168074360729e-05,
                "bub_pp%mu_g": 0.00017881922111898042,
                "bub_pp%gam_v": 1.33,
                "bub_pp%gam_g": 1.4,
                "bub_pp%M_v": 18.02,
                "bub_pp%M_g": 28.97,
                "bub_pp%k_v": 0.5583395141263873,
                "bub_pp%k_g": 0.7346421281308791,
                "bub_pp%R_v": 1334.8378710170155,
                "bub_pp%R_g": 830.2995663005393,
                "patch_icpp(1)%alpha_rho(1)": 0.96,
                "patch_icpp(1)%alpha(1)": 4e-02,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 4e-02,
                "patch_icpp(3)%alpha_rho(1)": 0.96,
                "patch_icpp(3)%alpha(1)": 4e-02,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%pres": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}, ppn=2))
        stack.pop()

        # Viscous with 2 MPI ranks
        stack.push(
            "MPI Consistency -> 3D -> Viscous",
            {
                **base_3d,
                "dt": 1e-11,
                "fluid_pp(1)%Re(1)": 0.0001,
                "viscous": "T",
                "patch_icpp(1)%vel(1)": 1.0,
                "patch_icpp(2)%vel(1)": 1.0,
                "patch_icpp(3)%vel(1)": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}, ppn=2))
        stack.pop()

        # Hypoelasticity with 2 MPI ranks
        stack.push(
            "MPI Consistency -> 3D -> Hypoelasticity",
            {
                **base_3d,
                "dt": 1e-06,
                "hypoelasticity": "T",
                "riemann_solver": 1,
                "fd_order": 4,
                "fluid_pp(1)%gamma": 0.3,
                "fluid_pp(1)%eos": "stiffened_gas",
                "fluid_pp(1)%pi_inf": 7.8e05,
                "fluid_pp(1)%G": 1.0e05,
                "patch_icpp(1)%pres": 1.0e06,
                "patch_icpp(1)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(2)%pres": 1.0e05,
                "patch_icpp(2)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(3)%pres": 5.0e05,
                "patch_icpp(3)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(1)%tau_e(1)": 0.0e00,
                "patch_icpp(2)%tau_e(1)": 0.0e00,
                "patch_icpp(3)%tau_e(1)": 0.0e00,
                "patch_icpp(1)%tau_e(2)": 0.0e00,
                "patch_icpp(1)%tau_e(3)": 0.0e00,
                "patch_icpp(2)%tau_e(2)": 0.0e00,
                "patch_icpp(2)%tau_e(3)": 0.0e00,
                "patch_icpp(3)%tau_e(2)": 0.0e00,
                "patch_icpp(3)%tau_e(3)": 0.0e00,
                "patch_icpp(1)%tau_e(4)": 0.0e00,
                "patch_icpp(1)%tau_e(5)": 0.0e00,
                "patch_icpp(1)%tau_e(6)": 0.0e00,
                "patch_icpp(2)%tau_e(4)": 0.0e00,
                "patch_icpp(2)%tau_e(5)": 0.0e00,
                "patch_icpp(2)%tau_e(6)": 0.0e00,
                "patch_icpp(3)%tau_e(4)": 0.0e00,
                "patch_icpp(3)%tau_e(5)": 0.0e00,
                "patch_icpp(3)%tau_e(6)": 0.0e00,
            },
        )
        cases.append(define_case_d(stack, "", {}, ppn=2))
        stack.pop()

    mpi_consistency_tests()

    def restart_roundtrip_tests():
        """Tests that verify save-restart roundtrip fidelity.

        Each test runs a straight simulation, then a restart from the
        midpoint. The restarted output is compared against the straight
        run output to verify restart I/O doesn't introduce drift.
        """
        # 1D restart
        stack.push(
            "Restart Roundtrip -> 1D",
            {
                "m": 299,
                "n": 0,
                "p": 0,
                "x_domain%beg": 0.0e00,
                "x_domain%end": 1.0e00,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "patch_icpp(1)%geometry": 1,
                "patch_icpp(2)%geometry": 1,
                "patch_icpp(3)%geometry": 1,
                "patch_icpp(1)%x_centroid": 0.05,
                "patch_icpp(1)%length_x": 0.1,
                "patch_icpp(2)%x_centroid": 0.45,
                "patch_icpp(2)%length_x": 0.7,
                "patch_icpp(3)%x_centroid": 0.9,
                "patch_icpp(3)%length_x": 0.2,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(3)%vel(1)": 0.0,
            },
        )
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        stack.pop()

        # 3D restart
        base_3d = {
            "m": 24,
            "n": 24,
            "p": 24,
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 1.0e00,
            "z_domain%beg": 0.0e00,
            "z_domain%end": 1.0e00,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
        }
        base_3d.update(make_3d_box_patches())
        stack.push("Restart Roundtrip -> 3D", base_3d)
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        stack.pop()

    restart_roundtrip_tests()

    def kernel_golden_tests():
        """Focused golden-value tests for specific physics kernels.

        Grid stretching in 3D: exercises non-uniform grid spacing in all
        three directions. Stretching interacts with WENO reconstruction
        and gradient calculations in direction-specific ways. Not covered
        by any dynamic test (only via examples at reduced resolution).

        MTHINC on stretched 2D grid: a circular bubble creates diagonal
        interface normals that have components in both x and y. With a
        non-uniform x grid, the reference-space normal weighting
        (Δx_j / (x_cc(j+1)-x_cc(j-1))) differs from 0.5, exercising the
        grid-spacing correction in s_compute_mthinc_normals.
        """
        base_3d = {
            "m": 24,
            "n": 24,
            "p": 24,
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "y_domain%beg": 0.0e00,
            "y_domain%end": 1.0e00,
            "z_domain%beg": 0.0e00,
            "z_domain%end": 1.0e00,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
        }
        base_3d.update(make_3d_box_patches())

        # 3D grid stretching in all directions.
        # The cosh-based stretching expands the domain beyond the original
        # bounds (e.g., [0,1] → ~[0,1.39] with a=2, x_a=0.3, x_b=0.7).
        # Patches must be enlarged to cover the stretched domain, otherwise
        # cells beyond the original bounds are uninitialized (zero density),
        # causing ICFL blowup.
        stack.push(
            "Kernel -> 3D -> Grid Stretching",
            {
                **base_3d,
                "stretch_x": "T",
                "a_x": 2.0,
                "x_a": 0.3,
                "x_b": 0.7,
                "loops_x": 1,
                "stretch_y": "T",
                "a_y": 2.0,
                "y_a": 0.3,
                "y_b": 0.7,
                "loops_y": 1,
                "stretch_z": "T",
                "a_z": 2.0,
                "z_a": 0.3,
                "z_b": 0.7,
                "loops_z": 1,
                # Enlarge x/y coverage for all patches (stretched domain reaches ~1.39)
                "patch_icpp(1)%x_centroid": 0.75,
                "patch_icpp(1)%length_x": 1.5,
                "patch_icpp(1)%y_centroid": 0.75,
                "patch_icpp(1)%length_y": 1.5,
                "patch_icpp(2)%x_centroid": 0.75,
                "patch_icpp(2)%length_x": 1.5,
                "patch_icpp(2)%y_centroid": 0.75,
                "patch_icpp(2)%length_y": 1.5,
                "patch_icpp(3)%x_centroid": 0.75,
                "patch_icpp(3)%length_x": 1.5,
                "patch_icpp(3)%y_centroid": 0.75,
                "patch_icpp(3)%length_y": 1.5,
                # Extend last z-patch to cover stretched z range
                "patch_icpp(3)%z_centroid": 1.15,
                "patch_icpp(3)%length_z": 0.7,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # 2D MTHINC on a stretched (non-uniform) x-grid.
        # A circular bubble creates diagonal interface normals (components in both
        # x and y), so the reference-space weighting Δx_j/(x_cc(j+1)-x_cc(j-1))
        # differs from 0.5 and changes the normalized normal on non-uniform grids.
        # Axis-aligned interfaces would not catch this because the unit normal is
        # the same regardless of the gradient scaling.
        eps = 1e-6
        stack.push(
            "Kernel -> 2D -> MTHINC -> Grid Stretching",
            {
                "m": 24,
                "n": 24,
                "p": 0,
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "num_patches": 2,
                "num_fluids": 2,
                "fluid_pp(2)%gamma": 2.5,
                "fluid_pp(2)%eos": "ideal_gas",
                "fluid_pp(2)%pi_inf": 0.0,
                # Patch 1: fluid 1 background rectangle; length covers stretched extent (~1.39).
                # vel(1)=0.5 provides advection so MTHINC reconstruction affects the solution.
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.75,
                "patch_icpp(1)%length_x": 1.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.5,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0 - eps,
                "patch_icpp(1)%alpha(1)": 1.0 - eps,
                "patch_icpp(1)%alpha_rho(2)": eps,
                "patch_icpp(1)%alpha(2)": eps,
                # Patch 2: fluid 2 circular bubble centered in the stretched region;
                # alter_patch(1)=T needed so patch 2 can overwrite patch 1 cells
                "patch_icpp(2)%geometry": 2,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(2)%radius": 0.2,
                "patch_icpp(2)%vel(1)": 0.5,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(2)%alpha_rho(1)": eps,
                "patch_icpp(2)%alpha(1)": eps,
                "patch_icpp(2)%alpha_rho(2)": 1.0 - eps,
                "patch_icpp(2)%alpha(2)": 1.0 - eps,
                # MTHINC
                "int_comp": 2,
                # x-stretching creates non-uniform cells at the bubble interface
                "stretch_x": "T",
                "a_x": 2.0,
                "x_a": 0.3,
                "x_b": 0.7,
                "loops_x": 1,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

    kernel_golden_tests()

    add_convergence_cases(cases)

    # Sanity Check 1
    if stack.size() != 0:
        raise common.MFCException("list_cases: stack isn't fully pop'ed")

    # Sanity Check 2
    uuids = [case.get_uuid() for case in cases]
    l1, l2 = len(uuids), len(set(uuids))
    if l1 != l2:
        raise common.MFCException(f"list_cases: uuids aren't unique ({l1} cases but {l2} unique uuids)")

    # Tag the always-run canary smoke set (see _CANARY_TRACES). Validate first so a renamed
    # or removed trace is a loud error, not a silently dropped canary.
    missing = _CANARY_TRACES - {case.trace for case in cases}
    if missing:
        raise common.MFCException(f"list_cases: canary trace(s) not found (renamed/removed?): {sorted(missing)}")
    for case in cases:
        if case.trace in _CANARY_TRACES:
            case.canary = True

    return cases
