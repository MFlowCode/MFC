import itertools
import os
import typing

from mfc import common

from ..state import ARG
from .case import CaseGeneratorStack, Nt, TestCaseBuilder, define_case_d, define_case_f, define_convergence_case
from .convergence import ConvergenceSpec, run_dt_sweep, run_h_sweep, run_sod_l1

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

    def alter_igr(amr_variant=False):
        stack.push("IGR", {"igr": "T", "alf_factor": 10, "num_igr_iters": 10, "elliptic_smoothing": "T", "elliptic_smoothing_iters": 10, "num_igr_warm_start_iters": 10})

        for order in [3, 5]:
            stack.push(f"igr_order={order}", {"igr_order": order})

            cases.append(define_case_d(stack, "Jacobi", {"igr_iter_solver": 1}))
            if order == 5:
                cases.append(define_case_d(stack, "Gauss Seidel", {"igr_iter_solver": 2}))
            # AMR (stage-1 restriction-only coupling): the fine block runs its own fixed-iteration
            # sigma solve, seeded and Dirichlet-bounded by the converged coarse sigma (frozen
            # ghost ring, per-iteration BC populate skipped); validated free-stream-exact, with
            # the AMR-vs-reference error at resolution scale (rho 1.3e-4 rel-L2) and a
            # truncation-level transverse seam artifact from the coarse/fine sigma jump
            if order == 3 and amr_variant:
                stack.push("AMR", {"amr": "T", "amr_block_beg(1)": 14, "amr_block_beg(2)": 12, "amr_block_end(1)": 33, "amr_block_end(2)": 27, "amr_regrid_int": 0, "igr_iter_solver": 1})
                cases.append(define_case_d(stack, "", {}))
                cases.append(define_case_d(stack, "dynamic regrid", {"amr_regrid_int": 5, "amr_tag_eps": 1.0e-2, "amr_buf": 2}))
                stack.pop()

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

    def alter_riemann_solvers(num_fluids):
        for riemann_solver in [1, 5, 2]:
            stack.push(f"riemann_solver={riemann_solver}", {"riemann_solver": riemann_solver})

            cases.append(define_case_d(stack, "mixture_err", {"mixture_err": "T"}))

            if riemann_solver in (1, 2):
                cases.append(define_case_d(stack, "avg_state=1", {"avg_state": 1}))
                cases.append(define_case_d(stack, "wave_speeds=2", {"wave_speeds": 2}))

                if riemann_solver == 2:
                    cases.append(define_case_d(stack, "model_eqns=3", {"model_eqns": 3}))

                if num_fluids == 2:
                    if riemann_solver == 2:
                        cases.append(define_case_d(stack, "alt_soundspeed", {"alt_soundspeed": "T"}))

                    cases.append(define_case_d(stack, "mpp_lim", {"mpp_lim": "T"}))

            stack.pop()

    def alter_low_Mach_correction():
        stack.push("", {"fluid_pp(1)%gamma": 0.16, "fluid_pp(1)%pi_inf": 3515.0, "dt": 1e-7})

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
            alter_low_Mach_correction()
            alter_ib(dimInfo)
            if len(dimInfo[0]) > 1:
                # AMR variants only on the targeted 2D 1-fluid inviscid base (one static + one
                # dynamic-regrid golden; the block indices are 2D)
                alter_igr(amr_variant=(len(dimInfo[0]) == 2 and num_fluids == 1))

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

            stack.push("", {"fluid_pp(1)%pi_inf": 351.5})
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

            reflective_params = {"bc_x%beg": -2, "bc_x%end": -2, "bc_y%beg": -2, "bc_y%end": -2}
            if len(dimInfo[0]) == 3:
                reflective_params.update({"bc_z%beg": -2, "bc_z%end": -2})

            if num_fluids == 1:
                cases.append(define_case_d(stack, "cont_damage", {"cont_damage": "T", "tau_star": 0.0, "cont_damage_s": 2.0, "alpha_bar": 1e-4}))
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
                            "fluid_pp(1)%pi_inf": 1.7409e09,
                            "fluid_pp(1)%cv": 1816,
                            "fluid_pp(1)%qv": -1167000,
                            "fluid_pp(1)%qvp": 0.0,
                            "fluid_pp(2)%gamma": 2.3266,
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
                            "fluid_pp(1)%pi_inf": 3515.0,
                            "fluid_pp(2)%gamma": 2.5,
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
                        stack.push("", {"acoustic(1)%support": 2})
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
                    # AMR with the bubble cloud EXCLUDED from the block (2D two-way, fixed dt):
                    # the block sits clear of the bubble (0.5, 0.5) and the acoustic support
                    # slab; EL alphas sum to the local liquid fraction, so they prolong WITHOUT
                    # the sum-to-one closure (which would corrupt the EL state - caught by the
                    # free-stream battery), and a per-stage guard keeps the cloud out of blocks
                    # KNOWN CI QUIRK: these two goldens fail with a post-detected NaN on the
                    # nvhpc 24.1/24.3 compat lanes ONLY (non-gating, continue-on-error; 24.5+
                    # green). Exhaustively unreproducible off GitHub's runners: the exact
                    # failing stack - NVHPC 24.3, -tp=px -Kieee, HPC-X MPI, the CI docker
                    # image itself (via apptainer) - passes on Phoenix, as do zen2/native
                    # builds. Suspected runner-hardware/virtualization interaction with old
                    # nvfortran codegen; revisit only if it starts failing on 24.5+.
                    if len(dimInfo[0]) == 2 and adap_dt == "F" and couplingMethod == 2:
                        stack.push("AMR", {"amr": "T", "amr_block_beg(1)": 7, "amr_block_end(1)": 13, "amr_block_beg(2)": 7, "amr_block_end(2)": 12, "amr_regrid_int": 0})
                        cases.append(define_case_d(stack, "", {}))
                        cases.append(define_case_d(stack, "dynamic regrid", {"amr_regrid_int": 5, "amr_tag_eps": 1.0e-3, "amr_buf": 2}))
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
                "2D_ibm_multiphase",
                "2D_acoustic_broadband",
                "1D_inert_shocktube",
                "1D_reactive_shocktube",
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
                "2D_zero_circ_vortex_analytical",
                "3D_TaylorGreenVortex_analytical",
                "3D_IGR_TaylorGreenVortex_nvidia",
                "2D_backward_facing_step",
                "2D_forward_facing_step",
                "1D_convergence",
                "3D_IGR_33jet",
                "1D_multispecies_diffusion",
                "2D_ibm_stl_MFCCharacter",
                "1D_qbmm",
                "2D_Thermal_Flatplate",  # formatted I/O field overflow on gfortran 12
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
                # Chaotic stiff collisional case: its step-50 fields diverge across
                # compilers/platforms, so the golden is not a portable regression target
                # (already GPU-marginal). Also exercises the upstream mibm central-diff IB
                # drag OOB (interior-only fd_coeff at boundary bodies) tracked in
                # MFlowCode/MFC#1633; fixed here in m_viscous but the golden stays non-portable.
                "3D_mibm_periodic_collision",
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

        # Chemistry AMR: a reactive H2/O2/AR shocktube with a static 2:1 fine block over the reaction zone.
        # Exercises the species sum/positivity prolongation closure, the per-block reaction on the fine level,
        # and species reflux. The ppn=2 variant places the block (coarse 16..31) across the rank seam, exercising
        # the fine halo exchange plus the temperature-ghost exchange the fine cons->prim Newton guess needs at
        # the seam (without it the widened conversion diverges to NaN).
        amr_chem_mods = {
            "m": 48,
            "t_step_start": 0,
            "t_step_stop": 20,
            "t_step_save": 20,
            "amr": "T",
            "amr_block_beg(1)": 16,
            "amr_block_end(1)": 31,
            "amr_regrid_int": 0,
        }
        for ppn, label in ((1, "Reactive Shocktube AMR"), (2, "Reactive Shocktube AMR -> 2 MPI Ranks")):
            cases.append(
                define_case_f(
                    f"1D -> Chemistry -> {label}",
                    "examples/1D_reactive_shocktube/case.py",
                    [],
                    ppn=ppn,
                    mods=amr_chem_mods,
                    override_tol=10 ** (-8),
                )
            )

        # Chemistry diffusion AMR: reactions + species mass diffusion with the same static block over the
        # reaction/diffusion zone. Exercises the flux_src reflux of the species (and energy) diffusion
        # fluxes into the coarse/fine registers - without it element mass/energy leak at the block boundary.
        cases.append(
            define_case_f(
                "1D -> Chemistry -> Reactive Shocktube AMR -> Species Diffusion",
                "examples/1D_reactive_shocktube/case.py",
                [],
                ppn=1,
                mods={**amr_chem_mods, "chem_params%diffusion": "T"},
                override_tol=10 ** (-8),
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

    foreach_example()

    chemistry_cases()

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

        # 3D active_box: localized central blast with a uniform ambient exterior so
        # the active-box initialization detects a strict subset of the domain
        # (corner cell (0,0,0) is ambient; blast occupies the central ~12 cells/dim).
        # The box grows by buff_size=4 cells/side each step, so on a 48^3 grid with
        # init box ~[14:33] it is still a strict subset after t_step_stop=3 grows
        # (-> ~[2:45]); the save at step 3 therefore pins a genuinely bounded state.
        # Requires single rank and the model_eqns=2 / WENO5 / HLLC / direct /
        # RK3 configuration that gates the optimization (all BASE_CFG defaults).
        stack.push(
            "Kernel -> 3D -> active_box",
            {
                "m": 47,
                "n": 47,
                "p": 47,
                "dt": 0.005,
                "t_step_stop": 3,
                "t_step_save": 3,
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
                "num_patches": 2,
                "num_fluids": 1,
                # Patch 1: uniform ambient that fills the whole domain.
                # Corner cell (0,0,0) samples this state as ab_ambient.
                "patch_icpp(1)%geometry": 9,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%z_centroid": 0.5,
                "patch_icpp(1)%length_z": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%vel(3)": 0.0,
                "patch_icpp(1)%alpha_rho(1)": 0.125,
                "patch_icpp(1)%pres": 0.1,
                "patch_icpp(1)%alpha(1)": 1.0,
                # Patch 2: high-pressure blast at center (~12 cells/dim out of 48).
                # alter_patch(1)=T overwrites the ambient in the blast region only.
                "patch_icpp(2)%geometry": 9,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.25,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(2)%length_y": 0.25,
                "patch_icpp(2)%z_centroid": 0.5,
                "patch_icpp(2)%length_z": 0.25,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(2)%vel(3)": 0.0,
                "patch_icpp(2)%alpha_rho(1)": 1.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(2)%alpha(1)": 1.0,
                "active_box": "T",
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # active_box + AMR (np=1 by active_box's own MPI gate): a 2D central blast with the
        # block strictly inside the initial active window. t_step_stop=10 keeps the window
        # partial for the whole run (it grows buff_size cells/step and self-disables at full
        # domain), so every step exercises the windowed coarse advance around a live fine
        # block. Validated: ab+AMR vs plain AMR agree to 9.8e-15 (the active_box round-off
        # spec) over 200 steps incl. the self-disable transition; the containment abort and
        # the regrid window-clamp are manually negative-tested.
        stack.push(
            "Kernel -> 2D -> active_box -> AMR",
            {
                "m": 127,
                "n": 127,
                "p": 0,
                "dt": 5.0e-5,
                "t_step_stop": 10,
                "t_step_save": 10,
                "num_patches": 2,
                "mixture_err": "F",
                "mapped_weno": "T",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%geometry": 2,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(2)%radius": 0.08,
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(2)%pres": 10.0,
                "patch_icpp(2)%alpha_rho(1)": 2.0,
                "patch_icpp(2)%alpha(1)": 1.0,
                "patch_icpp(2)%alter_patch(1)": "T",
                "active_box": "T",
                "amr": "T",
                "amr_block_beg(1)": 52,
                "amr_block_beg(2)": 52,
                "amr_block_end(1)": 75,
                "amr_block_end(2)": 75,
                "amr_regrid_int": 0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        cases.append(define_case_d(stack, "dynamic regrid", {"amr_regrid_int": 5, "amr_tag_eps": 0.02, "amr_buf": 3}))
        stack.pop()

    kernel_golden_tests()

    def amr_golden_tests():
        """Golden tests for the block-structured AMR module.

        Three minimal 1D cases built on the BASE_CFG Sod IC:
        (a) static block  – amr_regrid_int=0, no regrid, cheapest sanity check
        (b) dynamic regrid – amr_regrid_int=2, exercises the tag/rebuild path
        (c) subcycling    – amr_subcycle=T, exercises the dt/2 two-substep path
        (d) two-fluid     – num_fluids=2 + mpp_lim, regrid + subcycle: exercises the
            sum-preserving alpha prolongation and per-fluid reflux (SP9a)

        Grid: m=63 (indices 0..63); fine block beg=16, end=47 so that
        2*(47-16+1)-1 = 63 <= 63 satisfies the extent guard.
        """
        # Common 1D domain + Sod IC setup shared by all three AMR cases
        amr_1d_base = {
            "m": 63,
            "n": 0,
            "p": 0,
            "dt": 5.0e-4,
            "t_step_stop": 6,
            "t_step_save": 6,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # 1D geometry for the three BASE_CFG patches
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.05,
            "patch_icpp(1)%length_x": 0.1,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.45,
            "patch_icpp(2)%length_x": 0.7,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(3)%geometry": 1,
            "patch_icpp(3)%x_centroid": 0.9,
            "patch_icpp(3)%length_x": 0.2,
            "patch_icpp(3)%vel(1)": 0.0,
            # AMR: 2:1 fine block spanning coarse indices 16..47
            "amr": "T",
            "amr_block_beg(1)": 16,
            "amr_block_end(1)": 47,
        }

        # (a) static block
        stack.push("AMR -> 1D -> static block", {**amr_1d_base, "amr_regrid_int": 0})
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        stack.pop()

        # (b) dynamic regrid
        stack.push(
            "AMR -> 1D -> dynamic regrid",
            {
                **amr_1d_base,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
            },
        )
        # restart_check on the REGRIDDED layout: unlike the static block, a regridded block set
        # cannot be reconstructed from the ICs, so the roundtrip proves the restart file itself
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        # hybrid sensors on the fine level. eps 0.5 exceeds the Sod shock's Jameson phi (~0.3),
        # putting the shock cells under the sensor's control (central weights/flux where
        # phi < eps): a dead sensor moves the answer by ~5e-4. At physical eps (1e-2) the
        # combination is bitwise-identical to plain AMR by design (only constant/eps-dominated
        # cells go central, and those reconstruct identically under any convex weights).
        # override_tol 5e-5, sized by decision theory rather than drift-chasing: the harness
        # reports only the FIRST variable past tolerance, so each tighter bound just reveals
        # the next quantile of the same deterministic deviation field (intel-CPU and Cray-acc
        # produce IDENTICAL sensor-on results that differ from the gnu-generated golden by up
        # to ~6e-6 - compiler FP-contraction behavior, not noise). 5e-5 sits an order below
        # the ~5e-4 dead-sensor signal these goldens exist to catch and 8x above the worst
        # observed legitimate drift; if a platform ever drifts past 5e-5 the golden design
        # itself (static cross-compiler reference) stops being able to discriminate
        cases.append(define_case_d(stack, "hybrid_weno sensor", {"hybrid_weno": "T", "hybrid_weno_eps": 0.5}, override_tol=5.0e-5))
        cases.append(define_case_d(stack, "hybrid_riemann sensor", {"hybrid_riemann": "T", "hybrid_weno_eps": 0.5, "hybrid_smooth_flux": 2}, override_tol=5.0e-5))
        # 2 MPI ranks + parallel_io: the ONLY test that executes the MPI-IO AMR restart write/read
        # (EXSCAN offset arithmetic, per-rank-extents validation) and multi-rank dynamic regrid
        # (coarse-halo exchange before tagging, fine seam halo) - a rank-seam or restart-offset bug
        # is a silent wrong answer everywhere else in the suite
        stack.push("2 MPI Ranks", {"parallel_io": "T"})
        cases.append(define_case_d(stack, "", {}, ppn=2, restart_check=True, honor_io_keys=True))
        stack.pop()
        stack.pop()

        # (b') stretched grid + dynamic regrid: the ONLY test where the coarse grid is
        # nonuniform - exercises the exact parent-bisection ghost-shell coordinates and the
        # per-swap WENO coefficient recompute (amr_weno_coef_recompute armed at init).
        # stretch_x expands the domain beyond [0,1], so the end patches are widened to keep
        # the expanded cells covered; the fine block 16..47 straddles the uniform core
        # [x_a, x_b] so its ghost shells sit on nonuniform parents on both sides.
        stack.push(
            "AMR -> 1D -> stretched grid -> dynamic regrid",
            {
                **amr_1d_base,
                "stretch_x": "T",
                "a_x": 2.0,
                "x_a": 0.4,
                "x_b": 0.6,
                "loops_x": 1,
                "patch_icpp(1)%x_centroid": -1.95,
                "patch_icpp(1)%length_x": 4.1,
                "patch_icpp(3)%x_centroid": 2.9,
                "patch_icpp(3)%length_x": 4.2,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
            },
        )
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        # 2 MPI ranks: the ONLY case exercising the '- start_idx(d)' rank-offset terms of the
        # parent-bisection ghost formula on a grid where a wrong offset changes coordinates
        # (uniform spacing makes any parent index give the same value; the block spans the seam)
        cases.append(define_case_d(stack, "2 MPI Ranks", {}, ppn=2))
        stack.pop()

        # (b'') stretched in y, 2D: the y-direction parent-bisection formula (copy-pasted per
        # dim) reduces to the uniform formula on every other 2D/3D golden; here the block's
        # lower/upper y ghost shells sit on nonuniform parents (uniform core y in [0.35, 0.6])
        stack.push(
            "AMR -> 2D -> stretched grid y -> dynamic regrid",
            {
                **amr_1d_base,
                "n": 39,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "stretch_y": "T",
                "a_y": 2.0,
                "y_a": 0.35,
                "y_b": 0.6,
                "loops_y": 1,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(2)%geometry": 3,
                "patch_icpp(3)%geometry": 3,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(2)%y_centroid": 0.5,
                "patch_icpp(3)%y_centroid": 0.5,
                "patch_icpp(1)%length_y": 10.0,
                "patch_icpp(2)%length_y": 10.0,
                "patch_icpp(3)%length_y": 10.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(3)%vel(2)": 0.0,
                "amr_block_beg(2)": 10,
                "amr_block_end(2)": 25,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (c) subcycling
        stack.push(
            "AMR -> 1D -> subcycle",
            {
                **amr_1d_base,
                "amr_regrid_int": 0,
                "amr_subcycle": "T",
            },
        )
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        # body forces ride the subcycle: accel is evaluated at the coarse-step-frozen mytime on
        # fine substeps - the same per-step time freezing the coarse RK3 stages already apply, so
        # coarse and fine see one consistent forcing. Oscillatory + gravity per suite convention.
        stack.push("bodyforces", {"bf_x": "T", "k_x": 1, "w_x": 1, "p_x": 1, "g_x": 10})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (d) 3D static block — z-dimension coverage. 26^3 base grid (the
        # checker's WENO5 floor is 26 cells per axis) with the Sod-like slabs
        # stacked along z; 12^3 fine block at coarse indices 6..17 per axis
        # (fine extent 23 <= 25; >= buff_size=4 inside the domain).
        amr_3d_base = {
            "m": 25,
            "n": 25,
            "p": 25,
            "dt": 2.0e-3,
            "t_step_stop": 6,
            "t_step_save": 6,
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
            # 3D geometry: the three BASE_CFG states as full-x/y slabs along z
            **{
                f"patch_icpp({i})%{key}": val
                for i in (1, 2, 3)
                for key, val in (
                    ("geometry", 9),
                    ("x_centroid", 0.5),
                    ("length_x", 1.0),
                    ("y_centroid", 0.5),
                    ("length_y", 1.0),
                    ("vel(1)", 0.0),
                    ("vel(2)", 0.0),
                    ("vel(3)", 0.0),
                )
            },
            "patch_icpp(1)%z_centroid": 0.05,
            "patch_icpp(1)%length_z": 0.1,
            "patch_icpp(2)%z_centroid": 0.45,
            "patch_icpp(2)%length_z": 0.7,
            "patch_icpp(3)%z_centroid": 0.9,
            "patch_icpp(3)%length_z": 0.2,
            # AMR: 2:1 fine block spanning coarse indices 6..17 per axis
            "amr": "T",
            "amr_block_beg(1)": 6,
            "amr_block_beg(2)": 6,
            "amr_block_beg(3)": 6,
            "amr_block_end(1)": 17,
            "amr_block_end(2)": 17,
            "amr_block_end(3)": 17,
        }

        stack.push("AMR -> 3D -> static block", {**amr_3d_base, "amr_regrid_int": 0})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (d) two-fluid: material interface (density ratio 10) at x=0.5, inside the initial
        # block (cells 16..47); uniform p and u advect it under regrid + subcycle
        eps_a = 1.0e-6
        stack.push(
            "AMR -> 1D -> two-fluid",
            {
                **amr_1d_base,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "num_fluids": 2,
                "mpp_lim": "T",
                "fluid_pp(2)%gamma": 1.0e00 / (1.6e00 - 1.0e00),
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%cv": 0.0,
                "fluid_pp(2)%qv": 0.0,
                "fluid_pp(2)%qvp": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%pres": 1.0,
                "patch_icpp(1)%vel(1)": 0.5,
                "patch_icpp(2)%vel(1)": 0.5,
                "patch_icpp(3)%vel(1)": 0.5,
                "patch_icpp(2)%x_centroid": 0.3,
                "patch_icpp(2)%length_x": 0.4,
                "patch_icpp(3)%x_centroid": 0.75,
                "patch_icpp(3)%length_x": 0.5,
                "patch_icpp(1)%alpha_rho(1)": (1.0 - eps_a) * 1.0,
                "patch_icpp(1)%alpha_rho(2)": eps_a * 10.0,
                "patch_icpp(1)%alpha(1)": 1.0 - eps_a,
                "patch_icpp(1)%alpha(2)": eps_a,
                "patch_icpp(2)%alpha_rho(1)": (1.0 - eps_a) * 1.0,
                "patch_icpp(2)%alpha_rho(2)": eps_a * 10.0,
                "patch_icpp(2)%alpha(1)": 1.0 - eps_a,
                "patch_icpp(2)%alpha(2)": eps_a,
                "patch_icpp(3)%alpha_rho(1)": eps_a * 1.0,
                "patch_icpp(3)%alpha_rho(2)": (1.0 - eps_a) * 10.0,
                "patch_icpp(3)%alpha(1)": eps_a,
                "patch_icpp(3)%alpha(2)": 1.0 - eps_a,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        # THINC interface compression on the advecting interface: the sharpener reads the live
        # grid arrays (swapped per block) and its scratch spans idwbuff, so it is AMR-correct by
        # construction - this golden protects the reachable WENO+int_comp combo under regrid+subcycle
        stack.push("thinc", {"int_comp": 1})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        # 6-equation model on the same interface advection: the internal-energy equations ride the
        # generic conservative prolong/restrict/reflux, and the per-stage pressure relaxation
        # (cell-local) runs on the fine block mirroring the coarse stage order
        stack.push("6eq", {"model_eqns": 3})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (d2) hypoelasticity: the suite's 1D hypoelastic shock config (stiff water EOS + shear
        # modulus G) on a static 2:1 fine block over the wave region. Stress components prolong
        # via the generic conservative-linear path; the fine swap recomputes the spacing-dependent
        # FD coefficients the stress source uses (coarse coefficients would halve fine gradients).
        stack.push(
            "AMR -> 1D -> hypoelastic static block",
            {
                **amr_1d_base,
                "amr_regrid_int": 0,
                # stiff water EOS: c ~ 83; dt=5e-5 keeps the 2:1 fine block at CFL ~ 0.5
                "dt": 5.0e-5,
                "hypoelasticity": "T",
                "riemann_solver": 1,
                "fd_order": 4,
                "fluid_pp(1)%gamma": 0.3,
                "fluid_pp(1)%pi_inf": 7.8e05,
                "fluid_pp(1)%G": 1.0e05,
                "patch_icpp(1)%pres": 1.0e06,
                "patch_icpp(1)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(2)%pres": 1.0e05,
                "patch_icpp(2)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(3)%pres": 5.0e05,
                "patch_icpp(3)%alpha_rho(1)": 1000.0e00,
                "patch_icpp(1)%tau_e(1)": 0.0e-00,
                "patch_icpp(2)%tau_e(1)": 0.0e-00,
                "patch_icpp(3)%tau_e(1)": 0.0e-00,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        # continuum damage rides hypoelasticity: its source is cell-local (pointwise in the local
        # stress), so it is AMR-correct by construction - this golden protects the reachable combo
        stack.push("cont_damage", {"cont_damage": "T", "tau_star": 0.0, "cont_damage_s": 2.0, "alpha_bar": 1e-4})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (d3) acoustic source: a sine pulse emitted on the coarse grid (support 1 at x=0.1) with a
        # static fine block downstream (x in [0.44, 0.75]); the wave crosses the coarse/fine boundary
        # into the block during the run. The source acts on the coarse grid only - a support/block
        # overlap aborts at startup - and the fine advance skips it (coarse-index spatials).
        stack.push(
            "AMR -> 1D -> acoustic static block",
            {
                **amr_1d_base,
                "amr_regrid_int": 0,
                "amr_block_beg(1)": 28,
                "dt": 2.0e-3,
                "t_step_stop": 200,
                "t_step_save": 200,
                # uniform quiescent background (overrides the Sod-like patch states)
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(2)%alpha_rho(1)": 1.0,
                "patch_icpp(3)%alpha_rho(1)": 1.0,
                "acoustic_source": "T",
                "acoustic(1)%support": 1,
                "acoustic(1)%loc(1)": 0.1,
                "acoustic(1)%pulse": 1,
                "acoustic(1)%wavelength": 0.2,
            },
        )
        # override_tol 1e-8: intel -O3 drifts this case ~3e-10 rel past the 1e-12 default
        cases.append(define_case_d(stack, "", {}, override_tol=1.0e-8))
        # dynamic regrid chasing the emitted wave: the tagger fires on the travelling pulse, and
        # the regrid keeps its boxes clear of the source support (tags suppressed over it,
        # candidate boxes clipped) - the source region stays coarse while blocks track the wave
        stack.push("dynamic regrid", {"amr_regrid_int": 2, "amr_tag_eps": 0.01, "amr_buf": 2})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (e) viscous (SP11): single-fluid Sod with physical viscosity (Re=100), regrid + subcycle.
        # Exercises the viscous flux-register reflux (flux_src_n momentum/energy captured into the
        # same registers as the advective flux_n) so the c/f boundary sees matched total fluxes.
        stack.push(
            "AMR -> 1D -> viscous",
            {
                **amr_1d_base,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "viscous": "T",
                "weno_Re_flux": "T",
                "weno_avg": "T",
                "fluid_pp(1)%Re(1)": 100.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (f) multi-block (SP12a): three constant states with density interfaces at x=0.25 (cell 16)
        # and x=0.75 (cell 48) -- two features ~32 coarse cells apart (> buff_size + 2*amr_buf), so the
        # Berger-Rigoutsos clustering forms TWO blocks (one per interface) rather than one bounding box
        # spanning both plus the empty middle. Uniform pressure so the interfaces stay put; regrid on.
        # Exercises the per-slot advance + the single coarse-RHS flux-register capture filling both
        # blocks' registers (the whole SP12a capability).
        stack.push(
            "AMR -> 1D -> multi-block",
            {
                **amr_1d_base,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_max_blocks": 4,
                "patch_icpp(1)%x_centroid": 0.125,
                "patch_icpp(1)%length_x": 0.25,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.5,
                "patch_icpp(2)%alpha_rho(1)": 0.2,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%x_centroid": 0.875,
                "patch_icpp(3)%length_x": 0.25,
                "patch_icpp(3)%alpha_rho(1)": 1.0,
                "patch_icpp(3)%pres": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (g) Euler-Euler bubbles (SP13): monodisperse (nb=1) polytropic bubbles in a uniform
        # void-fraction liquid, with a left pressure slab (pres=2) launching a wave through the
        # block. Exercises the realizability-preserving radius-moment prolongation (floor at the
        # interior/regrid/ghost fills) and the flux-based bubble-moment reflux, under regrid +
        # subcycle. pb/mv are inert stubs (non-qbmm polytropic), so only the moment path runs.
        stack.push(
            "AMR -> 1D -> bubbles",
            {
                **amr_1d_base,
                "dt": 5.0e-5,  # stiff bubble liquid (pi_inf=3515): keeps coarse+subcycled fine ICFL < 1
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "bubbles_euler": "T",
                "bubble_model": 2,
                "polytropic": "T",
                "nb": 1,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%T0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%vd": 0.2404125083932959,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%gam_g": 1.4,
                "patch_icpp(1)%alpha_rho(1)": 0.96,
                "patch_icpp(1)%alpha(1)": 4e-02,
                "patch_icpp(1)%pres": 2.0,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 4e-02,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%alpha_rho(1)": 0.96,
                "patch_icpp(3)%alpha(1)": 4e-02,
                "patch_icpp(3)%pres": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (h) phase change (SP15): two-fluid (liquid water + vapor) pT-equilibrium relaxation
        # (relax_model=5) with a pressure/temperature-disequilibrium interface inside the block,
        # regrid + subcycle. Exercises the per-block relaxation on the fine solution BEFORE
        # restriction: s_amr_relax_fine equilibrates the fine cells (cell-local, mass/energy-
        # conserving) so the restricted coarse average is relax-consistent. Small dt keeps the
        # stiff water EOS (pi_inf~1.7e9) CFL-stable over the six captured steps.
        stack.push(
            "AMR -> 1D -> phase change",
            {
                **amr_1d_base,
                "dt": 1.0e-6,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "num_fluids": 2,
                "mpp_lim": "T",
                "relax": "T",
                "relax_model": 5,
                "palpha_eps": 1.0e-2,
                "ptgalpha_eps": 1.0e-2,
                "fluid_pp(1)%gamma": 0.7409,
                "fluid_pp(1)%pi_inf": 1.7409e09,
                "fluid_pp(1)%cv": 1816.0,
                "fluid_pp(1)%qv": -1167000.0,
                "fluid_pp(1)%qvp": 0.0,
                "fluid_pp(2)%gamma": 2.3266,
                "fluid_pp(2)%pi_inf": 0.0e00,
                "fluid_pp(2)%cv": 1040.0,
                "fluid_pp(2)%qv": 2030000.0,
                "fluid_pp(2)%qvp": -23400.0,
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
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (i) CROSS-FEATURE: viscous + two-fluid + multi-block + subcycle (SP11+SP9a+SP12a+SP6).
        # Two material interfaces (fluid1|fluid2|fluid1, total-density ratio 10) at x=0.25 (cell 16) and
        # x=0.75 (cell 48): the density-gradient tagger clusters TWO blocks (~32 coarse cells apart >
        # buff_size + 2*amr_buf). A velocity step across each interface drives a real viscous stress, so
        # both blocks' registers reflux the viscous momentum/energy AND the per-fluid species fluxes at
        # once, under regrid + subcycle. Conservation defect stays ~1e-13.
        stack.push(
            "AMR -> 1D -> viscous multifluid multiblock",
            {
                **amr_1d_base,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "amr_max_blocks": 4,
                "num_fluids": 2,
                "mpp_lim": "T",
                "viscous": "T",
                "weno_Re_flux": "T",
                "weno_avg": "T",
                "fluid_pp(1)%Re(1)": 100.0,
                "fluid_pp(2)%gamma": 1.0e00 / (1.6e00 - 1.0e00),
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%cv": 0.0,
                "fluid_pp(2)%qv": 0.0,
                "fluid_pp(2)%qvp": 0.0,
                "fluid_pp(2)%Re(1)": 100.0,
                # fluid1 | fluid2 | fluid1  => total-density interfaces at x=0.25 and x=0.75
                "patch_icpp(1)%x_centroid": 0.125,
                "patch_icpp(1)%length_x": 0.25,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%vel(1)": 0.5,
                "patch_icpp(1)%alpha_rho(1)": (1.0 - eps_a) * 1.0,
                "patch_icpp(1)%alpha_rho(2)": eps_a * 10.0,
                "patch_icpp(1)%alpha(1)": 1.0 - eps_a,
                "patch_icpp(1)%alpha(2)": eps_a,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.5,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(2)%vel(1)": 0.3,
                "patch_icpp(2)%alpha_rho(1)": eps_a * 1.0,
                "patch_icpp(2)%alpha_rho(2)": (1.0 - eps_a) * 10.0,
                "patch_icpp(2)%alpha(1)": eps_a,
                "patch_icpp(2)%alpha(2)": 1.0 - eps_a,
                "patch_icpp(3)%x_centroid": 0.875,
                "patch_icpp(3)%length_x": 0.25,
                "patch_icpp(3)%pres": 1.0,
                "patch_icpp(3)%vel(1)": 0.5,
                "patch_icpp(3)%alpha_rho(1)": (1.0 - eps_a) * 1.0,
                "patch_icpp(3)%alpha_rho(2)": eps_a * 10.0,
                "patch_icpp(3)%alpha(1)": 1.0 - eps_a,
                "patch_icpp(3)%alpha(2)": eps_a,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (j) CROSS-FEATURE: Euler-Euler bubbles + multi-block + subcycle + regrid (SP13+SP12a+SP6+SP5).
        # Polytropic bubbly liquid with a pressure slab at EACH end (pres=2 at the left/right quarters,
        # pres=1 in the middle): the two inward-running compression fronts create two separated density
        # features, so the tagger clusters up to two blocks. Exercises the realizability-preserving
        # radius-moment prolongation and the bubble-moment reflux across both blocks' registers, with the
        # moments read from the subcycle time-lerp ghosts. dt=5e-5 keeps coarse+subcycled fine ICFL < 1.
        stack.push(
            "AMR -> 1D -> bubbles multiblock",
            {
                **amr_1d_base,
                "dt": 5.0e-5,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "amr_max_blocks": 4,
                "bubbles_euler": "T",
                "bubble_model": 2,
                "polytropic": "T",
                "nb": 1,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%T0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%vd": 0.2404125083932959,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%gam_g": 1.4,
                "patch_icpp(1)%x_centroid": 0.125,
                "patch_icpp(1)%length_x": 0.25,
                "patch_icpp(1)%alpha_rho(1)": 0.96,
                "patch_icpp(1)%alpha(1)": 4e-02,
                "patch_icpp(1)%pres": 2.0,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%length_x": 0.5,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 4e-02,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%x_centroid": 0.875,
                "patch_icpp(3)%length_x": 0.25,
                "patch_icpp(3)%alpha_rho(1)": 0.96,
                "patch_icpp(3)%alpha(1)": 4e-02,
                "patch_icpp(3)%pres": 2.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (k) Euler-Euler bubbles, NON-POLYTROPIC + multi-bin (SP18): nb=3 R0 bins with
        # polytropic=F, so each bin carries FOUR conserved moments (radius nR, velocity nV,
        # partial pressure npb, vapor mass nmv) — 12 bubble moments total. A left pressure slab
        # (pres=2) drives a wave through the block under regrid + subcycle. Exercises the extended
        # realizability floor (positive moments nR/npb/nmv floored, signed nV free) across all bins
        # at the interior/regrid/ghost prolongation, and the flux-based reflux of the full moment
        # set. pb/mv live entirely in q_cons (non-qbmm), so no side-state advance is needed on the
        # fine level. dt=5e-5 keeps coarse+subcycled fine ICFL < 1.
        stack.push(
            "AMR -> 1D -> bubbles nonpolytropic",
            {
                **amr_1d_base,
                "dt": 5.0e-5,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "bubbles_euler": "T",
                "bubble_model": 2,
                "polytropic": "F",
                "nb": 3,
                "fluid_pp(1)%gamma": 0.16,
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
                "patch_icpp(1)%pres": 2.0,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 4e-02,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%alpha_rho(1)": 0.96,
                "patch_icpp(3)%alpha(1)": 4e-02,
                "patch_icpp(3)%pres": 1.0,
            },
        )
        # override_tol: the non-polytropic bubble source is a stiff thermal+mechanical bubble ODE integrated every
        # RHS substage in a stiff liquid (pi_inf=3515); it amplifies GPU/CPU float-reassociation differences into the
        # bubble thermal primitives (~1.2e-10 on OpenACC, ~1.6e-9 on OpenMP offload) - genuine stiff-source
        # ill-conditioning (fields stay bounded and smooth), well under its ~7e-10 model-level conservation defect. A
        # real regression in the O(1) bubble physics would be >>5e-9, so this still catches bugs.
        cases.append(define_case_d(stack, "", {}, override_tol=5 * 10 ** (-9)))
        stack.pop()

        # (l) Euler-Euler bubbles, QBMM + polytropic (SP19): nb=3 R0 bins, each carrying a bivariate
        # 6-moment set (m00,m10,m01,m20,m11,m02) inverted by CHyQMOM every RHS call. polytropic=T keeps
        # pb/mv inert (degenerate stubs), so no quadrature side-state is advanced on the fine level and
        # all bubble state lives in q_cons. A left pressure slab (pres=2) drives a wave through the block
        # under regrid + subcycle. Exercises the realizability-preserving prolongation: the whole bub
        # block is injected piecewise-constant so every fine/ghost child inherits the coarse cell's
        # realizable moment set (variance c20 = m20/m00 - (m10/m00)^2 > 0), keeping the CHyQMOM inversion
        # NaN-free; the moments still reflux/restrict on the standard conservative q_cons path.
        stack.push(
            "AMR -> 1D -> bubbles QBMM",
            {
                **amr_1d_base,
                "dt": 5.0e-5,
                "amr_regrid_int": 2,
                "amr_tag_eps": 0.1,
                "amr_buf": 2,
                "amr_subcycle": "T",
                "bubbles_euler": "T",
                "bubble_model": 2,
                "polytropic": "T",
                "polydisperse": "F",
                "thermal": 3,
                "qbmm": "T",
                "dist_type": 2,
                "poly_sigma": 0.3,
                "sigR": 0.1,
                "sigV": 0.1,
                "rhoRV": 0.0,
                "nb": 3,
                "fluid_pp(1)%gamma": 0.16,
                "fluid_pp(1)%pi_inf": 3515.0,
                "bub_pp%R0ref": 1.0,
                "bub_pp%p0ref": 1.0,
                "bub_pp%rho0ref": 1.0,
                "bub_pp%T0ref": 1.0,
                "bub_pp%ss": 0.07179866765358993,
                "bub_pp%pv": 0.02308216136195411,
                "bub_pp%vd": 0.2404125083932959,
                "bub_pp%mu_l": 0.009954269975623244,
                "bub_pp%gam_g": 1.4,
                "patch_icpp(1)%r0": 1,
                "patch_icpp(1)%v0": 0,
                "patch_icpp(1)%alpha_rho(1)": 0.96,
                "patch_icpp(1)%alpha(1)": 4e-02,
                "patch_icpp(1)%pres": 2.0,
                "patch_icpp(2)%alpha_rho(1)": 0.96,
                "patch_icpp(2)%alpha(1)": 4e-02,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%alpha_rho(1)": 0.96,
                "patch_icpp(3)%alpha(1)": 4e-02,
                "patch_icpp(3)%pres": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        # (l2) QBMM NON-POLYTROPIC on a STATIC block: polytropic=F activates the pb/mv quadrature
        # side-state (nnode x nb per cell), which the block now carries itself: prolonged
        # piecewise-constant (realizability, like the moments), advanced with the block's own rhs
        # scratch (the coarse rhs_pb/mv stay untouched at fine indices), and restricted back with
        # the moments. Static block, no subcycle (both gated for nonpoly). Same override_tol
        # rationale as the nonpolytropic bubbles case: stiff thermal+mechanical bubble ODE in a
        # stiff liquid amplifies float-reassociation noise; a real bug is orders larger.
        stack.push(
            "nonpolytropic",
            {
                "amr_regrid_int": 0,
                "amr_subcycle": "F",
                "polytropic": "F",
                "bub_pp%mu_v": 8.758168074360729e-05,
                "bub_pp%mu_g": 0.00017881922111898042,
                "bub_pp%gam_v": 1.33,
                "bub_pp%M_v": 18.02,
                "bub_pp%M_g": 28.97,
                "bub_pp%k_v": 0.5583395141263873,
                "bub_pp%k_g": 0.7346421281308791,
                "bub_pp%R_v": 1334.8378710170155,
                "bub_pp%R_g": 830.2995663005393,
            },
        )
        cases.append(define_case_d(stack, "", {}, override_tol=5 * 10 ** (-9)))
        # dynamic regrid + subcycle: the pb/mv side-state now regrids (stor bounce + re-prolong +
        # overlap copy) and subcycles (two-source ghost time-lerp) exactly like q_cons
        stack.push("regrid subcycle", {"amr_regrid_int": 2, "amr_tag_eps": 0.1, "amr_buf": 2, "amr_subcycle": "T"})
        cases.append(define_case_d(stack, "", {}, override_tol=5 * 10 ** (-9)))
        stack.pop()
        stack.pop()
        stack.pop()
        # (n) STATIC IMMERSED BOUNDARY (SP20): a fixed circular body resolved on a static fine block that
        # covers it. Each fine block carries its own fine-grid IB markers/ghost points computed from the
        # geometry; the fine advance applies the IB state correction on the block. buff_size is floored to
        # 10 by ib, so the 64x64 base and the 24-coarse-cell block (beg=20, end=43: fine extent 47 <= 63,
        # >= 10 cells inside the domain) satisfy the block guards. Exercises: per-block fine IB setup, the
        # fine-block correct-state, and the coarse/fine coupling around a body (non-conservative by IB
        # nature at the body; conservative reflux elsewhere).
        stack.push(
            "AMR -> 2D -> static IBM circle",
            {
                "m": 63,
                "n": 63,
                "p": 0,
                "dt": 1.0e-4,
                "t_step_stop": 10,
                "t_step_save": 10,
                "num_patches": 1,
                "mixture_err": "T",
                "mapped_weno": "T",
                "mp_weno": "T",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                # single uniform patch: quiescent flow drifting toward +x
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.1,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                # static circular body at the domain center, inside the fine block
                "ib": "T",
                "num_ibs": 1,
                "fd_order": 2,
                "viscous": "F",
                "patch_ib(1)%geometry": 2,
                "patch_ib(1)%x_centroid": 0.5,
                "patch_ib(1)%y_centroid": 0.5,
                "patch_ib(1)%radius": 0.1,
                "patch_ib(1)%slip": "F",
                # static fine block covering the body (2:1)
                "amr": "T",
                "amr_block_beg(1)": 20,
                "amr_block_beg(2)": 20,
                "amr_block_end(1)": 43,
                "amr_block_end(2)": 43,
                "amr_regrid_int": 0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        # dynamic regrid with a STATIC body: candidate boxes expand to contain the body + margin,
        # overlapping expansions merge, and the fine IB state rebuilds from geometry every regrid
        stack.push("dynamic regrid", {"amr_regrid_int": 2, "amr_tag_eps": 1.0e-4, "amr_buf": 2, "t_step_stop": 20, "t_step_save": 20})
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (n2) MULTI-BODY static IB AMR: TWO static circular bodies sharing one static 2:1 fine block, in
        # quiescent flow drifting +x. The fine-IB setup reuses the multi-body-capable core routines, so both
        # bodies are marked/ghosted on the fine block (validated: fine-block ghost points = 2x the single body).
        stack.push(
            "AMR -> 2D -> static IBM two circles",
            {
                "m": 63,
                "n": 63,
                "p": 0,
                "dt": 1.0e-4,
                "t_step_stop": 10,
                "t_step_save": 10,
                "num_patches": 1,
                "mixture_err": "T",
                "mapped_weno": "T",
                "mp_weno": "T",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.1,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                # two static circular bodies, both inside the fine block
                "ib": "T",
                "num_ibs": 2,
                "fd_order": 2,
                "viscous": "F",
                "patch_ib(1)%geometry": 2,
                "patch_ib(1)%x_centroid": 0.42,
                "patch_ib(1)%y_centroid": 0.5,
                "patch_ib(1)%radius": 0.06,
                "patch_ib(1)%slip": "F",
                "patch_ib(2)%geometry": 2,
                "patch_ib(2)%x_centroid": 0.58,
                "patch_ib(2)%y_centroid": 0.5,
                "patch_ib(2)%radius": 0.06,
                "patch_ib(2)%slip": "F",
                # single static 2:1 fine block covering both bodies
                "amr": "T",
                "amr_block_beg(1)": 20,
                "amr_block_beg(2)": 20,
                "amr_block_end(1)": 43,
                "amr_block_end(2)": 43,
                "amr_regrid_int": 0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

        # (o) PRESCRIBED-MOTION MOVING IMMERSED BOUNDARY (SP21): a single circular body translating at a
        # prescribed velocity (moving_ibm=1) through quiescent flow, resolved on a STATIC fine block that
        # contains its whole trajectory. Each fine RK substage rebuilds the block's IB markers/ghost points
        # at the body's sub-time position (the same linear time-interpolation the subcycle applies to the
        # fluid ghosts). Exercises the per-substep fine-IB recompute and subcycled body-time consistency;
        # force-driven motion (moving_ibm=2) stays gated under amr.
        stack.push(
            "AMR -> 2D -> moving IBM circle",
            {
                "m": 63,
                "n": 63,
                "p": 0,
                "dt": 1.0e-3,
                "t_step_stop": 20,
                "t_step_save": 20,
                "num_patches": 1,
                "mixture_err": "T",
                "mapped_weno": "T",
                "mp_weno": "T",
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "bc_y%beg": -3,
                "bc_y%end": -3,
                # quiescent uniform flow; the body's prescribed motion is the sole driver
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                # circular body at the domain center translating in +y (prescribed, moving_ibm=1)
                "ib": "T",
                "num_ibs": 1,
                "fd_order": 2,
                "viscous": "F",
                "patch_ib(1)%geometry": 2,
                "patch_ib(1)%x_centroid": 0.5,
                "patch_ib(1)%y_centroid": 0.5,
                "patch_ib(1)%radius": 0.1,
                "patch_ib(1)%slip": "F",
                "patch_ib(1)%moving_ibm": 1,
                "patch_ib(1)%vel(2)": 1.0,
                # static fine block containing the body's whole trajectory (2:1)
                "amr": "T",
                "amr_subcycle": "T",
                "amr_block_beg(1)": 20,
                "amr_block_beg(2)": 20,
                "amr_block_end(1)": 43,
                "amr_block_end(2)": 43,
                "amr_regrid_int": 0,
            },
        )
        cases.append(define_case_d(stack, "", {}, restart_check=True))
        # TWO prescribed-motion bodies: the per-substage fine-IB rebuild runs the multi-body core
        # for moving bodies too (allowed since the num_ibs gate lift, previously unvalidated)
        stack.push(
            "two bodies",
            {
                "num_ibs": 2,
                "patch_ib(1)%x_centroid": 0.40,
                "patch_ib(1)%radius": 0.06,
                "patch_ib(2)%geometry": 2,
                "patch_ib(2)%x_centroid": 0.60,
                "patch_ib(2)%y_centroid": 0.5,
                "patch_ib(2)%radius": 0.06,
                "patch_ib(2)%slip": "F",
                "patch_ib(2)%moving_ibm": 1,
                "patch_ib(2)%vel(2)": 1.0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        # MOVING body + DYNAMIC regrid: candidate boxes expand over the body's LIVE centroid and
        # the per-substage containment guard holds between regrids. Mach 0.25 keeps the start
        # transient compact (at Mach ~0.85 the tagged wake legitimately outgrows the per-rank
        # box cap - a named abort, not a supported regime at this domain size); validated
        # tracking: 200 steps / 29 regrids with the body rising 25 cells past its initial box,
        # error at/below the static-block control (rho L2 5.6e-3 vs 8.5e-3)
        stack.push(
            "dynamic regrid",
            {
                "patch_ib(1)%vel(2)": 0.3,
                "amr_regrid_int": 5,
                "amr_tag_eps": 0.05,
                "amr_buf": 3,
                "t_step_stop": 40,
                "t_step_save": 40,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()
        stack.pop()

        # (p) 2D AXISYMMETRIC: an off-axis pressure pulse drives genuinely radial flow (nonzero
        # geometric sources) with the static fine block's lower-r edge at the MINIMUM legal axis
        # distance (amr_block_beg(2) = buff_size) - the stiffest 1/r a block can see. The axis
        # half-width cell makes the coarse y-WENO coefficients per-cell, so this also exercises
        # the per-swap coefficient recompute (amr_weno_coef_recompute). Validated against a
        # no-AMR reference (diffs match the Cartesian control's resolution effects; r-weighted
        # mass drift 1.09e-6 vs the reference's 1.07e-6).
        stack.push(
            "AMR -> 2D -> axisymmetric",
            {
                "m": 63,
                "n": 63,
                "p": 0,
                "cyl_coord": "T",
                "dt": 2.0e-4,
                "t_step_stop": 40,
                "t_step_save": 40,
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "y_domain%beg": 0.0,
                "y_domain%end": 1.0,
                "bc_x%beg": -1,
                "bc_x%end": -1,
                "bc_y%beg": -2,
                "bc_y%end": -6,
                "num_patches": 2,
                "mixture_err": "T",
                "mapped_weno": "T",
                "mp_weno": "T",
                "patch_icpp(1)%geometry": 3,
                "patch_icpp(1)%x_centroid": 0.5,
                "patch_icpp(1)%length_x": 1.0,
                "patch_icpp(1)%y_centroid": 0.5,
                "patch_icpp(1)%length_y": 1.0,
                "patch_icpp(1)%vel(1)": 0.0,
                "patch_icpp(1)%vel(2)": 0.0,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": 1.0,
                "patch_icpp(1)%alpha(1)": 1.0,
                "patch_icpp(2)%geometry": 2,
                "patch_icpp(2)%x_centroid": 0.5,
                "patch_icpp(2)%y_centroid": 0.28,
                "patch_icpp(2)%radius": 0.1,
                "patch_icpp(2)%alter_patch(1)": "T",
                "patch_icpp(2)%vel(1)": 0.0,
                "patch_icpp(2)%vel(2)": 0.0,
                "patch_icpp(2)%pres": 5.0,
                "patch_icpp(2)%alpha_rho(1)": 1.0,
                "patch_icpp(2)%alpha(1)": 1.0,
                "amr": "T",
                "amr_block_beg(1)": 20,
                "amr_block_beg(2)": 4,
                "amr_block_end(1)": 43,
                "amr_block_end(2)": 27,
                "amr_regrid_int": 0,
            },
        )
        cases.append(define_case_d(stack, "", {}))
        stack.pop()

    amr_golden_tests()

    def hybrid_sensor_tests():
        """Golden tests for the hybrid WENO/Riemann smoothness sensors: a 1D Sod-type shock so the
        sensor genuinely partitions the domain (full nonlinear WENO/upwind flux at the
        discontinuities, central weights/flux in the smooth regions). These protect the sensor
        plumbing and the shared nonlinear-weight block in m_weno against silent divergence."""
        hybrid_1d_base = {
            "m": 63,
            "n": 0,
            "p": 0,
            "dt": 5.0e-4,
            "t_step_stop": 6,
            "t_step_save": 6,
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.05,
            "patch_icpp(1)%length_x": 0.1,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.45,
            "patch_icpp(2)%length_x": 0.7,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(3)%geometry": 1,
            "patch_icpp(3)%x_centroid": 0.9,
            "patch_icpp(3)%length_x": 0.2,
            "patch_icpp(3)%vel(1)": 0.0,
        }
        stack.push("Hybrid -> 1D -> WENO sensor", {**hybrid_1d_base, "hybrid_weno": "T", "hybrid_weno_eps": 1.0e-2})
        cases.append(define_case_d(stack, "", {}))
        # liveness golden: eps 0.5 > the shock's phi (~0.3) puts consequential cells under the
        # sensor (answer moves ~5e-4 vs plain WENO); the eps=1e-2 case is bitwise-identical to
        # plain by design (only constant/eps-dominated cells go central) so it protects
        # no-corruption but cannot detect a dead sensor
        cases.append(define_case_d(stack, "consequential eps", {"hybrid_weno_eps": 0.5}, override_tol=5.0e-5))
        stack.pop()
        stack.push(
            "Hybrid -> 1D -> Riemann sensor",
            {**hybrid_1d_base, "hybrid_riemann": "T", "hybrid_weno_eps": 1.0e-2, "hybrid_smooth_flux": 2},
        )
        cases.append(define_case_d(stack, "", {}))
        cases.append(define_case_d(stack, "consequential eps", {"hybrid_weno_eps": 0.5}, override_tol=5.0e-5))
        # the central smooth-flux (enum 1) is a distinct flux path from Rusanov (2) - cover both
        cases.append(define_case_d(stack, "central flux", {"hybrid_smooth_flux": 1}))
        stack.pop()

    hybrid_sensor_tests()

    def load_balance_tests():
        """Golden for the weighted init-time decomposition (load_balance): a two-fluid material
        interface at x=0.5 makes the alpha marginal asymmetric (fluid-1 volume fraction ~1 left,
        ~0 right), so the 2-rank weighted split genuinely differs from the equal split and
        s_apply_weighted_offsets re-decomposes. This is the only coverage of that path - the
        feature is a no-op at 1 rank by construction, and wrong weighted halo extents would
        corrupt the solution everywhere while looking plausible."""
        eps_lb = 1.0e-6
        stack.push(
            "LoadBalance -> 1D -> weighted split",
            {
                "m": 63,
                "n": 0,
                "p": 0,
                "dt": 5.0e-4,
                "t_step_stop": 6,
                "t_step_save": 6,
                "x_domain%beg": 0.0,
                "x_domain%end": 1.0,
                "bc_x%beg": -3,
                "bc_x%end": -3,
                "parallel_io": "T",
                "load_balance": "T",
                "num_fluids": 2,
                "fluid_pp(2)%gamma": 1.0e00 / (1.6e00 - 1.0e00),
                "fluid_pp(2)%pi_inf": 0.0,
                "fluid_pp(2)%cv": 0.0,
                "fluid_pp(2)%qv": 0.0,
                "fluid_pp(2)%qvp": 0.0,
                "patch_icpp(1)%geometry": 1,
                "patch_icpp(1)%x_centroid": 0.05,
                "patch_icpp(1)%length_x": 0.1,
                "patch_icpp(1)%vel(1)": 0.5,
                "patch_icpp(2)%geometry": 1,
                "patch_icpp(2)%x_centroid": 0.3,
                "patch_icpp(2)%length_x": 0.4,
                "patch_icpp(2)%vel(1)": 0.5,
                "patch_icpp(3)%geometry": 1,
                "patch_icpp(3)%x_centroid": 0.75,
                "patch_icpp(3)%length_x": 0.5,
                "patch_icpp(3)%vel(1)": 0.5,
                "patch_icpp(1)%pres": 1.0,
                "patch_icpp(2)%pres": 1.0,
                "patch_icpp(3)%pres": 1.0,
                "patch_icpp(1)%alpha_rho(1)": (1.0 - eps_lb) * 1.0,
                "patch_icpp(1)%alpha_rho(2)": eps_lb * 10.0,
                "patch_icpp(1)%alpha(1)": 1.0 - eps_lb,
                "patch_icpp(1)%alpha(2)": eps_lb,
                "patch_icpp(2)%alpha_rho(1)": (1.0 - eps_lb) * 1.0,
                "patch_icpp(2)%alpha_rho(2)": eps_lb * 10.0,
                "patch_icpp(2)%alpha(1)": 1.0 - eps_lb,
                "patch_icpp(2)%alpha(2)": eps_lb,
                "patch_icpp(3)%alpha_rho(1)": eps_lb * 1.0,
                "patch_icpp(3)%alpha_rho(2)": (1.0 - eps_lb) * 10.0,
                "patch_icpp(3)%alpha(1)": eps_lb,
                "patch_icpp(3)%alpha(2)": 1.0 - eps_lb,
            },
        )
        cases.append(define_case_d(stack, "", {}, ppn=2, honor_io_keys=True))
        stack.pop()

    load_balance_tests()

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
