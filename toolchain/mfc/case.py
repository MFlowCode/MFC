import copy
import dataclasses
import difflib
import json
import math
import re

import fastjsonschema

from . import common
from .analytic_expr import AnalyticExprError, fortranize_expr
from .params.eos_families import EOS_FAMILIES
from .printer import cons
from .run import case_dicts
from .state import ARG


def _suggest_similar_params(unknown_key: str, valid_keys: list, n: int = 3) -> list:
    """Find similar parameter names for typo suggestions."""
    return difflib.get_close_matches(unknown_key, valid_keys, n=n, cutoff=0.6)


QPVF_IDX_VARS = {
    "alpha_rho": "eqn_idx%cont%beg",
    "vel": "eqn_idx%mom%beg",
    "pres": "eqn_idx%E",
    "alpha": "eqn_idx%adv%beg",
    "tau_e": "eqn_idx%stress%beg",
    "Y": "eqn_idx%species%beg",
    "cf_val": "eqn_idx%c",
    "Bx": "eqn_idx%B%beg",
    "By": "eqn_idx%B%end-1",
    "Bz": "eqn_idx%B%end",
}

# fluid_pp(i)%eos values whose coefficients depend on the local state; derived so a new
# family never leaves any_state_dependent_eos stale (see EOS_FAMILIES).
EOS_STATE_DEPENDENT_VALUES = frozenset(f.value for f in EOS_FAMILIES if f.state_dependent)

MIBM_ANALYTIC_VARS = ["vel(1)", "vel(2)", "vel(3)", "angular_vel(1)", "angular_vel(2)", "angular_vel(3)"]
# "eqn_idx%B%end - 1" not "eqn_idx%B%beg + 1" must be used because 1D does not have Bx


@dataclasses.dataclass(init=False)
class Case:
    params: dict

    def __init__(self, params: dict) -> None:
        self.params = copy.deepcopy(params)
        for key, val in self.params.items():
            self.params[key] = self.__normalize_named_value(key, val)

    @staticmethod
    def __normalize_named_value(key: str, val):
        """Convert a named enumerated value (e.g. "hllc") to its integer code."""
        from .params.errors import constraint_error
        from .params.registry import REGISTRY

        if not isinstance(val, str) or common.is_number(val):
            return val
        param_def = REGISTRY.get_param_def(key)
        if param_def is None or not param_def.constraints:
            return val
        names = param_def.constraints.get("names")
        if not names:
            return val
        if val in names:
            return names[val]
        shown = ", ".join(f"{v} ({n})" for n, v in sorted(names.items(), key=lambda kv: kv[1]))
        raise common.MFCException(constraint_error(key, "choices", shown, val))

    def get_parameters(self) -> dict:
        return self.params

    def get_cell_count(self) -> int:
        return math.prod([max(1, int(self.params.get(dir, 0))) for dir in ["m", "n", "p"]])

    def has_parameter(self, key: str) -> bool:
        return key in self.params.keys()

    def gen_json_dict_str(self) -> str:
        return json.dumps(self.params, indent=4)

    def get_inp(self, _target) -> str:
        from . import build

        target = build.get_target(_target)

        cons.print(f"Generating [magenta]{target.name}.inp[/magenta]:")
        cons.indent()

        MASTER_KEYS = case_dicts.get_input_dict_keys(target.name)

        ignored = []

        # Create Fortran-style input file content string
        dict_str = ""
        for key, val in self.params.items():
            if key in MASTER_KEYS and key not in case_dicts.IGNORE:
                if self.__is_ic_analytical(key, val) or self.__is_mib_analytical(key, val):
                    dict_str += f"{key} = 0d0\n"
                    ignored.append(key)
                    continue

                if not isinstance(val, str) or len(val) == 1:
                    dict_str += f"{key} = {val}\n"
                else:
                    dict_str += f"{key} = '{val}'\n"
            else:
                ignored.append(key)

            if key not in case_dicts.ALL:
                suggestions = _suggest_similar_params(key, list(case_dicts.ALL.keys()))
                hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
                raise common.MFCException(f"Unknown parameter '{key}'.{hint}")

        cons.print(f"[yellow]INFO:[/yellow] Forwarded {len(self.params) - len(ignored)}/{len(self.params)} parameters.")
        cons.unindent()

        return f"&user_inputs\n{dict_str}&end/\n"

    def validate_params(self, origin_txt: str = None):
        """Validates parameters read from case file:
        1. Type checking via JSON schema
        2. Constraint validation (valid values, ranges)
        3. Dependency checking (required/recommended params)
        """
        # Type checking
        try:
            case_dicts.get_validator()(self.params)
        except fastjsonschema.JsonSchemaException as e:
            if origin_txt:
                raise common.MFCException(f"{origin_txt}: {e}")
            raise common.MFCException(f"{e}")

        # Constraint and dependency validation
        from .params.validate import validate_case

        errors, warnings = validate_case(self.params)

        # Show warnings (non-fatal)
        if warnings:
            cons.print()
            for w in warnings:
                cons.print(f"[yellow]Warning:[/yellow] {w}")

        # Raise errors (fatal)
        if errors:
            error_msg = "\n".join(f"  - {e}" for e in errors)
            if origin_txt:
                raise common.MFCException(f"{origin_txt}:\n{error_msg}")
            raise common.MFCException(f"Validation errors:\n{error_msg}")

    def __get_ndims(self) -> int:
        return 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)

    def __grid_geometry(self) -> int:
        # grid_geometry as m_global_parameters.fpp computes it (1 Cartesian, 2 axisymmetric, 3 cylindrical).
        if self.params.get("cyl_coord", "F") != "T":
            return 1
        return 2 if int(self.params.get("p", 0)) == 0 else 3

    def non_newtonian_fluids(self) -> list:
        """1-based indices of fluids with non_newtonian=T, the same fluids hllc_nn_fluids
        and the kernel's unrolled mixture sum enumerate."""
        num_fluids = int(self.params.get("num_fluids", 1))
        return [i for i in range(1, num_fluids + 1) if str(self.params.get(f"fluid_pp({i})%non_newtonian", "F")) == "T"]

    def any_state_dependent_eos(self) -> bool:
        """True when some fluid's EOS family carries density-dependent coefficients
        (EOS_STATE_DEPENDENT_VALUES from the EOS_FAMILIES registry)."""
        num_fluids = int(self.params.get("num_fluids", 1))
        return any(int(self.params.get(f"fluid_pp({f})%eos", 1)) in EOS_STATE_DEPENDENT_VALUES for f in range(1, num_fluids + 1))

    def fused_viscous_stress(self) -> bool:
        """The fused kernel assembles and differences the viscous stress. Matches
        f_use_fused_viscous_stress in m_rhs.fpp. Non-Newtonian fluids need a Cartesian
        grid, WENO reconstruction and Re(1) > 0, 3-D cylindrical grids need WENO
        reconstruction, and non-Cartesian grids need non-characteristic boundaries;
        other shapes use the host pass."""
        if self.params.get("viscous", "F") != "T":
            return False
        if self.__grid_geometry() == 3 and int(self.params.get("recon_type", 1)) != 1:
            return False
        if self.non_newtonian_fluids():
            if self.__grid_geometry() != 1:
                return False
            if int(self.params.get("recon_type", 1)) != 1:
                return False
            if any(float(self.params.get(f"fluid_pp({i})%Re(1)", 0)) <= 0 for i in self.non_newtonian_fluids()):
                return False
        if self.__grid_geometry() != 1 and any(-12 <= int(self.params.get(key, 0)) <= -5 for key in ("bc_x%beg", "bc_x%end", "bc_y%beg", "bc_y%end", "bc_z%beg", "bc_z%end")):
            return False
        return True

    def weno_re_flux_stencil(self) -> bool:
        """The weno_Re_flux shapes served with a plain WENO-stencil reconstruction
        (orders 3, 5 and 7; MUSCL and TENO7 keep the host families)."""
        return (
            self.fused_viscous_stress()
            and self.params.get("weno_Re_flux", "F") == "T"
            and int(self.params.get("recon_type", 1)) == 1
            and int(self.params.get("weno_order", 5)) in (3, 5, 7)
            and not (self.params.get("teno", "F") == "T" and int(self.params.get("weno_order", 5)) == 7)
        )

    def fused_grad_weno(self) -> bool:
        """The kernel WENO-reconstructs the weno_Re_flux face gradients in-column."""
        return self.weno_re_flux_stencil()

    def fused_vel_gradients(self) -> bool:
        """One raw-array kernel builds the weno_Re_flux cell-centred velocity
        gradients straight from q_prim."""
        return self.weno_re_flux_stencil()

    def fused_rhs_bubbles(self) -> bool:
        """Bubble transport shapes supported by the fused column."""
        qbmm = self.params.get("qbmm", "F") == "T"
        adv_n = self.params.get("adv_n", "F") == "T"
        bins = int(self.params.get("nb", 1))
        rows = bins * (6 if qbmm else 2 if self.params.get("polytropic", "T") == "T" else 4) + int(adv_n)
        return (
            self.params.get("bubbles_euler", "F") == "T"
            and int(self.params.get("model_eqns", 2)) == 2
            and int(self.params.get("num_fluids", 1)) == 1
            and self.params.get("mpp_lim", "F") != "T"
            and int(self.params.get("riemann_solver", 2)) == 2
            and bins >= 1
            and rows <= (18 if qbmm else 12)
            and (not qbmm or not adv_n)
        )

    def fused_rhs_hypoelastic(self) -> bool:
        """HLL stress transport in the finite-difference NC mode, with sources completed
        before fused RK. The interface-RHS mode uses the generic path, which exports the
        Riemann interface velocities."""
        return (
            self.params.get("hypoelasticity", "F") == "T"
            and self.params.get("hypo_hll_interface_rhs", "F") != "T"
            and int(self.params.get("riemann_solver", 2)) == 1
            and self.params.get("mpp_lim", "F") != "T"
            and int(self.params.get("low_Mach", 0)) == 0
            and (self.params.get("cyl_coord", "F") != "T" or int(self.params.get("p", 0)) == 0)
        )

    def fused_rhs_chemistry(self) -> bool:
        """Species transport with the existing thermochemical closure, in builds with a restructured solver lane."""
        return (
            (ARG("fused_rhs", False) or ARG("hllc_register_tuning", False))
            and self.params.get("chemistry", "F") == "T"
            and int(self.params.get("model_eqns", 2)) == 2
            and int(self.params.get("num_fluids", 1)) == 1
            and int(self.params.get("wave_speeds", 1)) == 1
            and self.params.get("mpp_lim", "F") != "T"
            and self.params.get("cyl_coord", "F") != "T"
            and not any(self.params.get(name, "F") == "T" for name in ("bubbles_euler", "hypoelasticity", "surface_tension", "mhd"))
        )

    def fused_rhs_mhd(self) -> bool:
        """Single-fluid HLL/HLLD magnetic transport."""
        solver = int(self.params.get("riemann_solver", 2))
        return (
            self.params.get("mhd", "F") == "T"
            and int(self.params.get("model_eqns", 2)) == 2
            and int(self.params.get("num_fluids", 1)) == 1
            and int(self.params.get("wave_speeds", 1)) == 1
            and int(self.params.get("low_Mach", 0)) == 0
            and self.params.get("mpp_lim", "F") != "T"
            and self.params.get("cyl_coord", "F") != "T"
            and not any(self.params.get(name, "F") == "T" for name in ("bubbles_euler", "hypoelasticity", "surface_tension", "chemistry", "cont_damage"))
            and (self.params.get("hyper_cleaning", "F") != "T" or int(self.params.get("n", 0)) > 0)
            and self.params.get("relativity", "F") != "T"
            and (solver == 1 or (solver == 4 and not any(self.params.get(name, "F") == "T" for name in ("viscous", "hyper_cleaning"))))
        )

    def fused_rhs_scalar_thermo(self) -> bool:
        """Fuse scalar thermodynamic reductions in HLL/LF/HLLC faces."""
        return self.fused_rhs_chemistry() and int(self.params.get("riemann_solver", 2)) in (1, 2, 5)

    def __is_ic_analytical(self, key: str, val: str) -> bool:
        """Is this initial condition analytical?
        More precisely, is this an arbitrary expression or a string representing a number?"""
        if common.is_number(val) or not isinstance(val, str):
            return False

        for array in QPVF_IDX_VARS:
            if re.match(rf"^patch_icpp\([0-9]+\)%{array}", key):
                return True

        return False

    def __is_mib_analytical(self, key: str, val: str) -> bool:
        """Is this initial condition analytical?
        More precisely, is this an arbitrary expression or a string representing a number?"""
        if common.is_number(val) or not isinstance(val, str):
            return False

        for variable in MIBM_ANALYTIC_VARS:
            if re.match(rf"^patch_ib\([0-9]+\)%{re.escape(variable)}", key):
                return True

        return False

    def __get_analytic_ic_fpp(self, print: bool) -> str:
        # generates the content of an FFP file that will hold the functions for
        # some initial condition
        DATA = {
            1: {"ptypes": [1, 15, 16], "sf_idx": "i, 0, 0"},
            2: {"ptypes": [2, 3, 4, 5, 6, 7, 13, 17, 18, 21], "sf_idx": "i, j, 0"},
            3: {"ptypes": [8, 9, 10, 11, 12, 14, 19, 21], "sf_idx": "i, j, k"},
        }[self.__get_ndims()]

        patches = {}

        # iterates over the parameters and checks if they are defined as an
        # analytical function. If so, append it to the `patches`` object
        for key, val in self.params.items():
            if not self.__is_ic_analytical(key, val):
                continue

            patch_id = re.search(r"[0-9]+", key).group(0)

            if patch_id not in patches:
                patches[patch_id] = []

            patches[patch_id].append((key, val))

        srcs = []

        # for each analytical patch that is required to be added, generate
        # the string that contains that function.
        for pid, items in patches.items():
            ptype = self.params[f"patch_icpp({pid})%geometry"]

            if ptype not in DATA["ptypes"]:
                raise common.MFCException(f"Patch #{pid} of type {ptype} cannot be analytically defined.")

            var_map = {
                "x": "x_cc(i)",
                "y": "y_cc(j)",
                "z": "z_cc(k)",
                "xc": f"patch_icpp({pid})%x_centroid",
                "yc": f"patch_icpp({pid})%y_centroid",
                "zc": f"patch_icpp({pid})%z_centroid",
                "lx": f"patch_icpp({pid})%length_x",
                "ly": f"patch_icpp({pid})%length_y",
                "lz": f"patch_icpp({pid})%length_z",
                "r": f"patch_icpp({pid})%radius",
                "eps": f"patch_icpp({pid})%epsilon",
                "beta": f"patch_icpp({pid})%beta",
                "tau_e": f"patch_icpp({pid})%tau_e",
                "radii": f"patch_icpp({pid})%radii",
            }

            lines = []
            # translate each analytic expression to Fortran and emit an assignment
            for attribute, expr in items:
                if print:
                    cons.print(f"* Codegen: {attribute} = {expr}")

                varname = re.findall(r"[a-zA-Z][a-zA-Z0-9_]*", attribute)[1]
                qpvf_idx = QPVF_IDX_VARS[varname][:]

                if len(re.findall(r"[0-9]+", attribute)) == 2:
                    idx = int(re.findall(r"[0-9]+", attribute)[1]) - 1
                    qpvf_idx = f"{qpvf_idx} + {idx}"

                lhs = f"q_prim_vf({qpvf_idx})%sf({DATA['sf_idx']})"
                try:
                    rhs = fortranize_expr(expr, var_map)
                except AnalyticExprError as err:
                    raise common.MFCException(f"{attribute}: {err}") from err

                lines.append(f"        {lhs} = {rhs}")

            # concatenates all of the analytic lines into a single string with
            # each element separated by new line characters. Then write those
            # new lines as a fully concatenated string with fortran syntax
            srcs.append(f"""\
    if (patch_id == {pid}) then
{f"{chr(10)}".join(lines)}
    end if\
""")

        content = f"""\
! This file was generated by MFC. It is only used when analytical patches are
! present in the input file. It is used to define the analytical patches with
! expressions that are evaluated at runtime from the input file.

#:def analytical()
{f"{chr(10)}{chr(10)}".join(srcs)}
#:enddef
"""
        return content

    # gets the analytic description of a moving IB's velocity and rotation rate
    def __get_analytic_mib_fpp(self, print: bool) -> str:
        # iterates over the parameters and checks if they are defined as an
        # analytical function. If so, append it to the `patches`` object
        ib_patches = {}

        for key, val in self.params.items():
            if not self.__is_mib_analytical(key, val):
                continue

            patch_id = re.search(r"[0-9]+", key).group(0)

            if patch_id not in ib_patches:
                ib_patches[patch_id] = []

            ib_patches[patch_id].append((key, val))

        srcs = []

        # for each analytical patch that is required to be added, generate
        # the string that contains that function.
        for pid, items in ib_patches.items():
            var_map = {
                "x": "x_cc(i)",
                "y": "y_cc(j)",
                "z": "z_cc(k)",
                "t": "mytime",
                "r": f"patch_ib({pid})%radius",
            }

            lines = []
            # translate each analytic expression to Fortran and emit an assignment
            for attribute, expr in items:
                if print:
                    cons.print(f"* Codegen: {attribute} = {expr}")

                lhs = attribute
                try:
                    rhs = fortranize_expr(expr, var_map)
                except AnalyticExprError as err:
                    raise common.MFCException(f"{attribute}: {err}") from err

                lines.append(f"        {lhs} = {rhs}")

            # concatenates all of the analytic lines into a single string with
            # each element separated by new line characters. Then write those
            # new lines as a fully concatenated string with fortran syntax
            srcs.append(f"""\
    if (gbl_id == {pid}) then
{f"{chr(10)}".join(lines)}
    end if\
""")

        content = f"""\
! This file was generated by MFC. It is only used when we analytically
! parameterize the velocity and rotation rate of a moving IB.

#:def mib_analytical()
gbl_id = patch_ib(i)%gbl_patch_id
{f"{chr(10)}{chr(10)}".join(srcs)}
#:enddef
"""
        return content

    def __get_sim_fpp(self, print: bool) -> str:
        if ARG("case_optimization"):
            if print:
                cons.print("Case optimization is enabled.")

            nterms = 1

            bubble_model = int(self.params.get("bubble_model", "-100"))

            if bubble_model == 2:
                nterms = 32
            elif bubble_model == 3:
                nterms = 7

            mapped_weno = 1 if self.params.get("mapped_weno", "F") == "T" else 0
            wenoz = 1 if self.params.get("wenoz", "F") == "T" else 0
            teno = 1 if self.params.get("teno", "F") == "T" else 0
            wenojs = 0 if (mapped_weno or wenoz or teno) else 1
            recon_type = self.params.get("recon_type", 1)

            # This fixes a bug on Frontier to do with allocating 0:0 arrays
            weno_order = int(self.params.get("weno_order", 0))
            if recon_type == 1:
                weno_polyn = int((weno_order - 1) / 2)
            else:
                weno_polyn = 1

            if self.params.get("igr", "F") == "T":
                weno_order = 5
                weno_polyn = 3

            if teno:
                weno_num_stencils = weno_order - 3
            else:
                weno_num_stencils = weno_polyn

            num_dims = 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)
            if self.params.get("mhd", "F") == "T":
                num_vels = 3
            else:
                num_vels = num_dims

            # Baking this lets the compiler drop the state-dependent EOS chain entirely. Left in the
            # call graph it costs registers, and so occupancy, in every kernel that can reach it.
            num_fluids_case = int(self.params.get("num_fluids", 1))
            any_state_dependent_eos = 1 if self.any_state_dependent_eos() else 0

            mhd = 1 if self.params.get("mhd", "F") == "T" else 0
            relativity = 1 if self.params.get("relativity", "F") == "T" else 0
            viscous = 1 if self.params.get("viscous", "F") == "T" else 0
            igr = 1 if self.params.get("igr", "F") == "T" else 0
            igr_pres_lim = 1 if self.params.get("igr_pres_lim", "F") == "T" else 0

            # Parameters the simulation sources additionally bake under case optimization.
            # Defaults mirror the Fortran namelist defaults each parameter reads.
            bubbles_euler = 1 if self.params.get("bubbles_euler", "F") == "T" else 0
            qbmm = 1 if self.params.get("qbmm", "F") == "T" else 0
            polytropic = 1 if self.params.get("polytropic", "T") == "T" else 0
            adv_n = 1 if self.params.get("adv_n", "F") == "T" else 0
            hypoelasticity = 1 if self.params.get("hypoelasticity", "F") == "T" else 0
            cyl_coord = 1 if self.params.get("cyl_coord", "F") == "T" else 0
            # Heat conduction is on when any fluid has k_therm > 0, as in m_global_parameters_common.
            heat_conduction = 1 if any(float(self.params.get(f"fluid_pp({i})%k_therm", 0.0)) > 0 for i in range(1, num_fluids_case + 1)) else 0
            surface_tension = 1 if self.params.get("surface_tension", "F") == "T" else 0
            hyper_cleaning = 1 if self.params.get("hyper_cleaning", "F") == "T" else 0
            cont_damage = 1 if self.params.get("cont_damage", "F") == "T" else 0
            qvs = [float(self.params.get(f"fluid_pp({i})%qv", 0.0)) for i in range(1, num_fluids_case + 1)]
            qv_zero = all(qv == 0.0 and math.copysign(1.0, qv) > 0.0 for qv in qvs)
            hllc_bcs = [int(self.params.get(f"bc_{axis}%{side}", 0)) for axis in "xyz" for side in ("beg", "end")]
            hllc_cbc = any(-12 <= bc <= -5 for bc in hllc_bcs)

            # Throw error if wenoz_q is required but not set
            out = f"""\
#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
#:set recon_type            = {recon_type}
#:set weno_order            = {weno_order}
#:set weno_polyn            = {weno_polyn}
#:set muscl_order           = {int(self.params.get("muscl_order", 1))}
#:set muscl_polyn           = {int(self.params.get("muscl_order", 1))}
#:set muscl_lim             = {int(self.params.get("muscl_lim", 1))}
#:set weno_num_stencils     = {weno_num_stencils}
#:set nb                    = {int(self.params.get("nb", 1))}
#:set num_dims              = {num_dims}
#:set num_vels              = {num_vels}
#:set nterms                = {nterms}
#:set num_fluids            = {int(self.params.get("num_fluids", 1))}
#:set wenojs                = {wenojs}
#:set mapped_weno           = {mapped_weno}
#:set wenoz                 = {wenoz}
#:set teno                  = {teno}
#:set wenoz_q               = {self.params.get("wenoz_q", -1)}
#:set mhd                   = {mhd}
#:set relativity            = {relativity}
#:set igr                   = {igr}
#:set igr_iter_solver       = {self.params.get("igr_iter_solver", 1)}
#:set igr_pres_lim          = {igr_pres_lim}
#:set igr_order             = {self.params.get("igr_order", 3)}
#:set viscous               = {viscous}
#:set any_state_dependent_eos = {any_state_dependent_eos}
#:set bubbles_euler         = {bubbles_euler}
#:set qbmm                  = {qbmm}
#:set polytropic            = {polytropic}
#:set adv_n                 = {adv_n}
#:set hypoelasticity        = {hypoelasticity}
#:set cyl_coord             = {cyl_coord}
#:set model_eqns            = {int(self.params.get("model_eqns", 2))}
#:set heat_conduction       = {heat_conduction}
#:set surface_tension       = {surface_tension}
#:set hyper_cleaning        = {hyper_cleaning}
#:set cont_damage           = {cont_damage}
#:set qv_zero               = {qv_zero}
#:set hllc_cbc              = {hllc_cbc}
#:set hllc_extrap           = {-4 in hllc_bcs}
#:set hllc_riemann          = {int(self.params.get("riemann_solver", 2))}
#:set hllc_hll_u_interface  = {self.params.get("hll_u_interface", "F") == "T"}
#:set hllc_wave_speeds      = {int(self.params.get("wave_speeds", 1))}
#:set hllc_low_mach         = {int(self.params.get("low_Mach", 0))}
#:set hllc_alt_soundspeed   = {self.params.get("alt_soundspeed", "F") == "T"}
#:set hllc_mpp_lim          = {self.params.get("mpp_lim", "F") == "T"}
#:set hllc_mp_weno          = {self.params.get("mp_weno", "F") == "T"}
#:set hllc_int_comp         = {int(self.params.get("int_comp", 0) or 0)}
#:set hllc_gamma_method     = {int(self.params.get("chem_params%gamma_method", 1))}
#:set hllc_bubbles          = {self.fused_rhs_bubbles()}
#:set hllc_hypoelastic      = {self.fused_rhs_hypoelastic()}
#:set hllc_chemistry        = {self.fused_rhs_chemistry()}
#:set hllc_mhd              = {self.fused_rhs_mhd()}
#:set hllc_scalar_thermo    = {self.fused_rhs_scalar_thermo()}
#:set hllc_viscous_stress   = {self.fused_viscous_stress()}
#:set hllc_grad_weno        = {self.fused_grad_weno()}
#:set hllc_vel_gradients    = {self.fused_vel_gradients()}
#:set hllc_weno_re_flux     = {self.params.get("weno_Re_flux", "F") == "T"}
#:set hllc_non_newtonian    = {any(str(self.params.get(f"fluid_pp({i})%non_newtonian", "F")) == "T" for i in range(1, int(self.params.get("num_fluids", 1)) + 1))}
#:set hllc_nn_fluids        = {[str(self.params.get(f"fluid_pp({i})%non_newtonian", "F")) == "T" for i in range(1, int(self.params.get("num_fluids", 1)) + 1)]}
#:set hllc_reactive_burn    = {self.params.get("reactive_burn", "F") == "T"}
"""

        else:
            out = ""

        out = out + self.__get_analytic_mib_fpp(print)

        # We need to also include the pre_processing includes so that common subroutines have access to the @:analytical function
        return out + f"\n{self.__get_pre_fpp(print)}"

    def __get_pre_fpp(self, print: bool) -> str:
        out = self.__get_analytic_ic_fpp(print)
        return out

    def get_fpp(self, target, print=True) -> str:
        from . import build

        def _prepend() -> str:
            return f"""\
#:set chemistry             = {self.params.get("chemistry", "F") == "T"}
"""

        def _default(_) -> str:
            return "! This file is purposefully empty."

        result = {
            "pre_process": self.__get_pre_fpp,
            "simulation": self.__get_sim_fpp,
        }.get(build.get_target(target).name, _default)(print)

        return _prepend() + result

    def __getitem__(self, key: str) -> str:
        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = self.__normalize_named_value(key, val)
