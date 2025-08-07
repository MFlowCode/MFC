import re, json, math, copy, dataclasses, fastjsonschema

from . import common
from . import build
from .printer import cons

from .state import ARG
from .run   import case_dicts

QPVF_IDX_VARS = {
    'alpha_rho': 'contxb', 'vel'  : 'momxb',         'pres': 'E_idx', 
    'alpha':     'advxb',  'tau_e': 'stress_idx%beg', 'Y':   'chemxb',
    'cf_val': 'c_idx', 'Bx': 'B_idx%beg', 'By': 'B_idx%end-1', 'Bz': 'B_idx%end',
}
# "B_idx%end - 1" not "B_idx%beg + 1" must be used because 1D does not have Bx

@dataclasses.dataclass(init=False)
class Case:
    params: dict

    def __init__(self, params: dict) -> None:
        self.params = copy.deepcopy(params)

    def get_parameters(self) -> dict:
        return self.params

    def get_cell_count(self) -> int:
        return math.prod([max(1, int(self.params.get(dir, 0))) for dir in ["m", "n", "p"]])

    def has_parameter(self, key: str)-> bool:
        return key in self.params.keys()

    def gen_json_dict_str(self) -> str:
        return json.dumps(self.params, indent=4)

    def get_inp(self, _target) -> str:
        target = build.get_target(_target)

        cons.print(f"Generating [magenta]{target.name}.inp[/magenta]:")
        cons.indent()

        MASTER_KEYS: list = case_dicts.get_input_dict_keys(target.name)

        ignored = []

        # Create Fortran-style input file content string
        dict_str = ""
        for key, val in self.params.items():
            if key in MASTER_KEYS and key not in case_dicts.IGNORE:
                if self.__is_ic_analytical(key, val):
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
                raise common.MFCException(f"MFCInputFile::dump: Case parameter '{key}' is not used by any MFC code. Please check your spelling or add it as a new parameter.")

        cons.print(f"[yellow]INFO:[/yellow] Forwarded {len(self.params)-len(ignored)}/{len(self.params)} parameters.")
        cons.unindent()

        return f"&user_inputs\n{dict_str}&end/\n"

    def validate_params(self, origin_txt: str = None):
        '''Typechecks parameters read from case file. If a parameter
        is assigned a vlaie of the wrong type, this method throws an exception
        highlighting the violating parameter and specifying what it expects.'''
        try:
            case_dicts.get_validator()(self.params)
        except fastjsonschema.JsonSchemaException as e:
            if origin_txt:
                raise common.MFCException(f"{origin_txt}: {e}")

            raise common.MFCException(f"{e}")

    def __get_ndims(self) -> int:
        return 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)

    def __is_ic_analytical(self, key: str, val: str) -> bool:
        '''Is this initial condition analytical?
        More precisely, is this an arbitrary expression or a string representing a number?'''
        if common.is_number(val) or not isinstance(val, str):
            return False

        for array in QPVF_IDX_VARS:
            if re.match(fr'^patch_icpp\([0-9]+\)%{array}', key):
                return True

        return False

    # pylint: disable=too-many-locals
    def __get_pre_fpp(self, print: bool) -> str:
        # generates the content of an FFP file that will hold the functions for
        # some initial condition
        DATA = {
            1: {'ptypes': [1, 15, 16],                         'sf_idx': 'i, 0, 0'},
            2: {'ptypes': [2,  3,  4,  5,  6,  7, 17, 18, 21], 'sf_idx': 'i, j, 0'},
            3: {'ptypes': [8,  9, 10, 11, 12, 13, 14, 19, 21], 'sf_idx': 'i, j, k'}
        }[self.__get_ndims()]

        patches = {}

        # iterates over the parameters and checks if they are defined as an
        # analytical function. If so, append it to the `patches`` object
        for key, val in self.params.items():
            if not self.__is_ic_analytical(key, val):
                continue

            patch_id = re.search(r'[0-9]+', key).group(0)

            if patch_id not in patches:
                patches[patch_id] = []

            patches[patch_id].append((key, val))

        srcs = []

        # for each analytical patch that is required to be added, generate
        # the string that contains that function.
        for pid, items in patches.items():
            ptype = self.params[f"patch_icpp({pid})%geometry"]

            if ptype not in DATA['ptypes']:
                raise common.MFCException(f"Patch #{pid} of type {ptype} cannot be analytically defined.")

            # function that defines how we will replace variable names with
            # values from the case file
            def rhs_replace(match):
                return {
                    'x': 'x_cc(i)', 'y': 'y_cc(j)', 'z': 'z_cc(k)',

                    'xc': f'patch_icpp({pid})%x_centroid', 'yc': f'patch_icpp({pid})%y_centroid', 'zc': f'patch_icpp({pid})%z_centroid',
                    'lx': f'patch_icpp({pid})%length_x',   'ly': f'patch_icpp({pid})%length_y',   'lz': f'patch_icpp({pid})%length_z',

                    'r':     f'patch_icpp({pid})%radius',  'eps':   f'patch_icpp({pid})%epsilon', 'beta':  f'patch_icpp({pid})%beta',
                    'tau_e': f'patch_icpp({pid})%tau_e',   'radii': f'patch_icpp({pid})%radii',

                    'e' : f'{math.e}', 'pi': f'{math.pi}',
                }.get(match.group(), match.group())

            lines = []
            # perform the replacement of strings for each analytic function
            # to generate some fortran string representing the code passed in
            for attribute, expr in items:
                if print:
                    cons.print(f"* Codegen: {attribute} = {expr}")

                varname  = re.findall(r"[a-zA-Z][a-zA-Z0-9_]*", attribute)[1]
                qpvf_idx = QPVF_IDX_VARS[varname][:]

                if len(re.findall(r"[0-9]+", attribute)) == 2:
                    idx = int(re.findall(r'[0-9]+', attribute)[1]) - 1
                    qpvf_idx = f"{qpvf_idx} + {idx}"

                lhs = f"q_prim_vf({qpvf_idx})%sf({DATA['sf_idx']})"
                rhs = re.sub(r"[a-zA-Z]+", rhs_replace, expr)

                lines.append(f"        {lhs} = {rhs}")

            # concatenates all of the analytic lines into a single string with
            # each element separated by new line characters. Then write those
            # new lines as a fully concatenated string with fortran syntax
            srcs.append(f"""\
    if (patch_id == {pid}) then
{f'{chr(10)}'.join(lines)}
    end if\
""")

        content = f"""\
! This file was generated by MFC. It is only used when analytical patches are
! present in the input file. It is used to define the analytical patches with
! expressions that are evaluated at runtime from the input file.

#:def analytical()
{f'{chr(10)}{chr(10)}'.join(srcs)}
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

            mapped_weno = 1 if self.params.get("mapped_weno", 'F') == 'T' else 0
            wenoz  = 1 if self.params.get("wenoz", 'F') == 'T' else 0
            teno   = 1 if self.params.get("teno", 'F') == 'T' else 0
            wenojs = 0 if (mapped_weno or wenoz or teno) else 1
            igr = 1 if self.params.get("igr", 'F') == 'T' else 0

            recon_type = self.params.get("recon_type", 1)

            # This fixes a bug on Frontier to do with allocating 0:0 arrays
            weno_order = int(self.params.get("weno_order", 0))
            if recon_type == 1:
                weno_polyn = int((weno_order - 1) / 2)
            else:
                weno_polyn = 1

            if self.params.get("igr", "F") == 'T':
                weno_order = 5
                weno_polyn = 3

            if teno:
                weno_num_stencils = weno_order - 3
            else:
                weno_num_stencils = weno_polyn

            num_dims = 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)
            if self.params.get("mhd", 'F') == 'T':
                num_vels = 3
            else:
                num_vels = num_dims

            mhd = 1 if self.params.get("mhd", 'F') == 'T' else 0
            relativity = 1 if self.params.get("relativity", 'F') == 'T' else 0
            viscous = 1 if self.params.get("viscous", 'F') == 'T' else 0
            igr = 1 if self.params.get("igr", 'F') == 'T' else 0
            igr_pres_lim = 1 if self.params.get("igr_pres_lim", 'F') == 'T' else 0

            # Throw error if wenoz_q is required but not set
            return f"""\
#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
#:set recon_type            = {recon_type}
#:set weno_order            = {weno_order}
#:set weno_polyn            = {weno_polyn}
#:set muscl_order           = {int(self.params.get("muscl_order", 0))}
#:set muscl_polyn           = {int(self.params.get("muscl_order", 0))}
#:set muscl_lim             = {int(self.params.get("muscl_lim", 1))}
#:set weno_num_stencils     = {weno_num_stencils}
#:set nb                    = {int(self.params.get("nb", 1))}
#:set num_dims              = {num_dims}
#:set num_vels              = {num_vels}
#:set nterms                = {nterms}
#:set num_fluids            = {int(self.params["num_fluids"])}
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
"""

        return """\
! This file is purposefully empty. It is only important for builds that make use
! of --case-optimization.
"""

    def get_fpp(self, target, print = True) -> str:
        def _prepend() -> str:
            return f"""\
#:set chemistry             = {self.params.get("chemistry", 'F') == 'T'}
"""

        def _default(_) -> str:
            return "! This file is purposefully empty."

        result = {
            "pre_process" : self.__get_pre_fpp,
            "simulation"  : self.__get_sim_fpp,
        }.get(build.get_target(target).name, _default)(print)

        return _prepend() + result

    def __getitem__(self, key: str) -> str:
        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = val
