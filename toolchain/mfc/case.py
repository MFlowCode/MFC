import re, json, math, copy, dataclasses, jsonschema
import jsonschema.exceptions

from . import common
from . import build
from .printer import cons

from .state import ARG
from .run   import case_dicts

QPVF_IDX_VARS = {
    'alpha_rho': 'contxb', 'vel'  : 'momxb',         'pres': 'E_idx', 
    'alpha':     'advxb',  'tau_e': 'stress_idx%beg', 'Y':   'chemxb',
    'cf_val': 'c_idx'
}

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
            if key in MASTER_KEYS:
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
            jsonschema.validate(self.params, case_dicts.SCHEMA)
        except jsonschema.exceptions.ValidationError as e:
            exception_txt = f"Case parameter '{e.path[0]}' is of the wrong type. Expected type: '{e.schema['type']}'. Got value '{e.instance}' of type '{type(e.instance).__name__}'."
            if origin_txt:
                exception_txt = f"Origin: {origin_txt}. {exception_txt}"

            raise common.MFCException(exception_txt)

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
        DATA = {
            1: {'ptypes': [1, 15, 16],                         'sf_idx': 'i, 0, 0'},
            2: {'ptypes': [2,  3,  4,  5,  6,  7, 17, 18, 21], 'sf_idx': 'i, j, 0'},
            3: {'ptypes': [8,  9, 10, 11, 12, 13, 14, 19, 21], 'sf_idx': 'i, j, k'}
        }[self.__get_ndims()]

        patches = {}

        for key, val in self.params.items():
            if not self.__is_ic_analytical(key, val):
                continue

            patch_id = re.search(r'[0-9]+', key).group(0)

            if patch_id not in patches:
                patches[patch_id] = []

            patches[patch_id].append((key, val))

        srcs = []

        for pid, items in patches.items():
            ptype = self.params[f"patch_icpp({pid})%geometry"]

            if ptype not in DATA['ptypes']:
                raise common.MFCException(f"Patch #{pid} of type {ptype} cannot be analytically defined.")

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

            return f"""\
#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
#:set weno_order            = {int(self.params["weno_order"])}
#:set weno_polyn            = {int((self.params["weno_order"] - 1) / 2)}
#:set nb                    = {int(self.params.get("nb", 1))}
#:set num_dims              = {1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)}
#:set nterms                = {nterms}
#:set num_fluids            = {int(self.params["num_fluids"])}
#:set wenojs                = {wenojs}
#:set mapped_weno           = {mapped_weno}
#:set wenoz                 = {wenoz}
#:set teno                  = {teno}
"""

        return """\
! This file is purposefully empty. It is only important for builds that make use
! of --case-optimization.
"""

    def get_fpp(self, target, print = True) -> str:
        def _default(_) -> str:
            return "! This file is purposefully empty."

        result = {
            "pre_process" : self.__get_pre_fpp,
            "simulation"  : self.__get_sim_fpp,
        }.get(build.get_target(target).name, _default)(print)

        return result

    def __getitem__(self, key: str) -> str:
        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = val
