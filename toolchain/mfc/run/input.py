import os, re, json, math, dataclasses

from ..printer import cons
from ..        import common, build
from ..state   import ARG, ARGS
from .         import case_dicts

QPVF_IDX_VARS = {
    'alpha_rho': 'contxb', 'vel'  : 'momxb',         'pres': 'E_idx', 
    'alpha':     'advxb',  'tau_e': 'stress_idx%beg'
}

@dataclasses.dataclass
class MFCInputFile:
    filename:     str
    case_dirpath: str
    case_dict:    dict

    def __get_ndims(self) -> int:
        return 1 + min(int(self.case_dict.get("n", 0)), 1) + min(int(self.case_dict.get("p", 0)), 1)

    def __is_ic_analytical(self, key: str, val: str) -> bool:
        if common.is_number(val) or not isinstance(val, str):
            return False

        for array in QPVF_IDX_VARS:
            if re.match(fr'^patch_icpp\([0-9]+\)%{array}', key):
                return True

        return False

    def generate_inp(self, target) -> None:
        target = build.get_target(target)

        cons.print(f"Generating [magenta]{target.name}.inp[/magenta]:")
        cons.indent()

        MASTER_KEYS: list = case_dicts.get_input_dict_keys(target.name)

        ignored = []

        # Create Fortran-style input file content string
        dict_str = ""
        for key, val in self.case_dict.items():
            if key in MASTER_KEYS:
                if self.__is_ic_analytical(key, val):
                    dict_str += f" {key} = 0d0\n"
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

        cons.print(f"[yellow]INFO:[/yellow] Forwarded {len(self.case_dict)-len(ignored)}/{len(self.case_dict)} parameters.")

        contents = f"&user_inputs\n{dict_str}&end/\n"

        # Save .inp input file
        common.file_write(f"{self.case_dirpath}/{target.name}.inp", contents)

        cons.unindent()

    def __save_fpp(self, target, contents: str) -> None:
        inc_dir = os.path.join(target.get_staging_dirpath(), "include", target.name)
        common.create_directory(inc_dir)

        fpp_path = os.path.join(inc_dir, "case.fpp")

        cons.print("Writing a (new) custom case.fpp file.")
        common.file_write(fpp_path, contents, True)

    # pylint: disable=too-many-locals
    def __get_pre_fpp(self, print: bool) -> str:
        DATA = {
            1: {'ptypes': [1, 15, 16],                         'sf_idx': 'i, 0, 0'},
            2: {'ptypes': [2,  3,  4,  5,  6,  7, 17, 18, 21], 'sf_idx': 'i, j, 0'},
            3: {'ptypes': [8,  9, 10, 11, 12, 13, 14, 19, 21], 'sf_idx': 'i, j, k'}
        }[self.__get_ndims()]

        patches = {}

        for key, val in self.case_dict.items():
            if not self.__is_ic_analytical(key, val):
                continue

            patch_id = re.search(r'[0-9]+', key).group(0)

            if patch_id not in patches:
                patches[patch_id] = []

            patches[patch_id].append((key, val))

        srcs = []

        for pid, items in patches.items():
            ptype = self.case_dict[f"patch_icpp({pid})%geometry"]

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

            bubble_model = int(self.case_dict.get("bubble_model", "-100"))

            if bubble_model == 2:
                nterms = 32
            elif bubble_model == 3:
                nterms = 7

            return f"""\
#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
#:set weno_order            = {int(self.case_dict["weno_order"])}
#:set weno_polyn            = {int((self.case_dict["weno_order"] - 1) / 2)}
#:set nb                    = {int(self.case_dict.get("nb", 1))}
#:set num_dims              = {1 + min(int(self.case_dict.get("n", 0)), 1) + min(int(self.case_dict.get("p", 0)), 1)}
#:set nterms                = {nterms}
"""

        return """\
! This file is purposefully empty. It is only important for builds that make use
! of --case-optimization.
"""

    def __get_post_fpp(self, _) -> str:
        return """\
! This file is purposefully empty for all post-process builds.
"""

    def get_fpp(self, target, print = True) -> str:
        def _default(_) -> str:
            return ""

        result = {
            "pre_process"  : self.__get_pre_fpp,
            "simulation"   : self.__get_sim_fpp,
            "post_process" : self.__get_post_fpp,
        }.get(build.get_target(target).name, _default)(print)

        return result

    def generate_fpp(self, target) -> None:
        if target.isDependency:
            return

        cons.print(f"Generating [magenta]case.fpp[/magenta].")
        cons.indent()

        self.__save_fpp(target, self.get_fpp(target))

        cons.unindent()

    # Generate case.fpp & [target.name].inp
    def generate(self, target) -> None:
        self.generate_inp(target)
        cons.print()
        self.generate_fpp(target)


# Load the input file
def load(empty_data: dict = None) -> MFCInputFile:
    if load.CACHED_MFCInputFile is not None:
        return load.CACHED_MFCInputFile

    if not ARG("input"):
        if empty_data is None:
            raise common.MFCException("Please provide an input file.")

        load.CACHED_MFCInputFile = MFCInputFile("empty.py", "empty.py", empty_data)
    else:
        filename: str = ARG("input").strip()

        cons.print(f"Acquiring [bold magenta]{filename}[/bold magenta]...")

        dirpath:    str  = os.path.abspath(os.path.dirname(filename))
        dictionary: dict = {}

        if not os.path.exists(filename):
            raise common.MFCException(f"Input file '{filename}' does not exist. Please check the path is valid.")

        if filename.endswith(".py"):
            (json_str, err) = common.get_py_program_output(filename, [json.dumps(ARGS())] + ARG("arguments"))

            if err != 0:
                raise common.MFCException(f"Input file {filename} terminated with a non-zero exit code. Please make sure running the file doesn't produce any errors.")
        elif filename.endswith(".json"):
            json_str = common.file_read(filename)
        else:
            raise common.MFCException("Unrecognized input file format. Only .py and .json files are supported. Please check the README and sample cases in the examples directory.")

        try:
            dictionary = json.loads(json_str)
        except Exception as exc:
            raise common.MFCException(f"Input file {filename} did not produce valid JSON. It should only print the case dictionary.\n\n{exc}\n")

        load.CACHED_MFCInputFile = MFCInputFile(filename, dirpath, dictionary)

    return load.CACHED_MFCInputFile


load.CACHED_MFCInputFile = None
