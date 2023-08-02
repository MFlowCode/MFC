import os, re, json, math, dataclasses

from ..printer import cons
from ..        import common
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
        
        for array in QPVF_IDX_VARS.keys():
            if re.match(fr'^patch_icpp\([0-9]+\)%{array}', key):
                return True
        
        return False

    def __generate_inp(self, target_name: str) -> None:
        cons.print(f"Generating [magenta]{target_name}.inp[/magenta].")
        cons.indent()
        
        MASTER_KEYS: list = case_dicts.get_input_dict_keys(target_name)

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
                continue
            else:
                ignored.append(key)

            if key not in case_dicts.ALL:
                raise common.MFCException(f"MFCInputFile::dump: Case parameter '{key}' is not used by any MFC code. Please check your spelling or add it as a new parameter.")

        cons.print(f"[yellow]INFO:[/yellow] Forwarded {len(self.case_dict)-len(ignored)}/{len(self.case_dict)} parameters.")
        
        contents = f"&user_inputs\n{dict_str}&end/\n"

        # Save .inp input file
        common.file_write(f"{self.case_dirpath}/{target_name}.inp", contents)
        
        cons.unindent()
        
    def __save_fpp(self, target_name: str, contents: str) -> None:
        filepath = os.path.join(os.getcwd(), "src", target_name, "include", "case.fpp")

        # Check if this case already has a case.fpp file.
        # If so, we don't need to generate a new one, which
        # would cause a partial and unnecessary rebuild.
        if os.path.exists(filepath):
            with open(filepath, "r") as f:
                lhs = [ l.strip() for l in f.read().splitlines() if not common.isspace(l) ]
                rhs = [ l.strip() for l in contents.splitlines() if not common.isspace(l) ]

                if lhs == rhs:
                    cons.print("[yellow]INFO:[/yellow] Existing case.fpp file is up to date.")
                    return

        cons.print("[yellow]INFO:[/yellow] Overwriting existing case.fpp file. This will cause a partial rebuild.")
        common.file_write(filepath, contents)            

    def __generate_pre_fpp(self) -> None:
        cons.print(f"Generating [magenta]pre_process/include/case.fpp[/magenta].")
        cons.indent()
        
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
                varname         = re.findall(r"[a-zA-Z][a-zA-Z0-9_]*", attribute)[1]
                qpvf_idx_var    = QPVF_IDX_VARS[varname]
                qpvf_idx_offset = ""
                
                if len(re.findall(r"[0-9]+", attribute)) == 2:
                    idx = int(re.findall(r'[0-9]+', attribute)[1]) - 1
                    if idx != 0:
                        qpvf_idx_offset = f" + {idx}"
                
                sf_idx = DATA['sf_idx']
                
                cons.print(f"[yellow]INFO:[/yellow] {self.__get_ndims()}D Analytical Patch #{pid}: Code generation for [magenta]{varname}[/magenta]...")
                
                lhs = f"q_prim_vf({qpvf_idx_var}{qpvf_idx_offset})%sf({sf_idx})"
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
        
        self.__save_fpp("pre_process", content)
        
        cons.unindent()
        
    def __generate_sim_fpp(self) -> None:
        cons.print(f"Generating [magenta]simulation/include/case.fpp[/magenta].")
        cons.indent()

        content = f"""\
! This file was generated by MFC. It is only used if the --case-optimization
! option is passed to ./mfc.sh run or test, enabling a GPU-oriented optimization
! that hard-codes certain case parameters from the input file.

#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
"""

        if ARG("case_optimization"):
            cons.print("[yellow]INFO:[/yellow] Case optimization is enabled.")
            
            nterms = -100

            bubble_model = int(self.case_dict.get("bubble_model", "-100"))

            if bubble_model == 2:
                nterms = 32
            elif bubble_model == 3:
                nterms = 7

            content = content + f"""
#:set weno_order = {int(self.case_dict["weno_order"])}
#:set weno_polyn = {int((self.case_dict["weno_order"] - 1) / 2)}
#:set nb         = {int(self.case_dict.get("nb", 1))}
#:set num_dims   = {1 + min(int(self.case_dict.get("n", 0)), 1) + min(int(self.case_dict.get("p", 0)), 1)}
#:set nterms     = {nterms}
"""
        else:
            cons.print("[yellow]INFO:[/yellow] Case optimization is disabled. Use --case-optimization to enable it.")

        self.__save_fpp("simulation", content)
        cons.unindent()

    def __generate_post_fpp(self) -> None:
        cons.print("Generating [magenta]post_process/include/case.fpp[/magenta].")
        cons.indent()
        cons.print("[yellow]INFO:[/yellow] No case.fpp file is generated for post_process.")
        cons.unindent()
        pass

    # Generate case.fpp & [target_name].inp
    def generate(self, target_name: str, bOnlyFPPs = False) -> None:
        if not bOnlyFPPs:
            self.__generate_inp(target_name)
        
        cons.print()
       
        def _default():
            cons.print(f"No additional input file generation needed for [bold magenta]{target_name}[/bold magenta].")

        {
            "pre_process"  : self.__generate_pre_fpp,
            "simulation"   : self.__generate_sim_fpp,
            "post_process" : self.__generate_post_fpp,
        }.get(target_name, _default)()
        

# Load the input file
def load() -> MFCInputFile:
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

    return MFCInputFile(filename, dirpath, dictionary)
