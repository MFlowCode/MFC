#!/usr/bin/env python3

import os
import re
import copy
import math
import colorama
import dataclasses

from pathlib import Path

import internal.common      as common
import internal.treeprint   as treeprint


@dataclasses.dataclass
class Case:
    name:       str
    parameters: list

    def __init__(self, data: dict) -> None:
        self.name       = data.get("name")
        self.parameters = {}
        
        for p in list(data.get("parameters").items()):
            self.parameters[p[0]] = p[1]

    def get_keys(self):
        keys = []
        for param in self.parameters:
            keys.append(param.name)
        
        return keys

    def has_parameter(self, key: str):
        return key in self.get_keys()

    def __getitem__(self, key: str) -> str:
        if key not in self.parameters:
            raise common.MFCException(f"Case: Parameter {key} does not exist.")
        
        return self.parameters[key]

    def __setitem__(self, key: str, val: str):
        self.parameters[key] = val

    def create_case_dict_str(self) -> str: 
        result: str = "{\n"

        for key,val in self.parameters.items():
            result = f'{result}\t"{key}": "{val}",\n'

        return result + "}"


@dataclasses.dataclass
class Test:
    case: Case

    def __init__(self, data: dict) -> None:
        self.case = data.get("case", {})


Tend = 0.25
Nt   = 500
mydt = Tend/(1.*Nt)

BASE_CASE = Case({
    "name": "Base Case",
    "parameters": {
        'case_dir'                     : '\'.\'',
        'run_time_info'                : 'F',
        'nodes'                        : 1,
        'ppn'                          : 1,
        'queue'                        : 'normal',
        'walltime'                     : '24:00:00',
        'mail_list'                    : '',
        'x_domain%beg'                 : 0.E+00,   
        'x_domain%end'                 : 1.E+00,   
        'y_domain%beg'                 : 0.E+00,   
        'y_domain%end'                 : 1.E+00,   
        'z_domain%beg'                 : 0.E+00,   
        'z_domain%end'                 : 1.E+00,   
        'm'                            : 0,       
        'n'                            : 0,        
        'p'                            : 0,        
        'dt'                           : mydt,     
        't_step_start'                 : 0,        
        't_step_stop'                  : int(Nt+1),
        't_step_save'                  : int(Nt),  
        'num_patches'                  : 2,     
        'model_eqns'                   : 2,     
        'alt_soundspeed'               : 'F',   
        'num_fluids'                   : 1,     
        'adv_alphan'                   : 'T',   
        'mpp_lim'                      : 'F',   
        'mixture_err'                  : 'F',   
        'time_stepper'                 : 3,     
        'weno_vars'                    : 2,     
        'weno_order'                   : 5,     
        'weno_eps'                     : 1.E-16,
        'mapped_weno'                  : 'T',
        'null_weights'                 : 'F',
        'mp_weno'                      : 'F',
        'riemann_solver'               : 2,
        'wave_speeds'                  : 1,
        'avg_state'                    : 2,
        'bc_x%beg'                     : -3, 
        'bc_x%end'                     : -3,
        'bc_y%beg'                     : -3, 
        'bc_y%end'                     : -3,
        'bc_z%beg'                     : -3, 
        'bc_z%end'                     : -3,
        'format'                       : 1,
        'precision'                    : 2,
        'prim_vars_wrt'                :'T',
        'parallel_io'                  :'F',
        'patch_icpp(1)%geometry'       : 1,
        'patch_icpp(1)%x_centroid'     : 0.25,
        'patch_icpp(1)%y_centroid'     : 0.5,
        'patch_icpp(1)%z_centroid'     : 0.5,
        'patch_icpp(1)%length_x'       : 0.5,
        'patch_icpp(1)%length_y'       : 1,
        'patch_icpp(1)%length_z'       : 1,
        'patch_icpp(1)%vel(1)'         : 0.0,
        'patch_icpp(1)%vel(2)'         : 0.0,
        'patch_icpp(1)%vel(3)'         : 0.0,
        'patch_icpp(1)%pres'           : 1.0,
        'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
        'patch_icpp(1)%alpha(1)'       : 1.,
        'patch_icpp(2)%geometry'       : 1,
        'patch_icpp(2)%x_centroid'     : 0.75,
        'patch_icpp(2)%y_centroid'     : 0.5,
        'patch_icpp(2)%z_centroid'     : 0.5,
        'patch_icpp(2)%length_x'       : 0.5,
        'patch_icpp(2)%length_y'       : 1,
        'patch_icpp(2)%length_z'       : 1,
        'patch_icpp(2)%vel(1)'         : 0.0,
        'patch_icpp(2)%vel(2)'         : 0.0,
        'patch_icpp(2)%vel(3)'         : 0.0,
        'patch_icpp(2)%pres'           : 0.1,
        'patch_icpp(2)%alpha_rho(1)'   : 0.125E+00,
        'patch_icpp(2)%alpha(1)'       : 1.,
        'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
        'fluid_pp(1)%pi_inf'           : 0.0,
    }
})

class MFCTest:
    def __init__(self, bootstrap):
        self.bootstrap = bootstrap

        # Aliases
        self.tree = self.bootstrap.tree
        self.args = self.bootstrap.args

    def test(self):
        self.tree.print(f"Testing mfc")
        self.tree.indent()

        self.text_id = 1
        self.test_acc_packed = ""

        if self.args["generate"]:
            common.delete_directory_recursive_safe(common.MFC_TESTDIR)
            common.create_directory_safe(common.MFC_TESTDIR)

        if not self.bootstrap.is_build_satisfied("mfc"):
            raise common.MFCException(f"Can't test mfc because its build isn't satisfied.")
        
        for dimId, dimParams in enumerate([ #{'patch_icpp(1)%geometry': 1, 'patch_icpp(2)%geometry': 1, 'm': 299},
                                            {'patch_icpp(1)%geometry': 3, 'patch_icpp(2)%geometry': 3, 'm': 49, 'n': 39},
                                            #{'patch_icpp(1)%geometry': 9, 'patch_icpp(2)%geometry': 9, 'm': 39, 'n': 29, 'p': 19}
                                            ]):
            for weno_order in [3, 5]:
                for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'mapped_weno': mapped_weno, 'mp_weno': mp_weno}})

                for riemann_solver in [1, 2]:
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'alt_soundspeed': 'T'}})
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mixture_err':    'T'}})
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mpp_lim':        'T'}})
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'avg_state':      1}})
                    self.handle_case(self.tree, {**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'wave_speeds':    2}})

                # TODO: num_comp

            for ppn in [2, 4]:
                self.handle_case(self.tree, {**dimParams, **{'ppn': ppn}})
    
        common.clear_line()
        self.tree.print(f"Tested. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})")
        self.tree.unindent()


    def get_case_dir_name(self, mods: dict):
        return "".join([f"{str(x[0]).split('_')[0][:4]}-{str(x[1])[:4]}_" for x in mods.items()])[:-1]

    def get_case_dir(self, mods: dict):
        return f"{common.MFC_TESTDIR}/{self.get_case_dir_name(mods)}"

    def create_case_dir(self, mods: dict):
        case     = copy.deepcopy(BASE_CASE)
        case_dir = self.get_case_dir(mods)

        for key, val in mods.items():
            case[key] = val

        content = f"""\
#!/usr/bin/env python3

from pathlib import Path
from sys     import path

path.insert(0, f"{{Path(__file__).parent.resolve()}}/../../src/master_scripts")

# Let Python find MFC's module
from m_python_proxy import f_execute_mfc_component

case_dict = {case.create_case_dict_str()}

f_execute_mfc_component('pre_process', case_dict, '..', 'serial')
f_execute_mfc_component('simulation',  case_dict, '..', 'serial')

"""

        os.makedirs(case_dir, exist_ok=True)

        f = open(f"{case_dir}/input.py", "w")
        f.write(content)
        f.close()

    def golden_file_compare_match(self, truth: str, candidate: str):
        for candidate_line in candidate.splitlines():
            if candidate_line == "":
                continue

            file_subpath: str = candidate_line.split(' ')[0]
            
            line_trusted: str = ""
            for l in truth.splitlines():
                if l.startswith(file_subpath):
                    line_trusted = l
                    break
            
            if len(line_trusted) == 0:
                continue
            
            numbers_cand  = [ float(x) for x in candidate_line.strip().split(' ')[1:] ]
            numbers_trust = [ float(x) for x in line_trusted.strip().split(' ')[1:]   ]

            # Different amount of spaces, means that there are more entires in one than in the other
            if len(numbers_cand) != len(numbers_trust):
                return (False, "Variable count didn't match.")

            # check values one by one
            for i in range(len(numbers_cand)):
                # FIXME: set abs_tol
                if not math.isclose(numbers_cand[i], numbers_trust[i], rel_tol=1e-13):
                    abs_delta    = abs(numbers_cand[i]-numbers_trust[i])
                    percent_diff = abs_delta/numbers_trust[i]
                    return (False, f"Error margin is too high for value at position #{i+1} in {file_subpath}: ~{round(percent_diff, 5)}% (~{round(abs_delta, 5)}).")

        # Both tests gave the same results within an acceptable tolerance
        return (True, "")

    def handle_case(self, tree: treeprint.TreePrinter, parameters: dict):
        global text_id

        self.create_case_dir(parameters)
        
        tree.print_progress(f"Running test #{self.text_id} - {self.get_case_dir_name(parameters)}", self.text_id, 26)

        def on_test_errror(msg: str = ""):
            common.clear_line()
            tree.print(f"Test #{self.text_id}: Failed! ({colorama.Fore.RED}FAILURE{colorama.Style.RESET_ALL})")
            if msg != "":
                tree.print(msg)
            tree.print(f"The test is available at: {self.get_case_dir(parameters)}/input.py")
            raise common.MFCException("Testing failed (view above).")

        common.execute_shell_command_safe(f"cd '{self.get_case_dir(parameters)}' && python3 input.py >> '../test.log' 2>&1", on_error=lambda: on_test_errror("MFC Execution Failed. Please refer to tests/test.log"))

        pack = self.pack_case_output(parameters)
        
        golden_filepath = f"{self.get_case_dir(parameters)}/golden.txt"

        if self.args["generate"]:
            common.delete_file_safe(golden_filepath)
            with open(golden_filepath, "w") as f:
                f.write(pack)
        
        if not os.path.isfile(golden_filepath):
            common.clear_line()
            tree.print(f"Test #{self.text_id}: Golden file doesn't exist! ({colorama.Fore.RED}FAILURE{colorama.Style.RESET_ALL})")
            tree.print(f"To generate golden files, use the '-g' flag.")
            on_test_errror()

        with open(golden_filepath, "r") as f:                
            bSuccess, errorMsg = self.golden_file_compare_match(f.read(), pack)
            if not bSuccess:
                on_test_errror(errorMsg)

        self.text_id+=1

    def pack_case_output(self, params: dict):
        result: str = ""

        case_dir = self.get_case_dir(params)
        D_dir    = f"{case_dir}/D/"

        for filepath in list(Path(D_dir).rglob("*.dat")):
            short_filepath = str(filepath).replace(f'{case_dir}/', '')
            with open(filepath, "r") as file:
                result += f"{short_filepath} " + re.sub(r' +', ' ', file.read().replace('\n', ' ')).strip() + '\n'
        
        return result
