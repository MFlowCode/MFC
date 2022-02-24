#!/usr/bin/env python3

import os
import copy
import hashlib
import colorama
import random
import dataclasses

import internal.bootstrap   as boot
import internal.common      as common
import internal.treeprint   as treeprint
import internal.configfiles as configfiles


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

    def create_case_dict(self) -> str: 
        result: str = "{\n"

        for key,val in self.parameters.items():
            result = f'{result}"{key}": "{val}",\n'

        return result + "}"


@dataclasses.dataclass
class Test:
    case: Case

    def __init__(self, data: dict) -> None:
        self.case = data.get("case", {})


Nx   = 299
dx   = 1./(1.*(Nx+1))
Tend = 0.25
Nt   = 2500
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
        'm'                            : Nx,       
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
        'format'                       : 1,
        'precision'                    : 1,
        'prim_vars_wrt'                :'T',
        'parallel_io'                  :'F',
        'patch_icpp(1)%geometry'       : 1,
        'patch_icpp(1)%x_centroid'     : 0.25,
        'patch_icpp(1)%length_x'       : 0.5,
        'patch_icpp(1)%vel(1)'         : 0.0,
        'patch_icpp(1)%pres'           : 1.0,
        'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
        'patch_icpp(1)%alpha(1)'       : 1.,
        'patch_icpp(2)%geometry'       : 1,
        'patch_icpp(2)%x_centroid'     : 0.75,
        'patch_icpp(2)%length_x'       : 0.5,
        'patch_icpp(2)%vel(1)'         : 0.0,
        'patch_icpp(2)%pres'           : 0.1,
        'patch_icpp(2)%alpha_rho(1)'   : 0.125E+00,
        'patch_icpp(2)%alpha(1)'       : 1.,
        'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
        'fluid_pp(1)%pi_inf'           : 0.0,
    }
})


class MFCTest(configfiles.ConfigFileBase):
    def __init__(self, bootstrap):
        super().__init__(common.MFC_TEST_FILEPATH, {
            "base": dataclasses.asdict(BASE_CASE)
        })

        self.base      = self.tree_get("base", {})
        self.tests     = [ Test(e) for e in self.tree_get("tests", []) ]
        self.bootstrap = bootstrap

        # Aliases
        self.tree = self.bootstrap.tree

    def test(self):
        self.tree.print(f"Testing mfc")
        self.tree.indent()

        if not self.bootstrap.is_build_satisfied("mfc", ignoreIfSource=True):
            raise common.MFCException(f"Can't test mfc because its build isn't satisfied.")

        # TODO: 1d, 2d, 3d

        for weno_order in [3, 5]:
            for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                handle_case(self.tree, {'weno_order': weno_order, 'mapped_weno': mapped_weno, 'mp_weno': mp_weno})

            for riemann_solver in [1, 2]:
                handle_case(self.tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'alt_soundspeed': 'T'})
                handle_case(self.tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mixture_err':    'T'})
                handle_case(self.tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mpp_lim':        'T'})
                handle_case(self.tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'avg_state':      1})
                handle_case(self.tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'wave_speeds':    2})

                # TODO: num_comp

            # TODO: mpi_rank

        self.tree.print(f"Tested. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})")
        self.tree.unindent()


def get_case_dir(mods: dict):
    case_dir = f"{common.MFC_TESTDIR}/"

    return case_dir + hashlib.md5(str(mods).encode()).hexdigest()

def create_case_dir(mods: dict):
    case     = copy.deepcopy(BASE_CASE)
    case_dir = get_case_dir(mods)

    common.delete_directory_recursive_safe(case_dir)

    for key, val in mods.items():
        case[key] = val

    content = f"""\
#!/usr/bin/env python3

import math

from pathlib import Path
from os      import chdir
from os.path import dirname
from sys     import argv, path

path.insert(0, f"{{Path(__file__).parent.resolve()}}/../../src/master_scripts")

from m_python_proxy import f_execute_mfc_component    

case_dict = {case.create_case_dict()}

f_execute_mfc_component('pre_process', case_dict, '..', 'serial')
f_execute_mfc_component('simulation',  case_dict, '..', 'serial')

"""

    os.makedirs(case_dir, exist_ok=True)

    f = open(f"{case_dir}/input.py", "w")
    f.write(content)
    f.close()


text_id = 1

def handle_case(tree: treeprint.TreePrinter, parameters: dict):
    global text_id

    create_case_dir(parameters)
    
    tree.print_progress(f"Running test #{text_id}", text_id, 26)

    common.execute_shell_command_safe(f"cd '{get_case_dir(parameters)}' && python3 input.py >> '../test.log' 2>&1")

    text_id+=1
