#!/usr/bin/env python3

import copy
import json
import colorama

from os import makedirs, mkdir, rmdir, system

import internal.common    as common
import internal.treeprint as treeprint

Nx   = 299
dx   = 1./(1.*(Nx+1))
Tend = 0.25
Nt   = 2500
mydt = Tend/(1.*Nt)

BASE_CASE = {
    # Logistics ================================================
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'F',
    'nodes'                        : 1,
    'ppn'                          : 1,
    'queue'                        : 'normal',
    'walltime'                     : '24:00:00',
    'mail_list'                    : '',
    # ==========================================================
    # Computational Domain Parameters ==========================
    'x_domain%beg'                 : 0.E+00,   
    'x_domain%end'                 : 1.E+00,   
    'm'                            : Nx,       
    'n'                            : 0,        
    'p'                            : 0,        
    'dt'                           : mydt,     
    't_step_start'                 : 0,        
    't_step_stop'                  : int(Nt+1),
    't_step_save'                  : int(Nt),  
    # ==========================================================
    # Simulation Algorithm Parameters ==========================
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
    # ==========================================================
    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 1,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    # ==========================================================
    # Patch 1 L ================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.25,
    'patch_icpp(1)%length_x'       : 0.5,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%pres'           : 1.0,
    'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # 'patch_icpp(1)%alpha_rho(2)'   : 0,
    # 'patch_icpp(1)%alpha(2)'       : 0.   
    # ==========================================================
    # Patch 2 R ================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 0.75,
    'patch_icpp(2)%length_x'       : 0.5,
    'patch_icpp(2)%vel(1)'         : 0.0,
    'patch_icpp(2)%pres'           : 0.1,
    'patch_icpp(2)%alpha_rho(1)'   : 0.125E+00,
    'patch_icpp(2)%alpha(1)'       : 1.,
    # 'patch_icpp(2)%alpha_rho(2)'   : 0.125E+00,
    # 'patch_icpp(2)%alpha(2)'       : 1.   
    # ==========================================================
    # Fluids Physical Parameters ===============================
    # Surrounding liquid
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # 'fluid_pp(2)%gamma'            : 1.E+00/(1.4-1.E+00),
    # 'fluid_pp(2)%pi_inf'           : 0.0,
    # ==========================================================
}

def get_case_dir(mods: dict):
    case_dir = f"{common.MFC_TESTDIR}/"

    for key, val in mods.items():
        if key not in BASE_CASE:
            raise RuntimeError(f"Case key {key} does not exist in the BASE_CASE.")

        case_dir+=str(key).upper()[0] + str(val)[0]

    return case_dir

def create_case(mods: dict):
    case     = copy.deepcopy(BASE_CASE)
    case_dir = get_case_dir(mods)

    for key, val in mods.items():
        if key not in BASE_CASE:
            raise RuntimeError(f"Case key {key} does not exist in the BASE_CASE.")

        case[key]=val

    content = f"""\
#!/usr/bin/env python3

import math

from pathlib import Path
from os      import chdir
from os.path import dirname
from sys     import argv, path

path.insert(0, f"{{Path(__file__).parent.resolve()}}/../../src/master_scripts")

from m_python_proxy import f_execute_mfc_component    

case_dict = {str(json.dumps(case, indent=4, sort_keys=True))}

f_execute_mfc_component('pre_process', case_dict, '..', 'serial')
f_execute_mfc_component('simulation',  case_dict, '..', 'serial')

"""

    makedirs(case_dir, exist_ok=True)

    f = open(f"{case_dir}/input.py", "w")
    f.write(content)
    f.close()


text_id = 1

def handle_case(tree: treeprint.TreePrinter, mods: dict):
    global text_id

    create_case(mods)
    
    tree.print_progress(f"Running test #{text_id}", text_id, 26)

    common.execute_shell_command_safe(f"cd '{get_case_dir(mods)}' && python3 input.py >> '../test.log' 2>&1")

    text_id+=1

def test(tree: treeprint.TreePrinter):
    tree.print("Testing MFC")
    tree.indent()

    # TODO: 1d, 2d, 3d

    for weno_order in [3, 5]:
        for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
            handle_case(tree, {'weno_order': weno_order, 'mapped_weno': mapped_weno, 'mp_weno': mp_weno})

        for riemann_solver in [1, 2]:
            handle_case(tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'alt_soundspeed': 'T'})
            handle_case(tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mixture_err':    'T'})
            handle_case(tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mpp_lim':        'T'})
            handle_case(tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'avg_state':      1})
            handle_case(tree, {'weno_order': weno_order, 'riemann_solver': riemann_solver, 'wave_speeds':    2})

            # TODO: num_comp

        # TODO: mpi_rank

    tree.print(f"Tested. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})")
    tree.unindent()
