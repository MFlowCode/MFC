import os, glob, hashlib, binascii, subprocess, itertools, dataclasses, shutil

from typing import List, Set, Union, Callable, Optional

from ..      import case, common
from ..state import ARG
from ..run   import input
from ..build import MFCTarget, get_target

Tend = 0.25
Nt   = 50
mydt = 0.0005

BASE_CFG = {
    'run_time_info'                : 'T',
    'm'                            : 0,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : int(Nt),
    't_step_save'                  : int(Nt),
    'num_patches'                  : 3,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
    'recon_type'                   : 1,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'format'                       : 1,
    'precision'                    : 2,

    'patch_icpp(1)%pres'           : 1.0,
    'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
    'patch_icpp(1)%alpha(1)'       : 1.,

    'patch_icpp(2)%pres'           : 0.5,
    'patch_icpp(2)%alpha_rho(1)'   : 0.5,
    'patch_icpp(2)%alpha(1)'       : 1.,

    'patch_icpp(3)%pres'           : 0.1,
    'patch_icpp(3)%alpha_rho(1)'   : 0.125,
    'patch_icpp(3)%alpha(1)'       : 1.,

    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    'fluid_pp(1)%cv'               : 0.0,
    'fluid_pp(1)%qv'               : 0.0,
    'fluid_pp(1)%qvp'              : 0.0,   
    'bubbles_euler'                 : 'F',
    'Ca'                            : 0.9769178386380458,
    'Web'                           : 13.927835051546392,
    'Re_inv'                        : 0.009954269975623245,
    'pref'                          : 101325.0,
    'rhoref'                        : 1000.0,
    'bubble_model'                  :  3,
    'polytropic'                    : 'T',
    'polydisperse'                  : 'F',
    'thermal'                       :  3,
    'R0ref'                         : 1e-05,
    'patch_icpp(1)%r0'              :  1,
    'patch_icpp(1)%v0'              :  0,
    'patch_icpp(2)%r0'              :  1,
    'patch_icpp(2)%v0'              :  0,
    'patch_icpp(3)%r0'              :  1,
    'patch_icpp(3)%v0'              :  0,

    'qbmm'                          : 'F',
    'dist_type'                     : 2,
    'poly_sigma'                    : 0.3,
    'sigR'                          : 0.1,
    'sigV'                          : 0.1,
    'rhoRV'                         : 0.0,

    'acoustic_source'                   : 'F',
    'num_source'                        : 1,
    'acoustic(1)%loc(1)'                : 0.5,
    'acoustic(1)%mag'                   : 0.2,
    'acoustic(1)%length'                : 0.25,
    'acoustic(1)%dir'                   : 1.0,
    'acoustic(1)%npulse'                : 1,
    'acoustic(1)%pulse'                 : 1,
    'rdma_mpi'                          : 'F',

    'bubbles_lagrange'                 : 'F',
    'lag_params%nBubs_glb'             : 1,
    'lag_params%solver_approach'       : 0,
    'lag_params%cluster_type'          : 2,
    'lag_params%pressure_corrector'    : 'F',
    'lag_params%smooth_type'           : 1,
    'lag_params%epsilonb'              : 1.0,
    'lag_params%heatTransfer_model'    : 'F',
    'lag_params%massTransfer_model'    : 'F',
    'lag_params%valmaxvoid'            : 0.9,
    'lag_params%c0'                    : 10.1,
    'lag_params%rho0'                  : 1000.,
    'lag_params%T0'                    : 298.,
    'lag_params%x0'                    : 1.,
    'lag_params%diffcoefvap'           : 2.5e-5,
    'lag_params%Thost'                 : 298.,
}

def trace_to_uuid(trace: str) -> str:
    return hex(binascii.crc32(hashlib.sha1(str(trace).encode()).digest())).upper()[2:].zfill(8)

@dataclasses.dataclass(init=False)
class TestCase(case.Case):
    ppn:          int
    trace:        str
    override_tol: Optional[float] = None

    def __init__(self, trace: str, mods: dict, ppn: int = None, override_tol: float = None) -> None:
        self.trace        = trace
        self.ppn          = ppn or 1
        self.override_tol = override_tol
        super().__init__({**BASE_CFG.copy(), **mods})

    def run(self, targets: List[Union[str, MFCTarget]], gpus: Set[int]) -> subprocess.CompletedProcess:
        if gpus is not None and len(gpus) != 0:
            gpus_select = ["--gpus"] + [str(_) for _ in gpus]
        else:
            gpus_select = []

        filepath          = f'{self.get_dirpath()}/case.py'
        tasks             = ["-n", str(self.ppn)]
        jobs              = ["-j", str(ARG("jobs"))] if ARG("case_optimization") else []
        case_optimization = ["--case-optimization"]  if ARG("case_optimization") else []

        if self.params.get("bubbles_lagrange", 'F') == 'T':
            input_bubbles_lagrange(self)

        mfc_script = ".\\mfc.bat" if os.name == 'nt' else "./mfc.sh"

        target_names = [ get_target(t).name for t in targets ]

        command = [
            mfc_script, "run", filepath, "--no-build", *tasks, *case_optimization,
            *jobs, "-t", *target_names, *gpus_select, *ARG("--")
        ]

        return common.system(command, print_cmd=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def get_trace(self) -> str:
        return self.trace

    def get_uuid(self) -> str:
        return trace_to_uuid(self.trace)

    def get_dirpath(self):
        return os.path.join(common.MFC_TEST_DIR, self.get_uuid())

    def get_filepath(self):
        filepath = os.path.join(self.get_dirpath(), "case.py")
        if os.name == 'nt':
            return filepath.replace('\\', '\\\\')
        return filepath

    def delete_output(self):
        dirpath = self.get_dirpath()

        exts = ["*.inp", "*.1", "*.dat", "*.inf", "*.sh", "*.txt"]
        for f in list(itertools.chain.from_iterable(glob.glob(os.path.join(dirpath, ext)) for ext in exts)):
            if "golden" in f:
                continue

            common.delete_file(f)

        common.delete_directory(os.path.join(dirpath, "D"))
        common.delete_directory(os.path.join(dirpath, "p_all"))
        common.delete_directory(os.path.join(dirpath, "silo_hdf5"))
        common.delete_directory(os.path.join(dirpath, "restart_data"))
        if self.params.get("bubbles_lagrange", 'F') == 'T':
            common.delete_directory(os.path.join(dirpath, "input"))
            common.delete_directory(os.path.join(dirpath, "lag_bubbles_post_process"))

        for f in ["pack", "pre_process", "simulation", "post_process"]:
            common.delete_file(os.path.join(dirpath, f"{f}.txt"))

    def create_directory(self):
        dirpath = self.get_dirpath()

        common.create_directory(dirpath)

        common.file_write(self.get_filepath(), f"""\
#!/usr/bin/env python3
#
# {self.get_filepath()}:
# {self.trace}

import json
import argparse

parser = argparse.ArgumentParser(
    prog="{self.get_filepath()}",
    description="{self.get_filepath()}: {self.trace}",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{{}}', metavar="DICT",
                    help="MFC's toolchain's internal state.")

ARGS = vars(parser.parse_args())

case = {self.gen_json_dict_str()}
mods = {{}}

if "post_process" in ARGS["mfc"]["targets"]:
    mods = {{
        'parallel_io'  : 'T', 'cons_vars_wrt'   : 'T',
        'prim_vars_wrt': 'T', 'alpha_rho_wrt(1)': 'T',
        'rho_wrt'      : 'T', 'mom_wrt(1)'      : 'T',
        'vel_wrt(1)'   : 'T', 'E_wrt'           : 'T',
        'pres_wrt'     : 'T', 'alpha_wrt(1)'    : 'T',
        'gamma_wrt'    : 'T', 'heat_ratio_wrt'  : 'T',
        'pi_inf_wrt'   : 'T', 'pres_inf_wrt'    : 'T',
        'c_wrt'        : 'T',
    }}

    if case['p'] != 0:
        mods['fd_order']  = 1
        mods['omega_wrt(1)'] = 'T'
        mods['omega_wrt(2)'] = 'T'
        mods['omega_wrt(3)'] = 'T'
else:
    mods['parallel_io']   = 'F'
    mods['prim_vars_wrt'] = 'F'

print(json.dumps({{**case, **mods}}))
""")

    def __str__(self) -> str:
        return f"tests/[bold magenta]{self.get_uuid()}[/bold magenta]: {self.trace}"

    def to_input_file(self) -> input.MFCInputFile:
        return input.MFCInputFile(
            os.path.basename(self.get_filepath()),
            self.get_dirpath(),
            self.get_parameters())

    # pylint: disable=too-many-return-statements
    def compute_tolerance(self) -> float:
        if self.override_tol:
            return self.override_tol

        tolerance = 1e-12  # Default
        single = ARG("single")

        if "Example" in self.trace.split(" -> "):
            tolerance = 1e-3
        elif "Cylindrical" in self.trace.split(" -> "):
            tolerance = 1e-9
        elif self.params.get("hypoelasticity", 'F') == 'T':
            tolerance = 1e-7
        elif self.params.get("mixlayer_perturb", 'F') == 'T':
            tolerance = 1e-7
        elif any(self.params.get(key, 'F') == 'T' for key in ['relax', 'ib', 'qbmm', 'bubbles_euler', 'bubbles_lagrange']):
            tolerance = 1e-10
        elif self.params.get("low_Mach") in [1, 2]:
            tolerance = 1e-10
        elif self.params.get("acoustic_source", 'F') == 'T':
            if self.params.get("acoustic(1)%pulse") == 3:  # Square wave
                return 1e-1 if single else 1e-5
            tolerance = 3e-12
        elif self.params.get("weno_order") == 7:
            tolerance = 1e-9
        elif self.params.get("mhd", 'F') == 'T':
            tolerance = 1e-8

        return 1e8 * tolerance if single else tolerance

@dataclasses.dataclass
class TestCaseBuilder:
    trace:        str
    mods:         dict
    path:         str
    args:         List[str]
    ppn:          int
    functor:      Optional[Callable]
    override_tol: Optional[float] = None

    def get_uuid(self) -> str:
        return trace_to_uuid(self.trace)

    def to_case(self) -> TestCase:
        dictionary = {}
        if self.path:
            dictionary.update(input.load(self.path, self.args, do_print=False).params)

            for key, value in dictionary.items():
                if not isinstance(value, str):
                    continue

                for path in [value, os.path.join(os.path.dirname(self.path), value)]:
                    path = os.path.abspath(path)
                    if os.path.exists(path):
                        dictionary[key] = path

        if self.mods:
            dictionary.update(self.mods)

        if self.functor:
            self.functor(dictionary)

        return TestCase(self.trace, dictionary, self.ppn, self.override_tol)


@dataclasses.dataclass
class CaseGeneratorStack:
    trace: list # list of strs
    mods:  list # list of dicts

    def __init__(self) -> None:
        self.trace, self.mods = [], []

    def size(self) -> int:
        return len(self.trace)

    def push(self, trace: str, mods: dict) -> None:
        self.trace.append(trace)
        self.mods.append(mods)

    def pop(self) -> None:
        return (self.mods.pop(), self.trace.pop())


# pylint: disable=too-many-arguments, too-many-positional-arguments
def define_case_f(trace: str, path: str, args: List[str] = None, ppn: int = None, mods: dict = None, functor: Callable = None, override_tol: float = None) -> TestCaseBuilder:
    return TestCaseBuilder(trace, mods or {}, path, args or [], ppn or 1, functor, override_tol)


# pylint: disable=too-many-arguments, too-many-positional-arguments
def define_case_d(stack: CaseGeneratorStack, newTrace: str, newMods: dict, ppn: int = None, functor: Callable = None, override_tol: float = None) -> TestCaseBuilder:
    mods: dict = {}

    for mod in stack.mods:
        mods.update(mod)

    mods.update(newMods)

    if isinstance(newTrace, str):
        newTrace = [newTrace]

    traces: list = []
    for trace in stack.trace[:] + newTrace:
        if not common.isspace(trace):
            traces.append(trace)

    return TestCaseBuilder(' -> '.join(traces), mods, None, None, ppn or 1, functor, override_tol)

def input_bubbles_lagrange(self):
    if "lagrange_bubblescreen" in self.trace:
        copy_input_lagrange(f'/3D_lagrange_bubblescreen',f'{self.get_dirpath()}')
    elif "lagrange_shbubcollapse" in self.trace:
        copy_input_lagrange(f'/3D_lagrange_shbubcollapse',f'{self.get_dirpath()}')
    else:
        create_input_lagrange(f'{self.get_dirpath()}')

def create_input_lagrange(path_test):
    folder_path_lagrange = path_test + '/input'
    file_path_lagrange = folder_path_lagrange + '/lag_bubbles.dat'
    if not os.path.exists(folder_path_lagrange):
        os.mkdir(folder_path_lagrange)

    with open(file_path_lagrange, "w") as file:
        file.write('0.5\t0.5\t0.5\t0.0\t0.0\t0.0\t8.0e-03\t0.0')

def copy_input_lagrange(path_example_input, path_test):
    folder_path_dest = path_test + '/input/'
    fite_path_dest = folder_path_dest + 'lag_bubbles.dat'
    file_path_src = common.MFC_EXAMPLE_DIRPATH + path_example_input + '/input/lag_bubbles.dat'
    if not os.path.exists(folder_path_dest):
        os.mkdir(folder_path_dest)

    shutil.copyfile(file_path_src, fite_path_dest)
