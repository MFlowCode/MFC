import os, glob, typing, hashlib, binascii, subprocess, itertools, dataclasses

from ..      import case, common
from ..state import ARG
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
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 3,
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
    'prim_vars_wrt'                :'F',
    'parallel_io'                  :'F',

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
    'bubbles'                       : 'F',
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
    'R0_type'                       : 1,
    'sigR'                          : 0.1,
    'sigV'                          : 0.1,
    'rhoRV'                         : 0.0,

    'Monopole'                      : 'F',
    'num_mono'                      : 1,
    'Mono(1)%loc(1)'                : 0.5,
    'Mono(1)%mag'                   : 1.0,
    'Mono(1)%length'                : 0.25,
    'Mono(1)%dir'                   : 1.0,
    'Mono(1)%npulse'                : 1,
    'Mono(1)%pulse'                 : 1,
    'cu_mpi'                        :'F',
}

@dataclasses.dataclass(init=False)
class TestCase(case.Case):
    ppn:   int
    trace: str
    opt:   bool

    def __init__(self, trace: str, mods: dict, ppn: int = None, opt: bool = None) -> None:
        self.trace = trace
        self.ppn   = ppn or 1
        self.opt   = opt or False
        super().__init__({**BASE_CFG.copy(), **mods})

    def run(self, targets: typing.List[typing.Union[str, MFCTarget]], gpus: typing.Set[int]) -> subprocess.CompletedProcess:
        if gpus is not None and len(gpus) != 0:
            gpus_select = ["--gpus"] + [str(_) for _ in gpus]
        else:
            gpus_select = []

        filepath          = f'{self.get_dirpath()}/case.py'
        tasks             = ["-n", str(self.ppn)]
        jobs              = ["-j", str(ARG("jobs"))] if ARG("case_optimization") else []
        case_optimization = ["--case-optimization"] if ARG("case_optimization") or self.opt else ["--no-build"]

        mfc_script = ".\\mfc.bat" if os.name == 'nt' else "./mfc.sh"

        target_names = [ get_target(t).name for t in targets ]

        command = [
            mfc_script, "run", filepath, *tasks, *case_optimization,
            *jobs, "-t", *target_names, *gpus_select, *ARG("--")
        ]

        return common.system(command, print_cmd=False, text=True, capture_output=True)

    def get_uuid(self) -> str:
        return hex(binascii.crc32(hashlib.sha1(str(self.trace).encode()).digest())).upper()[2:].zfill(8)

    def get_dirpath(self):
        return os.path.join(common.MFC_TESTDIR, self.get_uuid())

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

        for f in ["pack", "pre_process", "simulation", "post_process"]:
            common.delete_file(os.path.join(dirpath, f"{f}.txt"))

    def create_directory(self):
        dirpath = self.get_dirpath()

        common.create_directory(dirpath)

        common.file_write(f"{dirpath}/case.py", f"""\
#!/usr/bin/env python3
#
# tests/{self.get_uuid()}/case.py:
# {self.trace}

import json
import argparse

parser = argparse.ArgumentParser(
    prog="tests/{self.get_uuid()}/case.py",
    description="tests/{self.get_uuid()}/case.py: {self.trace}",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
parser.add_argument("dict", type=str, metavar="DICT", help=argparse.SUPPRESS)

ARGS = vars(parser.parse_args())

ARGS["dict"] = json.loads(ARGS["dict"])

case = {self.gen_json_dict_str()}
mods = {{}}

if "post_process" in ARGS["dict"]["targets"]:
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

print(json.dumps({{**case, **mods}}))
""")

    def __str__(self) -> str:
        return f"tests/[bold magenta]{self.get_uuid()}[/bold magenta]: {self.trace}"

    def compute_tolerance(self) -> float:
        if self.params.get("qbmm", 'F') == 'T':
            return 1e-10

        if self.params.get("bubbles", 'F') == 'T':
            return 1e-10

        if self.params.get("hypoelasticity", 'F') == 'T':
            return 1e-7

        if self.params.get("relax", 'F') == 'T':
            return 1e-10

        if self.params.get("ib", 'F') == 'T':
            return 1e-10

        return 1e-12


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


def create_case(stack: CaseGeneratorStack, newTrace: str, newMods: dict, ppn: int = None) -> TestCase:
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

    return TestCase(' -> '.join(traces), mods, ppn)
