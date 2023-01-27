import os, hashlib, binascii, subprocess, dataclasses

from ..      import case, common
from ..state import ARG

Tend = 0.25
Nt   = 50
mydt = 0.0005

BASE_CFG = {
    'case_dir'                     : '\'.\'',
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
    'weno_vars'                    : 2,
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
    'prim_vars_wrt'                :'T',
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

    def __init__(self, trace: str, mods: dict, ppn: int = None) -> None:
        self.trace = trace
        self.ppn   = ppn if ppn is not None else 1
        super().__init__({**BASE_CFG.copy(), **mods})

    def run(self) -> subprocess.CompletedProcess:
        filepath          = f'"{self.get_dirpath()}/case.py"'
        tasks             = f"-n {self.ppn}"
        jobs              = f"-j {ARG('jobs')}"    if ARG("case_optimization")  else ""
        binary_option     = f"-b {ARG('binary')}"  if ARG("binary") is not None else ""
        case_optimization =  "--case-optimization" if ARG("case_optimization")  else "--no-build"
        
        mfc_script = ".\mfc.bat" if os.name == 'nt' else "./mfc.sh"
                
        command: str = f'''\
{mfc_script} run {filepath} {tasks} {binary_option} {case_optimization} \
{jobs} -t pre_process simulation 2>&1\
'''

        return subprocess.run(command, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True,
                              shell=True)

    def get_uuid(self) -> str:
        return hex(binascii.crc32(hashlib.sha1(str(self.trace).encode()).digest())).upper()[2:].zfill(8)

    def get_dirpath(self):
        return os.path.join(common.MFC_TESTDIR, self.get_uuid())

    def create_directory(self):
        dirpath = self.get_dirpath()

        content = f"""\
#!/usr/bin/env python3

import json

print(json.dumps({self.gen_json_dict_str()}))
"""

        common.create_directory(dirpath)

        common.file_write(f"{dirpath}/case.py", content)

    def __str__(self) -> str:
        return f"tests/[bold magenta]{self.get_uuid()}[/bold magenta]: {self.trace}"


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

    for dict in stack.mods:
        mods.update(dict)
    mods.update(newMods)

    if isinstance(newTrace, str):
        newTrace = [newTrace]

    traces: list = []
    for trace in stack.trace[:] + newTrace:
        if not common.isspace(trace):
            traces.append(trace)

    return TestCase(' -> '.join(traces), mods, ppn)
