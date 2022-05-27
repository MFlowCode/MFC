
import json
import hashlib
import binascii
import subprocess
import dataclasses

import common

Tend = 0.25
Nt   = 50
mydt = 0.0005

BASE_CFG = {
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'F',
    'ppn'                          : 1,
    'queue'                        : 'normal',
    'walltime'                     : '24:00:00',
    'mail_list'                    : '',
    'm'                            : 0,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : int(Nt+1),
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
    'nnode'                         : 4,
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

@dataclasses.dataclass
class Case:
    trace:  str
    params: dict

    def __init__(self, trace: str, mods: dict) -> None:
        self.trace  = trace
        self.params = {**BASE_CFG.copy(), **mods}

    def run(self, args: dict) -> subprocess.CompletedProcess:
        binary_option = ""
        if args["binary"] is not None:
            binary_option = f"-b {args['binary']}"

        command: str = f'''\
./mfc.sh run "{self.get_dirpath()}/case.py" -m "{args["mode"]}" -n {self["ppn"]} \
-t pre_process simulation {binary_option} 2>&1\
'''

        return subprocess.run(command, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True,
                              shell=True)

    def get_keys(self) -> str:
        return self.params.keys()

    def has_parameter(self, key: str)-> bool:
        return key in self.get_keys()

    def gen_json_dict_str(self) -> str:
        return json.dumps(self.params, indent=4)

    def get_uuid(self) -> str:
        return hex(binascii.crc32(hashlib.sha1(str(self.trace).encode()).digest())).upper()[2:].zfill(8)

    def get_dirpath(self):
        return f"{common.MFC_TESTDIR}/{self.get_uuid()}"

    def create_directory(self):
        dirpath = self.get_dirpath()

        content = f"""\
#!/usr/bin/env python3

import json

print(json.dumps({self.gen_json_dict_str()}))

"""

        common.create_directory(dirpath)

        common.file_write(f"{dirpath}/case.py", content)

    def __getitem__(self, key: str) -> str:
        if key not in self.params:
            raise common.MFCException(f"Case {self.trace}: Parameter {key} does not exist.")

        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = val


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


def create_case(stack: CaseGeneratorStack, newTrace, newMods) -> Case:
    mods: dict = {}

    for dict in stack.mods:
        mods.update(dict)
    mods.update(newMods)

    traces: list = []
    for trace in stack.trace[:] + [newTrace]:
        if not common.isspace(trace):
            traces.append(trace)

    return Case(' -> '.join(traces), mods)
