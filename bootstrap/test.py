#!/usr/bin/env python3

import common

import os
import re
import copy
import time
import signal
import binascii
import threading
import subprocess
import dataclasses

from pathlib import Path

import rich, rich.progress

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


@dataclasses.dataclass
class TestThreadHolder:
    thread: threading.Thread
    ppn:    int


Tend = 0.25
Nt   = 50
mydt = 0.0005

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
})

class TestCaseConfiguration:
    parameters: dict = {}
    traceback:  str  = ""

    def __init__(self, parameters: list, traceback: list) -> None:
        self.parameters = {}
        for e in parameters:
            self.parameters.update(e)
        self.traceback  = ' -> '.join(traceback)

class MFCTest:
    def __init__(self, mfc):
        self.mfc = mfc

    def get_test_params(self):
        tests = []

        traceback  = []
        parameters = []

        for dimInfo in [ (["x"],           {'m': 299, 'n': 0,  'p': 0},  {"geometry": 1}),
                         (["x", "y"],      {'m': 49,  'n': 39, 'p': 0},  {"geometry": 3}),
                         (["x", "y", "z"], {'m': 24,  'n': 24, 'p': 24}, {"geometry": 9}) ]:
            dimParams = {**dimInfo[1]}

            for dimCmp in dimInfo[0]:
                dimParams[f"{dimCmp}_domain%beg"] = 0.E+00
                dimParams[f"{dimCmp}_domain%end"] = 1.E+00

            for patchID in range(1, 3+1):
                dimParams[f"patch_icpp({patchID})%geometry"] = dimInfo[2].get("geometry")

                if "z" in dimInfo[0]:
                    dimParams[f"patch_icpp({1})%z_centroid"] = 0.05
                    dimParams[f"patch_icpp({1})%length_z"]   = 0.1

                    dimParams[f"patch_icpp({2})%z_centroid"] = 0.45
                    dimParams[f"patch_icpp({2})%length_z"]   = 0.7

                    dimParams[f"patch_icpp({3})%z_centroid"] = 0.9
                    dimParams[f"patch_icpp({3})%length_z"]   = 0.2


                    dimParams[f"patch_icpp({patchID})%y_centroid"] = 0.5
                    dimParams[f"patch_icpp({patchID})%length_y"]   = 1
                    dimParams[f"patch_icpp({patchID})%x_centroid"] = 0.5
                    dimParams[f"patch_icpp({patchID})%length_x"]   = 1

                elif "y" in dimInfo[0]:
                    dimParams[f"patch_icpp({1})%y_centroid"] = 0.05
                    dimParams[f"patch_icpp({1})%length_y"]   = 0.1

                    dimParams[f"patch_icpp({2})%y_centroid"] = 0.45
                    dimParams[f"patch_icpp({2})%length_y"]   = 0.7

                    dimParams[f"patch_icpp({3})%y_centroid"] = 0.9
                    dimParams[f"patch_icpp({3})%length_y"]   = 0.2


                    dimParams[f"patch_icpp({patchID})%x_centroid"] = 0.5
                    dimParams[f"patch_icpp({patchID})%length_x"]   = 1
                else:
                    dimParams[f"patch_icpp({1})%x_centroid"] = 0.05
                    dimParams[f"patch_icpp({1})%length_x"]   = 0.1

                    dimParams[f"patch_icpp({2})%x_centroid"] = 0.45
                    dimParams[f"patch_icpp({2})%length_x"]   = 0.7

                    dimParams[f"patch_icpp({3})%x_centroid"] = 0.9
                    dimParams[f"patch_icpp({3})%length_x"]   = 0.2

                if "x" in dimInfo[0]:
                    dimParams[f"patch_icpp({patchID})%vel(1)"]     = 0.0

                if "y" in dimInfo[0]:
                    dimParams[f"patch_icpp({patchID})%vel(2)"]     = 0.0

                if "z" in dimInfo[0]:
                    dimParams[f"patch_icpp({patchID})%vel(3)"]     = 0.0

            traceback.append (f"{len(dimInfo[0])}D (m={dimInfo[1].get('m')},n={dimInfo[1].get('n')},p={dimInfo[1].get('p')})")
            parameters.append(dimParams)

            for bc in [ -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12]:
                params = {}
                for dimCmp in dimInfo[0]:
                    params = {**params, **{f'bc_{dimCmp}%beg': bc, f'bc_{dimCmp}%end': bc}}

                trace = f"bc={bc}"
                tests.append(TestCaseConfiguration(parameters + [params], traceback + [trace]))

                if bc == -3:
                    parameters.append(params)
                    traceback.append(trace)

            for weno_order in [3, 5]:
                traceback.append (f"weno_order={weno_order}")
                parameters.append({'weno_order': weno_order})
                for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                    traceback.append (f"(mapped_weno={mapped_weno},mp_weno={mp_weno})")
                    parameters.append({'mapped_weno': mapped_weno, 'mp_weno': mp_weno})
                    if not (mp_weno == 'T' and weno_order != 5):
                        tests.append(TestCaseConfiguration(parameters, traceback))
                    traceback.pop()
                    parameters.pop()
                traceback.pop()
                parameters.pop()

            for num_fluids in [1, 2]:
                traceback.append(f"num_fluids={num_fluids}")
                parameters.append({"num_fluids": num_fluids})

                if num_fluids == 2:
                    parameters.append({'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0,'patch_icpp(1)%alpha_rho(1)': 0.81, 'patch_icpp(1)%alpha(1)': 0.9, 'patch_icpp(1)%alpha_rho(2)': 0.19, 'patch_icpp(1)%alpha(2)': 0.1, 'patch_icpp(2)%alpha_rho(1)': 0.25, 'patch_icpp(2)%alpha(1)': 0.5, 'patch_icpp(2)%alpha_rho(2)': 0.25, 'patch_icpp(2)%alpha(2)': 0.5, 'patch_icpp(3)%alpha_rho(1)': 0.08, 'patch_icpp(3)%alpha(1)': 0.2, 'patch_icpp(3)%alpha_rho(2)': 0.0225, 'patch_icpp(3)%alpha(2)': 0.8,})

                for riemann_solver in [1, 2]:
                    traceback.append(f"riemann_solver={riemann_solver}")
                    parameters.append({'riemann_solver': riemann_solver})

                    tests.append(TestCaseConfiguration(parameters + [{'mixture_err': 'T'}], traceback + ['mixture_err=T']))
                    tests.append(TestCaseConfiguration(parameters + [{'avg_state':   '1'}], traceback + ['avg_state=1']))
                    tests.append(TestCaseConfiguration(parameters + [{'wave_speeds': '2'}], traceback + ['wave_speeds=2']))

                    if num_fluids == 2:
                        if riemann_solver == 2:
                            tests.append(TestCaseConfiguration(parameters + [{'alt_soundspeed': 'T'}], traceback + ['alt_soundspeed=T']))

                        tests.append(TestCaseConfiguration(parameters + [{'avg_state':     1}], traceback + ['avg_state=1']))
                        tests.append(TestCaseConfiguration(parameters + [{'wave_speeds':   2}], traceback + ['wave_speeds=2']))
                        tests.append(TestCaseConfiguration(parameters + [{'mpp_lim':     'T'}], traceback + ['mpp_lim=T']))

                    traceback.pop()
                    parameters.pop()

                if num_fluids == 2:
                    parameters.pop()

                traceback.pop()
                parameters.pop()

            if len(dimInfo[0]) == 3:
                tests.append(TestCaseConfiguration(parameters + [{'ppn': 2, 'm': 29, 'n': 29, 'p': 49}], traceback + [f'ppn=2,m=29,n=29,p=49']))
            else:
                tests.append(TestCaseConfiguration(parameters + [{'ppn': 2}], traceback + [f'ppn=2']))

            if len(dimInfo[0]) == 1:
                parameters.append({'dt': 1e-07})
            elif len(dimInfo[0]) == 2:
                parameters.append({'dt': 1e-06})
            else:
                parameters.append({'dt':1e-06})

            if len(dimInfo[0]) > 0:
                traceback.append(f"bubbles={'T'}")
                parameters.append({"bubbles": 'T'})


                parameters.append({'nb' : 3,  'fluid_pp(1)%gamma' : 0.16, 'fluid_pp(1)%pi_inf': 3515.0, 'fluid_pp(2)%gamma': 2.5, 'fluid_pp(2)%pi_inf': 0.0, 'fluid_pp(1)%mul0' : 0.001002, 'fluid_pp(1)%ss' : 0.07275,'fluid_pp(1)%pv' : 2338.8,'fluid_pp(1)%gamma_v' : 1.33,'fluid_pp(1)%M_v' : 18.02,'fluid_pp(1)%mu_v' : 8.816e-06,'fluid_pp(1)%k_v' : 0.019426,'fluid_pp(2)%gamma_v' : 1.4,'fluid_pp(2)%M_v' : 28.97,'fluid_pp(2)%mu_v' : 1.8e-05, 'fluid_pp(2)%k_v' : 0.02556, 'patch_icpp(1)%alpha_rho(1)': 0.999999999999, 'patch_icpp(1)%alpha(1)': 1e-12, 'patch_icpp(2)%alpha_rho(1)': 0.96, 'patch_icpp(2)%alpha(1)': 4e-02,  'patch_icpp(3)%alpha_rho(1)': 0.999999999999, 'patch_icpp(3)%alpha(1)': 1e-12, 'patch_icpp(1)%pres': 1.0, 'patch_icpp(2)%pres': 1.0, 'patch_icpp(3)%pres': 1.0 })

                traceback.append(f"Monopole={'T'}")
                parameters.append({"Monopole": 'T'})

                if len(dimInfo[0]) >= 2:
                    parameters.append({'Mono(1)%loc(2)': 0.5})


                if len(dimInfo[0]) >= 3:
                    parameters.append({'Mono(1)%loc(3)': 0.5, 'Mono(1)%support': 3})

                for polytropic in ['T', 'F']:

                    for bubble_model in [3, 2]:


                        traceback.append(f"polytropic={polytropic}")
                        parameters.append({'polytropic' : polytropic})

                        traceback.append(f"bubble_model={bubble_model}")
                        parameters.append({'bubble_model' : bubble_model})

                        if not (polytropic == 'F' and bubble_model == 3):
                            tests.append(TestCaseConfiguration(parameters, traceback))

                        traceback.pop()
                        parameters.pop()

                        traceback.pop()
                        parameters.pop()


                parameters.append({'polytropic' : 'T'})
                parameters.append({'bubble_model' : 2})

                parameters.append({'nb':1})
                traceback.append(f"nb={'1'}")
                tests.append(TestCaseConfiguration(parameters, traceback))

                traceback.pop()
                parameters.pop()

                traceback.append(f"qbmm={'T'}")
                parameters.append({'qbmm': 'T'})

                tests.append(TestCaseConfiguration(parameters, traceback))

                parameters.append({'bubble_model' : 3})

                tests.append(TestCaseConfiguration(parameters, traceback))

                traceback.pop()
                for i in range(4):
                    parameters.pop()

                if len(dimInfo[0]) >= 2:
                    parameters.pop()


                if len(dimInfo[0]) >= 3:
                    parameters.pop()

                parameters.pop()
                traceback.pop()

                parameters.pop()

                parameters.pop()
                traceback.pop()

            parameters.pop()

            for i in range(2):
                traceback.pop()
                parameters.pop()

        return tests

    def filter_tests(self, tests: list):
        if len(self.mfc.args["only"]) > 0:
            for i, test in enumerate(tests[:]):
                doKeep = False
                for o in self.mfc.args["only"]:
                    if str(o).isnumeric():
                        testID = i+1
                        if testID == int(o):
                            doKeep = True
                            break

                    testHash = self.get_case_dir_name(test.parameters)
                    if str(o) == testHash:
                        doKeep = True
                        break

                if not doKeep:
                    tests.remove(test)
        
        return tests

    def test(self):
        if self.mfc.args["generate"]:
            common.delete_directory_recursive(common.MFC_TESTDIR)
            common.create_directory(common.MFC_TESTDIR)

        for target in ["pre_process", "simulation"]:
            if not self.mfc.build.is_built(target):
                rich.print(f"> {target} needs (re)building...")
                self.mfc.build.build_target(f"{target}", "> > ")

        tests = self.filter_tests(self.get_test_params())

        nAvailableThreads = self.mfc.args["jobs"]
        threads           = []

        for i, test in enumerate(rich.progress.track(tests, "Queue Tests [bold blue]mfc[/bold blue]...")):
            test: TestCaseConfiguration
            
            # Use the correct test ID if --only is selected
            testID = i+1
            if len(self.mfc.args["only"]):
                testID = self.mfc.args["only"][i]

            ppn = self.get_case_from_mods(test.parameters)["ppn"]

            # Wait until there are threads available
            while nAvailableThreads < ppn:
                # This is important if "-j 1" is used (the default) since there
                # are test cases that require ppn=2
                if ppn > self.mfc.args["jobs"] and nAvailableThreads > 0:
                    break

                # Keep track of threads that are done
                for threadID, threadHolder in enumerate(threads):
                    threadHolder: TestThreadHolder

                    if not threadHolder.thread.is_alive():
                        nAvailableThreads += threadHolder.ppn
                        del threads[threadID]
                        break
                
                # Do not overwhelm this core with this loop
                time.sleep(0.2)

            # nMaxAvailableThreads >= ppn

            rich.print(f"Submit [T{len(threads)}] #{i}/{len(tests)} - {test.traceback}")
            thread = threading.Thread(target=self.handle_case, args=(testID, test))
            thread.start()
            threads.append(TestThreadHolder(thread, ppn))
            nAvailableThreads -= ppn

        # Join remaining threads
        while len(threads) != 0:
            for threadID, threadHolder in enumerate(threads):
                threadHolder: TestThreadHolder

                if not threadHolder.thread.is_alive():
                    del threads[threadID]
                    break
            
            time.sleep(0.2)

        rich.print(f"> Tested [bold green]✓[/bold green]")

    def get_case_dir_name(self, mods: dict):
        return hex(binascii.crc32(str(mods.items()).encode()))[2:]

    def get_case_dir(self, mods: dict):
        return f"{common.MFC_TESTDIR}/{self.get_case_dir_name(mods)}"

    def get_case_from_mods(self, mods: dict):
        case = copy.deepcopy(BASE_CASE)

        for key, val in mods.items():
            case[key] = val

        return case

    def create_case_dir(self, mods: dict):
        case     = self.get_case_from_mods(mods)
        case_dir = self.get_case_dir(mods)

        content = f"""\
#!/usr/bin/env python3

import json

print(json.dumps({case.create_case_dict_str()}))

"""

        common.create_directory(case_dir)

        common.file_write(f"{case_dir}/case.py", content)

    def golden_file_compare_match(self, truth: str, candidate: str, tol):
        if truth.count('\n') != candidate.count('\n'):
            return (False, "Line count didn't match.")

        if "NaN" in truth:
            return (False, "NaN in golden file")

        if "NaN" in candidate:
            return (False, "NaN in packed file")

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
                abs_delta = abs(numbers_cand[i]-numbers_trust[i])
                rel_diff  = abs(abs_delta/numbers_trust[i]) if numbers_trust[i] != 0 else 0
                if    (abs_delta > tol and rel_diff > tol):
                    percent_diff = rel_diff*100
                    return (False, f"Error margin is too high for the value #{i+1} in {file_subpath}: ~{round(percent_diff, 5)}% (~{round(abs_delta, 5)}).")

        # Both tests gave the same results within an acceptable tolerance
        return (True, "")

    def get_test_summary(self, mods: dict):
        return "".join([f"{str(x[0]).split('_')[0][:4]}-{str(x[1])[:4]}_" for x in mods.items()])[:-1]

    def handle_case(self, testID, test: TestCaseConfiguration):
        try:
            self.create_case_dir(test.parameters)
            if('qbmm' in test.parameters):
                tol = 1e-7
            elif('bubbles' in test.parameters):
                tol = 1e-10
            else:
                tol = 1e-12

            def on_test_errror(msg: str = "", term_out: str = ""):
                common.clear_line()
                rich.print(f"> Test #{testID}: Failed! [bold green]✓[/bold green]")
                if msg != "":
                    rich.print(f"> {msg}")

                common.file_write(f"{common.MFC_TESTDIR}/failed_test.txt", f"""\
(1/3) Test #{testID}:
- Test ID:  {testID}
- Summary:  {test.traceback}
- Location: {self.get_case_dir(test.parameters)}
- Error:    {msg}
- When:     {common.get_datetime_str()}

(2/3) Test case:
{self.get_case_from_mods(test.parameters).create_case_dict_str()}

(3/3) Terminal output:
{term_out}
""")

                rich.print(f"> Please read {common.MFC_TESTDIR}/failed_test.txt for more information.")
                raise common.MFCException("Testing failed (view above).")

            cmd = subprocess.run(f'./mfc.sh run "{self.get_case_dir(test.parameters)}/case.py" -m "{self.mfc.args["mode"]}" -c {self.get_case_from_mods(test.parameters).parameters["ppn"]} -t pre_process simulation 2>&1',
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, shell=True)
            common.file_write(f"{self.get_case_dir(test.parameters)}/out.txt", cmd.stdout)

            if cmd.returncode != 0:
                on_test_errror("MFC Execution Failed.", cmd.stdout)

            pack = self.pack_case_output(test.parameters)
            common.file_write(f"{self.get_case_dir(test.parameters)}/pack.txt", pack)

            golden_filepath = f"{self.get_case_dir(test.parameters)}/golden.txt"

            if self.mfc.args["generate"]:
                common.delete_file(golden_filepath)
                common.file_write(golden_filepath, pack)

            if not os.path.isfile(golden_filepath):
                common.clear_line()
                on_test_errror("Golden file doesn't exist! To generate golden files, use the '-g' flag.", cmd.stdout)

            golden_file_content = common.file_read(golden_filepath)
            bSuccess, errorMsg  = self.golden_file_compare_match(golden_file_content, pack, tol)
            if not bSuccess:
                on_test_errror(errorMsg, cmd.stdout)
        except BaseException as exc:
            print(exc)
            
            # Exit
            os.kill(os.getpid(), signal.SIGTERM)

    def pack_case_output(self, params: dict):
        result: str = ""

        case_dir = self.get_case_dir(params)
        D_dir    = f"{case_dir}/D/"

        for filepath in list(Path(D_dir).rglob("*.dat")):
            file_content   = common.file_read(filepath)
            short_filepath = str(filepath).replace(f'{case_dir}/', '')

            result += f"{short_filepath} " + re.sub(r' +', ' ', file_content.replace('\n', ' ')).strip() + '\n'

        return result
