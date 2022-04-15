import os, json, rich, typing, dataclasses

import common

PRE_PROCESS = ['case_dir', 'old_grid', 'old_ic', 't_step_old', 'm', 'n', 'p',
               'cyl_coord', 'model_eqns', 'num_fluids', 'adv_alphan', 'mpp_lim',
               'weno_order', 'precision', 'parallel_io', 'perturb_flow',
               'perturb_flow_fluid', 'perturb_sph', 'perturb_sph_fluid',
               'fluid_rho', 'hypoelasticity', 'num_patches', 'Ca', 'Web',
               'Re_inv', 'pref', 'rhoref', 'bubbles' , 'polytropic',
               'polydisperse', 'poly_sigma', 'thermal', 'nb', 'R0ref', 'qbmm',
               'dist_type', 'R0_type', 'nnode', 'sigR', 'sigV', 'rhoRV']

for cmp in ["x", "y", "z"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS.append(f"{cmp}_{prepend}")
    
    for append in ["stretch", "a", "loops"]:
        PRE_PROCESS.append(f"{append}_{cmp}")
    
    PRE_PROCESS.append(f"bc_{cmp}%beg")
    PRE_PROCESS.append(f"bc_{cmp}%end")

for f_id in range(1, 10+1):
    PRE_PROCESS.append(f'fluid_rho({f_id})')

    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G"]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{attribute}")

for p_id in range(1, 10+1):
    for attribute in ["geometry", "radius", "radii", "epsilon", "beta",
                      "normal", "smoothen", "smooth_patch_id", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf", "r0", "v0", "p0", "m0"]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{attribute}")

    for cmp_id, cmp in enumerate(["x", "y", "z"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_icpp({p_id})%{cmp}_centroid')
        PRE_PROCESS.append(f'patch_icpp({p_id})%length_{cmp}')

        for append in ["radii", "normal", "vel"]:
            PRE_PROCESS.append(f'patch_icpp({p_id})%{append}({cmp_id})')

    for arho_id in range(1, 10+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha({arho_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha_rho({arho_id})')
    
    for taue_id in range(1, 6+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%tau_e({arho_id})')

SIMULATION = ['case_dir', 'run_time_info', 't_step_old', 't_tol', 'debug', 'm',
              'n', 'p', 'cyl_coord', 'dt', 't_step_start', 't_step_stop',
              't_step_save', 'model_eqns', 'num_fluids', 'adv_alphan',
              'mpp_lim', 'time_stepper', 'weno_vars', 'weno_order', 'weno_eps',
              'char_decomp', 'mapped_weno', 'mp_weno', 'weno_avg',
              'weno_Re_flux', 'riemann_solver', 'wave_speeds', 'avg_state',
              'commute_err', 'split_err', 'alt_crv', 'alt_soundspeed',
              'regularization', 'reg_eps', 'null_weights', 'mixture_err',
              'tvd_riemann_flux', 'tvd_rhs_flux', 'tvd_wave_speeds', 'flux_lim',
              'We_riemann_flux', 'We_rhs_flux', 'We_src', 'We_wave_speeds',
              'lsq_deriv', 'parallel_io', 'precision', 'hypoelasticity',
              'fd_order' , 'com_wrt', 'num_probes', 'probe_wrt', 'cb_wrt',
              'threshold_mf', 'moment_order', 'pref', 'rhoref', 'polydisperse',
              'poly_sigma', 'bubbles', 'bubble_model', 'polytropic', 'thermal',
              'R0ref', 'Ca', 'Web', 'Re_inv', 'nb', 'Monopole', 'num_mono',
              'qbmm', 'R0_type', 'nnode', 'integral_wrt', 'num_integrals']

for cmp in ["x", "y", "z"]:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')

for wrt_id in range(1,10+1):
    SIMULATION.append(f'com_wrt({wrt_id})')
    SIMULATION.append(f'cb_wrt({wrt_id})')

    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,3+1):
    for cmp in ["x", "y", "z"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for mf_id in range(1,5+1):
    SIMULATION.append(f'threshold_mf({mf_id})')

for order_id in range(1,5+1):
    SIMULATION.append(f'moment_order({order_id})')

for f_id in range(1,10+1):
    for attribute in ["gamma", "pi_inf", "mul0", "ss", "pv", "gamma_v", "M_v",
                      "mu_v", "k_v", "G"]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for mono_id in range(1,4+1):
        for attribute in ["mag", "length", "dir", "npulse", "pulse", "support",
                          "delay"]:
            SIMULATION.append(f"Mono({mono_id})%{attribute}")
        
        for cmp_id in range(1,3+1):
            SIMULATION.append(f"Mono({mono_id})%loc({cmp_id})")

    for int_id in range(1,5+1):
        for cmp in ["x", "y", "z"]:
            SIMULATION.append(f"integral({int_id})%{cmp}min")
            SIMULATION.append(f"integral({int_id})%{cmp}max")

    for r_id in range(1,10+1):
        # FIXME: Add all fluid_pp({f_id})%...(...)
        # I don't get how some indices are selected

        #print("FIX THIS")
        pass

POST_PROCESS = ['case_dir', 'cyl_coord', 'm', 'n', 'p', 't_step_start',
                't_step_stop', 't_step_save', 'model_eqns', 'num_fluids',
                'adv_alphan', 'mpp_lim', 'weno_order', 'alt_soundspeed',
                'mixture_err', 'parallel_io', 'hypoelasticity',
                'polydisperse', 'poly_sigma', 'polytropic', 'thermal',
                'pref', 'Ca', 'Web', 'Re_inv', 'rhoref', 'bubbles',
                'R0ref', 'nb', 'format', 'precision', 'coarsen_silo',
                'fourier_decomp', 'fourier_modes%beg',
                'fourier_modes%end', 'alpha_rho_wrt', 'rho_wrt',
                'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt',
                'pres_wrt', 'alpha_wrt', 'kappa_wrt', 'gamma_wrt',
                'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt',
                'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt',
                'schlieren_wrt', 'schlieren_alpha', 'fd_order']

for cmp_id in range(1,3+1):
    cmp = ["x", "y", "z"][cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg')
    POST_PROCESS.append(f'bc_{cmp}%end')

    for attribute in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{attribute}({cmp_id})')

for fl_id in range(1,10+1):
    for append in ["schlieren_alpha", "alpha_rho_wrt", "alpha_wrt", "kappa_wrt"]:
        POST_PROCESS.append(f'{append}({fl_id})')

    for attribute in ["gamma", "pi_inf", "ss", "pv", "gamma_v", "M_v", "mu_v", "k_v", "G", "mul0"]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{attribute}")

def get_input_dict_keys(target_name: str) -> list:
    if target_name == "pre_process":  return PRE_PROCESS.copy()
    if target_name == "simulation":   return SIMULATION.copy()
    if target_name == "post_process": return POST_PROCESS.copy()
    
    raise common.MFCException(f"[INPUT DICTS] Target {target_name} doesn't have an input dict.")


@dataclasses.dataclass
class QueueSystem:
    name: str

    def is_active(self) -> bool:
        raise common.MFCException("is_active: not implemented.")

    def gen_batch_header(self, args: dict, target: str) -> str:
        raise common.MFCException("gen_batch_header: not implemented.")


class PBSSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("PBS")
    
    def is_active(self) -> bool:
        return True
        #FIXME:
        #return 0 == os.system(f"qsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, target: str) -> str:
        return f"""\
#PBS -l nodes={args["nodes"]}:ppn={args["tasks_per_node"]}
#PBS -l walltime={args["walltime"]}
#PBS -q {args["partition"]}
#PBS -N {target}\
# """

    def gen_batch_cmd(self):
        pass


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")
    
    def is_active(self) -> bool:
        return True
        #FIXME:
        #return 0 == os.system(f"bsub -h > /dev/null 2>&1")
    
    def gen_batch_header(self, args: dict, target: str) -> str:
        return f"""\
#SBATCH --job-name="{target}"
#SBATCH --time={args["walltime"]}
#SBATCH --ntasks={0}
#SBATCH --nodes={args["nodes"]}
#SBATCH --ntasks-per-node={args["tasks_per_node"]}
#SBATCH --cpus-per-task={1}
#SBATCH --partition={args["partition"]}
#SBATCH --account={args["account"]}
"""

    def gen_batch_cmd(self):
        pass


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")
    
    def is_active(self) -> bool:
        return True
        #FIXME:
        #return 0 == os.system(f"sbatch -h > /dev/null 2>&1")
    
    def gen_batch_header(self, args: dict, target: str) -> str:
        pass

    def gen_batch_cmd(self):
        pass



QUEUE_SYSTEMS = [ PBSSystem(), LSFSystem(), SLURMSystem() ]


@dataclasses.dataclass
class Engine:
    name: str
    slug: str

    def run(self, mfc, target_name: str) -> None:
        raise common.MFCException(f"MFCEngine::run: not implemented for {self.name}.")


class SerialEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Serial", "serial")

    def run(self, mfc, target_name: str) -> None:
        self.mfc = mfc

        date = f"> > [bold cyan][{common.get_datetime_str()}][/bold cyan]"
        rich.print(f"{date} Running...")
        common.execute_shell_command(self.mfc.run.get_exec_cmd(target_name))

        rich.print(f"> > Done [bold green]✓[/bold green]")

class ParallelEngine(Engine):
    def __init__(self) -> None:
        super().__init__("Parallel", "parallel")

    def run(self, mfc, target_name: str) -> None:
        self.mfc = mfc

        queue_sys = self.mfc.run.detect_queue_system()

        self.create_batch_file(queue_sys, target_name)

        self.execute_batch_file(queue_sys)
    
    def create_batch_file(self, system: QueueSystem, target_name: str):
        case_dirpath = self.mfc.run.get_case_dirpath()

        BATCH_CONTENT: str = f"""\
#!/bin/sh -l
{system.gen_batch_header(self.mfc.args, target_name)}

echo "================================================="
echo "| Starting job #{target_name}"
echo "| - Start-date: `date +%D`"
echo "| - Start-time: `date +%T`"
echo "================================================="

t_start=$(date +%s)

{self.mfc.run.get_exec_cmd(target_name)}

code=$?

status_msg="SUCCESS"
if [ "$code" -ne "0" ]; then
    status_msg="FAILURE"
fi

t_stop=$(date +%s)

echo "================================================="
echo "| Finished job {target_name}: $status_msg"
echo "| - End-date: `date +%D`"
echo "| - End-time: `date +%T`"
echo "| - Total-time: $(expr $t_stop - $t_start)s"
echo "================================================="

exit $code
"""

        common.file_write(f"{case_dirpath}/{target_name}.sh", BATCH_CONTENT)

    def execute_batch_file(self, system: QueueSystem):
        if system == "PBS":
            # handle
            pass
        else:
            raise common.MFCException(f"Running batch file for {system.name} is not supported.")


ENGINES = [ SerialEngine(), ParallelEngine() ]


class MFCRun:
    def __init__(self, mfc):
        self.mfc = mfc

    def get_case_dict(self):
        case:  dict = {}
        input: str  = self.mfc.args["input"].strip()

        rich.print(f"> > Fetching case dir from {input}...")

        if input.endswith(".py"):
            (output, err) = common.get_py_program_output(input)

            if err != 0:
                rich.print(f"> > Input file {input} terminated with a non-zero exit code. View the output bellow: [bold red]❌[/bold red]")
                for line in output.splitlines():
                    rich.print(line)

                raise common.MFCException(f"> > Input file {input} terminated with a non-zero exit code. View above.")

            case = json.loads(output)
        else:
            rich.print(f"> > Unrecognized input file format for '{input}'. Please check the extension. [bold red]✘[/bold red]")
            raise common.MFCException("Unrecognized input file format.")
        
        return case

    def create_input_file(self, target_name: str, case_dict: dict):
        MASTER_KEYS: list = get_input_dict_keys(target_name)

        # Create Fortran-style input file content string
        dict_str = ""
        for key,val in case_dict.items():
            if key in MASTER_KEYS:
                dict_str += f"{key} = {val}\n"
        
        contents = f"&user_inputs\n{dict_str}&end/"

        # Save .inp input file
        dirpath  = os.path.abspath(os.path.dirname(self.mfc.args["input"]))
        filename = f"{target_name}.inp"
        common.file_write(f"{dirpath}/{filename}", contents)

    def detect_queue_system(self) -> str:
        for system in QUEUE_SYSTEMS:
            if system.is_active():
                rich.print(f"> > Detected the [bold magenta]{system.name}[/bold magenta] queueing system.")
                return system
        
        raise common.MFCException(f"Failed to detect a queueing system.")

    def get_binpath(self, target: str) -> str:
        return f'{self.mfc.build.get_build_path(target)}/bin/{target}'

    def get_ld(self) -> str:
        return f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{common.MFC_SUBDIR}/common/build/lib"'

    def get_case_dirpath(self) -> str:
        return os.path.abspath(os.path.dirname(self.mfc.args["input"]))

    def get_exec_cmd(self, target_name: str):
        bin  = self.get_binpath(target_name)
        
        cd   = f'cd {self.get_case_dirpath()}'
        ld   = self.get_ld()

        if os.system("jsrun -h > /dev/null 2>&1") == 0:
            return f'{cd} && {ld} jsrun -r{self.mfc.args["tasks_per_node"]} -g1 -c1 -a1 "{bin}"'
        else:
            return f'{cd} && {ld} mpiexec -np {self.mfc.args["tasks_per_node"]} "{bin}"'

    def run(self):
        targets = self.mfc.args["targets"]
        if targets[0] == "mfc":
            targets = self.mfc.conf.get_dependency_names(targets[0], recursive=False)

        rich.print(f"""\
[bold][u]Run:[/u][/bold]
> Targets       (-t)  {self.mfc.args['targets']}
> Engine        (-e)  {self.mfc.args['engine']}
> Config        (-cc) {self.mfc.args['compiler_configuration']}
> Input         (-i)  {self.mfc.args['input']}
> Nodes         (-N)  {self.mfc.args['nodes']}
> Tasks (/node) (-n)  {self.mfc.args['tasks_per_node']}
""")

        for target_name in targets:
            rich.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")
            
            #FIXME:
            #if not self.mfc.build.is_build_satisfied(target_name):
            #    rich.print(f"> > Target {target_name} needs (re)building...")
            #    self.mfc.build.build_target(target_name)

            self.create_input_file(target_name, self.get_case_dict())

            engine_slug = self.mfc.args["engine"]
            engine: Engine = None
            for candidate in ENGINES:
                candidate: Engine

                if candidate.slug == engine_slug:
                    engine = candidate
                    break
            
            if engine == None:
                raise common.MFCException(f"Unsupported engine {engine_slug}.")

            engine.run(self.mfc, target_name)
