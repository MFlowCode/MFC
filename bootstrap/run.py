import os, json, rich, typing, dataclasses
from queue import Queue

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
        raise common.MFCException("QueueSystem::is_active: not implemented.")

    def gen_batch_header(self, args: dict, target_name: str) -> str:
        raise common.MFCException("QueueSystem::gen_batch_header: not implemented.")

    def gen_submit_cmd(self, filename: str) -> None:
        raise common.MFCException("QueueSystem::gen_submit_cmd: not implemented.")


class PBSSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("PBS")

    def is_active(self) -> bool:
        return 0 == os.system(f"qsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#PBS -A {args["account"]}
#PBS -l nodes={args["nodes"]}:ppn={args["cpus_per_node"]}
#PBS -l walltime={args["walltime"]}
#PBS -q {args["partition"]}
#PBS -N {job_name}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"qsub {filename}"


class LSFSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("LSF")

    def is_active(self) -> bool:
        return 0 == os.system(f"bsub -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#BSUB -P {args["account"]}
#BSUB -W {args["walltime"][:-3]}
#BSUB -J {job_name}
#BSUB -nnodes {args["nodes"]}
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"bsub {filename}"


class SLURMSystem(QueueSystem):
    def __init__(self) -> None:
        super().__init__("SLURM")

    def is_active(self) -> bool:
        return 0 == os.system(f"sbatch -h > /dev/null 2>&1")

    def gen_batch_header(self, args: dict, job_name: str) -> str:
        return f"""\
#SBATCH --job-name="{job_name}"
#SBATCH --time={args["walltime"]}
#SBATCH --nodes={args["nodes"]}
#SBATCH --ntasks-per-node={args["cpus_per_node"]}
#SBATCH --cpus-per-task={1}
#SBATCH --partition={args["partition"]}
#SBATCH --account={args["account"]}
#SBATCH --mail-user={args["email"]}
#SBATCH --mail-type=BEGIN, END, FAIL
"""

    def gen_submit_cmd(self, filename: str) -> None:
        return f"sbatch {filename}"


QUEUE_SYSTEMS = [ LSFSystem(), SLURMSystem(), PBSSystem() ]


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

        self.execute_batch_file(queue_sys, target_name)

        self.remove_batch_file(queue_sys, target_name)

    def get_batch_filepath(self, system: QueueSystem, target_name: str):
        case_dirpath = self.mfc.run.get_case_dirpath()

        return os.path.abspath(f"{case_dirpath}/{target_name}.sh")

    def create_batch_file(self, system: QueueSystem, target_name: str):
        job_name = self.mfc.run.get_job_name(target_name)

        BATCH_CONTENT: str = f"""\
#!/usr/bin/env bash
{system.gen_batch_header(self.mfc.args, job_name)}

RED="\\u001b[31m";   CYAN="\\u001b[36m";   BLUE="\\u001b[34m";    WHITE="\\u001b[37m"
GREEN="\\u001b[32m"; YELLOW="\\u001b[33m"; MAGENTA="\\u001b[35m"; COLOR_RESET="\\033[m"

TABLE_FORMAT_LINE="| - %-14s %-25s - %-14s %-25s |\\n"
TABLE_HEADER="/‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\ \\n"
TABLE_FOOTER="\_______________________________________________________________________________________/\\n"
TABLE_TITLE_FORMAT="| %8s $MAGENTA%-51s$COLOR_RESET                          |\\n"
TABLE_CONTENT=$(cat <<-END
$(printf "$TABLE_FORMAT_LINE" "Start-time:"    "$(date +%T)"                       "Start-date:"    "$(date +%T)")
$(printf "$TABLE_FORMAT_LINE" "Partition:"     "{self.mfc.args["partition"]}"      "Walltime:"      "{self.mfc.args["walltime"]}")
$(printf "$TABLE_FORMAT_LINE" "Account:"       "{self.mfc.args["account"]}"        "Nodes:"         "{self.mfc.args["nodes"]}")
$(printf "$TABLE_FORMAT_LINE" "CPUs (/node):"  "{self.mfc.args["cpus_per_node"]}"  "GPUs (/node):"  "{self.mfc.args["gpus_per_node"]}")
$(printf "$TABLE_FORMAT_LINE" "Input File:"    "{self.mfc.args["input"]}"          "Engine"         "{self.mfc.args["engine"]}")
$(printf "$TABLE_FORMAT_LINE" "Queue System:"  "{system.name}"                     "Mode:"          "{self.mfc.args["mode"]}")
$(printf "$TABLE_FORMAT_LINE" "Email:"         "{self.mfc.args["email"]}"          "Job Name:"      "{job_name}")\\n
END
)

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Starting" "{job_name}"
printf "$TABLE_CONTENT"
printf "$TABLE_FOOTER"

t_start=$(date +%s)

echo ""

{self.mfc.run.get_exec_cmd(target_name)}

echo ""

code=$?

status_msg="SUCCESS"
if [ "$code" -ne "0" ]; then
    status_msg="FAILURE"
fi

t_stop="$(date +%s)"

printf "$TABLE_HEADER"
printf "$TABLE_TITLE_FORMAT" "Finished" "{job_name}"
printf "$TABLE_FORMAT_LINE" "Total-time:"  "$(expr $t_stop - $t_start)s"  "Status:"    "$status_msg"
printf "$TABLE_CONTENT"
printf "$TABLE_FORMAT_LINE" "End-time:"    "$(date +%T)"                  "End-date:"  "$(date +%T)"
printf "$TABLE_FOOTER"

printf "\\nThank you for using MFC !\\n\\n"

exit $code
"""

        filepath = self.get_batch_filepath(system, target_name)

        common.file_write(filepath, BATCH_CONTENT)

    def remove_batch_file(self, system: QueueSystem, target_name: str):
        os.remove(self.get_batch_filepath(system, target_name))

    def execute_batch_file(self, system: QueueSystem, target_name: str):
        if 0 != os.system(system.gen_submit_cmd(self.get_batch_filepath(system, target_name))):
            raise common.MFCException(f"Running batch file for {system.name} failed.")


ENGINES = [ SerialEngine(), ParallelEngine() ]


class MFCRun:
    def __init__(self, mfc):
        self.mfc = mfc

    def get_case_dict(self):
        case:  dict = {}
        input: str  = self.mfc.args["input"].strip()

        rich.print(f"> > Fetching case dictionary from {input}...")

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

    def get_job_name(self, target_name: str):
        return f'MFC-{str(self.mfc.args["name"])}-{target_name}'

    def get_input_filepath(self, target_name: str):
        dirpath  = os.path.abspath(os.path.dirname(self.mfc.args["input"]))
        filename = f"{target_name}.inp"

        return f"{dirpath}/{filename}"

    def create_input_file(self, target_name: str, case_dict: dict):
        MASTER_KEYS: list = get_input_dict_keys(target_name)

        # Create Fortran-style input file content string
        dict_str = ""
        for key,val in case_dict.items():
            if key in MASTER_KEYS:
                dict_str += f"{key} = {val}\n"

        contents = f"&user_inputs\n{dict_str}&end/"

        # Save .inp input file
        common.file_write(self.get_input_filepath(target_name), contents)

    def remove_input_file(self, target_name: str):
        os.remove(self.get_input_filepath(target_name))

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
        bin = self.get_binpath(target_name)

        cd = f'cd "{self.get_case_dirpath()}"'
        ld = self.get_ld()

        if os.system("jsrun -h > /dev/null 2>&1") == 0:
            # ORNL Summit: https://docs.olcf.ornl.gov/systems/summit_user_guide.html?highlight=lsf#launching-a-job-with-jsrun
                        
            if int(self.mfc.args["cpus_per_node"]) != int(self.mfc.args["gpus_per_node"]) \
               and int(self.mfc.args["gpus_per_node"]) != 0:
               raise common.MFCException("JSRUN: Conflicting job execution parameters. If using GPUs, CPUs per node and GPUs per node must match.")

            # One resource set per CPU(Core)/GPU pair.
            rs=int(self.mfc.args["cpus_per_node"])
            cpus_per_rs=1
            gpus_per_rs=min(int(self.mfc.args["gpus_per_node"]), 1)
            tasks_per_rs=1

            options = f'--smpiargs="-gpu" --nrs{rs} --cpu_per_rs{cpus_per_rs} --gpu_per_rs{gpus_per_rs} --tasks_per_rs{tasks_per_rs}'

            return f'{cd} && {ld} jsrun {options} "{bin}"'
        elif os.system("srun -h > /dev/null 2>&1") == 0:
            raise common.MFCException("srun not implemented")
#            return f'{cd} && {ld} srun --nodes={self.mfc.args["nodes"]} -n  "{bin}"'
        else:
            return f'{cd} && {ld} mpiexec -N {self.mfc.args["nodes"]} -n {self.mfc.args["cpus_per_node"]} "{bin}"'

    def run(self):
        targets = self.mfc.args["targets"]
        if targets[0] == "mfc":
            targets = self.mfc.conf.get_dependency_names(targets[0], recursive=False)

        rich.print(f"""\
[bold][u]Run:[/u][/bold]
> Targets       (-t)  {self.mfc.args['targets']}
> Engine        (-e)  {self.mfc.args['engine']}
> Mode          (-m)  {self.mfc.args['mode']}
> Input         (-i)  {self.mfc.args['input']}
> Nodes         (-N)  {self.mfc.args['nodes']}
> CPUs (/node)  (-n)  {self.mfc.args['cpus_per_node']}
> GPUs (/node)  (-g)  {self.mfc.args["gpus_per_node"]}
> Walltime      (-w)  {self.mfc.args["walltime"]}
> Partition     (-p)  {self.mfc.args["partition"]}
> Account       (-a)  {self.mfc.args["account"]}
> Email         (??)  {self.mfc.args["email"]}
""")

        for target_name in targets:
            rich.print(f"> Running [bold magenta]{target_name}[/bold magenta]:")

            if not self.mfc.build.is_built(target_name):
                rich.print(f"> > Target {target_name} needs (re)building...")
                self.mfc.build.build_target(target_name)

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

            self.remove_input_file(target_name)
