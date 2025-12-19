import os, json, glob, typing, dataclasses

import pyrometheus as pyro
import cantera     as ct

from ..printer import cons
from ..        import common, build
from ..state   import ARGS, ARG, gpuConfigOptions
from ..case    import Case
from ..        import case_validator

@dataclasses.dataclass(init=False)
class MFCInputFile(Case):
    filename: str
    dirpath:  str

    def __init__(self, filename: str, dirpath: str, params: dict) -> None:
        super().__init__(params)
        self.filename = filename
        self.dirpath  = dirpath

    def generate_inp(self, target) -> None:
        target = build.get_target(target)

        # Save .inp input file
        common.file_write(f"{self.dirpath}/{target.name}.inp", self.get_inp(target))

    def __save_fpp(self, target, contents: str) -> None:
        inc_dir = os.path.join(target.get_staging_dirpath(self), "include", target.name)
        common.create_directory(inc_dir)

        fpp_path = os.path.join(inc_dir, "case.fpp")

        cons.print("Writing a (new) custom case.fpp file.")
        common.file_write(fpp_path, contents, True)

    def get_cantera_solution(self) -> ct.Solution:
        if self.params.get("chemistry", 'F') == 'T':
            cantera_file = self.params["cantera_file"]

            candidates = [
                cantera_file,
                os.path.join(self.dirpath, cantera_file),
                os.path.join(common.MFC_MECHANISMS_DIR, cantera_file),
            ]
        else:
            # If Chemistry is turned off, we return a default (dummy) solution
            # that will not be used in the simulation, so that MFC can still
            # be compiled.
            candidates = ["h2o2.yaml"]

        for candidate in candidates:
            try:
                return ct.Solution(candidate)
            except Exception:
                continue

        raise common.MFCException(f"Cantera file '{cantera_file}' not found. Searched: {', '.join(candidates)}.")

    def generate_fpp(self, target) -> None:
        if target.isDependency:
            return

        cons.print(f"Generating [magenta]case.fpp[/magenta].")
        cons.indent()

        # Case FPP file
        self.__save_fpp(target, self.get_fpp(target))

        # (Thermo)Chemistry source file
        modules_dir = os.path.join(target.get_staging_dirpath(self), "modules", target.name)
        common.create_directory(modules_dir)

        # Determine the real type based on the single precision flag
        real_type = 'real(sp)' if (ARG('single') or ARG('mixed')) else 'real(dp)'

        if ARG("gpu") == gpuConfigOptions.MP.value:
            directive_str = 'mp'
        elif ARG("gpu") == gpuConfigOptions.ACC.value:
            directive_str = 'acc'
        else:
            directive_str = None

        # Write the generated Fortran code to the m_thermochem.f90 file with the chosen precision
        common.file_write(
            os.path.join(modules_dir, "m_thermochem.f90"),
            pyro.FortranCodeGenerator().generate(
                "m_thermochem",
                self.get_cantera_solution(),
                pyro.CodeGenerationOptions(scalar_type = real_type, directive_offload = directive_str)
            ),
            True
        )

        cons.unindent()


    def validate_constraints(self, target) -> None:
        """Validate case parameter constraints for a given target stage"""
        target_obj = build.get_target(target)
        stage = target_obj.name

        try:
            case_validator.validate_case_constraints(self.params, stage)
        except case_validator.CaseConstraintError as e:
            raise common.MFCException(f"Case validation failed for {stage}:\n{e}") from e

    # Generate case.fpp & [target.name].inp
    def generate(self, target) -> None:
        # Validate constraints before generating input files
        self.validate_constraints(target)
        self.generate_inp(target)
        cons.print()
        self.generate_fpp(target)

    def clean(self, _targets) -> None:
        targets = [build.get_target(target) for target in _targets]

        files = set()
        dirs  = set()

        files = set([
            "equations.dat", "run_time.inf", "time_data.dat",
            "io_time_data.dat", "fort.1", "pre_time_data.dat"
        ] + [f"{target.name}.inp" for target in targets])

        if build.PRE_PROCESS in targets:
            files = files | set(glob.glob(os.path.join(self.dirpath, "D", "*.000000.dat")))
            dirs  = dirs  | set(glob.glob(os.path.join(self.dirpath, "p_all", "p*", "0")))

        if build.SIMULATION in targets:
            restarts = set(glob.glob(os.path.join(self.dirpath, "restart_data", "*.dat")))
            restarts = restarts - set(glob.glob(os.path.join(self.dirpath, "restart_data", "lustre_0.dat")))
            restarts = restarts - set(glob.glob(os.path.join(self.dirpath, "restart_data", "lustre_*_cb.dat")))

            Ds = set(glob.glob(os.path.join(self.dirpath, "D", "*.dat")))
            Ds = Ds - set(glob.glob(os.path.join(self.dirpath, "D", "*.000000.dat")))

            files = files | restarts
            files = files | Ds

        if build.POST_PROCESS in targets:
            dirs.add("silo_hdf5")

        for relfile in files:
            if not os.path.isfile(relfile):
                relfile = os.path.join(self.dirpath, relfile)
            common.delete_file(relfile)

        for reldir in dirs:
            if not os.path.isdir(reldir):
                reldir = os.path.join(self.dirpath, reldir)
            common.delete_directory(reldir)


# Load the input file
def load(filepath: str = None, args: typing.List[str] = None, empty_data: dict = None, do_print: bool = True) -> MFCInputFile:
    if not filepath:
        if empty_data is None:
            raise common.MFCException("Please provide an input file.")

        input_file = MFCInputFile("empty.py", "empty.py", empty_data)
        input_file.validate_params()
        return input_file

    filename: str = filepath.strip()

    if do_print:
        cons.print(f"Acquiring [bold magenta]{filename}[/bold magenta]...")

    dirpath:    str  = os.path.abspath(os.path.dirname(filename))
    dictionary: dict = {}

    if not os.path.exists(filename):
        raise common.MFCException(f"Input file '{filename}' does not exist. Please check the path is valid.")

    if filename.endswith(".py"):
        (json_str, err) = common.get_py_program_output(filename, ["--mfc", json.dumps(ARGS())] + (args or []))

        if err != 0:
            raise common.MFCException(f"Input file {filename} terminated with a non-zero exit code. Please make sure running the file doesn't produce any errors.")
    elif filename.endswith(".json"):
        json_str = common.file_read(filename)
    else:
        raise common.MFCException("Unrecognized input file format. Only .py and .json files are supported. Please check the README and sample cases in the examples directory.")

    try:
        dictionary = json.loads(json_str)
    except Exception as exc:
        raise common.MFCException(f"Input file {filename} did not produce valid JSON. It should only print the case dictionary.\n\n{exc}\n")

    input_file = MFCInputFile(filename, dirpath, dictionary)
    input_file.validate_params(f"Input file {filename}")
    return input_file


load.CACHED_MFCInputFile = None
