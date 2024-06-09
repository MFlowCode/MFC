import os, json, glob, typing, dataclasses

from ..printer import cons
from ..        import common, build
from ..state   import ARGS
from ..case    import Case

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

    def generate_fpp(self, target) -> None:
        if target.isDependency:
            return

        cons.print(f"Generating [magenta]case.fpp[/magenta].")
        cons.indent()

        # Case FPP file
        self.__save_fpp(target, self.get_fpp(target))

        cons.unindent()

    # Generate case.fpp & [target.name].inp
    def generate(self, target) -> None:
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
def load(filepath: str = None, args: typing.List[str] = None, empty_data: dict = None) -> MFCInputFile:
    if not filepath:
        if empty_data is None:
            raise common.MFCException("Please provide an input file.")

        input_file = MFCInputFile("empty.py", "empty.py", empty_data)
        input_file.validate_params()
        return input_file

    filename: str = filepath.strip()

    cons.print(f"Acquiring [bold magenta]{filename}[/bold magenta]...")

    dirpath:    str  = os.path.abspath(os.path.dirname(filename))
    dictionary: dict = {}

    if not os.path.exists(filename):
        raise common.MFCException(f"Input file '{filename}' does not exist. Please check the path is valid.")

    if filename.endswith(".py"):
        (json_str, err) = common.get_py_program_output(filename, [json.dumps(ARGS())] + (args or []))

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
