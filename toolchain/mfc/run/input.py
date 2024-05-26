import os, json, typing, dataclasses, jsonschema, traceback

from ..printer import cons
from ..        import common, build
from ..state   import ARGS
from ..case    import Case
from ..run     import case_dicts

@dataclasses.dataclass(init=False)
class MFCInputFile(Case):
    filename: str
    dirpath:  str

    def __init__(self, filename: str, dirpath: str, params: dict) -> None:
        super().__init__(params)
        # Typecheck parameters
        jsonschema.validate(self.params, case_dicts.SCHEMA)

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


    def clean(self, targets) -> None:
        for relfile in [
            "equations.dat", "run_time.inf", "time_data.dat",
            "io_time_data.dat", "fort.1"
        ] + [f"{build.get_target(target).name}.inp" for target in targets]:
            common.delete_file(os.path.join(self.dirpath, relfile))

        for reldir in ["D", "p_all", "silo_hdf5", "viz"]:
            common.delete_directory(os.path.join(self.dirpath, reldir))


# Load the input file
def load(filepath: str = None, args: typing.List[str] = None, empty_data: dict = None) -> MFCInputFile:
    if not filepath:
        if empty_data is None:
            raise common.MFCException("Please provide an input file.")

        return MFCInputFile("empty.py", "empty.py", empty_data)

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

    return MFCInputFile(filename, dirpath, dictionary)


load.CACHED_MFCInputFile = None
