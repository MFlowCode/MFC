
import os
import json
import dataclasses

import rich

import common

import run.case_dicts as case_dicts

@dataclasses.dataclass
class MFCInputFile:
    filename:     str
    case_dirpath: str
    case_dict:    dict

    # Generate .inp input file.
    def dump(self,  target_name: str) -> None:
        MASTER_KEYS: list = case_dicts.get_input_dict_keys(target_name)

        # Create Fortran-style input file content string
        dict_str = ""
        for key,val in self.case_dict.items():
            if key in MASTER_KEYS:
                dict_str += f"{key} = {val}\n"

        contents = f"&user_inputs\n{dict_str}&end/\n"

        # Save .inp input file
        common.file_write(f"{self.case_dirpath}/{target_name}.inp", contents)


# Load from Python input file
def load(filename: str) -> MFCInputFile:
    dirpath:    str  = os.path.abspath(os.path.dirname(filename))
    dictionary: dict = {}

    rich.print(f"> > Fetching case dictionary from {input}...")

    if not filename.endswith(".py"):
        raise common.MFCException("Unrecognized input file format. Only .py files are supported.")

    if not os.path.exists(filename):
        raise common.MFCException(f"Input file '{filename}' does not exist.")

    (output, err) = common.get_py_program_output(filename)

    if err != 0:
        raise common.MFCException(f"Input file {filename} terminated with a non-zero exit code.")

    try:
        dictionary = json.loads(output)
    except Exception as exc:
        raise common.MFCException(f"Input file {filename} did not produce valid JSON. It should only print the case dictionary.\n\n{exc}\n")

    return MFCInputFile(filename, dirpath, dictionary)
