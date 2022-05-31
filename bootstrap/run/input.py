
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
                continue
                
            if key not in case_dicts.PRE_PROCESS  and \
               key not in case_dicts.POST_PROCESS and \
               key not in case_dicts.SIMULATION:
                raise common.MFCException(f"MFCInputFile::dump: Case parameter '{key}' is not used by any MFC code. Please check your spelling or add it as a new parameter.")

        contents = f"&user_inputs\n{dict_str}&end/\n"

        # Save .inp input file
        common.file_write(f"{self.case_dirpath}/{target_name}.inp", contents)


# Load from Python input file
def load(filename: str) -> MFCInputFile:
    dirpath:    str  = os.path.abspath(os.path.dirname(filename))
    dictionary: dict = {}

    rich.print(f"> > Fetching case dictionary from {input}...")

    if not filename.endswith(".py"):
        raise common.MFCException("Unrecognized input file format. Only .py files are supported. Please check the README and sample cases in the samples directory.")

    if not os.path.exists(filename):
        raise common.MFCException(f"Input file '{filename}' does not exist. Please check the path is valid.")

    (output, err) = common.get_py_program_output(filename)

    if err != 0:
        raise common.MFCException(f"Input file {filename} terminated with a non-zero exit code. Please make sure running the file doesn't produce any errors.")

    try:
        dictionary = json.loads(output)
    except Exception as exc:
        raise common.MFCException(f"Input file {filename} did not produce valid JSON. It should only print the case dictionary.\n\n{exc}\n")

    return MFCInputFile(filename, dirpath, dictionary)
