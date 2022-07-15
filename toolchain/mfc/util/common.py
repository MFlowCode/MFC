import os
import yaml
import enum
import typing
import shutil
import subprocess


MFC_ROOTDIR         = os.path.normpath(f"{os.path.dirname(os.path.realpath(__file__))}/../../..")
MFC_TESTDIR         = os.path.abspath(f"{MFC_ROOTDIR}/tests")
MFC_SUBDIR          = os.path.abspath(f"{MFC_ROOTDIR}/build")
MFC_DEV_FILEPATH    = os.path.abspath(f"{MFC_ROOTDIR}/toolchain/mfc.dev.yaml")
MFC_USER_FILEPATH   = os.path.abspath(f"{MFC_ROOTDIR}/mfc.user.yaml")
MFC_LOCK_FILEPATH   = os.path.abspath(f"{MFC_SUBDIR}/mfc.lock.yaml")

MFC_LOGO = f"""\
     ___            ___          ___
    /__/\          /  /\        /  /\\
   |  |::\        /  /:/_      /  /:/
   |  |:|:\      /  /:/ /\    /  /:/
 __|__|:|\:\    /  /:/ /:/   /  /:/  ___
/__/::::| \:\  /__/:/ /:/   /__/:/  /  /\\
\  \:\~~\__\/  \  \:\/:/    \  \:\ /  /:/
 \  \:\         \  \::/      \  \:\  /:/
  \  \:\         \  \:\       \  \:\/:/
   \  \:\         \  \:\       \  \::/
    \__\/          \__\/        \__\/
"""


class MFCException(Exception):
    pass


def system(command: str, no_exception: bool = False, exception_text=None, on_error=lambda: None) -> int:
    status = os.system(command)

    if status != 0:
        on_error()

        if not(no_exception):
            if exception_text is None:
                exception_text = f'Failed to execute command "{command}".'

            raise MFCException(exception_text)

    return status


def file_write(filepath: str, content: str):
    try:
        with open(filepath, "w") as f:
            f.write(content)
    except IOError as exc:
        raise MFCException(f'Failed to write to "{filepath}": {exc}')


def file_read(filepath: str):
    try:
        with open(filepath, "r") as f:
            return f.read()
    except IOError as exc:
        raise MFCException(f'Failed to read from "{filepath}": {exc}')


def file_load_yaml(filepath: str):
    try:
        with open(filepath, "r") as f:
            return yaml.safe_load(f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to load YAML from "{filepath}": {exc}')


def file_dump_yaml(filepath: str, data) -> None:
    try:
        with open(filepath, "w") as f:
            yaml.dump(data, f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to dump YAML to "{filepath}": {exc}.')


def delete_file(filepath: str) -> None:
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file(filepath: str) -> None:
    if not os.path.exists(filepath):
        try:
            open(filepath, "w").close()
        except IOError as exc:
            raise MFCException(f"Failed to create file {filepath}: {exc}")


def create_directory(dirpath: str) -> None:
    os.makedirs(dirpath, exist_ok=True)


def delete_directory(dirpath: str) -> None:
    if os.path.isdir(dirpath):
        shutil.rmtree(dirpath)


def get_py_program_output(filepath: str):
    dirpath:  str = os.path.abspath (os.path.dirname(filepath))
    filename: str = os.path.basename(filepath)

    command: str = f"cd {dirpath} && python3 {filename} 2>&1"

    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

    return (proc.stdout, proc.returncode)


def isspace(s: str) -> bool:
    if s == None:
        return True

    return len(s.strip()) == 0


def does_command_exist(s: str) -> bool:
    # If system is (not) Posix compliant
    if shutil.which("command") is None:
        return shutil.which(s) is not None

    # command -v is useful because it allows for finding much more
    # than with "which". Checkout "man command"
    return 0 == os.system(f"command -v {s} > /dev/null 2>&1")


def format_list_to_string(arr: list, empty = "nothing"):
    if len(arr) == 0:
        return empty

    if len(arr) == 1:
        return arr[0]

    if len(arr) == 2:
        return f"{arr[0]} and {arr[1]}"

    lhs = ', '.join(arr[:-1])
    rhs = f", and {arr[-1]}"

    return lhs + rhs


def find(predicate, arr: list):
    for index, item in enumerate(arr):
        if predicate(index, item):
            return index, item

    return None, None


def quit(sig):
    os.kill(os.getpid(), sig)


def does_system_use_modules() -> bool:
    """
    Returns True if the system uses modules.
    """

    return does_command_exist("module")


def get_loaded_modules() -> typing.List[str]:
    """
    Returns a list of loaded modules.
    """

    return subprocess.getoutput("module -t list").splitlines()


def enumint(x: typing.Union[int, enum.Enum]) -> int:
    return x.value if isinstance(x, enum.Enum) else x


def enumeq(lhs: typing.Union[enum.Enum, int], rhs: typing.Union[enum.Enum, int]) -> bool:
    return enumint(lhs) == enumint(rhs)


def enumne(lhs: typing.Union[enum.Enum, int], rhs: typing.Union[enum.Enum, int]) -> bool:
    return not enumeq(lhs, rhs)

