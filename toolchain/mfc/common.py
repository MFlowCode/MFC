import os, yaml, typing, shutil, subprocess, dataclasses

from .printer import cons

from os.path import abspath, normpath, dirname, realpath


MFC_ROOTDIR       = normpath(f"{dirname(realpath(__file__))}/../..")
MFC_TESTDIR       = abspath(f"{MFC_ROOTDIR}/tests")
MFC_SUBDIR        = abspath(f"{MFC_ROOTDIR}/build")
MFC_DEV_FILEPATH  = abspath(f"{MFC_ROOTDIR}/toolchain/mfc.dev.yaml")
MFC_USER_FILEPATH = abspath(f"{MFC_ROOTDIR}/defaults.yaml")
MFC_LOCK_FILEPATH = abspath(f"{MFC_SUBDIR}/lock.yaml")

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


def system(command: typing.List[str], no_exception: bool = False, exception_text=None, on_error=lambda: None, cwd=None, stdout=None, stderr=None) -> int:
    cmd = [ str(x) for x in command if not isspace(str(x)) ]

    if stdout != subprocess.DEVNULL:
        cons.print(no_indent=True)
    
    cons.print(f"$ {' '.join(cmd)}")
    
    if stdout != subprocess.DEVNULL:
        cons.print(no_indent=True)

    r = subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=stderr)

    if r.returncode != 0:
        on_error()

        if not(no_exception):
            if exception_text is None:
                exception_text = f'Failed to execute command "{command}".'

            raise MFCException(exception_text)

    return r.returncode


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
    dirpath  = os.path.abspath (os.path.dirname(filepath))
    filename = os.path.basename(filepath)

    proc = subprocess.Popen(["python3", filename], cwd=dirpath,
                            stdout=subprocess.PIPE)

    return (proc.communicate()[0], proc.returncode)


def isspace(s: str) -> bool:
    """
    Returns whether a string, s, is empty, whitespace, or None.
    """

    if s == None:
        return True

    return len(s.strip()) == 0


def does_command_exist(s: str) -> bool:
    """
    Returns whether a program/alias is accessible.
    """

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
    """
    Returns the first element in array that satisfies the predicate.
    """

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

