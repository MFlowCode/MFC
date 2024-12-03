import os, yaml, typing, shutil, subprocess

from os.path import join, abspath, normpath, dirname, realpath

from .printer import cons


MFC_ROOT_DIR       = abspath(normpath(f"{dirname(realpath(__file__))}/../.."))
MFC_TEST_DIR       = abspath(join(MFC_ROOT_DIR, "tests"))
MFC_BUILD_DIR      = abspath(join(MFC_ROOT_DIR, "build"))
MFC_TOOLCHAIN_DIR  = abspath(join(MFC_ROOT_DIR, "toolchain"))
MFC_EXAMPLE_DIRPATH = abspath(join(MFC_ROOT_DIR, "examples"))
MFC_LOCK_FILEPATH  = abspath(join(MFC_BUILD_DIR, "lock.yaml"))
MFC_TEMPLATE_DIR   = abspath(join(MFC_TOOLCHAIN_DIR, "templates"))
MFC_BENCH_FILEPATH = abspath(join(MFC_TOOLCHAIN_DIR, "bench.yaml"))
MFC_MECHANISMS_DIR = abspath(join(MFC_TOOLCHAIN_DIR, "mechanisms"))

MFC_LOGO = """\
     .=++*:          -+*+=.
    :+   -*-        ==   =* .
  :*+      ==      ++    .+-
 :*##-.....:*+   .#%+++=--+=:::.
 -=-++-======#=--**+++==+*++=::-:.
.:++=----------====+*= ==..:%.....
 .:-=++++===--==+=-+=   +.  :=
 +#=::::::::=%=. -+:    =+   *:
.*=-=*=..    :=+*+:      -...--
"""


class MFCException(Exception):
    pass

def system(command: typing.List[str], print_cmd = None, **kwargs) -> subprocess.CompletedProcess:
    cmd = [ str(x) for x in command if not isspace(str(x)) ]

    if print_cmd in [True, None]:
        cons.print(f"$ {' '.join(cmd)}")
        cons.print(no_indent=True)

    return subprocess.run(cmd, **kwargs, check=False)


def file_write(filepath: str, content: str, if_different: bool = False):
    try:
        if if_different and os.path.isfile(filepath):
            with open(filepath, "r") as f:
                if f.read() == content:
                    return

        with open(filepath, "w") as f:
            f.write(content)
    except IOError as exc:
        raise MFCException(f'Failed to write to "{filepath}": {exc}') from exc


def file_read(filepath: str):
    try:
        with open(filepath, "r") as f:
            return f.read()
    except IOError as exc:
        raise MFCException(f'Failed to read from "{filepath}": {exc}') from exc


def file_load_yaml(filepath: str):
    try:
        with open(filepath, "r") as f:
            return yaml.safe_load(f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to load YAML from "{filepath}": {exc}') from exc


def file_dump_yaml(filepath: str, data) -> None:
    try:
        with open(filepath, "w") as f:
            yaml.dump(data, f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to dump YAML to "{filepath}": {exc}.') from exc


def delete_file(filepath: str) -> None:
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file(filepath: str) -> None:
    if not os.path.exists(filepath):
        try:
            with open(filepath, "w"):
                pass
        except IOError as exc:
            raise MFCException(f"Failed to create file {filepath}: {exc}") from exc


def create_directory(dirpath: str) -> None:
    os.makedirs(dirpath, exist_ok=True)


def delete_directory(dirpath: str) -> None:
    if os.path.isdir(dirpath):
        shutil.rmtree(dirpath)


def get_program_output(arguments: typing.List[str] = None, cwd=None):
    with subprocess.Popen([ str(_) for _ in arguments ] or [], cwd=cwd, stdout=subprocess.PIPE) as proc:
        return (proc.communicate()[0].decode(), proc.returncode)


def get_py_program_output(filepath: str, arguments: typing.List[str] = None):
    dirpath  = os.path.abspath (os.path.dirname(filepath))
    filename = os.path.basename(filepath)

    return get_program_output(["python3", filename] + arguments, cwd=dirpath)


def isspace(s: str) -> bool:
    """
    Returns whether a string, s, is empty, whitespace, or None.
    """

    if s is None:
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
    return 0 == os.system(f"command -v '{s}' > /dev/null 2>&1")


def format_list_to_string(arr: list, item_style=None, empty=None):
    if empty is None:
        empty = "nothing"

    pre, post = "", ""
    if item_style is not None:
        pre  = f"[{item_style}]"
        post = f"[/{item_style}]"

    if len(arr) == 0:
        return f"{pre}{empty}{post}"

    if len(arr) == 1:
        return f"{pre}{arr[0]}{post}"

    if len(arr) == 2:
        return f"{pre}{arr[0]}{post} and {pre}{arr[1]}{post}"

    lhs = ', '.join([ f"{pre}{e}{post}" for e in arr[:-1]])
    rhs = f", and {pre}{arr[-1]}{post}"

    return lhs + rhs


def find(predicate, arr: list):
    """
    Returns the first element in array that satisfies the predicate.
    """

    for index, item in enumerate(arr):
        if predicate(index, item):
            return index, item

    return None, None


# pylint: disable=redefined-builtin
def quit(sig):
    os.kill(os.getpid(), sig)


def does_system_use_modules() -> bool:
    """
    Returns True if the system uses modules.
    """

    return does_command_exist("module")


def is_number(x: str) -> bool:
    if x is None:
        return False

    if isinstance(x, (int, float)):
        return True

    try:
        float(x)
        return True
    except ValueError:
        return False


def get_cpuinfo():
    if does_command_exist("lscpu"):
        # Linux
        with subprocess.Popen(['lscpu'], stdout=subprocess.PIPE, universal_newlines=True) as proc:
            output = f"From lscpu\n{proc.communicate()[0]}"
    elif does_command_exist("sysctl"):
        # MacOS
        with subprocess.Popen(['sysctl', '-a'], stdout=subprocess.PIPE) as proc1:
            with subprocess.Popen(['grep', 'machdep.cpu'], stdin=proc1.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  universal_newlines=True) as proc2:
                proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
                output = f"From sysctl -a \n{proc2.communicate()[0]}"
    else:
        output = "No CPU info found"

    return f"CPU Info:\n{output}"

def generate_git_tagline() -> str:
    if not does_command_exist("git"):
        return "Could not find git"

    rev    = system(["git", "rev-parse",                 "HEAD"], print_cmd=False, stdout=subprocess.PIPE).stdout.decode().strip()
    branch = system(["git", "rev-parse", "--abbrev-ref", "HEAD"], print_cmd=False, stdout=subprocess.PIPE).stdout.decode().strip()
    dirty  = "dirty" if system(["git", "diff", "--quiet"], print_cmd=False).returncode != 0 else "clean"

    return f"{rev} on {branch} ({dirty})"
