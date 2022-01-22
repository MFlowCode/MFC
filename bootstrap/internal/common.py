import os
import re
import sys
import yaml     # *: PyYAML package
import shutil
import colorama # *: Colorama package
import tarfile


MFC_ROOTDIR       = os.path.normpath(f"{os.path.dirname(os.path.realpath(__file__))}/../..")
MFC_SUBDIR        = f"{MFC_ROOTDIR}/.mfc"
MFC_CONF_FILEPATH = f"{MFC_ROOTDIR}/mfc.conf.yaml"
MFC_LOCK_FILEPATH = f"{MFC_SUBDIR}/mfc.lock.yaml"


class MFCException(Exception):
    pass


def execute_shell_command_safe(command: str, no_exception: bool = False) -> int:
    status = os.system(command)

    if status != 0 and not(no_exception):
        raise MFCException(f'Failed to execute command "{command}".')

    return status


def clear_line() -> None:
    sys.stdout.write("\033[K")


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


# TODO: Better solution
def uncompress_archive_to(archive_filepath: str, destination: str) -> None:
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    src = "/".join(destination.rstrip("/").split("/")[:-1])+"/"+archive_filepath.split("/")[-1].replace(".tar.gz", "")

    os.rename(src, destination)


def delete_file_safe(filepath: str) -> None:
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file_safe(filepath: str) -> None:
    if not os.path.exists(filepath):
        open(filepath, "w").close()


def delete_directory_recursive_safe(directory_path: str) -> None:
    if os.path.isdir(directory_path):
        shutil.rmtree(directory_path)


def create_directory_safe(directory_path: str) -> None:
    os.makedirs(directory_path, exist_ok=True)


def center_ansi_escaped_text(message: str) -> str:
    nCols = shutil.get_terminal_size((80, 20)).columns

    to_escape = [re.escape(colorama.Style.RESET_ALL)]
    for key in dir(colorama.Fore):
        if not callable(getattr(colorama.Fore, key)) and not key.startswith("__"):
            to_escape.append(re.escape(getattr(colorama.Fore, key)))

    longest_string_len = max([len(re.compile("|".join(to_escape), flags=re.DOTALL).sub("", line)) for line in message.splitlines()])

    padding = " "*((nCols - longest_string_len) // 2)

    return "\n".join([f'{padding}{line}{padding}' for line in message.splitlines()])


def clear_print(message, end='\n') -> None:
    clear_line()
    print(message, end=end)


def update_symlink(at: str, to: str) -> None:
    if os.path.islink(at):
        os.remove(at)

    # Use a relative symlink so that users can
    # move and rename MFC's root folder
    os.symlink(os.path.relpath(to, MFC_SUBDIR), at)
