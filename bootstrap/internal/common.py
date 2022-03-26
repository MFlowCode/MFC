import os
import sys
import yaml       # *: PyYAML package
import shutil
import tarfile
import subprocess

from datetime import datetime


MFC_ROOTDIR         = os.path.normpath(f"{os.path.dirname(os.path.realpath(__file__))}/../..")
MFC_TESTDIR         = os.path.abspath(f"{MFC_ROOTDIR}/tests")
MFC_SUBDIR          = os.path.abspath(f"{MFC_ROOTDIR}/build")
MFC_DEV_FILEPATH    = os.path.abspath(f"{MFC_ROOTDIR}/bootstrap/mfc.dev.yaml")
MFC_USER_FILEPATH   = os.path.abspath(f"{MFC_ROOTDIR}/mfc.user.yaml")
MFC_LOCK_FILEPATH   = os.path.abspath(f"{MFC_SUBDIR}/mfc.lock.yaml")


class MFCException(Exception):
    pass


def execute_shell_command(command: str, no_exception: bool = False, exception_text=None, on_error=lambda: None) -> int:
    status = os.system(command)

    if status != 0:
        on_error()

        if not(no_exception):
            if exception_text is None:
                exception_text = f'Failed to execute command "{command}".'

            raise MFCException(exception_text)

    return status

def get_datetime_str() -> str:
    return datetime.now().strftime('%d/%m/%Y %H:%M:%S')

def clear_line() -> None:
    sys.stdout.write("\033[K")


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


# TODO: Better solution
def uncompress_archive_to(archive_filepath: str, destination: str) -> None:
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    src = "/".join(destination.rstrip("/").split("/")[:-1])+"/"+archive_filepath.split("/")[-1].replace(".tar.gz", "")

    os.rename(src, destination)


def delete_file(filepath: str) -> None:
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file(filepath: str) -> None:
    if not os.path.exists(filepath):
        try:
            open(filepath, "w").close()
        except IOError as exc:
            raise MFCException(f"Failed to create file {filepath}: {exc}")

def delete_directory_recursive(directory_path: str) -> None:
    if os.path.isdir(directory_path):
        shutil.rmtree(directory_path)


def create_directory(directory_path: str) -> None:
    os.makedirs(directory_path, exist_ok=True)


def clear_print(message, end='\n') -> None:
    clear_line()
    print(message, end=end)


def update_symlink(at: str, to: str) -> None:
    if os.path.islink(at):
        os.remove(at)

    # Use a relative symlink so that users can
    # move and rename MFC's root folder
    os.symlink(os.path.relpath(to, MFC_SUBDIR), at)

def get_py_program_output(filepath: str):
    dirpath:  str = os.path.abspath (os.path.dirname(filepath))
    filename: str = os.path.basename(filepath)

    command: str = f"cd {dirpath} && python3 {filename} 2>&1"

    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

    return (proc.stdout, proc.returncode)
