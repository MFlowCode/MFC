import os, dataclasses

from .        import build, state, common
from .state   import MFCConfig
from .printer import cons


MFC_LOCK_CURRENT_VERSION: int = 3


@dataclasses.dataclass
class MFCLockData:
    config:  MFCConfig
    version: int


data: MFCLockData = None


def init():
    global data
    
    if not os.path.exists(common.MFC_LOCK_FILEPATH):
        config = MFCConfig()
        data   = MFCLockData(config, MFC_LOCK_CURRENT_VERSION)
        state.gCFG = config

        common.create_file(common.MFC_LOCK_FILEPATH)
        write()
    else:
        load()


def load():
    global data
    
    d = common.file_load_yaml(common.MFC_LOCK_FILEPATH)
    
    # 0 is the default version in order to accommodate versions of mfc.sh
    # prior to the introduction of the "version" attribute to the lock file.

    if d["version"] < MFC_LOCK_CURRENT_VERSION:
        raise common.MFCException(f"""\
There has been a breaking change to the MFC build system. Please delete your \
build/ directory and run MFC again. (v{d["version"]} -> v{MFC_LOCK_CURRENT_VERSION}).\
""")

    config = MFCConfig.from_dict(d["config"])
    data   = MFCLockData(config, d["version"])
    state.gCFG = config


def write():
    global data
    
    common.file_dump_yaml(common.MFC_LOCK_FILEPATH, dataclasses.asdict(data))


def switch(to: MFCConfig):
    global data
    
    if to == data.config:
        return

    cons.print(f"[bold yellow]Switching from {data.config} to {to}[/bold yellow]")
    cons.print("")

    data.config = to
    state.gCFG  = to
    write()

    for target_name in build.get_mfc_target_names():
        dirpath = build.get_build_dirpath(build.get_target(target_name))
        cons.print(f"[bold red]Removing {os.path.relpath(dirpath)}[/bold red]")
        common.delete_directory(dirpath)
