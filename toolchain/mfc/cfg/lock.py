from ..util     import common
from ..cfg.user import MFCUser

import os
import dataclasses

MFC_LOCK_CURRENT_VERSION: int = 2

@dataclasses.dataclass
class MFCLock:
    mode:    str
    mpi:     bool
    version: int

    def __init__(self, user: MFCUser):
        if not os.path.exists(common.MFC_LOCK_FILEPATH):
            common.create_file(common.MFC_LOCK_FILEPATH)
            common.file_dump_yaml(common.MFC_LOCK_FILEPATH, {
                "mpi":     True,
                "mode":    user.modes[0].name,
                "version": MFC_LOCK_CURRENT_VERSION
            })

        data: dict = common.file_load_yaml(common.MFC_LOCK_FILEPATH)

        self.version = int(data.get("version", "0"))

        # 0 is the default version in order to accommodate versions of mfc.sh
        # prior to the introduction of the "version" attribute to the lock file.

        if self.version < MFC_LOCK_CURRENT_VERSION:
            raise common.MFCException(f"""\
There has been a breaking change to the MFC build system. Please delete your \
build/ directory and run MFC again. (v{self.version} -> v{MFC_LOCK_CURRENT_VERSION}).\
""")

        self.mode    = data["mode"]
        self.mpi     = data['mpi']


    def write(self):
        common.file_dump_yaml(common.MFC_LOCK_FILEPATH, dataclasses.asdict(self))

