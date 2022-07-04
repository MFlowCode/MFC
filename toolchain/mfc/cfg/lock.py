import mfc.util.common as common
from mfc.cfg.user import MFCUser

import os
import dataclasses

@dataclasses.dataclass
class MFCLock:
    mode: str

    def __init__(self, user: MFCUser):
        if not os.path.exists(common.MFC_LOCK_FILEPATH):
            common.create_file(common.MFC_LOCK_FILEPATH)
            common.file_dump_yaml(common.MFC_LOCK_FILEPATH, {
                "mode": user.modes[0].name
            })
        
        data: dict = common.file_load_yaml(common.MFC_LOCK_FILEPATH)
        self.mode = data["mode"]

    def write(self):
        common.file_dump_yaml(common.MFC_LOCK_FILEPATH, dataclasses.asdict(self))
