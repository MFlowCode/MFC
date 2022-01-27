import os

import internal.common      as common
import internal.configfiles as configfiles

class MFCLock(configfiles.ConfigFileBase):
    def __init__(self):
        super().__init__(common.MFC_LOCK_FILEPATH, {"targets": []})

        self.data = common.file_load_yaml(common.MFC_LOCK_FILEPATH)

    def was_target_built(self, name: str, restrict_cc: str = None):
        matches = self.get_target_matches(name, restrict_cc)

    def get_target_matches(self, name: str, restrict_cc: str = None):
        def peek_filter(e: dict):
            if e["name"] != name:
                return False

            if restrict_cc is None:
                return True

            return e.get("compiler_configuration", None) == restrict_cc

        return list(filter(peek_filter, self.data["targets"]))

    def does_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) > 0

    def does_unique_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) == 1

    def get_target(self, name: str, restrict_cc: str = None):
        matches = self.get_target_matches(name, restrict_cc)

        if len(matches) == 0:
            raise common.MFCException(f'Failed to retrieve dependency "{name}" with restrict_cc="{restrict_cc}".')

        if len(matches) > 1:
            raise common.MFCException(f'More than one dependency to choose from for "{name}" with restrict_cc="{restrict_cc}".')

        return matches[0]
