import os

import internal.common as common

class MFCLock:
    def __init__(self):
        if not os.path.exists(common.MFC_LOCK_FILEPATH):
            with open(common.MFC_LOCK_FILEPATH, 'w') as f:
                f.write("targets: []")

        self.data = common.file_load_yaml(common.MFC_LOCK_FILEPATH)

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise common.MFCException(f'MFCLock: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]

    def get_target_matches(self, name: str, restrict_cc: str = None):
        def peek_filter(e: dict):
            if e["name"] != name:
                return False

            if restrict_cc is None:
                return True

            return e.get("compiler_configuration", None) == restrict_cc

        return list(filter(peek_filter, self["targets"]))

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

    def save(self):
        common.file_dump_yaml(common.MFC_LOCK_FILEPATH, self.data)
