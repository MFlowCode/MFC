import os

import internal.common     as common
import internal.objecttree as objecttree

class ConfigFileBase(objecttree.ObjectTree):
    def __init__(self, filepath: str, emptyData = {}, noexist_ok: bool = True):
        self.filepath = filepath

        if not os.path.exists(self.filepath):
            if not noexist_ok:
                raise common.MFCException(f'ConfigFile: File "{self.filepath}" doesn\'t exist.')

            common.file_dump_yaml(self.filepath, emptyData)

        super().__init__(common.file_load_yaml(self.filepath))

    def save(self):
        common.file_dump_yaml(self.filepath, self.data)

    def __setitem__(self, key: str, value):
        self.data[key] = value

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise common.MFCException(f'ConfigFile: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]
