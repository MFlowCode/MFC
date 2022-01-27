import copy

import internal.common as common

# TODO: NOT VERY DRY CODE

class ObjectTree:
    def __init__(self, data={}):
        self.data = data

    def exists(self, key_sequence, curr=None):
        if curr is None:
            curr = self.data

        if isinstance(key_sequence, str):
            key_sequence = [key_sequence]

        key = key_sequence[0]
        if key not in curr:
            return False

        if len(key_sequence) == 1:
            return True

        return self.exists(key_sequence[1:], curr[key])

    def set(self, key_sequence, value, curr=None):
        if curr is None:
            curr = self.data

        if isinstance(key_sequence, str):
            key_sequence = [key_sequence]

        key = key_sequence[0]
        if key not in curr:
            raise common.MFCException(f'ObjectTree: Key "{key}" doesn\'t exist.')

        if len(key_sequence) == 1:
            curr[key] = value
            return

        return self.set(key_sequence[1:], value, curr[key])

    def __setitem__(self, key_sequence, value):
        self.set(key_sequence, value)

    def get(self, key_sequence, default=None, curr=None):
        if curr is None:
            curr = self.data

        if isinstance(key_sequence, str):
            key_sequence = [key_sequence]

        key = key_sequence[0]

        if key not in curr:
            if default is None:
                raise common.MFCException(f'ObjectTree: Key "{key}" doesn\'t exist.')

            return default

        if len(key_sequence) == 1:
            return curr[key]

        return self.get(key_sequence[1:], default, curr[key])

    def __getitem__(self, key_sequence, default=None):
        return self.get(key_sequence, default)
