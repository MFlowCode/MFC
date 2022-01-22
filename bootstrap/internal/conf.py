import internal.common as common

class MFCConf:
    def __init__(self):
        self.data = common.file_load_yaml(common.MFC_CONF_FILEPATH)

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise common.MFCException(f'MFCConf: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]

    def get_target_matches(self, name: str):
        return list(filter(lambda x: x["name"] == name, self["targets"]))

    def does_target_exist(self, name: str):
        return len(self.get_target_matches(name)) > 0

    def does_unique_target_exist(self, name: str):
        return len(self.get_target_matches(name)) == 1

    def get_target(self, name: str):
        matches = self.get_target_matches(name)

        if len(matches) == 0:
            raise common.MFCException(f'Failed to retrieve dependency "{name}".')

        if len(matches) > 1:
            raise common.MFCException(f'More than one dependency to choose from for "{name}".')

        return matches[0]

    def get_dependency_names(self, name: str, recursive=False, visited: list = None):
        result: list = []

        if visited == None:
            visited = []

        if name not in visited:
            visited.append(name)

            desc = self.get_target(name)

            for dependency_name in desc.get("depends", []):
                result.append(dependency_name)

                if recursive:
                    result  += self.get_dependency_names(dependency_name, recursive=recursive, visited=visited)
                    visited += result

        return list(set(result))
