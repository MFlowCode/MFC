import internal.common      as common
import internal.configfiles as configfiles

class MFCConf(configfiles.ConfigFileBase):
    def __init__(self):
        super().__init__(common.MFC_CONF_FILEPATH, noexist_ok=False)

    def get_configuration(self, name: str):
        for configuration in self["configurations"]:
            if configuration["name"] == name:
                return configuration

        return common.MFCException(f'MFCConf: Configuration "{name}" doesn\'t exist')

    def get_target_configuration_name(self, name: str, default: str):
        target = self.get_target(name)

        if "common_configuration" in target:
            return target["common_configuration"]

        return default

    def get_target_configuration_folder_name(self, name: str, default: str):
        if "common_configuration" in self.get_target(name):
            return "common"

        return default

    def get_target_configuration(self, name: str, default: str):
        return self.get_configuration(self.get_target_configuration_name(name, default))

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
