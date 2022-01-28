import traceback

import internal.common      as common
import internal.configfiles as configfiles

from typing      import Any
from dataclasses import dataclass

@dataclass
class Compilers:
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data):
        self.c       = data["c"]
        self.cpp     = data["cpp"]
        self.fortran = data["fortran"]

@dataclass
class Configuration:
    name:  str
    flags: Compilers

    def __init__(self, data):
        self.name  = data["name"]
        self.flags = Compilers(data["flags"])

@dataclass
class Target_Download:
    link:    str
    version: str

    def __init__(self, data):
        self.link    = data["link"]
        self.version = data["version"]

@dataclass
class Target_Clone:
    git:  str
    hash: str

    def __init__(self, data):
        self.git  = data["git"]
        self.hash = data["hash"]

@dataclass
class Target_Source:
    source: str

    def __init__(self, data):
        self.source = data["source"]

@dataclass
class Target_Collection:
    def __init__(self, data):
        pass

@dataclass
class Target_Fetch:
    method: str
    params: Any

    def __init__(self, data):
        self.method = data["method"]

        if self.method == "download":
            self.params = Target_Download(data["params"])
        elif self.method == "clone":
            self.params = Target_Clone(data["params"])
        elif self.method == "source":
            self.params = Target_Source(data["params"])
        elif self.method == "collection":
            self.params = Target_Collection(data["params"])
        else:
            raise MFCException(f"[mfc.conf.yaml]: '{target}' - Unrecognized fetch method '{method}'.")


@dataclass
class Target:
    name:  str
    build: list
    test:  list
    clean: list
    fetch: Target_Fetch
    common_configuration: str

    def __init__(self, data):
        self.name    = data["name"]
        self.build   = data["build"]
        self.depends = data.get("depends", [])
        self.test    = data.get("test",    [])
        self.clean   = data.get("clean",   [])
        self.fetch   = Target_Fetch(data["fetch"])
        self.common_configuration = data.get("common_configuration", None)

class MFCConf:
    def __init__(self):
        data = configfiles.ConfigFileBase(common.MFC_CONF_FILEPATH, noexist_ok=False)

        self.compilers      = Compilers(data.tree_get("compilers"))
        self.configurations = [ Configuration(e) for e in data.tree_get("configurations") ]
        self.targets        = [ Target       (e) for e in data.tree_get("targets")        ]

    def get_configuration(self, name: str) -> Configuration:
        for configuration in self.configurations:
            if configuration.name == name:
                return configuration

        raise common.MFCException(f'MFCConf: Configuration "{name}" doesn\'t exist')

    def get_target_configuration_name(self, name: str, default: str) -> str:
        target = self.get_target(name)

        if target.common_configuration is not None:
            return target.common_configuration

        return default

    def get_target_configuration_folder_name(self, name: str, default: str) -> str:
        if self.get_target(name).common_configuration is not None:
            return "common"

        return default

    def get_target_configuration(self, name: str, default: str) -> Configuration:
        return self.get_configuration(self.get_target_configuration_name(name, default))

    def get_target_matches(self, name: str) -> list:
        return list(filter(lambda x: x.name == name, self.targets))

    def does_target_exist(self, name: str) -> bool:
        return len(self.get_target_matches(name)) > 0

    def does_unique_target_exist(self, name: str) -> bool:
        return len(self.get_target_matches(name)) == 1

    def get_target(self, name: str) -> Target:
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

            for dependency_name in desc.depends:
                result.append(dependency_name)

                if recursive:
                    result  += self.get_dependency_names(dependency_name, recursive=recursive, visited=visited)
                    visited += result

        return list(set(result))
