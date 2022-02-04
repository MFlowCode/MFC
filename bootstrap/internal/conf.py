import internal.common      as common
import internal.configfiles as configfiles

from typing      import Any
from dataclasses import dataclass

@dataclass
class CompilerVersion:
    name:    str
    is_used: str
    fetch:   str
    minimum: str

    def __init__(self, data) -> None:
        self.name    = data.get("name")
        self.is_used = data.get("is_used")
        self.fetch   = data.get("fetch")
        self.minimum = data.get("minimum")

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
            raise common.MFCException(f"[mfc.conf.yaml]: Unrecognized fetch method '{self.method}'.")

@dataclass
class Target_Build_Options:
    def __init__(self, data) -> None:
        pass

@dataclass
class Target:
    name:  str
    build: str
    test:  list
    clean: list
    fetch: Target_Fetch
    common_configuration: str

    def __init__(self, data):
        self.name    = data["name"]
        self.build   = data.get("build",   [])
        self.depends = data.get("depends", [])
        self.test    = data.get("test",    [])
        self.clean   = data.get("clean",   [])
        self.fetch   = Target_Fetch(data["fetch"])
        self.common_configuration = data.get("common_configuration", None)

class MFCConf:
    def __init__(self):
        data = configfiles.ConfigFileBase(common.MFC_DEV_FILEPATH, noexist_ok=False)

        self.compiler_verions = [ CompilerVersion(e) for e in data.tree_get("compiler_versions") ]
        self.targets          = [ Target(e)          for e in data.tree_get("targets")           ]

    def get_target_configuration_name(self, name: str, default: str) -> str:
        target = self.get_target(name)

        if target.common_configuration is not None:
            return target.common_configuration

        return default

    def get_target_configuration_folder_name(self, name: str, default: str) -> str:
        if self.get_target(name).common_configuration is not None:
            return "common"

        return default
    
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
