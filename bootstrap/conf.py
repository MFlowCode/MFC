import os, typing, dataclasses

import common


@dataclasses.dataclass
class Compiler:
    name:        str
    is_used_cmd: str
    get_version: str
    min_version: str
    supported:   bool

    def __init__(self, data: dict) -> None:
        self.name        = data.get("name")
        self.is_used_cmd = data.get("is_used_cmd")
        self.get_version = data.get("get_version", "")
        self.min_version = data.get("min_version",  0)
        self.supported   = data.get("supported")

@dataclasses.dataclass
class Target_Download:
    link:    str
    version: str

    def __init__(self, data: dict):
        self.link    = data["link"]
        self.version = data["version"]

@dataclasses.dataclass
class Target_Clone:
    git:  str
    hash: str

    def __init__(self, data: dict):
        self.git  = data["git"]
        self.hash = data["hash"]

@dataclasses.dataclass
class Target_Source:
    source: str
    check:  str

    def __init__(self, data: dict):
        self.source = data.get("source")
        self.check  = data.get("check")

@dataclasses.dataclass
class Target_Collection:
    def __init__(self, data):
        pass

@dataclasses.dataclass
class Target_Fetch:
    method: str
    params: typing.Any

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

@dataclasses.dataclass
class Target_Build_Options:
    def __init__(self, data) -> None:
        pass

@dataclasses.dataclass
class Target:
    name:  str
    build: str
    test:  list
    clean: list
    fetch: Target_Fetch
    common_mode: str

    def __init__(self, data):
        self.name    = str.lower(data["name"])
        self.build   = data.get("build",   [])
        self.depends = [str.lower(dep) for dep in data.get("depends", [])]
        self.test    = data.get("test",    [])
        self.clean   = data.get("clean",   [])
        self.fetch   = Target_Fetch(data["fetch"])
        self.common_mode = data.get("common_mode", None)

class MFCConf:
    def __init__(self, mfc):
        data = common.file_load_yaml(common.MFC_DEV_FILEPATH)

        self.compilers = [ Compiler(e) for e in data["compilers"] ]
        self.targets   = [ Target  (e) for e in data["targets"]   ]
        self.mfc       = mfc

    def is_target_common(self, name: str) -> bool:
        return self.get_target(name).common_mode is not None

    def get_desired_target_mode_name(self, name: str) -> str:
        target = self.get_target(name)

        if self.is_target_common(name):
            return target.common_mode

        return self.mfc.args["mode"]

    def get_target_mode_folder_name(self, name: str, default: str) -> str:
        if self.is_target_common(name):
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
