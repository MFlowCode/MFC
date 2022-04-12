import os, dataclasses

import common

@dataclasses.dataclass
class General:
    configuration: str

    def __init__(self, data: dict):
        self.configuration = data.get("configuration")


@dataclasses.dataclass
class Compilers:
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data: dict) -> None:
        self.c       = data["c"]
        self.cpp     = data["cpp"]
        self.fortran = data["fortran"]


@dataclasses.dataclass
class Configuration:
    name:    str
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data: dict) -> None:
        self.name    = data["name"]
        self.c       = data.get("c",       "")
        self.cpp     = data.get("cpp",     "")
        self.fortran = data.get("fortran", "")


@dataclasses.dataclass
class Build:
    def __init__(self, data: dict) -> None:
        self.threads   = data.get("threads")
        self.compilers = Compilers(data.get("compilers"))


@dataclasses.dataclass
class Run:
    nodes:          int
    partition:      str
    tasks_per_node: int
    gpus_per_node:  int
    walltime:       str
    account:        str

    def __init__(self, data: dict) -> None:
        self.nodes          = int(data.get("nodes"))
        self.partition      = data.get("partition")
        self.tasks_per_node = int(data.get("tasks-per-node"))
        self.gpus_per_node  = int(data.get("gpus-per-node"))
        self.walltime       = data.get("walltime")
        self.account        = data.get("account")


class MFCUser:
    def __init__(self) -> None:
        if not os.path.exists(common.MFC_USER_FILEPATH):
            common.create_file(common.MFC_USER_FILEPATH)
        
        data: dict = common.file_load_yaml(common.MFC_USER_FILEPATH)

        self.general        = General(data["general"])
        self.build          = Build  (data["build"])
        self.run            = Run    (data["run"])
        self.configurations = [ Configuration(e) for e in data["configurations"] ]

    def get_configuration(self, name: str) -> Configuration:
        for configuration in self.configurations:
            if configuration.name == name:
                return configuration

        raise common.MFCException(f'MFCConf: Configuration "{name}" doesn\'t exist')
