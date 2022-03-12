import internal.common      as common
import internal.configfiles as configfiles

from typing      import Any
from dataclasses import dataclass

@dataclass
class General:
    configuration: str

    def __init__(self, data: dict):
        self.configuration = data.get("configuration")

@dataclass
class Compilers:
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data: dict) -> None:
        self.c       = data["c"]
        self.cpp     = data["cpp"]
        self.fortran = data["fortran"]

@dataclass
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

@dataclass
class Build:
    def __init__(self, data: dict) -> None:
        self.threads   = data.get("threads")
        self.compilers = Compilers(data.get("compilers"))

@dataclass
class Run:
    nodes:          int
    partition:      str
    tasks_per_node: int
    gpus_per_node:  int
    walltime:       str

    def __init__(self, data: dict) -> None:
        self.nodes          = int(data.get("nodes"))
        self.partition      = data.get("partition")
        self.tasks_per_node = int(data.get("tasks-per-node"))
        self.gpus_per_node  = int(data.get("gpus-per-node"))
        self.walltime       = data.get("walltime")

class MFCUser:
    def __init__(self) -> None:
        data = configfiles.ConfigFileBase(common.MFC_USER_FILEPATH, noexist_ok=False)

        self.general        = General(data.tree_get("general"))
        self.build          = Build  (data.tree_get("build"))
        self.run            = Run    (data.tree_get("run"))
        self.configurations = [ Configuration(e) for e in data.tree_get("configurations") ]

    def get_configuration(self, name: str) -> Configuration:
        for configuration in self.configurations:
            if configuration.name == name:
                return configuration

        raise common.MFCException(f'MFCConf: Configuration "{name}" doesn\'t exist')

