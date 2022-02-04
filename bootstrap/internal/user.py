import internal.common      as common
import internal.configfiles as configfiles

from typing      import Any
from dataclasses import dataclass

@dataclass
class Defaults:
    threads:       str
    configuration: str

    def __init__(self, data) -> None:
        self.threads       = data.get("threads", 1)
        self.configuration = data.get("configuration")
    
@dataclass
class Compilers:
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data) -> None:
        self.c       = data["c"]
        self.cpp     = data["cpp"]
        self.fortran = data["fortran"]

@dataclass
class Configuration:
    name:    str
    c:       str
    cpp:     str
    fortran: str

    def __init__(self, data):
        self.name    = data["name"]
        self.c       = data.get("c",       "")
        self.cpp     = data.get("cpp",     "")
        self.fortran = data.get("fortran", "")

class MFCUser:
    def __init__(self) -> None:
        data = configfiles.ConfigFileBase(common.MFC_USER_FILEPATH, noexist_ok=False)

        self.defaults       = Defaults(data.tree_get("defaults"))
        self.compilers      = Compilers(data.tree_get("compilers"))
        self.configurations = [ Configuration(e) for e in data.tree_get("configurations") ]

    def get_configuration(self, name: str) -> Configuration:
        for configuration in self.configurations:
            if configuration.name == name:
                return configuration

        raise common.MFCException(f'MFCConf: Configuration "{name}" doesn\'t exist')

