import os
import typing
import dataclasses

from ..util.common import MFC_USER_FILEPATH, MFCException, create_file, \
                            file_load_yaml


@dataclasses.dataclass
class Mode:
    name:  str
    type:  str
    flags: typing.List[str]

    def __init__(self, data: dict) -> None:
        self.name  = data["name"]
        self.type  = data["type"]
        self.flags = data.get("flags", [])


@dataclasses.dataclass
class Build:
    threads: int

    def __init__(self, data: dict) -> None:
        self.threads = data.get("threads", 1)


@dataclasses.dataclass
class Run:
    nodes:          int
    partition:      str
    cpus_per_node:  int
    gpus_per_node:  int
    walltime:       str
    account:        str
    email:          str
    name:           str
    flags:          list

    def __init__(self, data: dict) -> None:
        self.nodes         = int(data.get("nodes",         ""))
        self.partition     =     data.get("partition",     "")
        self.cpus_per_node = int(data.get("cpus-per-node", ""))
        self.gpus_per_node = int(data.get("gpus-per-node", ""))
        self.walltime      =     data.get("walltime",      "")
        self.account       =     data.get("account",       "")
        self.email         =     data.get("email",         "")
        self.name          =     data.get("name",          "")
        self.flags         =     data.get("flags",         [])


@dataclasses.dataclass
class MFCUser:
    modes: list
    build: Build
    run:   Run

    def __init__(self) -> None:
        if not os.path.exists(MFC_USER_FILEPATH):
            create_file(MFC_USER_FILEPATH)
        
        data: dict = file_load_yaml(MFC_USER_FILEPATH)

        self.build = Build(data["build"])
        self.run   = Run  (data["run"])
        self.modes = [ Mode(e) for e in data["modes"] ]

    def get_mode(self, name: str) -> Mode:
        for mode in self.modes:
            if mode.name == name:
                return mode

        raise MFCException(f'MFCConf: Mode "{mode}" doesn\'t exist')
