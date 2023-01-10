import typing, dataclasses


@dataclasses.dataclass
class MFCConfig:
    mpi:   bool = True
    gpu:   bool = False
    debug: bool = False

    def from_dict(d: dict):
        r = MFCConfig()

        for key in d:
            setattr(r, key, d[key])

        return r

    def __eq__(self, other) -> bool:
        for field in dataclasses.fields(self):
            if getattr(self, field.name) != getattr(other, field.name):
                return False

        return True

    def __str__(self) -> str:
        m = { False: "No", True: "Yes" }
        r = ' & '.join([ f"{field.name}={m[getattr(self, field.name)]}" for field in dataclasses.fields(self) ])

        return r


gCFG: MFCConfig = MFCConfig()
gARG: dict      = {}


def ARG(arg: str) -> typing.Any:
    global gARG
    return gARG[arg]

def ARGS() -> dict:
    global gARG
    return gARG

def CFG() -> MFCConfig:
    global gCFG
    return gCFG
