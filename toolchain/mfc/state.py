import typing, dataclasses
from enum import Enum, unique

@unique
class gpuConfigOptions(Enum):
    NONE = 'no'
    ACC = 'acc'
    MP = 'mp'

@dataclasses.dataclass
class MFCConfig:
    # pylint: disable=too-many-instance-attributes
    mpi:       bool = True
    gpu:     str = gpuConfigOptions.NONE.value
    debug:     bool = False
    gcov:      bool = False
    unified:   bool = False
    single:    bool = False
    mixed:   bool = False
    fastmath: bool = False

    @staticmethod
    def from_dict(d: dict):
        """ Create a MFCConfig object from a dictionary with the same keys
            as the fields of MFCConfig """
        r = MFCConfig()

        for field in dataclasses.fields(MFCConfig):
            setattr(r, field.name, d[field.name])

        return r

    def items(self) -> typing.Iterable[typing.Tuple[str, typing.Any]]:
        return dataclasses.asdict(self).items()

    def make_options(self) -> typing.List[str]:
        """ Returns a list of options that could be passed to mfc.sh again.
            Example: --no-debug --mpi --no-gpu --no-gcov --no-unified"""
        options = []
        for k, v in self.items():
            if k == 'gpu':
                options.append(f"--{v}-{k}")
            else:
                options.append(f"--{'no-' if not v else ''}{k}")
        return options

    def make_slug(self) -> str:
        """ Sort the items by key, then join them with underscores. This uniquely 
            identifies the configuration. Example: no-debug_no-gpu_no_mpi_no-gcov """
        options = []
        for k, v in sorted(self.items(), key=lambda x: x[0]):
            if k == 'gpu':
                options.append(f"--{v}-{k}")
            else:
                options.append(f"--{'no-' if not v else ''}{k}")
        return '_'.join(options)

    def __eq__(self, other) -> bool:
        """ Check if two MFCConfig objects are equal, field by field. """
        for field in dataclasses.fields(self):
            if getattr(self, field.name) != getattr(other, field.name):
                return False

        return True

    def __str__(self) -> str:
        """ Returns a string like "mpi=No & gpu=No & debug=No & gcov=No & unified=No" """
        strings = []
        for k,v in self.items():
            if isinstance(v, bool):
                strings.append(f"{k}={'Yes' if v else 'No'}")
            elif isinstance(v, str):
                strings.append(f"{k}={v.capitalize()}")
            elif isinstance(v, int):
                strings.append(f"{k}={v}")
            else:
                strings.append(f"{k}={v.__str__()}")

        return ' & '.join(strings)


gCFG: MFCConfig = MFCConfig()
gARG: dict      = {"rdma_mpi": False}

def ARG(arg: str, dflt = None) -> typing.Any:
    # pylint: disable=global-variable-not-assigned
    global gARG
    if arg in gARG:
        return gARG[arg]
    if dflt is not None:
        return dflt

    raise KeyError(f"{arg} is not an argument.")

def ARGS() -> dict:
    # pylint: disable=global-variable-not-assigned
    global gARG
    return gARG

def CFG() -> MFCConfig:
    # pylint: disable=global-variable-not-assigned
    global gCFG
    return gCFG
