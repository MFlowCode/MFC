import typing, dataclasses


@dataclasses.dataclass
class MFCConfig:
    mpi:       bool = True
    gpu:       bool = False
    debug:     bool = False
    gcov:      bool = False
    unified:   bool = False
    single:    bool = False
    fastmath : bool = False

    @staticmethod
    def from_dict(d: dict):
        """ Create a MFCConfig object from a dictionary with the same keys
            as the fields of MFCConfig """
        r = MFCConfig()

        for field in dataclasses.fields(MFCConfig):
            setattr(r, field.name, d[field.name])

        return r

    def items(self) -> typing.List[typing.Tuple[str, bool]]:
        return dataclasses.asdict(self).items()

    def make_options(self) -> typing.List[str]:
        """ Returns a list of options that could be passed to mfc.sh again.
            Example: --no-debug --mpi --no-gpu --no-gcov --no-unified"""
        return [ f"--{'no-' if not v else ''}{k}" for k, v in self.items() ]

    def make_slug(self) -> str:
        """ Sort the items by key, then join them with underscores. This uniquely 
            identifies the configuration. Example: no-debug_no-gpu_no_mpi_no-gcov """
        return '_'.join([ f"{'no-' if not v else ''}{k}" for k, v in sorted(self.items(), key=lambda x: x[0]) ])

    def __eq__(self, other) -> bool:
        """ Check if two MFCConfig objects are equal, field by field. """
        for field in dataclasses.fields(self):
            if getattr(self, field.name) != getattr(other, field.name):
                return False

        return True

    def __str__(self) -> str:
        """ Returns a string like "mpi=No & gpu=No & debug=No & gcov=No & unified=No" """

        return ' & '.join([ f"{k}={'Yes' if v else 'No'}" for k, v in self.items() ])


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
