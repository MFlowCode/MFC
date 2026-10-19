import dataclasses
import typing
from enum import Enum, unique


@unique
class gpuConfigOptions(Enum):
    NONE = "no"
    ACC = "acc"
    MP = "mp"


# Opt-in build features, each mapped to the CMake option MFC_<NAME>.
OPT_IN_FEATURES = (
    "nvshmem",
    "fused_rhs",
    "halo_single_round",
    "async_queue",
    "cuda_graphs",
    "direct_riemann_flux",
    "cbc_no_host_sync",
    "cbc_parallel_scatter",
    "accurate_grind_time",
    "halo_index_cache",
    "igr_shallow_halo",
    "igr_tuned_sweeps",
    "muscl_direct_loads",
    "lagrangian_device_cell_list",
    "hllc_register_tuning",
    "weno_direct_loads",
    "active_box",
    "persistent_bc",
    "fused_bubble_source",
    "qbmm_device_inversion",
)


@dataclasses.dataclass
class MFCConfig:
    mpi: bool = True
    gpu: str = gpuConfigOptions.NONE.value
    debug: bool = False
    reldebug: bool = False
    gcov: bool = False
    unified: bool = False
    single: bool = False
    mixed: bool = False
    fastmath: bool = False
    nvshmem: bool = False
    fused_rhs: bool = False
    halo_single_round: bool = False
    async_queue: bool = False
    cuda_graphs: bool = False
    direct_riemann_flux: bool = False
    cbc_no_host_sync: bool = False
    cbc_parallel_scatter: bool = False
    accurate_grind_time: bool = False
    halo_index_cache: bool = False
    igr_shallow_halo: bool = False
    igr_tuned_sweeps: bool = False
    muscl_direct_loads: bool = False
    lagrangian_device_cell_list: bool = False
    hllc_register_tuning: bool = False
    weno_direct_loads: bool = False
    active_box: bool = False
    persistent_bc: bool = False
    fused_bubble_source: bool = False
    qbmm_device_inversion: bool = False

    def __hash__(self):
        return hash(tuple(getattr(self, f.name) for f in dataclasses.fields(self)))

    @staticmethod
    def from_dict(d: dict):
        """Create a MFCConfig object from a dictionary with the same keys
        as the fields of MFCConfig"""
        r = MFCConfig()

        for field in dataclasses.fields(MFCConfig):
            setattr(r, field.name, d.get(field.name, getattr(r, field.name)))

        return r

    def items(self) -> typing.Iterable[typing.Tuple[str, typing.Any]]:
        return dataclasses.asdict(self).items()

    def make_options(self) -> typing.List[str]:
        """Returns a list of options that could be passed to mfc.sh again.
        Example: --no-debug --mpi --no-gpu --no-gcov --no-unified"""
        options = []
        for k, v in self.items():
            cli_k = k.replace("_", "-")
            if k == "gpu":
                options.append(f"--{v}-{cli_k}")
            else:
                options.append(f"--{'no-' if not v else ''}{cli_k}")
        return options

    def make_slug(self) -> str:
        """Sort the items by key, then join them with underscores. This uniquely
        identifies the configuration. Example: no-debug_no-gpu_no_mpi_no-gcov"""
        options = []
        for k, v in sorted(self.items(), key=lambda x: x[0]):
            # A disabled opt-in leaves the slug unchanged.
            if k in OPT_IN_FEATURES and not v:
                continue
            cli_k = k.replace("_", "-")
            if k == "gpu":
                options.append(f"--{v}-{cli_k}")
            else:
                options.append(f"--{'no-' if not v else ''}{cli_k}")
        return "_".join(options)

    def __eq__(self, other) -> bool:
        """Check if two MFCConfig objects are equal, field by field."""
        if not isinstance(other, MFCConfig):
            return NotImplemented
        for field in dataclasses.fields(self):
            if getattr(self, field.name) != getattr(other, field.name):
                return False

        return True

    def __str__(self) -> str:
        """Returns a string like "mpi=No & gpu=No & debug=No & gcov=No & unified=No" """
        strings = []
        for k, v in self.items():
            if isinstance(v, bool):
                strings.append(f"{k}={'Yes' if v else 'No'}")
            elif isinstance(v, str):
                strings.append(f"{k}={v.capitalize()}")
            elif isinstance(v, int):
                strings.append(f"{k}={v}")
            else:
                strings.append(f"{k}={v.__str__()}")

        return " & ".join(strings)


gCFG: MFCConfig = MFCConfig()
gARG: dict = {"rdma_mpi": False}


def ARG(arg: str, dflt=None) -> typing.Any:
    global gARG  # noqa: PLW0603
    if arg in gARG:
        return gARG[arg]
    if dflt is not None:
        return dflt

    raise KeyError(f"{arg} is not an argument.")


def ARGS() -> dict:
    global gARG  # noqa: PLW0603
    return gARG


def CFG() -> MFCConfig:
    global gCFG  # noqa: PLW0603
    return gCFG
