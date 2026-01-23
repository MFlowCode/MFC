import typing, os.path

from ..printer import cons
from ..state import ARG, ARGS
from . import pack as _pack
from . import errors
from . import tol as packtol
from ..common import MFCException

def load(packpath: str) -> _pack.Pack:
    return _pack.load(packpath)

def pack(casepath: str, packpath: str = None) -> typing.Tuple[_pack.Pack, str]:
    if packpath is None:
        packpath = casepath

    p, err = _pack.compile(casepath)

    if err is not None:
        return None, err

    p.save(packpath)

    return p, None

def compare(lhs: str = None, rhs: str = None, tol: packtol.Tolerance = None) -> typing.Tuple[errors.Error, str]:
    if isinstance(lhs, str):
        lhs = load(lhs)
    if isinstance(rhs, str):
        rhs = load(rhs)

    return packtol.compare(lhs, rhs, tol)

def packer():
    if ARG("packer") == "pack":
        if not os.path.isdir(ARG("input")):
            ARGS()["input"] = os.path.dirname(ARG("input"))

        out_dir = os.path.sep.join([ARG("input"), ARG("output")]) if ARG("output") is not None else None
        _, err = pack(ARG("input"), out_dir)
        if err is not None:
            raise MFCException(err)
    elif ARG("packer") == "compare":
        cons.print(
            f"Comparing [magenta]{os.path.relpath(ARG('input1'))}[/magenta] to "
            f"[magenta]{os.path.relpath(ARG('input2'))}[/magenta]:"
        )

        cons.indent()
        cons.print()

        tol = packtol.Tolerance(ARG("abstol"), ARG("reltol"))

        err, msg = compare(ARG("input1"), ARG("input2"), tol)
        if msg is not None:
            cons.print(f"[bold red]ERROR[/bold red]: The two packs are not within tolerance ({tol}).")
            cons.print(msg)
        else:
            cons.print(f"[bold green]OK[/bold green]: The two packs are within tolerance ({tol}).")
            cons.print(f"Average error: {err}.")

        cons.print()
        cons.unindent()
    else:
        raise MFCException(f"Unknown packer command: {ARG('packer')}")
