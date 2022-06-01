import enum
import json
import typing
import dataclasses


class Coordinates(enum.Enum):
    CARTESIAN   = 0
    CYLINDRICAL = 1
    # TODO: ...


class PatchGeometry(enum.Enum):
    LINE_SEGMENT  = 1
    D1_ANALYTICAL = 2
    # TODO: ...


@dataclasses.dataclass
class Fluid:
    gamma:   float      = None
    pi_inf:  float      = None
    mul0:    typing.Any = None
    ss:      typing.Any = None
    gamma_v: typing.Any = None
    M_v:     typing.Any = None
    mu_v:    typing.Any = None
    k_v:     typing.Any = None
    G:       typing.Any = None


@dataclasses.dataclass
class Grid:
    m: int
    n: int = dataclasses.field(default=0)
    p: int = dataclasses.field(default=0)
    coordinates: Coordinates =  dataclasses.field(default=Coordinates.CARTESIAN)


@dataclasses.dataclass
class Domain:
    begin: float
    end:   float


@dataclasses.dataclass
class Timing:
    start: int
    stop:  int
    save:  int
    dt:    float


@dataclasses.dataclass
class Point:
    x: float = dataclasses.field(default=0)
    y: float = dataclasses.field(default=0)
    z: float = dataclasses.field(default=0)


class Vector(Point):
    pass


@dataclasses.dataclass
class Patch:
    geometry:        PatchGeometry = None
    radius:          typing.Any    = None
    radii:           typing.Any    = None
    epsilon:         typing.Any    = None
    beta:            typing.Any    = None
    normal:          typing.Any    = None
    smoothen:        bool          = None
    smooth_patch_id: typing.Any    = None
    alpha_rho:       typing.Any    = None
    smooth_coeff:    typing.Any    = None
    rho:             typing.Any    = None
    vel:             typing.Any    = None
    pres:            typing.Any    = None
    alpha:           typing.Any    = None
    gamma:           typing.Any    = None
    pi_inf:          typing.Any    = None
    r0:              typing.Any    = None
    v0:              typing.Any    = None
    p0:              typing.Any    = None
    m0:              typing.Any    = None
    centroid:        Point         = None
    length:          Vector        = None
    radii:           Vector        = None
    normal:          Vector        = None
    velocity:        Vector        = None

    # TODO: Add for arho_id in range(1, 10+1):

    alter_patch:     typing.Any = None


@dataclasses.dataclass
class Case:
    grid:   Grid
    timing: Timing
    fluids: typing.List[Fluid]

    def to_json(self):
        def handle_enums(data):
            def convert_value(obj):
                if isinstance(obj, enum.Enum):
                    return obj.value
                return obj

            return dict((key, convert_value(value)) for key, value in data)

        return dataclasses.asdict(self, dict_factory=handle_enums)
