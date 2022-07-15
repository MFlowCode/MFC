import enum
import json
import typing
import dataclasses

from mfc.util.common import MFCException, enumeq

# === (RE)DEFINE FORTRAN CONSTANTS === #

DFLT_REAL, DFLT_INT = -1e6, -100

# === DEFINE ALL ENUMERATIONS === #


class PatchGeometry(enum.Enum):
    D1_LINE_SEGMENT =  1; D2_CIRCLE      =  2; D2_RECTANGLE =  3
    D2_SWEEP_LINE   =  4; D2_ELLIPSE     =  5; D2_VORTEX    =  6
    D2_ANALYTICAL   =  7; D3_SPHERE      =  8; D3_CUBOID    =  9
    D3_CYLINDER     = 10; D3_SWEEP_PLANE = 11; D3_ELLIPSOID = 12
    D3_ANALYTICAL   = 13


class FluxLimiter(enum.Enum):
    MINMOD = 1; MC         = 2; OSPRE    = 3; SUPERBEE = 4
    SWEBY  = 5; VAN_ALBADA = 6; VAN_LEER = 7


class MulticomponentModel(enum.Enum):
    GAMMA_PI_INF = 1; EQUATION_5 = 2; EQUATION_6 = 3


class BubbleModel(enum.Enum):
    GILMORE = 1; KELLER_MIKSIS = 2; RAYLEIGH_PLESSET = 3


class ThermalModel(enum.Enum):
    ADIABATIC = 1; ISOTHERMAL = 2; TRANSFER = 3


class BoundaryCondition(enum.Enum):
    PERIODIC                           =  -1
    REFLECTIVE                         =  -2
    GHOST_CELL_EXTRAPOLATION           =  -3
    RIEMANN_EXTRAPOLATION              =  -4
    SLIP_WALL                          =  -5
    NON_REFLECTING_SUBSONIC_BUFFER     =  -6
    NON_REFLECTING_SUBSONIC_INFLOW     =  -7
    NON_REFLECTING_SUBSONIC_OUTFLOW    =  -8
    FORCE_FREE_SUBSONIC_OUTFLOW        =  -9
    CONSTANT_PRESSURE_SUBSONIC_OUTFLOW = -10
    SUPERSONIC_INFLOW                  = -11
    SUPERSONIC_OUTFLOW                 = -12


class DatabaseFormat(enum.Enum):
    SILO_HDF5 = 1; BINARY = 2


class FloatingPrecision(enum.Enum):
    SINGLE = 1; DOUBLE = 2


class RiemannSolver(enum.Enum):
    HLL = 1; HLLC = 2; EXACT = 3


class WaveSpeedEstimation(enum.Enum):
    DIRECT = 1; PRESSURE_VELOCITY = 2


class AverageStateEvaluation(enum.Enum):
    ROE_MEAN = 1; ARITHMETIC_MEAN = 2


class TimeStepper(enum.Enum):
    RUNGE_KUTTA_1 = 1; RUNGE_KUTTA_2 = 2; RUNGE_KUTTA_3 = 3; RUNGE_KUTTA_4 = 4
    RUNGE_KUTTA_5 = 5


class WenoVariables(enum.Enum):
    CONSERVATIVE = 1; PRIMITIVE = 2


class AcousticWaveForm(enum.Enum):
    SINE = 1; GAUSSIAN = 2; SQUARE = 3


class AcousticSpacialSupport(enum.Enum):
    D1 = 1; D2_FINITE_WIDTH = 2; D3_FINITE_LINE_PATCH = 3


class BubbleDistribution(enum.Enum):
    BINORMAL = 1; LOGNORMAL_NORMAL = 2


# === DEFINE UTILITY FUNCTIONS === #


def dflt(factory: typing.Callable):
    return dataclasses.field(default_factory=factory)


def serialize(e: typing.List) -> str:
    if callable(getattr(e, "serialize", None)):
        return e.serialize()
    
    if isinstance(e, enum.Enum):
        return f"{e.value}"

    # We add "d0" to each "double", which we store in Python as "float",
    # because F90 considers any floating point value without "d0" to be single-
    # precision. Thus, we must specify "d0" ourselves to avoid eronious implicit
    # conversions, from "double" to "float", and to "double" again, thereby losing
    # precision.
    if isinstance(e, float):
        return f"{e:.16f}d0"

    if isinstance(e, bool):
        return f".{e}."

    return f"{e}"


def serialize_array_constructor(f90_type: str, arr: list) -> str:
    if len(arr) == 0:
        return f"[{f90_type} ::]"

    return f"(/ {f', '.join([ serialize(e) for e in arr ])} /)"


# === DEFINE UTILITY CLASSES === #


@dataclasses.dataclass
class Point:
    x: float = dflt(lambda: 0)
    y: float = dflt(lambda: 0)
    z: float = dflt(lambda: 0)


class Vector(Point):
    pass


# ===    DEFINE CLASSES    === #
# === COMPUTATIONAL DOMAIN === #


@dataclasses.dataclass
class AxisDomain:
    begin: float = dflt(lambda: 0)
    end:   float = dflt(lambda: 0)

    def serialize(self) -> str:
        return f"bounds_info({float(self.begin)}, {float(self.end)})"


@dataclasses.dataclass
class AxisStretch(AxisDomain):
    stretch:   bool   = dflt(lambda: False)
    rate:      float  = dflt(lambda: 0.0)
    loops:     int    = dflt(lambda: 1)
    begin_pos: float  = dflt(lambda: DFLT_REAL)
    begin_neg: float  = dflt(lambda: DFLT_REAL)


@dataclasses.dataclass
class Stretch:
    x: AxisStretch = dflt(lambda: AxisStretch())
    y: AxisStretch = dflt(lambda: AxisStretch())
    z: AxisStretch = dflt(lambda: AxisStretch())


@dataclasses.dataclass
class Cells:
    x: int
    y: int = dflt(lambda: 0)
    z: int = dflt(lambda: 0)


@dataclasses.dataclass
class SpacialDomain:
    x: AxisDomain
    y: AxisDomain = dflt(lambda: AxisDomain(begin=DFLT_REAL, end=DFLT_REAL))
    z: AxisDomain = dflt(lambda: AxisDomain(begin=DFLT_REAL, end=DFLT_REAL))


@dataclasses.dataclass
class Time:
    end:   int
    save:  int
    dt:    float
    begin: int   = dflt(lambda: 0)


@dataclasses.dataclass
class ComputationalDomain:
    cells:        Cells
    domain:       SpacialDomain
    time:         Time
    cyl_coord:    bool    = dflt(lambda: False)
    stretch:      Stretch = dflt(lambda: Stretch())


# ===        FLUIDS        === #


@dataclasses.dataclass
class Fluid:
    gamma:   float              = dflt(lambda: DFLT_REAL)
    pi_inf:  float              = dflt(lambda: DFLT_REAL)
    Re:      typing.List[float] = dflt(lambda: [DFLT_REAL, DFLT_REAL])
    mul0:    float              = dflt(lambda: DFLT_REAL)
    ss:      float              = dflt(lambda: DFLT_REAL)
    pv:      float              = dflt(lambda: DFLT_REAL)
    gamma_v: float              = dflt(lambda: DFLT_REAL)
    M_v:     float              = dflt(lambda: DFLT_REAL)
    mu_v:    float              = dflt(lambda: DFLT_REAL)
    k_v:     float              = dflt(lambda: DFLT_REAL)
    G:       float              = dflt(lambda: DFLT_REAL)

    def serialize(self) -> str:
        members = []
        for field in dataclasses.fields(self):
            val = getattr(self, field.name)

            if isinstance(val, list):
                py_type = typing.get_args(field.type)[0]

                if py_type == float:
                    f90_type = "real(kind(0d0))"

                members.append(serialize_array_constructor(f90_type, list(map(str, val))))
            else:
                members.append(serialize(val))
            
        return f"physical_parameters({', '.join(members)})"


@dataclasses.dataclass
class Fluids:
    count:  int                = dflt(lambda: 0)
    fluids: typing.List[Fluid] = dflt(lambda: [])

    def __getitem__(self, id) -> Fluid:
        return self.fluids[id]


# === SIMULATION ALGORITHM === #


@dataclasses.dataclass
class WenoParameters:
    order:     int
    variables: WenoVariables
    mapped:    bool  = dflt(lambda: False)
    mp:        bool  = dflt(lambda: False)
    flat:      bool  = dflt(lambda: True)
    epsilon:   float = dflt(lambda: 1e-16)
    average:   bool  = dflt(lambda: False)
    Re_flux:   bool  = dflt(lambda: False)


@dataclasses.dataclass
class AxisBoundaryCondition:
    begin: BoundaryCondition = DFLT_INT
    end:   BoundaryCondition = DFLT_INT
    # Note: Using "DFLT_REAL" instead of "DFLT_INT" doesn't make sense, but it is
    #       what is used by the fortran code.


@dataclasses.dataclass
class BoundaryConditions:
    x: AxisBoundaryCondition
    y: AxisBoundaryCondition = dflt(lambda: AxisBoundaryCondition())
    z: AxisBoundaryCondition = dflt(lambda: AxisBoundaryCondition())


@dataclasses.dataclass
class SimulationAlgorithm:
    model:              MulticomponentModel
    weno:               WenoParameters
    boundary:           BoundaryConditions
    time_stepper:       TimeStepper
    wave_speeds:        WaveSpeedEstimation
    avg_state:          AverageStateEvaluation
    riemann_solver:     RiemannSolver
    riemann_flat:       bool        = dflt(lambda: True)
    char_decomp:        bool        = dflt(lambda: False)
    commute_err:        bool        = dflt(lambda: False)
    split_err:          bool        = dflt(lambda: False)
    mixture_err:        bool        = dflt(lambda: False)
    mpp_lim:            bool        = dflt(lambda: False)
    adv_alphan:         bool        = dflt(lambda: False)
    hypoelasticity:     bool        = dflt(lambda: False)
    null_weights:       bool        = dflt(lambda: False)
    alt_crv:            bool        = dflt(lambda: False)
    alt_soundspeed:     bool        = dflt(lambda: False)
    regularization:     bool        = dflt(lambda: False)
    reg_eps:            float       = dflt(lambda: DFLT_REAL)
    tvd_riemann_flux:   bool        = dflt(lambda: False)
    tvd_rhs_flux:       bool        = dflt(lambda: False)
    tvd_wave_speeds:    bool        = dflt(lambda: False)
    flux_lim:           FluxLimiter = dflt(lambda: DFLT_INT)
    We_riemann_flux:    bool        = dflt(lambda: False)
    We_rhs_flux:        bool        = dflt(lambda: False)
    We_src:             bool        = dflt(lambda: False)
    We_wave_speeds:     bool        = dflt(lambda: False)
    lsq_deriv:          bool        = dflt(lambda: False)
    hypoelasticity:     bool        = dflt(lambda: False)
    perturb_flow:       bool        = dflt(lambda: False)
    perturb_flow_fluid: int         = dflt(lambda: DFLT_INT)
    perturb_sph:        bool        = dflt(lambda: False)
    perturb_sph_fluid:  int         = dflt(lambda: DFLT_INT)
    fluid_rho:          float       = dflt(lambda: DFLT_REAL)
    rhoref:             float       = dflt(lambda: DFLT_REAL) 
    pref:               float       = dflt(lambda: DFLT_REAL)



# ===   PATCH PARAMETERS   === #



@dataclasses.dataclass
class Patch:
    geometry:        PatchGeometry      = dflt(lambda: DFLT_INT)
    centroid:        Point              = dflt(lambda: Vector(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    length:          Vector             = dflt(lambda: Vector(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    radius:          float              = dflt(lambda: DFLT_REAL)
    radii:           Vector             = dflt(lambda: Vector(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    epsilon:         float              = dflt(lambda: DFLT_REAL)
    beta:            float              = dflt(lambda: DFLT_REAL)
    normal:          Vector             = dflt(lambda: Vector(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    alter:           typing.List[bool]  = dflt(lambda: [])
    smoothen:        bool               = dflt(lambda: False)
    smooth_patch_id: int                = dflt(lambda: DFLT_INT)
    smooth_coeff:    float              = dflt(lambda: DFLT_REAL)
    alpha_rho:       typing.List[float] = dflt(lambda: [])
    rho:             float              = dflt(lambda: DFLT_REAL)
    velocity:        Vector             = dflt(lambda: Vector(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    pressure:        float              = dflt(lambda: DFLT_REAL)
    alpha:           typing.List[float] = dflt(lambda: [])
    gamma:           float              = dflt(lambda: DFLT_REAL)
    pi_inf:          float              = dflt(lambda: DFLT_REAL)
    r0:              float              = dflt(lambda: DFLT_REAL)
    v0:              float              = dflt(lambda: DFLT_REAL)
    p0:              float              = dflt(lambda: DFLT_REAL)
    m0:              float              = dflt(lambda: DFLT_REAL)


    def serialize(self) -> str:
        members = []
        for field in dataclasses.fields(self):
            val = getattr(self, field.name)

            if field.name in ["radii", "normal", "velocity"]:
                members.append(f"(/ {float(val.x)}, {float(val.y)}, {float(val.z)} /)")
            elif field.name == "length":
                members.append(f"{serialize(val.x)}, {serialize(val.y)}, {serialize(val.z)}")
            elif isinstance(val, Point) or isinstance(val, Vector):
                members.append(f"{float(val.x)}, {float(val.y)}, {float(val.z)}")
            elif isinstance(val, list):
                py_type = typing.get_args(field.type)[0]

                if py_type == bool:
                    f90_type = "logical"
                elif py_type == float:
                    f90_type = "real(kind(0d0))"

                members.append(serialize_array_constructor(f90_type, val))
            else:
                members.append(serialize(val))
            
        return f"ic_patch_parameters({', '.join(members)})"


# ===  DATABASE STRUCTURE  === #


@dataclasses.dataclass
class AxisMinMax:
    min: float = dflt(lambda: DFLT_REAL)
    max: float = dflt(lambda: DFLT_REAL)


@dataclasses.dataclass
class Integral:
    x: AxisMinMax = dflt(lambda: AxisMinMax())
    y: AxisMinMax = dflt(lambda: AxisMinMax())
    z: AxisMinMax = dflt(lambda: AxisMinMax())

    def serialize(self) -> str:
        return f"integral_parameters({float(self.x.min)}, {float(self.x.max)}, {float(self.y.min)}, {float(self.y.max)}, {float(self.z.min)}, {float(self.z.max)})"


@dataclasses.dataclass
class Probe(Vector):
    def serialize(self) -> str:
        return f"probe_parameters({float(self.x)}, {float(self.y)}, {float(self.z)})"


@dataclasses.dataclass
class DatabaseWrite:
    alpha_rho:    bool = dflt(lambda: False)
    rho:          bool = dflt(lambda: False)
    mom:          bool = dflt(lambda: False)
    velocity:     bool = dflt(lambda: False)
    flux:         bool = dflt(lambda: False)
    E:            bool = dflt(lambda: False)
    pressure:     bool = dflt(lambda: False)
    alpha:        bool = dflt(lambda: False)
    gamma:        bool = dflt(lambda: False)
    heat_ratio:   bool = dflt(lambda: False)
    pi_inf:       bool = dflt(lambda: False)
    pressure_inf: bool = dflt(lambda: False)
    prim_vars:    bool = dflt(lambda: False)
    cons_vars:    bool = dflt(lambda: False)
    c:            bool = dflt(lambda: False)
    omega:        bool = dflt(lambda: False)
    schlieren:    bool = dflt(lambda: False)
    probe:        bool = dflt(lambda: False)
    integral:     bool = dflt(lambda: False)
    #TODO: Some (i) params


@dataclasses.dataclass
class DatabseStructure:
    precision:      FloatingPrecision = dflt(lambda: FloatingPrecision.DOUBLE)
    fd_order:       int               = dflt(lambda: DFLT_INT)
    write:          DatabaseWrite     = dflt(lambda: DatabaseWrite())
    alt_soundspeed: bool              = dflt(lambda: False)
    parallel_io:    bool              = dflt(lambda: False)
    coarsen_silo:   bool              = dflt(lambda: False)
    format:         DatabaseFormat    = dflt(lambda: DatabaseFormat.SILO_HDF5)
    probes:         typing.List[Probe]    = dflt(lambda: [])
    integrals:      typing.List[Integral] = dflt(lambda: [])


# ===       BUBBLES       === #


@dataclasses.dataclass
class Bubbles:
    bubbles:      bool               = dflt(lambda: False)
    model:        BubbleModel        = dflt(lambda: BubbleModel.GILMORE)
    polytropic:   bool               = dflt(lambda: True)
    thermal:      ThermalModel       = dflt(lambda: DFLT_INT)
    R0ref:        float              = dflt(lambda: DFLT_REAL)
    number:       int                = dflt(lambda: DFLT_INT)
    cavitation:   float              = dflt(lambda: DFLT_REAL)
    weber:        float              = dflt(lambda: DFLT_REAL)
    Re_inv:       float              = dflt(lambda: DFLT_REAL)
    mu_10:        float              = dflt(lambda: DFLT_REAL)
    ss:           float              = dflt(lambda: DFLT_REAL)
    pv:           float              = dflt(lambda: DFLT_REAL)
    gamma_v:      float              = dflt(lambda: DFLT_REAL)
    M_v:          float              = dflt(lambda: DFLT_REAL)
    mu_v:         float              = dflt(lambda: DFLT_REAL)
    k_v:          float              = dflt(lambda: DFLT_REAL)
    qbmm:         bool               = dflt(lambda: False)
    polydisperse: bool               = dflt(lambda: False)
    nnode:        int                = dflt(lambda: 1)
    sigR:         float              = dflt(lambda: DFLT_REAL)
    sigV:         float              = dflt(lambda: DFLT_REAL)
    rhoRV:        float              = dflt(lambda: 0)
    poly_sigma:   float              = dflt(lambda: DFLT_REAL)
    distribution: BubbleDistribution = dflt(lambda: DFLT_INT)
    R0_type:      int                = dflt(lambda: DFLT_INT)

# === LOGISTICS === #


@dataclasses.dataclass
class Logistics:
    case_dir:      str  = dflt(lambda: '.')
    run_time_info: bool = dflt(lambda: False)
    cu_mpi:        bool = dflt(lambda: False)
    cu_tensor:     bool = dflt(lambda: False)
    debug:         bool = dflt(lambda: False)
    old_grid:      bool = dflt(lambda: False)
    old_ic:        bool = dflt(lambda: False)
    t_step_old:    int  = dflt(lambda: DFLT_INT)


# === === #


@dataclasses.dataclass
class Monopole:    
    location:  Point                  = dflt(lambda: Point(x=DFLT_REAL,y=DFLT_REAL,z=DFLT_REAL))
    magnitude: float                  = dflt(lambda: DFLT_REAL)
    length:    float                  = dflt(lambda: DFLT_REAL)
    npulse:    float                  = dflt(lambda: 1.0)
    direction: float                  = dflt(lambda: 1.0)
    delay:     float                  = dflt(lambda: DFLT_REAL)
    pulse:     AcousticWaveForm       = dflt(lambda: AcousticWaveForm.SINE)
    support:   AcousticSpacialSupport = dflt(lambda: AcousticSpacialSupport.D1)


    def serialize(self) -> str:
        members = []
        for field in dataclasses.fields(self):
            val = getattr(self, field.name)

            if isinstance(val, Point) or isinstance(val, Vector):
                members.append(f"(/ {float(val.x)}, {float(val.y)}, {float(val.z)} /)")
            else:
                members.append(serialize(val))
            
        return f"mono_parameters({', '.join(members)})"


@dataclasses.dataclass
class AcousticParameters:
    monopole:  bool = dflt(lambda: False)
    monopoles: typing.List[Monopole] = dflt(lambda: [])


# === DEFINE A CASE === #


@dataclasses.dataclass
class Case:
    logistics:  Logistics
    domain:     ComputationalDomain
    fluids:     Fluids
    patches:    typing.List[Patch]
    bubbles:    Bubbles
    database:   DatabseStructure
    algorithm:  SimulationAlgorithm
    acoustic:   AcousticParameters


    def __repr__(self) -> str:
        return json.dumps(self.to_json(), indent=4)
    
    def set(self, path: str, value: typing.Any):
        expr = f"self.{path}={value}"

        # We may need to "allocate" another fluid
        match = "self.fluids["
        if expr.startswith(match):
            fluid_idx = int(expr[len(match):].split(']')[0])

            if fluid_idx + 1 > len(self.fluids.fluids):
                self.fluids.fluids.append(Fluid())

        exec(expr)

    def to_json(self) -> dict:
        def handle_enums(data):
            def convert_value(obj):
                if isinstance(obj, enum.Enum):
                    return obj.value

                return obj

            return dict((key, convert_value(value)) for key, value in data)

        # We autogenerate here some parameters to make the Fortran file cleaner.
        # Generating these parameters using Fypp would lead to a lot fo code
        # duplication as well as very long and unreadable python one-liners.
        # It also has the benefit of consolidating all input parameter logic
        # in one place.

        for id, patch in enumerate(self.patches):
            patch.smooth_patch_id = id + 1

        # Compute grid_geometry
        grid_geometry = None
        if self.domain.cyl_coord is False:
            grid_geometry = 1
        elif self.domain.cyl_coord and self.domain.cells.z == 0:
            grid_geometry = 2
        else:
            grid_geometry = 3

        nterms = DFLT_INT
        # Not set for BubbleModel.GILMORE
        if enumeq(self.bubbles.model, BubbleModel.KELLER_MIKSIS):
            nterms = 12
        elif enumeq(self.bubbles.model, BubbleModel.RAYLEIGH_PLESSET):
            nterms = 6

        # num_fluids
        if enumeq(self.algorithm.model, 1) and len(self.fluids) != 1:
            raise MFCException("Invalid combination.")

        num_fluids       = self.fluids.count
        num_fluids_alloc = max(self.fluids.count, len(self.fluids.fluids)) + 1 # +1 is a workaround for the bellow:

        # Workaround for https://github.com/henryleberre/MFC/blob/f6620bba0f248b729491dd3b0b1dd070ddd63c2d/src/common/m_global_parameters.fpp#L511-L515
        # if num_fluids == 1 and self.bubbles.bubbles and self.algorithm.model == MulticomponentModel.EQUATION_5:
        #     num_fluids_alloc += 1

        # patches & alter_patches
        patches = self.patches.copy()
        for patch in patches:
            patch.alter     = [True] + patch.alter     + [False]    *(len(self.patches)-len(patch.alter))
            patch.alpha     =          patch.alpha     + [DFLT_REAL]*(num_fluids-len(patch.alpha))
            patch.alpha_rho =          patch.alpha_rho + [DFLT_REAL]*(num_fluids-len(patch.alpha_rho))

        return {**dataclasses.asdict(self, dict_factory=handle_enums), **{
            "autogen": {
                "num_fluids":       num_fluids,
                "num_fluids_alloc": num_fluids_alloc,
                "grid_geometry":    grid_geometry,
                "nterms":           nterms,
                "x_domain":         serialize(self.domain.domain.x),
                "y_domain":         serialize(self.domain.domain.y),
                "z_domain":         serialize(self.domain.domain.z),
                "patch_icpp":       serialize_array_constructor("ic_patch_parameters", patches),
                "fluid_pp":         serialize_array_constructor("physical_parameters", self.fluids.fluids + [Fluid()]),
                "mono":             serialize_array_constructor("mono_parameters",     self.acoustic.monopoles),
                "probe":            serialize_array_constructor("probe_parameters",    self.database.probes),
                "integral":         serialize_array_constructor("integral_parameters", self.database.integrals),
            }
        }}
