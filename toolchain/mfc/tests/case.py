import os
import copy
import typing
import hashlib
import binascii
import subprocess
import dataclasses


import common

from mfc.case import *


Tend = 0.25
Nt   = 50
mydt = 0.0005


BASE_CASE = Case(
    logistics=Logistics(),
    domain=ComputationalDomain(
        cells=Cells(x=0),
        domain=SpacialDomain(
            x=AxisDomain(begin=DFLT_REAL, end=DFLT_REAL)
        ),
        time=Time(dt=mydt, end=int(Nt), save=int(Nt)),
    ),
    database=DatabseStructure(write=DatabaseWrite(prim_vars=True)),
    algorithm=SimulationAlgorithm(
        pref=101325.0,
        rhoref=1000.0,
        adv_alphan=True,
        boundary=BoundaryConditions(x=AxisBoundaryCondition()),
        weno=WenoParameters(
            order=5,
            epsilon=1.E-16,            
            variables=WenoVariables.PRIMITIVE,
        ),
        riemann_solver=RiemannSolver.HLLC,
        model=MulticomponentModel.EQUATION_5,
        wave_speeds=WaveSpeedEstimation.DIRECT,
        time_stepper=TimeStepper.RUNGE_KUTTA_3,
        avg_state=AverageStateEvaluation.ARITHMETIC_MEAN
    ),
    patches=[
        Patch(
            r0=1,
            v0=0,
            pressure=1.0,
            alpha_rho=[1.E+00],
            alpha=[1.]
        ),
        Patch(
            r0=1,
            v0=0,
            pressure=0.5,
            alpha_rho=[0.5],
            alpha=[1.]
        ),
        Patch(
            r0=1,
            v0=0,
            pressure=0.1,
            alpha_rho=[0.125],
            alpha=[1.]
        ),
    ],
    fluids=Fluids(1, [
        Fluid(
            gamma=1.E+00/(1.4-1.E+00),
            pi_inf=0.0
        )
    ]),
    bubbles=Bubbles(
        R0ref=1e-05,
        weber=13.927835051546392,
        Re_inv=0.009954269975623245,
        cavitation=0.9769178386380458,
        thermal=ThermalModel.TRANSFER,
        model=BubbleModel.RAYLEIGH_PLESSET,
        distribution=BubbleDistribution.LOGNORMAL_NORMAL,
        poly_sigma=0.3,
        R0_type=1,
        nnode=4,
        sigR=0.1,
        sigV=0.1,
    ),
    acoustic=AcousticParameters(
        monopoles=[
            Monopole(
                location=Point(x=0.5),
                magnitude=1.0,
                length=0.25
            )
        ]
    )
)


@dataclasses.dataclass
class TestCaseVariation:
    path:  str
    value: typing.Any


class TCV(TestCaseVariation):
    pass


@dataclasses.dataclass
class TestCaseVariations:
    trace:      str
    variations: typing.List[TCV]

    def add(self, variation: TCV):
        self.variations.append(variation)

    def adds(self, variations: typing.List[TCV]):
        for v in variations:
            self.add(v)


class TCVS(TestCaseVariations):
    pass


@dataclasses.dataclass(init=False)
class TestCase:
    trace: str
    case:  Case
    ppn:   int

    def __init__(self, stack: typing.List[TCVS], ppn = None) -> None:
        self.trace = stack.gen_trace()
        self.ppn   = ppn if ppn is not None else 1
        self.case  = copy.deepcopy(BASE_CASE)

        for tcv in stack.gen_variations():
            self.case.set(tcv.path, tcv.value)
        
    def run(self, args: dict) -> subprocess.CompletedProcess:
        command: str = f'''\
./mfc.sh run "{self.get_dirpath()}/case.json" -m "{args["mode"]}" -n {self.ppn} \
-t pre_process simulation -b {args['binary']} -j {args["jobs"]} 2>&1\
'''

        return subprocess.run(command, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True,
                              shell=True)

    def get_keys(self) -> str:
        return self.params.keys()

    def has_parameter(self, key: str)-> bool:
        return key in self.get_keys()

    def get_uuid(self) -> str:
        return hex(binascii.crc32(hashlib.sha1(str(self.trace).encode()).digest())).upper()[2:].zfill(8)

    def get_dirpath(self):
        return os.sep.join([common.MFC_TESTDIR, self.get_uuid()])

    def create_directory(self):
        dirpath = self.get_dirpath()

        common.create_directory(dirpath)

        common.file_write(f"{dirpath}/case.json", str(self.case))

    def __getitem__(self, key: str) -> str:
        if key not in self.params:
            raise common.MFCException(f"Case {self.trace}: Parameter {key} does not exist.")

        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = val


@dataclasses.dataclass
class CaseGeneratorStack:
    stack: typing.List[TCVS]

    def __init__(self) -> None:
        self.stack = []

    def size(self) -> int:
        return len(self.stack)

    def push(self, tcvs: TCVS) -> None:
        self.stack.append(tcvs)

    def pop(self) -> None:
        return self.stack.pop()
    
    def gen_trace(self) -> str:     
        traces = []
        
        for trace in [ tcvs.trace for tcvs in self.stack ]:
            if not common.isspace(trace):
                traces.append(trace)
        
        return ' -> '.join(traces)

    def gen_variations(self) -> TCV:
        result: TCV = []

        for tcvs in self.stack:
            for variations in tcvs.variations:
                result.append(variations)

        return result


def create_case(stack: CaseGeneratorStack,
                tcvs:  TCVS = None,
                ppn:   int  = None) -> TestCase:
    newStack = copy.deepcopy(stack)
    if tcvs is not None:
        newStack.push(tcvs)

    return TestCase(newStack, ppn)
