#!/usr/bin/env python3

from mfc.case import *

# Numerical setup ==========================
Nx, Ny   =            100,             70
dx, dy   = 1./(1.*(Nx+1)), 1./(1.*(Ny+1))
Tend, Nt =           0.03,            100
mydt     = Tend/(1.*Nt)

# ==============================================================================

print(Case(
    logistics=Logistics(
        run_time_info=True
    ),
    domain=ComputationalDomain(
        cells=Cells(x=Nx, y=Ny),
        domain=SpacialDomain(
            x=AxisDomain(begin=0.E+00, end=1.E+00),
            y=AxisDomain(begin=0.E+00, end=1.E+00)
        ),
        time=Time(dt=mydt, end=int(Nt), save=int(Nt)),
    ),
    algorithm=SimulationAlgorithm(
        model=MulticomponentModel.EQUATION_5,
        adv_alphan=True,
        mpp_lim=False,
        mixture_err=False,
        time_stepper=TimeStepper.RUNGE_KUTTA_1,
        weno=WenoParameters(
            variables=WenoVariables.PRIMITIVE,
            order=5,
            epsilon=1.E-16,
            mapped=True,
            mp=False
        ),
        null_weights=False,
        riemann_solver=RiemannSolver.HLLC,
        wave_speeds=WaveSpeedEstimation.DIRECT,
        avg_state=AverageStateEvaluation.ARITHMETIC_MEAN,
        boundary=BoundaryConditions(
            x=AxisBoundaryCondition(
                begin=BoundaryCondition.PERIODIC,
                end=BoundaryCondition.PERIODIC
            ),
            y=AxisBoundaryCondition(
                begin=BoundaryCondition.PERIODIC,
                end=BoundaryCondition.PERIODIC
            )
        ),
    ),
    fluids=[
        Fluid(
            gamma=1.E+00/(1.4-1.E+00),
            pi_inf=0.0
        )
    ],
    patches=[
        Patch(
            geometry=7,
            centroid=Point(x=0.5, y=0.5),
            length=Vector(x=1.0, y=1.0),
            velocity=Vector(0.05, 0.05),
            pressure=1.1,
            alpha_rho=[1.E+00],
            alpha=[1.]
        )
    ],
    acoustic=AcousticParameters(),
    database=DatabseStructure(
        parallel_io=False,
        alt_soundspeed=False,
        format=DatabaseFormat.SILO_HDF5,
        precision=FloatingPrecision.DOUBLE,
        write=DatabaseWrite(
            prim_vars=True
        )
    ),
    bubbles=Bubbles()
))
