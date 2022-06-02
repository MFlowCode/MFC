#!/usr/bin/env python3


import math, json

from mfc.case import *

x0      = 10.E-06
p0      = 101325.
rho0    = 1.E+03
c0      = math.sqrt( p0/rho0 )
patm    = 1.

# water props ==================================================================
n_tait  = 7.1
B_tait  = 306.E+06 / p0
mul0    = -10.002E-03     #viscosity
ss      = 0.07275       #surface tension
pv      = 2.3388E+03    #vapor pressure

gamma_v = 1.33
M_v     = 18.02
mu_v    = 0.8816E-05
k_v     = 0.019426

#air props
gamma_n = 1.4
M_n     = 28.97
mu_n    = 1.8E-05
k_n     = 0.02556

#air props
# gamma_gas = gamma_n
gamma_gas = 1.4

#reference bubble size
R0ref   = 10.E-06

pa      = 0.1 * 1.E+06 / 101325.

#Characteristic velocity
uu = math.sqrt( p0/rho0 )
#Cavitation number
Ca = 1.
# Ca = (p0 - pv)/(rho0*(uu**2.))
#Weber number
We = rho0*(uu**2.)*R0ref/ss
#Inv. bubble Reynolds number
Re_inv = mul0/(rho0*uu*R0ref)

# IC setup =====================================================================
vf0     = 1.E-4
n0      = vf0/(math.pi*4.E+00/3.E+00)

cact    = 1475.
t0      = x0/c0

nbubbles = 1
myr0    = R0ref

cfl     = 0.1
Nx      = 30
Ldomain = 20.E-03
L       = Ldomain/x0
dx      = L/float(Nx)
dt      = cfl*dx*c0/cact
Lpulse  = 0.3*Ldomain
Tpulse  = Lpulse/cact
Tfinal  = 0.25*10.*Tpulse*c0/x0
Nt      = int(Tfinal/dt)

dt = dt * 0.1

Nfiles = 20.
Nout   = int(math.ceil(Nt/Nfiles))
Nt     = int(Nout*Nfiles)

# ==============================================================================

print(Case(
    logistics=Logistics(
        run_time_info=True
    ),
    domain=ComputationalDomain(
        cells=Cells(x=Nx),
        domain=SpacialDomain(
            x=AxisDomain(begin=-10.E-03/x0, end=10.E-03/x0)
        ),
        time=Time(dt=dt, end=30000, save=1000),
    ),
    fluids=[
        Fluid(
            gamma=1.E+00/(n_tait-1.E+00),
            pi_inf=n_tait*B_tait/(n_tait-1.)
        ),
#       Fluid(
#           gamma=1./(gamma_gas-1.),
#           pi_inf=0.0E+00
#       )
    ],
    patches=[
        Patch(
            smooth_patch_id=1,
            geometry=PatchGeometry.D1_LINE_SEGMENT,
            centroid=Point(x=0),
            length=Vector(x=20.E-03/x0),
            velocity=Vector(x=0.0),
            pressure=1,
            r0=1,
            v0=0,
            alpha_rho=[(1.-vf0)*1.E+03/rho0],
            alpha=[vf0],
        ),
    ],
    bubbles=Bubbles(
        bubbles=True,
        model=BubbleModel.RAYLEIGH_PLESSET,
        number=1,
        nnode=4,
        qbmm=True,
        thermal=ThermalModel.TRANSFER,
        polytropic=True,
        R0ref=myr0,
        cavitation=Ca,
        Re_inv=Re_inv,
        sigR=0.1,
        sigV=0.1,
        rhoRV=0.0,
        distribution=BubbleDistribution.BINORMAL
    ),
    database=DatabseStructure(
        precision=FloatingPrecision.DOUBLE,
        write=DatabaseWrite(
            prim_vars=True,
            probe=True
        ),
        fd_order=1,
        parallel_io=True,
        format=DatabaseFormat.SILO_HDF5,
        alt_soundspeed=False
    ),
    algorithm=SimulationAlgorithm(
        model=MulticomponentModel.EQUATION_5,
        weno=WenoParameters(
            order=3,
            mapped=True,
            mp=False,
            variables=WenoVariables.PRIMITIVE
        ),
        boundary=BoundaryConditions(
            x=AxisBoundaryCondition(
                begin=BoundaryCondition.PERIODIC,
                end=BoundaryCondition.PERIODIC
            )
        ),
        riemann_solver=RiemannSolver.HLLC,
        time_stepper=TimeStepper.RUNGE_KUTTA_1,
        wave_speeds=WaveSpeedEstimation.DIRECT,
        avg_state=AverageStateEvaluation.ARITHMETIC_MEAN,
        mixture_err=False,
        adv_alphan=True,
        mpp_lim=False,
        null_weights=False,
        pref=p0,
        rhoref=rho0
    ),
    acoustic=AcousticParameters()
))

"""
# Configuring case dictionary
print(json.dumps({
    # Formatted Database Files Structure Parameters =============================
    'num_probes'                   : 1,
    'probe(1)%x'                   : 0.,
    # ===========================================================================

}))

# ==============================================================================
"""
