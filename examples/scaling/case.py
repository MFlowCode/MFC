#!/usr/bin/env python3
import sys, json, math, typing, argparse

parser = argparse.ArgumentParser(
    prog="scaling_and_perf",
    description="Weak- and strong-scaling and performance benchmark case.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "--mfc",
    type=json.loads,
    default="{}",
    metavar="DICT",
    help="MFC's toolchain's internal state.",
)
parser.add_argument(
    "-s",
    "--scaling",
    type=str,
    metavar="SCALING",
    choices=["weak", "strong"],
    help="Whether weak- or strong-scaling is being exercised.",
)
parser.add_argument(
    "-m",
    "--memory",
    type=int,
    metavar="MEMORY",
    help="Weak scaling: memory per rank in GB. Strong scaling: global memory in GB. Used to determine cell count.",
)
parser.add_argument("--rdma_mpi", metavar="RDMA", type=str, choices=["T", "F"], default="F", help="Enable RDMA-aware MPI optimizations.")
parser.add_argument("--n-steps", metavar="N", type=int, default=20, help="Number of time steps to simulate.")
parser.add_argument("--n-save", metavar="NS", type=int, default=20, help="Number of time steps between saves.")
args = parser.parse_args()

if args.scaling is None:
    parser.print_help()
    sys.exit(1)

# approx The number of cells per GB of memory. The exact value is not important.
cpg = 8000000 / 16.0

# Number of ranks.
nranks = args.mfc["nodes"] * args.mfc["tasks_per_node"]


# This subroutine finds three factors of n that are as close to each other as possible.
def closest_three_factors(n):
    best_triplet = None
    min_range = float("inf")

    # Iterate over possible first factor a
    for factor_one in range(1, int(n ** (1 / 3)) + 2):  # factor_one should be around the cube root of n
        if n % factor_one == 0:
            n1 = n // factor_one  # Remaining part

            # Iterate over possible second factor b
            for factor_two in range(factor_one, int(math.sqrt(n1)) + 2):  # factor_two should be around sqrt of n1
                if n1 % factor_two == 0:
                    factor_three = n1 // factor_two  # Third factor

                    triplet_range = factor_three - factor_one  # Spread of the numbers
                    if triplet_range < min_range:
                        min_range = triplet_range
                        best_triplet = (factor_one, factor_two, factor_three)

    return best_triplet


def nxyz_from_ncells_weak(ncells: float) -> typing.Tuple[int, int, int]:
    s = math.floor(ncells ** (1 / 3))
    ND = closest_three_factors(nranks)
    if any(N < 4 for N in ND) and nranks > 64:
        raise RuntimeError(f"Cannot represent {nranks} ranks with at least 4 partitions in each direction.")
    N1 = ND[0] * s - 1
    N2 = ND[1] * s - 1
    N3 = ND[2] * s - 1
    L1 = ND[0]
    L2 = ND[1]
    L3 = ND[2]
    return N1, N2, N3, L1, L2, L3


def nxyz_from_ncells_strong(ncells: float) -> typing.Tuple[int, int, int]:
    s = round(ncells ** (1 / 3))
    L1 = 4
    L2 = 4
    L3 = 4
    return s, s, s, L1, L2, L3


if args.scaling == "weak":
    Nx, Ny, Nz, Lx, Ly, Lz = nxyz_from_ncells_weak(cpg * args.memory)
else:
    Nx, Ny, Nz, Lx, Ly, Lz = nxyz_from_ncells_strong(cpg * args.memory)

# Atmospheric pressure - Pa (used as reference value)
patm = 101325

# Initial Droplet Diameter / Reference length - m
D0 = 1.0e-3

# initial shock distance from the y axis. Note that the droplet center is located at y = 0. Thus, the distance from the shock to
# the droplet is about D0/8
ISD = 5.0 / 8 * D0

## pre-shock properties - AIR
p0a = patm  # pressure - Pa
rho0a = 1.204  # density - kg/m3
gama = 1.40  # gamma
pia = 0  # pi infinity - Pa
c_a = math.sqrt(gama * (p0a + pia) / rho0a)  # speed of sound - M/s

## Droplet - WATER
rho0w = 1000  # density - kg/m3
p0w = patm  # pressure - Pa
gamw = 6.12  # gamma
piw = 3.43e08  # pi infty - Pa
c_w = math.sqrt(gamw * (p0w + piw) / rho0w)  # speed of sound - m/s

# Shock Mach number of interest
Min = 2.4

# Pos to pre shock ratios - AIR
psOp0a = (Min**2 - 1) * 2 * gama / (gama + 1) + 1  # pressure
rhosOrho0a = (1 + (gama + 1) / (gama - 1) * psOp0a) / ((gama + 1) / (gama - 1) + psOp0a)  # density
ss = Min * c_a  # shock speed of sound - m/s

# post-shock conditions - AIR
ps = psOp0a * p0a  # pressure - Pa
rhos = rhosOrho0a * rho0a  # density - kg / m3
c_s = math.sqrt(gama * (ps + pia) / rhos)  # post shock speed of sound - m/s
vel = c_a / gama * (psOp0a - 1.0) * p0a / (p0a + pia) / Min  # velocity at the post shock - m/s

# Domain extents
xb = -Lx * D0 / 2
xe = Lx * D0 / 2
yb = -Ly * D0 / 2
ye = Ly * D0 / 2
zb = -Lz * D0 / 2
ze = Lz * D0 / 2

# Calculating time step
dx = (xe - xb) / Nx
cfl = 0.05
dt = dx * cfl / ss

# Configuring case args.mfcionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            "rdma_mpi": args.rdma_mpi,
            # Computational Domain Parameters
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "y_domain%beg": yb,
            "y_domain%end": ye,
            "z_domain%beg": zb,
            "z_domain%end": ze,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": args.n_steps,
            "t_step_save": args.n_save,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "file_per_process": "T",
            # I will use 1 for WATER properties, and 2 for AIR properties
            # Patch 1: Background (AIR - 2)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": (xb + xe) / 2,
            "patch_icpp(1)%y_centroid": (yb + ye) / 2,
            "patch_icpp(1)%z_centroid": (zb + ze) / 2,
            "patch_icpp(1)%length_x": (xe - xb),
            "patch_icpp(1)%length_y": (ye - yb),
            "patch_icpp(1)%length_z": (ze - zb),
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%vel(3)": 0.0e00,
            "patch_icpp(1)%pres": p0a,
            "patch_icpp(1)%alpha_rho(1)": 0.0e00,
            "patch_icpp(1)%alpha_rho(2)": rho0a,
            "patch_icpp(1)%alpha(1)": 0.0e00,
            "patch_icpp(1)%alpha(2)": 1.0e00,
            # Patch 2: Shocked state (AIR - 2)
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -ISD - (xe - xb) / 2,
            "patch_icpp(2)%y_centroid": (yb + ye) / 2,
            "patch_icpp(2)%z_centroid": (zb + ze) / 2,
            "patch_icpp(2)%length_x": (xe - xb),
            "patch_icpp(2)%length_y": (ye - yb),
            "patch_icpp(2)%length_z": (ze - zb),
            "patch_icpp(2)%vel(1)": vel,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%vel(3)": 0.0e00,
            "patch_icpp(2)%pres": ps,
            "patch_icpp(2)%alpha_rho(1)": 0.0e00,
            "patch_icpp(2)%alpha_rho(2)": rhos,
            "patch_icpp(2)%alpha(1)": 0.0e00,
            "patch_icpp(2)%alpha(2)": 1.0e00,
            # Patch 3: Droplet (WATER - 1)
            "patch_icpp(3)%geometry": 8,
            "patch_icpp(3)%x_centroid": 0.0e00,
            "patch_icpp(3)%y_centroid": 0.0e00,
            "patch_icpp(3)%z_centroid": 0.0e00,
            "patch_icpp(3)%radius": D0 / 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": 0.0e00,
            "patch_icpp(3)%vel(2)": 0.0e00,
            "patch_icpp(3)%vel(3)": 0.0e00,
            "patch_icpp(3)%pres": p0w,
            "patch_icpp(3)%alpha_rho(1)": rho0w,
            "patch_icpp(3)%alpha_rho(2)": 0.0e00,
            "patch_icpp(3)%alpha(1)": 1.0e00,
            "patch_icpp(3)%alpha(2)": 0.0e00,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gamw - 1),
            "fluid_pp(1)%pi_inf": gamw * piw / (gamw - 1),
            "fluid_pp(2)%gamma": 1.0e00 / (gama - 1),
            "fluid_pp(2)%pi_inf": gama * pia / (gama - 1),
        }
    )
)
