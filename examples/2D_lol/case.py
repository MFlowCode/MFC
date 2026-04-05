#!/usr/bin/env python3
import json
import math
import cantera as ct
pA = 101325
rhoA = 1.29
gam = 1.4
c_l = math.sqrt(1.4 * pA / rhoA)

pS = 2.458 * pA
velS = 0.601 * c_l
rhoS = 1.862 * rhoA

leng = 1e-2
Ny = 700
Nx = Ny
dx = leng / Nx

time_end = 15 * leng / velS
cfl = 0.5

dt = cfl * dx / c_l
Nt = int(time_end / dt)

eps = 1e-5

ctfile = "h2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPY = 1125, ct.one_atm, "O2:0.21,N2:0.79"
# Configuring case dictionary
case =      {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -leng / 2.0,
            "x_domain%end": 3 * leng / 2,
            "y_domain%beg": -leng,
            "y_domain%end": leng,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "cfl_adap_dt": "T",
            "t_stop": time_end,
            "t_save": time_end / 50,
            "n_start": 0,
            "cfl_target": 0.8,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "viscous": "T",
        #    "elliptic_smoothing": "T",
        #    "elliptic_smoothing_iters": 50,
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
                        "bc_x%Twall_in" : 500.0,
                                             "bc_x%isothermal_in": "T",
                                             "bc_x%isothermal_out": "T",
                        "bc_x%Twall_out" : 500.0,
                                   "bc_y%Twall_in" : 500.0,
                                             "bc_y%isothermal_in": "T",
                                             "bc_y%isothermal_out": "T",
                        "bc_y%Twall_out" : 500.0,
            "num_bc_patches": 4,
            "patch_bc(1)%dir": 1,
            "patch_bc(1)%loc": -1,
            "patch_bc(1)%geometry": 1,
            "patch_bc(1)%type": -17,
            "patch_bc(1)%centroid(2)": 0.0,
            "patch_bc(1)%length(2)": leng / 4,
                        "patch_bc(1)%dir": 1,
            "patch_bc(2)%loc": 1,
            "patch_bc(2)%geometry": 1,
            "patch_bc(2)%type": -17,
            "patch_bc(2)%centroid(2)": 0.0,
            "patch_bc(2)%length(2)": leng / 4,
                        "patch_bc(3)%dir": 2,
            "patch_bc(3)%loc": -1,
            "patch_bc(3)%geometry": 1,
            "patch_bc(3)%type": -3,
            "patch_bc(3)%centroid(1)": leng/2,
            "patch_bc(3)%length(1)": leng / 3,
                        "patch_bc(4)%dir": 2,
            "patch_bc(4)%loc": 1,
            "patch_bc(4)%geometry": 1,
            "patch_bc(4)%type": -3,
            "patch_bc(4)%centroid(1)": leng/2,
            "patch_bc(4)%length(1)": leng / 3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
              "chemistry": "T" ,
      "chem_params%diffusion": "T",
      "chem_params%reactions": "F",
      "cantera_file": ctfile,
       "chem_wrt_T": "T",
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * leng,
            "patch_icpp(1)%length_y": 10 * leng,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": pA,
            "patch_icpp(1)%alpha_rho(1)": rhoA,
            "patch_icpp(1)%alpha(1)": 1.0 ,
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -leng / 2,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": leng / 4,
            "patch_icpp(2)%length_y": leng / 3,
            "patch_icpp(2)%vel(1)": velS,
            "patch_icpp(2)%vel(2)": 0,
            "patch_icpp(2)%pres": pS,
            "patch_icpp(2)%alpha_rho(1)": rhoS,
            "patch_icpp(2)%alpha(1)": 1,

                        "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": 3*leng / 2,
            "patch_icpp(3)%y_centroid": 0.0,
            "patch_icpp(3)%length_x": leng / 4,
            "patch_icpp(3)%length_y": leng / 3,
            "patch_icpp(3)%vel(1)": -velS,
            "patch_icpp(3)%vel(2)": 0,
            "patch_icpp(3)%pres": pS,
            "patch_icpp(3)%alpha_rho(1)": rhoS,
            "patch_icpp(3)%alpha(1)": 1,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 100000,
        }
for i in range(len(sol_L.Y)):
           case[f"chem_wrt_Y({i + 1})"] = "T"
           case[f"patch_icpp(1)%Y({i+1})"] = sol_L.Y[i]
           case[f"patch_icpp(2)%Y({i+1})"] = sol_L.Y[i]
           case[f"patch_icpp(3)%Y({i+1})"] = sol_L.Y[i]
if __name__ == "__main__":
     print(json.dumps(case))
