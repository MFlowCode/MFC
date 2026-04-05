#!/usr/bin/env python3
import json
import math
import cantera as ct

Lx = 0.04347
Ly = 0.01466+3.048/1000
Cav_Lx = 43.47/1000-21.99/1000-14.19/1000
Cav_Ly = 3.048/1000
Cav_CenterX = (Lx - 21.99/1000-14.15/1000)/2
Cav_CenterY = 3.048/2000 
Cav_Center2X = 2*Cav_CenterX+14.15/1000+21.99/2000
Cav_22 = Cav_Lx+14.15/2000.0

ctfile = "h2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPY = 1125, ct.one_atm, "O2:0.21,N2:0.79"
# Configuring case dictionary
case =     {    
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.04347,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "m": 699,
            "n": 699,
            "p": 0,
            "dt": 4.0e-08,
            "t_step_start": 0,
                        "t_step_stop": 75000,
            "t_step_save": 4500,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "mp_weno": "F",
            # 'recon_type'                   : 1,
            "weno_order": 5,
            "weno_eps": 1e-16,
            #'muscl_order'                  : 2,
            #'muscl_lim'                    : 1,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -7,
            "bc_x%end": -3,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%isothermal_in": "T",
            "bc_y%Twall_in" : 600.0,
            "bc_y%isothermal_out": "T",
            "bc_y%Twall_out" : 850.0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
                # "ib": "T",
         #   "num_ibs": 1,
            "omega_wrt(3)": "T",
          "fd_order": 2,
              "chemistry": "T" ,
      "chem_params%diffusion": "T",
      "chem_params%reactions": "F",
      "cantera_file": ctfile,
       "chem_wrt_T": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 291,
            "patch_icpp(1)%x_centroid": Lx/2,
            "patch_icpp(1)%y_centroid": Ly/2,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y":Ly,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%pres": 101325,
            "patch_icpp(1)%alpha_rho(1)": 1.00,
            "patch_icpp(1)%alpha(1)": 1,
            # Patch 1: IBM
 #            "patch_ib(1)%geometry": 3,
 #           "patch_ib(1)%x_centroid": 0,
 #           "patch_ib(1)%y_centroid": 0,
 #           "patch_ib(1)%length_x": 2*Cav_Lx,
 #           "patch_ib(1)%length_y": 2*Cav_Ly+Ly/5,
 #           "patch_ib(1)%slip": "F",

     #       # Patch 2: IBM
#            "patch_ib(2)%geometry": 3,
#            "patch_ib(2)%x_centroid": Lx,
#            "patch_ib(2)%y_centroid": 0,
#            "patch_ib(2)%length_x": 2*21.99/1000,
#            "patch_ib(2)%length_y": 2*Cav_Ly+Ly/5,
#            "patch_ib(2)%slip": "F",

            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0e00,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 100000,
        }
for i in range(len(sol_L.Y)):
           case[f"chem_wrt_Y({i + 1})"] = "T"
           case[f"patch_icpp(1)%Y({i+1})"] = sol_L.Y[i]
   
if __name__ == "__main__":
     print(json.dumps(case))
