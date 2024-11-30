#!/usr/bin/env python3
import math
import json

# Numerical setup
Nx      = 201  # Number of grid points in x
Ny      = 201  # Number of grid points in y
dx      = 1./(1.*(Nx+1))  # Grid spacing in x
dy      = 1./(1.*(Ny+1))  # Grid spacing in y

Tend    = 64E-06  # End time
Nt      = 2000 # 2000  # Number of time steps
mydt    = Tend/(1.*Nt)  # Time step size

# Configuring case dictionary
print(json.dumps({
        # Logistics ================================================                 
        'run_time_info'                : 'F',                       
        # ==========================================================

        # Computational Domain Parameters ==========================
        'x_domain%beg'                 : 0.E+00,  # x start                    
        'x_domain%end'                 : 2.E+00,  # x end                    
        'y_domain%beg'                 : 0.E+00,  # y start
        'y_domain%end'                 : 1.E+00,  # y end 
        'm'                            : Nx,  # Number of grid points in x direction                       
        'n'                            : Ny,  # Number of grid points in y direction                       
        'p'                            : 0,  # Number of grid points in z (for 3D, change this)                      
        'dt'                           : 1e-6,  # Time step size                      
        't_step_start'                 : 0,  # Start time                         
        't_step_stop'                  : Nt,  # End time                         
        't_step_save'                  : 500,  # Save frequency
        # ==========================================================
                                                                            
        # Simulation Algorithm Parameters ==========================
        'num_patches'                  : 1,  # Two patches                        
        'model_eqns'                   : 2,  # Number of model equations                       
        'alt_soundspeed'               : 'F',                      
        'num_fluids'                   : 2,        
        'low_Mach'                     : 0,
        'mpp_lim'                      : 'F',                      
        # ' mixture_err'                  : 'F',                      
        'time_stepper'                 : 3,                                               
        'weno_order'                   : 5,                        
        'weno_eps'                     : 1.E-16,
        'weno_Re_flux'                 : 'F',  
        'weno_avg'                     : 'F',
        'mapped_weno'                  : 'F',                     
        'null_weights'                 : 'F',                      
        'mp_weno'                      : 'F',                      
        'riemann_solver'               : 1,                        
        'wave_speeds'                  : 1,                        
        'avg_state'                    : 2,                        
        'bc_x%beg'                     : -3,                       
        'bc_x%end'                     : -3,                       
        'bc_y%beg'                     : -3,  # Boundary conditions for y direction
        'bc_y%end'                     : -3,
        'num_bc_patches'               : 1,
        'patch_bc(1)%type'             : -17,
        'patch_bc(1)%dir'              : 1,
        'patch_bc(1)%loc'              : -1,
        'patch_bc(1)%geometry'         : 1,
        'patch_bc(1)%centroid(1)'      : 0, 
        'patch_bc(1)%centroid(2)'      : 0.5,
        'patch_bc(1)%length(2)'        : 0.26,
        'patch_bc(1)%vel(1)'           : 10,
        'patch_bc(1)%vel(2)'           : 0,
        # ==========================================================

        # Turning on IB ================================
        'ib'               : 'T',
        'num_ibs'          : 2,    

        # ==========================================================

        # Formatted Database Files Structure Parameters ============
        'format'                       : 1,                        
        'precision'                    : 2,                        
        'prim_vars_wrt'                :'T',                       
        'parallel_io'                  :'T',                       
        # ==========================================================

        # Patch 1 (background flow) ===================
        'patch_icpp(1)%geometry'       : 3,  # 2D geometry                
        'patch_icpp(1)%x_centroid'     : 1.0,  # x-center                    
        'patch_icpp(1)%y_centroid'     : 0.5,  # y-center
        'patch_icpp(1)%length_x'       : 2.0,   # x-length                    
        'patch_icpp(1)%length_y'       : 1.0,   # y-length
        'patch_icpp(1)%vel(1)'         : 0.0,
        'patch_icpp(1)%vel(2)'         : 0.0,   # y-velocity    '100*sin(3*x*pi)'
        'patch_icpp(1)%pres'           : 1.E5,  # Pressure                    
        'patch_icpp(1)%alpha_rho(1)'   : 1000,  # Density    
        'patch_icpp(1)%alpha_rho(2)'   : 0.,        
        'patch_icpp(1)%alpha(1)'       : 1,                       
        'patch_icpp(1)%alpha(2)'       : 0.,                          
        'patch_icpp(1)%tau_e(1)'       : 0.0,             


        # ==========================================================

        # Patch 2 (hypo material in the center) ================
        'patch_ib(1)%geometry'       : 3,  # 2D geometry      
        # 'patch_ib(1)%hcid'           : 201,                    
        'patch_ib(1)%x_centroid'     : 0.5,  # x-center                    
        'patch_ib(1)%y_centroid'     : 0.65,  # y-center                    
        'patch_ib(1)%length_x'       : 1.0,   # x-length                   
        'patch_ib(1)%length_y'       : 0.04,   # y-length       
        'patch_ib(1)%slip'           : 'T',

        # ==========================================================

        # Patch 3 (hypo material in the center) ================
        'patch_ib(2)%geometry'       : 3,  # 2D geometry      
        # 'patch_ib(1)%hcid'           : 201,                    
        'patch_ib(2)%x_centroid'     : 0.5,  # x-center                    
        'patch_ib(2)%y_centroid'     : 0.35,  # y-center                    
        'patch_ib(2)%length_x'       : 1.0,   # x-length                   
        'patch_ib(2)%length_y'       : 0.04,   # y-length       
        'patch_ib(2)%slip'           : 'T',

        # Fluids Physical Parameters ===============================
        'fluid_pp(1)%gamma'            : 1.E+00/(6.12E+00-1.E+00),   
        'fluid_pp(1)%pi_inf'           : 6.12E+00*3.43E+08/(6.12E+00 - 1.E+00), 
        # 'fluid_pp(1)%G'                : 0,        
        'fluid_pp(2)%gamma'            : 1.E+00/(1.3E+00-1.E+00),   
        'fluid_pp(2)%pi_inf'           : 1.3E+00*2.E+08/(1.3E+00 - 1.E+00),                        
        # 'fluid_pp(2)%G'                : 2.7E+05/(2.E+00*(1.E+00 + 0.4E+00)),     
        'fluid_pp(2)%G'                : 1.E7,                
        # ==========================================================
}))
