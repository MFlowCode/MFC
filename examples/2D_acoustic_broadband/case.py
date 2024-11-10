#!/usr/bin/env python3

import json, math

print(json.dumps({
# Logistics ================================================================
'run_time_info' : 'T',
# ==========================================================================

# Computational Domain Parameters ==========================================
'x_domain%beg'                 : 0,
'x_domain%end'                 : 0.3,
'y_domain%beg'                 : 0,
'y_domain%end'                 : 0.1,
'm'                            : 299,
'n'                            : 99,
'p'                            : 0,
'dt'                           : 5e-7,
't_step_start'                 : 0,
't_step_stop'                  : 50000,
't_step_save'                  : 500,
# ==========================================================================

# Simulation Algorithm Parameters ==========================================
'num_patches'                  : 1,
'model_eqns'                   : 2,
'alt_soundspeed'               : 'F',
'num_fluids'                   : 1,
'mpp_lim'                      : 'F',
'mixture_err'                  : 'F',
'time_stepper'                 : 3,
'weno_order'                   : 5,
'weno_eps'                     : 1.E-16,
'teno'                         : 'T',
'teno_CT'                      : 1E-8,
'null_weights'                 : 'F',
'mp_weno'                      : 'F',
'riemann_solver'               : 2,
'wave_speeds'                  : 1,
'avg_state'                    : 2,
'bc_x%beg'                     : -6,
'bc_x%end'                     : -6,
'bc_y%beg'                     : -6,
'bc_y%end'                     : -6,
# ==========================================================================

# Formatted Database Files Structure Parameters ============================
'format'                       : 1,
'precision'                    : 2,
'prim_vars_wrt'                :'T',
'parallel_io'                  :'T',
'probe_wrt'                    :'T',
'fd_order'                     : 2,
'num_probes'                   : 1,
'probe(1)%x'                   : 0.13,
'probe(1)%y'                   : 0.05,

# ==========================================================================

# Patch 1 Liquid ===========================================================
'patch_icpp(1)%geometry'       : 3,
'patch_icpp(1)%x_centroid'     : 0.15,
'patch_icpp(1)%y_centroid'     : 0.05,
'patch_icpp(1)%length_x'       : 0.3,
'patch_icpp(1)%length_y'       : 0.1,
'patch_icpp(1)%vel(1)'         : 0.0,
'patch_icpp(1)%vel(2)'         : 0.0,
'patch_icpp(1)%pres'           : 1E+05,
'patch_icpp(1)%alpha_rho(1)'   : 1.19,
'patch_icpp(1)%alpha(1)'       : 1.0,
# ==========================================================================

# Acoustic source ==========================================================
'acoustic_source'              : 'T',
'num_source'                   : 1,
'acoustic(1)%support'          : 2,
'acoustic(1)%loc(1)'           : 0.1,
'acoustic(1)%loc(2)'           : 0.05,
'acoustic(1)%dir'              : math.pi * 0,
'acoustic(1)%length'           : 0.1,
'acoustic(1)%pulse'            : 4,
'acoustic(1)%npulse'           : 1000000,
'acoustic(1)%mag'              : 10.,
# ==========================================================================

# Fluids Physical Parameters ===============================================
'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),
'fluid_pp(1)%pi_inf'           : 0,
# ==========================================================================
}))

# ==============================================================================
