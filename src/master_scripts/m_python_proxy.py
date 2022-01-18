#!/usr/bin/env python3

## @brief This module contains the parameters and functions necessary to
##              establish a workable proxy between Python scripting and the MFC.
##              Included in this module are the dictionary definitions for the
##              user inputs to the pre-process, simulation, and post-process
##              components of the MFC, as well as functions which interconnect
##              these procedures' execution. Low-level access to the portable
##              batch system (PBS) is also included through additional dictionary
##              definitions and provides the MFC with parallel run capabilities.

# Dependencies =================================================================

# Command used to query the path of the current working directory
from os import getcwd

# Command used to query the name of the current working directory
from os.path import basename

# Output piping parameter and command-line execution function, respectively
from subprocess import PIPE, Popen

# Command used to exit the script
from sys import exit

import sys

import pathlib

# ==============================================================================


# Pre-process Dictionary =======================================================
pre_process_dict =                                                             \
    {                                                                          \
                    'case_dir'                      : None,                    \
                    'old_grid'                      : None,                    \
                    'old_ic'                        : None,                    \
                    't_step_old'                    : None,                    \
                    'm'                             : None,                    \
                    'n'                             : None,                    \
                    'p'                             : None,                    \
                    'x_domain%beg'                  : None,                    \
                    'x_domain%end'                  : None,                    \
                    'y_domain%beg'                  : None,                    \
                    'y_domain%end'                  : None,                    \
                    'z_domain%beg'                  : None,                    \
                    'z_domain%end'                  : None,                    \
                    'stretch_x'                     : None,                    \
                    'stretch_y'                     : None,                    \
                    'stretch_z'                     : None,                    \
                    'a_x'                           : None,                    \
                    'a_y'                           : None,                    \
                    'a_z'                           : None,                    \
                    'loops_x'                       : None,                    \
                    'loops_y'                       : None,                    \
                    'loops_z'                       : None,                    \
                    'x_a'                           : None,                    \
                    'y_a'                           : None,                    \
                    'z_a'                           : None,                    \
                    'x_b'                           : None,                    \
                    'y_b'                           : None,                    \
                    'z_b'                           : None,                    \
                    'cyl_coord'                     : None,                    \
                    'model_eqns'                    : None,                    \
                    'num_fluids'                    : None,                    \
                    'adv_alphan'                    : None,                    \
                    'mpp_lim'                       : None,                    \
                    'weno_order'                    : None,                    \
                    'precision'                     : None,                    \
                    'parallel_io'                   : None,                    \
                    'perturb_flow'                  : None,                    \
                    'perturb_flow_fluid'            : None,                    \
                    'perturb_sph'                   : None,                    \
                    'perturb_sph_fluid'             : None,                    \
                    'fluid_rho'                     : None,                    \
                    'fluid_rho(1)'                  : None,                    \
                    'fluid_rho(2)'                  : None,                    \
                    'fluid_rho(3)'                  : None,                    \
                    'fluid_rho(4)'                  : None,                    \
                    'fluid_rho(5)'                  : None,                    \
                    'fluid_rho(6)'                  : None,                    \
                    'fluid_rho(7)'                  : None,                    \
                    'fluid_rho(8)'                  : None,                    \
                    'fluid_rho(9)'                  : None,                    \
                    'fluid_rho(10)'                 : None,                    \
                    'bc_x%beg'                      : None,                    \
                    'bc_x%end'                      : None,                    \
                    'bc_y%beg'                      : None,                    \
                    'bc_y%end'                      : None,                    \
                    'bc_z%beg'                      : None,                    \
                    'bc_z%end'                      : None,                    \
                    'hypoelasticity'                : None,                    \
                    'num_patches'                   : None,                    \
                    'patch_icpp(1)%geometry'        : None,                    \
                    'patch_icpp(1)%x_centroid'      : None,                    \
                    'patch_icpp(1)%y_centroid'      : None,                    \
                    'patch_icpp(1)%z_centroid'      : None,                    \
                    'patch_icpp(1)%length_x'        : None,                    \
                    'patch_icpp(1)%length_y'        : None,                    \
                    'patch_icpp(1)%length_z'        : None,                    \
                    'patch_icpp(1)%radius'          : None,                    \
                    'patch_icpp(1)%radii'           : None,                    \
                    'patch_icpp(1)%radii(1)'        : None,                    \
                    'patch_icpp(1)%radii(2)'        : None,                    \
                    'patch_icpp(1)%radii(3)'        : None,                    \
                    'patch_icpp(1)%epsilon'         : None,                    \
                    'patch_icpp(1)%beta'            : None,                    \
                    'patch_icpp(1)%normal'          : None,                    \
                    'patch_icpp(1)%normal(1)'       : None,                    \
                    'patch_icpp(1)%normal(2)'       : None,                    \
                    'patch_icpp(1)%normal(3)'       : None,                    \
                    'patch_icpp(1)%smoothen'        : None,                    \
                    'patch_icpp(1)%smooth_patch_id' : None,                    \
                    'patch_icpp(1)%smooth_coeff'    : None,                    \
                    'patch_icpp(1)%alpha_rho'       : None,                    \
                    'patch_icpp(1)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(1)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(1)%rho'             : None,                    \
                    'patch_icpp(1)%vel'             : None,                    \
                    'patch_icpp(1)%vel(1)'          : None,                    \
                    'patch_icpp(1)%vel(2)'          : None,                    \
                    'patch_icpp(1)%vel(3)'          : None,                    \
                    'patch_icpp(1)%pres'            : None,                    \
                    'patch_icpp(1)%alpha'           : None,                    \
                    'patch_icpp(1)%alpha(1)'        : None,                    \
                    'patch_icpp(1)%alpha(2)'        : None,                    \
                    'patch_icpp(1)%alpha(3)'        : None,                    \
                    'patch_icpp(1)%alpha(4)'        : None,                    \
                    'patch_icpp(1)%alpha(5)'        : None,                    \
                    'patch_icpp(1)%alpha(6)'        : None,                    \
                    'patch_icpp(1)%alpha(7)'        : None,                    \
                    'patch_icpp(1)%alpha(8)'        : None,                    \
                    'patch_icpp(1)%alpha(9)'        : None,                    \
                    'patch_icpp(1)%alpha(10)'       : None,                    \
                    'patch_icpp(1)%gamma'           : None,                    \
                    'patch_icpp(1)%pi_inf'          : None,                    \
                    'patch_icpp(1)%tau_e(1)'        : None,                    \
                    'patch_icpp(1)%tau_e(2)'        : None,                    \
                    'patch_icpp(1)%tau_e(3)'        : None,                    \
                    'patch_icpp(1)%tau_e(4)'        : None,                    \
                    'patch_icpp(1)%tau_e(5)'        : None,                    \
                    'patch_icpp(1)%tau_e(6)'        : None,                    \
                    'patch_icpp(2)%geometry'        : None,                    \
                    'patch_icpp(2)%x_centroid'      : None,                    \
                    'patch_icpp(2)%y_centroid'      : None,                    \
                    'patch_icpp(2)%z_centroid'      : None,                    \
                    'patch_icpp(2)%length_x'        : None,                    \
                    'patch_icpp(2)%length_y'        : None,                    \
                    'patch_icpp(2)%length_z'        : None,                    \
                    'patch_icpp(2)%radius'          : None,                    \
                    'patch_icpp(2)%radii'           : None,                    \
                    'patch_icpp(2)%radii(1)'        : None,                    \
                    'patch_icpp(2)%radii(2)'        : None,                    \
                    'patch_icpp(2)%radii(3)'        : None,                    \
                    'patch_icpp(2)%epsilon'         : None,                    \
                    'patch_icpp(2)%beta'            : None,                    \
                    'patch_icpp(2)%normal'          : None,                    \
                    'patch_icpp(2)%normal(1)'       : None,                    \
                    'patch_icpp(2)%normal(2)'       : None,                    \
                    'patch_icpp(2)%normal(3)'       : None,                    \
                    'patch_icpp(2)%alter_patch'     : None,                    \
                    'patch_icpp(2)%alter_patch(1)'  : None,                    \
                    'patch_icpp(2)%smoothen'        : None,                    \
                    'patch_icpp(2)%smooth_patch_id' : None,                    \
                    'patch_icpp(2)%smooth_coeff'    : None,                    \
                    'patch_icpp(2)%alpha_rho'       : None,                    \
                    'patch_icpp(2)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(2)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(2)%rho'             : None,                    \
                    'patch_icpp(2)%vel'             : None,                    \
                    'patch_icpp(2)%vel(1)'          : None,                    \
                    'patch_icpp(2)%vel(2)'          : None,                    \
                    'patch_icpp(2)%vel(3)'          : None,                    \
                    'patch_icpp(2)%pres'            : None,                    \
                    'patch_icpp(2)%alpha'           : None,                    \
                    'patch_icpp(2)%alpha(1)'        : None,                    \
                    'patch_icpp(2)%alpha(2)'        : None,                    \
                    'patch_icpp(2)%alpha(3)'        : None,                    \
                    'patch_icpp(2)%alpha(4)'        : None,                    \
                    'patch_icpp(2)%alpha(5)'        : None,                    \
                    'patch_icpp(2)%alpha(6)'        : None,                    \
                    'patch_icpp(2)%alpha(7)'        : None,                    \
                    'patch_icpp(2)%alpha(8)'        : None,                    \
                    'patch_icpp(2)%alpha(9)'        : None,                    \
                    'patch_icpp(2)%alpha(10)'       : None,                    \
                    'patch_icpp(2)%gamma'           : None,                    \
                    'patch_icpp(2)%pi_inf'          : None,                    \
                    'patch_icpp(2)%tau_e(1)'        : None,                    \
                    'patch_icpp(2)%tau_e(2)'        : None,                    \
                    'patch_icpp(2)%tau_e(3)'        : None,                    \
                    'patch_icpp(2)%tau_e(4)'        : None,                    \
                    'patch_icpp(2)%tau_e(5)'        : None,                    \
                    'patch_icpp(2)%tau_e(6)'        : None,                    \
                    'patch_icpp(3)%geometry'        : None,                    \
                    'patch_icpp(3)%x_centroid'      : None,                    \
                    'patch_icpp(3)%y_centroid'      : None,                    \
                    'patch_icpp(3)%z_centroid'      : None,                    \
                    'patch_icpp(3)%length_x'        : None,                    \
                    'patch_icpp(3)%length_y'        : None,                    \
                    'patch_icpp(3)%length_z'        : None,                    \
                    'patch_icpp(3)%radius'          : None,                    \
                    'patch_icpp(3)%radii'           : None,                    \
                    'patch_icpp(3)%radii(1)'        : None,                    \
                    'patch_icpp(3)%radii(2)'        : None,                    \
                    'patch_icpp(3)%radii(3)'        : None,                    \
                    'patch_icpp(3)%epsilon'         : None,                    \
                    'patch_icpp(3)%beta'            : None,                    \
                    'patch_icpp(3)%normal'          : None,                    \
                    'patch_icpp(3)%normal(1)'       : None,                    \
                    'patch_icpp(3)%normal(2)'       : None,                    \
                    'patch_icpp(3)%normal(3)'       : None,                    \
                    'patch_icpp(3)%alter_patch'     : None,                    \
                    'patch_icpp(3)%alter_patch(1)'  : None,                    \
                    'patch_icpp(3)%alter_patch(2)'  : None,                    \
                    'patch_icpp(3)%smoothen'        : None,                    \
                    'patch_icpp(3)%smooth_patch_id' : None,                    \
                    'patch_icpp(3)%smooth_coeff'    : None,                    \
                    'patch_icpp(3)%alpha_rho'       : None,                    \
                    'patch_icpp(3)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(3)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(3)%rho'             : None,                    \
                    'patch_icpp(3)%vel'             : None,                    \
                    'patch_icpp(3)%vel(1)'          : None,                    \
                    'patch_icpp(3)%vel(2)'          : None,                    \
                    'patch_icpp(3)%vel(3)'          : None,                    \
                    'patch_icpp(3)%pres'            : None,                    \
                    'patch_icpp(3)%alpha'           : None,                    \
                    'patch_icpp(3)%alpha(1)'        : None,                    \
                    'patch_icpp(3)%alpha(2)'        : None,                    \
                    'patch_icpp(3)%alpha(3)'        : None,                    \
                    'patch_icpp(3)%alpha(4)'        : None,                    \
                    'patch_icpp(3)%alpha(5)'        : None,                    \
                    'patch_icpp(3)%alpha(6)'        : None,                    \
                    'patch_icpp(3)%alpha(7)'        : None,                    \
                    'patch_icpp(3)%alpha(8)'        : None,                    \
                    'patch_icpp(3)%alpha(9)'        : None,                    \
                    'patch_icpp(3)%alpha(10)'       : None,                    \
                    'patch_icpp(3)%gamma'           : None,                    \
                    'patch_icpp(3)%pi_inf'          : None,                    \
                    'patch_icpp(3)%tau_e(1)'        : None,                    \
                    'patch_icpp(3)%tau_e(2)'        : None,                    \
                    'patch_icpp(3)%tau_e(3)'        : None,                    \
                    'patch_icpp(3)%tau_e(4)'        : None,                    \
                    'patch_icpp(3)%tau_e(5)'        : None,                    \
                    'patch_icpp(3)%tau_e(6)'        : None,                    \
                    'patch_icpp(4)%geometry'        : None,                    \
                    'patch_icpp(4)%x_centroid'      : None,                    \
                    'patch_icpp(4)%y_centroid'      : None,                    \
                    'patch_icpp(4)%z_centroid'      : None,                    \
                    'patch_icpp(4)%length_x'        : None,                    \
                    'patch_icpp(4)%length_y'        : None,                    \
                    'patch_icpp(4)%length_z'        : None,                    \
                    'patch_icpp(4)%radius'          : None,                    \
                    'patch_icpp(4)%radii'           : None,                    \
                    'patch_icpp(4)%radii(1)'        : None,                    \
                    'patch_icpp(4)%radii(2)'        : None,                    \
                    'patch_icpp(4)%radii(3)'        : None,                    \
                    'patch_icpp(4)%epsilon'         : None,                    \
                    'patch_icpp(4)%beta'            : None,                    \
                    'patch_icpp(4)%normal'          : None,                    \
                    'patch_icpp(4)%normal(1)'       : None,                    \
                    'patch_icpp(4)%normal(2)'       : None,                    \
                    'patch_icpp(4)%normal(3)'       : None,                    \
                    'patch_icpp(4)%alter_patch'     : None,                    \
                    'patch_icpp(4)%alter_patch(1)'  : None,                    \
                    'patch_icpp(4)%alter_patch(2)'  : None,                    \
                    'patch_icpp(4)%alter_patch(3)'  : None,                    \
                    'patch_icpp(4)%smoothen'        : None,                    \
                    'patch_icpp(4)%smooth_patch_id' : None,                    \
                    'patch_icpp(4)%smooth_coeff'    : None,                    \
                    'patch_icpp(4)%alpha_rho'       : None,                    \
                    'patch_icpp(4)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(4)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(4)%rho'             : None,                    \
                    'patch_icpp(4)%vel'             : None,                    \
                    'patch_icpp(4)%vel(1)'          : None,                    \
                    'patch_icpp(4)%vel(2)'          : None,                    \
                    'patch_icpp(4)%vel(3)'          : None,                    \
                    'patch_icpp(4)%pres'            : None,                    \
                    'patch_icpp(4)%alpha'           : None,                    \
                    'patch_icpp(4)%alpha(1)'        : None,                    \
                    'patch_icpp(4)%alpha(2)'        : None,                    \
                    'patch_icpp(4)%alpha(3)'        : None,                    \
                    'patch_icpp(4)%alpha(4)'        : None,                    \
                    'patch_icpp(4)%alpha(5)'        : None,                    \
                    'patch_icpp(4)%alpha(6)'        : None,                    \
                    'patch_icpp(4)%alpha(7)'        : None,                    \
                    'patch_icpp(4)%alpha(8)'        : None,                    \
                    'patch_icpp(4)%alpha(9)'        : None,                    \
                    'patch_icpp(4)%alpha(10)'       : None,                    \
                    'patch_icpp(4)%gamma'           : None,                    \
                    'patch_icpp(4)%pi_inf'          : None,                    \
                    'patch_icpp(4)%tau_e(1)'        : None,                    \
                    'patch_icpp(4)%tau_e(2)'        : None,                    \
                    'patch_icpp(4)%tau_e(3)'        : None,                    \
                    'patch_icpp(4)%tau_e(4)'        : None,                    \
                    'patch_icpp(4)%tau_e(5)'        : None,                    \
                    'patch_icpp(4)%tau_e(6)'        : None,                    \
                    'patch_icpp(5)%geometry'        : None,                    \
                    'patch_icpp(5)%x_centroid'      : None,                    \
                    'patch_icpp(5)%y_centroid'      : None,                    \
                    'patch_icpp(5)%z_centroid'      : None,                    \
                    'patch_icpp(5)%length_x'        : None,                    \
                    'patch_icpp(5)%length_y'        : None,                    \
                    'patch_icpp(5)%length_z'        : None,                    \
                    'patch_icpp(5)%radius'          : None,                    \
                    'patch_icpp(5)%radii'           : None,                    \
                    'patch_icpp(5)%radii(1)'        : None,                    \
                    'patch_icpp(5)%radii(2)'        : None,                    \
                    'patch_icpp(5)%radii(3)'        : None,                    \
                    'patch_icpp(5)%epsilon'         : None,                    \
                    'patch_icpp(5)%beta'            : None,                    \
                    'patch_icpp(5)%normal'          : None,                    \
                    'patch_icpp(5)%normal(1)'       : None,                    \
                    'patch_icpp(5)%normal(2)'       : None,                    \
                    'patch_icpp(5)%normal(3)'       : None,                    \
                    'patch_icpp(5)%alter_patch'     : None,                    \
                    'patch_icpp(5)%alter_patch(1)'  : None,                    \
                    'patch_icpp(5)%alter_patch(2)'  : None,                    \
                    'patch_icpp(5)%alter_patch(3)'  : None,                    \
                    'patch_icpp(5)%alter_patch(4)'  : None,                    \
                    'patch_icpp(5)%smoothen'        : None,                    \
                    'patch_icpp(5)%smooth_patch_id' : None,                    \
                    'patch_icpp(5)%smooth_coeff'    : None,                    \
                    'patch_icpp(5)%alpha_rho'       : None,                    \
                    'patch_icpp(5)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(5)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(5)%rho'             : None,                    \
                    'patch_icpp(5)%vel'             : None,                    \
                    'patch_icpp(5)%vel(1)'          : None,                    \
                    'patch_icpp(5)%vel(2)'          : None,                    \
                    'patch_icpp(5)%vel(3)'          : None,                    \
                    'patch_icpp(5)%pres'            : None,                    \
                    'patch_icpp(5)%alpha'           : None,                    \
                    'patch_icpp(5)%alpha(1)'        : None,                    \
                    'patch_icpp(5)%alpha(2)'        : None,                    \
                    'patch_icpp(5)%alpha(3)'        : None,                    \
                    'patch_icpp(5)%alpha(4)'        : None,                    \
                    'patch_icpp(5)%alpha(5)'        : None,                    \
                    'patch_icpp(5)%alpha(6)'        : None,                    \
                    'patch_icpp(5)%alpha(7)'        : None,                    \
                    'patch_icpp(5)%alpha(8)'        : None,                    \
                    'patch_icpp(5)%alpha(9)'        : None,                    \
                    'patch_icpp(5)%alpha(10)'       : None,                    \
                    'patch_icpp(5)%gamma'           : None,                    \
                    'patch_icpp(5)%pi_inf'          : None,                    \
                    'patch_icpp(5)%tau_e(1)'        : None,                    \
                    'patch_icpp(5)%tau_e(2)'        : None,                    \
                    'patch_icpp(5)%tau_e(3)'        : None,                    \
                    'patch_icpp(5)%tau_e(4)'        : None,                    \
                    'patch_icpp(5)%tau_e(5)'        : None,                    \
                    'patch_icpp(5)%tau_e(6)'        : None,                    \
                    'patch_icpp(6)%geometry'        : None,                    \
                    'patch_icpp(6)%x_centroid'      : None,                    \
                    'patch_icpp(6)%y_centroid'      : None,                    \
                    'patch_icpp(6)%z_centroid'      : None,                    \
                    'patch_icpp(6)%length_x'        : None,                    \
                    'patch_icpp(6)%length_y'        : None,                    \
                    'patch_icpp(6)%length_z'        : None,                    \
                    'patch_icpp(6)%radius'          : None,                    \
                    'patch_icpp(6)%radii'           : None,                    \
                    'patch_icpp(6)%radii(1)'        : None,                    \
                    'patch_icpp(6)%radii(2)'        : None,                    \
                    'patch_icpp(6)%radii(3)'        : None,                    \
                    'patch_icpp(6)%epsilon'         : None,                    \
                    'patch_icpp(6)%beta'            : None,                    \
                    'patch_icpp(6)%normal'          : None,                    \
                    'patch_icpp(6)%normal(1)'       : None,                    \
                    'patch_icpp(6)%normal(2)'       : None,                    \
                    'patch_icpp(6)%normal(3)'       : None,                    \
                    'patch_icpp(6)%alter_patch'     : None,                    \
                    'patch_icpp(6)%alter_patch(1)'  : None,                    \
                    'patch_icpp(6)%alter_patch(2)'  : None,                    \
                    'patch_icpp(6)%alter_patch(3)'  : None,                    \
                    'patch_icpp(6)%alter_patch(4)'  : None,                    \
                    'patch_icpp(6)%alter_patch(5)'  : None,                    \
                    'patch_icpp(6)%smoothen'        : None,                    \
                    'patch_icpp(6)%smooth_patch_id' : None,                    \
                    'patch_icpp(6)%smooth_coeff'    : None,                    \
                    'patch_icpp(6)%alpha_rho'       : None,                    \
                    'patch_icpp(6)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(6)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(6)%rho'             : None,                    \
                    'patch_icpp(6)%vel'             : None,                    \
                    'patch_icpp(6)%vel(1)'          : None,                    \
                    'patch_icpp(6)%vel(2)'          : None,                    \
                    'patch_icpp(6)%vel(3)'          : None,                    \
                    'patch_icpp(6)%pres'            : None,                    \
                    'patch_icpp(6)%alpha'           : None,                    \
                    'patch_icpp(6)%alpha(1)'        : None,                    \
                    'patch_icpp(6)%alpha(2)'        : None,                    \
                    'patch_icpp(6)%alpha(3)'        : None,                    \
                    'patch_icpp(6)%alpha(4)'        : None,                    \
                    'patch_icpp(6)%alpha(5)'        : None,                    \
                    'patch_icpp(6)%alpha(6)'        : None,                    \
                    'patch_icpp(6)%alpha(7)'        : None,                    \
                    'patch_icpp(6)%alpha(8)'        : None,                    \
                    'patch_icpp(6)%alpha(9)'        : None,                    \
                    'patch_icpp(6)%alpha(10)'       : None,                    \
                    'patch_icpp(6)%gamma'           : None,                    \
                    'patch_icpp(6)%pi_inf'          : None,                    \
                    'patch_icpp(6)%tau_e(1)'        : None,                    \
                    'patch_icpp(6)%tau_e(2)'        : None,                    \
                    'patch_icpp(6)%tau_e(3)'        : None,                    \
                    'patch_icpp(6)%tau_e(4)'        : None,                    \
                    'patch_icpp(6)%tau_e(5)'        : None,                    \
                    'patch_icpp(6)%tau_e(6)'        : None,                    \
                    'patch_icpp(7)%geometry'        : None,                    \
                    'patch_icpp(7)%x_centroid'      : None,                    \
                    'patch_icpp(7)%y_centroid'      : None,                    \
                    'patch_icpp(7)%z_centroid'      : None,                    \
                    'patch_icpp(7)%length_x'        : None,                    \
                    'patch_icpp(7)%length_y'        : None,                    \
                    'patch_icpp(7)%length_z'        : None,                    \
                    'patch_icpp(7)%radius'          : None,                    \
                    'patch_icpp(7)%radii'           : None,                    \
                    'patch_icpp(7)%radii(1)'        : None,                    \
                    'patch_icpp(7)%radii(2)'        : None,                    \
                    'patch_icpp(7)%radii(3)'        : None,                    \
                    'patch_icpp(7)%epsilon'         : None,                    \
                    'patch_icpp(7)%beta'            : None,                    \
                    'patch_icpp(7)%normal'          : None,                    \
                    'patch_icpp(7)%normal(1)'       : None,                    \
                    'patch_icpp(7)%normal(2)'       : None,                    \
                    'patch_icpp(7)%normal(3)'       : None,                    \
                    'patch_icpp(7)%alter_patch'     : None,                    \
                    'patch_icpp(7)%alter_patch(1)'  : None,                    \
                    'patch_icpp(7)%alter_patch(2)'  : None,                    \
                    'patch_icpp(7)%alter_patch(3)'  : None,                    \
                    'patch_icpp(7)%alter_patch(4)'  : None,                    \
                    'patch_icpp(7)%alter_patch(5)'  : None,                    \
                    'patch_icpp(7)%alter_patch(6)'  : None,                    \
                    'patch_icpp(7)%smoothen'        : None,                    \
                    'patch_icpp(7)%smooth_patch_id' : None,                    \
                    'patch_icpp(7)%smooth_coeff'    : None,                    \
                    'patch_icpp(7)%alpha_rho'       : None,                    \
                    'patch_icpp(7)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(7)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(7)%rho'             : None,                    \
                    'patch_icpp(7)%vel'             : None,                    \
                    'patch_icpp(7)%vel(1)'          : None,                    \
                    'patch_icpp(7)%vel(2)'          : None,                    \
                    'patch_icpp(7)%vel(3)'          : None,                    \
                    'patch_icpp(7)%pres'            : None,                    \
                    'patch_icpp(7)%alpha'           : None,                    \
                    'patch_icpp(7)%alpha(1)'        : None,                    \
                    'patch_icpp(7)%alpha(2)'        : None,                    \
                    'patch_icpp(7)%alpha(3)'        : None,                    \
                    'patch_icpp(7)%alpha(4)'        : None,                    \
                    'patch_icpp(7)%alpha(5)'        : None,                    \
                    'patch_icpp(7)%alpha(6)'        : None,                    \
                    'patch_icpp(7)%alpha(7)'        : None,                    \
                    'patch_icpp(7)%alpha(8)'        : None,                    \
                    'patch_icpp(7)%alpha(9)'        : None,                    \
                    'patch_icpp(7)%alpha(10)'       : None,                    \
                    'patch_icpp(7)%gamma'           : None,                    \
                    'patch_icpp(7)%pi_inf'          : None,                    \
                    'patch_icpp(7)%tau_e(1)'        : None,                    \
                    'patch_icpp(7)%tau_e(2)'        : None,                    \
                    'patch_icpp(7)%tau_e(3)'        : None,                    \
                    'patch_icpp(7)%tau_e(4)'        : None,                    \
                    'patch_icpp(7)%tau_e(5)'        : None,                    \
                    'patch_icpp(7)%tau_e(6)'        : None,                    \
                    'patch_icpp(8)%geometry'        : None,                    \
                    'patch_icpp(8)%x_centroid'      : None,                    \
                    'patch_icpp(8)%y_centroid'      : None,                    \
                    'patch_icpp(8)%z_centroid'      : None,                    \
                    'patch_icpp(8)%length_x'        : None,                    \
                    'patch_icpp(8)%length_y'        : None,                    \
                    'patch_icpp(8)%length_z'        : None,                    \
                    'patch_icpp(8)%radius'          : None,                    \
                    'patch_icpp(8)%radii'           : None,                    \
                    'patch_icpp(8)%radii(1)'        : None,                    \
                    'patch_icpp(8)%radii(2)'        : None,                    \
                    'patch_icpp(8)%radii(3)'        : None,                    \
                    'patch_icpp(8)%epsilon'         : None,                    \
                    'patch_icpp(8)%beta'            : None,                    \
                    'patch_icpp(8)%normal'          : None,                    \
                    'patch_icpp(8)%normal(1)'       : None,                    \
                    'patch_icpp(8)%normal(2)'       : None,                    \
                    'patch_icpp(8)%normal(3)'       : None,                    \
                    'patch_icpp(8)%alter_patch'     : None,                    \
                    'patch_icpp(8)%alter_patch(1)'  : None,                    \
                    'patch_icpp(8)%alter_patch(2)'  : None,                    \
                    'patch_icpp(8)%alter_patch(3)'  : None,                    \
                    'patch_icpp(8)%alter_patch(4)'  : None,                    \
                    'patch_icpp(8)%alter_patch(5)'  : None,                    \
                    'patch_icpp(8)%alter_patch(6)'  : None,                    \
                    'patch_icpp(8)%alter_patch(7)'  : None,                    \
                    'patch_icpp(8)%smoothen'        : None,                    \
                    'patch_icpp(8)%smooth_patch_id' : None,                    \
                    'patch_icpp(8)%smooth_coeff'    : None,                    \
                    'patch_icpp(8)%alpha_rho'       : None,                    \
                    'patch_icpp(8)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(8)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(8)%rho'             : None,                    \
                    'patch_icpp(8)%vel'             : None,                    \
                    'patch_icpp(8)%vel(1)'          : None,                    \
                    'patch_icpp(8)%vel(2)'          : None,                    \
                    'patch_icpp(8)%vel(3)'          : None,                    \
                    'patch_icpp(8)%pres'            : None,                    \
                    'patch_icpp(8)%alpha'           : None,                    \
                    'patch_icpp(8)%alpha(1)'        : None,                    \
                    'patch_icpp(8)%alpha(2)'        : None,                    \
                    'patch_icpp(8)%alpha(3)'        : None,                    \
                    'patch_icpp(8)%alpha(4)'        : None,                    \
                    'patch_icpp(8)%alpha(5)'        : None,                    \
                    'patch_icpp(8)%alpha(6)'        : None,                    \
                    'patch_icpp(8)%alpha(7)'        : None,                    \
                    'patch_icpp(8)%alpha(8)'        : None,                    \
                    'patch_icpp(8)%alpha(9)'        : None,                    \
                    'patch_icpp(8)%alpha(10)'       : None,                    \
                    'patch_icpp(8)%gamma'           : None,                    \
                    'patch_icpp(8)%pi_inf'          : None,                    \
                    'patch_icpp(8)%tau_e(1)'        : None,                    \
                    'patch_icpp(8)%tau_e(2)'        : None,                    \
                    'patch_icpp(8)%tau_e(3)'        : None,                    \
                    'patch_icpp(8)%tau_e(4)'        : None,                    \
                    'patch_icpp(8)%tau_e(5)'        : None,                    \
                    'patch_icpp(8)%tau_e(6)'        : None,                    \
                    'patch_icpp(9)%geometry'        : None,                    \
                    'patch_icpp(9)%x_centroid'      : None,                    \
                    'patch_icpp(9)%y_centroid'      : None,                    \
                    'patch_icpp(9)%z_centroid'      : None,                    \
                    'patch_icpp(9)%length_x'        : None,                    \
                    'patch_icpp(9)%length_y'        : None,                    \
                    'patch_icpp(9)%length_z'        : None,                    \
                    'patch_icpp(9)%radius'          : None,                    \
                    'patch_icpp(9)%radii'           : None,                    \
                    'patch_icpp(9)%radii(1)'        : None,                    \
                    'patch_icpp(9)%radii(2)'        : None,                    \
                    'patch_icpp(9)%radii(3)'        : None,                    \
                    'patch_icpp(9)%epsilon'         : None,                    \
                    'patch_icpp(9)%beta'            : None,                    \
                    'patch_icpp(9)%normal'          : None,                    \
                    'patch_icpp(9)%normal(1)'       : None,                    \
                    'patch_icpp(9)%normal(2)'       : None,                    \
                    'patch_icpp(9)%normal(3)'       : None,                    \
                    'patch_icpp(9)%alter_patch'     : None,                    \
                    'patch_icpp(9)%alter_patch(1)'  : None,                    \
                    'patch_icpp(9)%alter_patch(2)'  : None,                    \
                    'patch_icpp(9)%alter_patch(3)'  : None,                    \
                    'patch_icpp(9)%alter_patch(4)'  : None,                    \
                    'patch_icpp(9)%alter_patch(5)'  : None,                    \
                    'patch_icpp(9)%alter_patch(6)'  : None,                    \
                    'patch_icpp(9)%alter_patch(7)'  : None,                    \
                    'patch_icpp(9)%alter_patch(8)'  : None,                    \
                    'patch_icpp(9)%smoothen'        : None,                    \
                    'patch_icpp(9)%smooth_patch_id' : None,                    \
                    'patch_icpp(9)%smooth_coeff'    : None,                    \
                    'patch_icpp(9)%alpha_rho'       : None,                    \
                    'patch_icpp(9)%alpha_rho(1)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(2)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(3)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(4)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(5)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(6)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(7)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(8)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(9)'    : None,                    \
                    'patch_icpp(9)%alpha_rho(10)'   : None,                    \
                    'patch_icpp(9)%rho'             : None,                    \
                    'patch_icpp(9)%vel'             : None,                    \
                    'patch_icpp(9)%vel(1)'          : None,                    \
                    'patch_icpp(9)%vel(2)'          : None,                    \
                    'patch_icpp(9)%vel(3)'          : None,                    \
                    'patch_icpp(9)%pres'            : None,                    \
                    'patch_icpp(9)%alpha'           : None,                    \
                    'patch_icpp(9)%alpha(1)'        : None,                    \
                    'patch_icpp(9)%alpha(2)'        : None,                    \
                    'patch_icpp(9)%alpha(3)'        : None,                    \
                    'patch_icpp(9)%alpha(4)'        : None,                    \
                    'patch_icpp(9)%alpha(5)'        : None,                    \
                    'patch_icpp(9)%alpha(6)'        : None,                    \
                    'patch_icpp(9)%alpha(7)'        : None,                    \
                    'patch_icpp(9)%alpha(8)'        : None,                    \
                    'patch_icpp(9)%alpha(9)'        : None,                    \
                    'patch_icpp(9)%alpha(10)'       : None,                    \
                    'patch_icpp(9)%gamma'           : None,                    \
                    'patch_icpp(9)%pi_inf'          : None,                    \
                    'patch_icpp(9)%tau_e(1)'        : None,                    \
                    'patch_icpp(9)%tau_e(2)'        : None,                    \
                    'patch_icpp(9)%tau_e(3)'        : None,                    \
                    'patch_icpp(9)%tau_e(4)'        : None,                    \
                    'patch_icpp(9)%tau_e(5)'        : None,                    \
                    'patch_icpp(9)%tau_e(6)'        : None,                    \
                    'patch_icpp(10)%geometry'       : None,                    \
                    'patch_icpp(10)%x_centroid'     : None,                    \
                    'patch_icpp(10)%y_centroid'     : None,                    \
                    'patch_icpp(10)%z_centroid'     : None,                    \
                    'patch_icpp(10)%length_x'       : None,                    \
                    'patch_icpp(10)%length_y'       : None,                    \
                    'patch_icpp(10)%length_z'       : None,                    \
                    'patch_icpp(10)%radius'         : None,                    \
                    'patch_icpp(10)%radii'          : None,                    \
                    'patch_icpp(10)%radii(1)'       : None,                    \
                    'patch_icpp(10)%radii(2)'       : None,                    \
                    'patch_icpp(10)%radii(3)'       : None,                    \
                    'patch_icpp(10)%epsilon'        : None,                    \
                    'patch_icpp(10)%beta'           : None,                    \
                    'patch_icpp(10)%normal'         : None,                    \
                    'patch_icpp(10)%normal(1)'      : None,                    \
                    'patch_icpp(10)%normal(2)'      : None,                    \
                    'patch_icpp(10)%normal(3)'      : None,                    \
                    'patch_icpp(10)%alter_patch'    : None,                    \
                    'patch_icpp(10)%alter_patch(1)' : None,                    \
                    'patch_icpp(10)%alter_patch(2)' : None,                    \
                    'patch_icpp(10)%alter_patch(3)' : None,                    \
                    'patch_icpp(10)%alter_patch(4)' : None,                    \
                    'patch_icpp(10)%alter_patch(5)' : None,                    \
                    'patch_icpp(10)%alter_patch(6)' : None,                    \
                    'patch_icpp(10)%alter_patch(7)' : None,                    \
                    'patch_icpp(10)%alter_patch(8)' : None,                    \
                    'patch_icpp(10)%alter_patch(9)' : None,                    \
                    'patch_icpp(10)%smoothen'       : None,                    \
                    'patch_icpp(10)%smooth_patch_id': None,                    \
                    'patch_icpp(10)%smooth_coeff'   : None,                    \
                    'patch_icpp(10)%alpha_rho'      : None,                    \
                    'patch_icpp(10)%alpha_rho(1)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(2)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(3)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(4)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(5)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(6)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(7)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(8)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(9)'   : None,                    \
                    'patch_icpp(10)%alpha_rho(10)'  : None,                    \
                    'patch_icpp(10)%rho'            : None,                    \
                    'patch_icpp(10)%vel'            : None,                    \
                    'patch_icpp(10)%vel(1)'         : None,                    \
                    'patch_icpp(10)%vel(2)'         : None,                    \
                    'patch_icpp(10)%vel(3)'         : None,                    \
                    'patch_icpp(10)%pres'           : None,                    \
                    'patch_icpp(10)%alpha'          : None,                    \
                    'patch_icpp(10)%alpha(1)'       : None,                    \
                    'patch_icpp(10)%alpha(2)'       : None,                    \
                    'patch_icpp(10)%alpha(3)'       : None,                    \
                    'patch_icpp(10)%alpha(4)'       : None,                    \
                    'patch_icpp(10)%alpha(5)'       : None,                    \
                    'patch_icpp(10)%alpha(6)'       : None,                    \
                    'patch_icpp(10)%alpha(7)'       : None,                    \
                    'patch_icpp(10)%alpha(8)'       : None,                    \
                    'patch_icpp(10)%alpha(9)'       : None,                    \
                    'patch_icpp(10)%alpha(10)'      : None,                    \
                    'patch_icpp(10)%gamma'          : None,                    \
                    'patch_icpp(10)%pi_inf'         : None,                    \
                    'patch_icpp(10)%tau_e(1)'       : None,                    \
                    'patch_icpp(10)%tau_e(2)'       : None,                    \
                    'patch_icpp(10)%tau_e(3)'       : None,                    \
                    'patch_icpp(10)%tau_e(4)'       : None,                    \
                    'patch_icpp(10)%tau_e(5)'       : None,                    \
                    'patch_icpp(10)%tau_e(6)'       : None,                    \
                    'fluid_pp(1)%gamma'             : None,                    \
                    'fluid_pp(1)%pi_inf'            : None,                    \
                    'fluid_pp(2)%gamma'             : None,                    \
                    'fluid_pp(2)%pi_inf'            : None,                    \
                    'fluid_pp(3)%gamma'             : None,                    \
                    'fluid_pp(3)%pi_inf'            : None,                    \
                    'fluid_pp(4)%gamma'             : None,                    \
                    'fluid_pp(4)%pi_inf'            : None,                    \
                    'fluid_pp(5)%gamma'             : None,                    \
                    'fluid_pp(5)%pi_inf'            : None,                    \
                    'fluid_pp(6)%gamma'             : None,                    \
                    'fluid_pp(6)%pi_inf'            : None,                    \
                    'fluid_pp(7)%gamma'             : None,                    \
                    'fluid_pp(7)%pi_inf'            : None,                    \
                    'fluid_pp(8)%gamma'             : None,                    \
                    'fluid_pp(8)%pi_inf'            : None,                    \
                    'fluid_pp(9)%gamma'             : None,                    \
                    'fluid_pp(9)%pi_inf'            : None,                    \
                    'fluid_pp(10)%gamma'            : None,                    \
                    'fluid_pp(10)%pi_inf'           : None,                     \
                    'fluid_pp(1)%mul0'              : None,                    \
                    'fluid_pp(1)%ss'                : None,                     \
                    'fluid_pp(1)%pv'                : None,                     \
                    'fluid_pp(1)%gamma_v'           : None,                     \
                    'fluid_pp(1)%M_v'               : None,                     \
                    'fluid_pp(1)%mu_v'              : None,                     \
                    'fluid_pp(1)%k_v'               : None,                     \
                    'fluid_pp(2)%mul0'              : None,                    \
                    'fluid_pp(2)%ss'                : None,                     \
                    'fluid_pp(2)%pv'                : None,                     \
                    'fluid_pp(2)%gamma_v'           : None,                     \
                    'fluid_pp(2)%M_v'               : None,                     \
                    'fluid_pp(2)%mu_v'              : None,                     \
                    'fluid_pp(2)%k_v'               : None,                     \
                    'fluid_pp(3)%mul0'              : None,                    \
                    'fluid_pp(3)%ss'                : None,                     \
                    'fluid_pp(3)%pv'                : None,                     \
                    'fluid_pp(3)%gamma_v'           : None,                     \
                    'fluid_pp(3)%M_v'               : None,                     \
                    'fluid_pp(3)%mu_v'              : None,                     \
                    'fluid_pp(3)%k_v'               : None,                     \
                    'fluid_pp(4)%mul0'              : None,                    \
                    'fluid_pp(4)%ss'                : None,                     \
                    'fluid_pp(4)%pv'                : None,                     \
                    'fluid_pp(4)%gamma_v'           : None,                     \
                    'fluid_pp(4)%M_v'               : None,                     \
                    'fluid_pp(4)%mu_v'              : None,                     \
                    'fluid_pp(4)%k_v'               : None,                     \
                    'fluid_pp(5)%mul0'              : None,                    \
                    'fluid_pp(5)%ss'                : None,                     \
                    'fluid_pp(5)%pv'                : None,                     \
                    'fluid_pp(5)%gamma_v'           : None,                     \
                    'fluid_pp(5)%M_v'               : None,                     \
                    'fluid_pp(5)%mu_v'              : None,                     \
                    'fluid_pp(5)%k_v'               : None,                     \
                    'fluid_pp(6)%mul0'              : None,                    \
                    'fluid_pp(6)%ss'                : None,                     \
                    'fluid_pp(6)%pv'                : None,                     \
                    'fluid_pp(6)%gamma_v'           : None,                     \
                    'fluid_pp(6)%M_v'               : None,                     \
                    'fluid_pp(6)%mu_v'              : None,                     \
                    'fluid_pp(6)%k_v'               : None,                     \
                    'fluid_pp(7)%mul0'              : None,                    \
                    'fluid_pp(7)%ss'                : None,                     \
                    'fluid_pp(7)%pv'                : None,                     \
                    'fluid_pp(7)%gamma_v'           : None,                     \
                    'fluid_pp(7)%M_v'               : None,                     \
                    'fluid_pp(7)%mu_v'              : None,                     \
                    'fluid_pp(7)%k_v'               : None,                     \
                    'fluid_pp(8)%mul0'              : None,                    \
                    'fluid_pp(8)%ss'                : None,                     \
                    'fluid_pp(8)%pv'                : None,                     \
                    'fluid_pp(8)%gamma_v'           : None,                     \
                    'fluid_pp(8)%M_v'               : None,                     \
                    'fluid_pp(8)%mu_v'              : None,                     \
                    'fluid_pp(8)%k_v'               : None,                     \
                    'fluid_pp(9)%mul0'              : None,                    \
                    'fluid_pp(9)%ss'                : None,                     \
                    'fluid_pp(9)%pv'                : None,                     \
                    'fluid_pp(9)%gamma_v'           : None,                     \
                    'fluid_pp(9)%M_v'               : None,                     \
                    'fluid_pp(9)%mu_v'              : None,                     \
                    'fluid_pp(9)%k_v'               : None,                     \
                    'fluid_pp(10)%mul0'             : None,                    \
                    'fluid_pp(10)%ss'               : None,                     \
                    'fluid_pp(10)%pv'               : None,                     \
                    'fluid_pp(10)%gamma_v'          : None,                     \
                    'fluid_pp(10)%M_v'              : None,                     \
                    'fluid_pp(10)%mu_v'             : None,                     \
                    'fluid_pp(10)%k_v'              : None,                     \
                    'fluid_pp(1)%G'                 : None,                     \
                    'fluid_pp(2)%G'                 : None,                     \
                    'fluid_pp(3)%G'                 : None,                     \
                    'fluid_pp(4)%G'                 : None,                     \
                    'fluid_pp(5)%G'                 : None,                     \
                    'fluid_pp(6)%G'                 : None,                     \
                    'fluid_pp(7)%G'                 : None,                     \
                    'fluid_pp(8)%G'                 : None,                     \
                    'fluid_pp(9)%G'                 : None,                     \
                    'fluid_pp(10)%G'                : None,                     \
                    'Ca'                            : None,                     \
                    'Web'                           : None,                     \
                    'Re_inv'                        : None,                     \
                    'pref'                          : None,                     \
                    'rhoref'                        : None,                     \
                    'bubbles'                       : None,                     \
                    'polytropic'                    : None,                     \
                    'polydisperse'                  : None,                     \
                    'poly_sigma'                    : None,                     \
                    'thermal'                       : None,                     \
                    'nb'                            : None,                     \
                    'R0ref'                         : None,                     \
                    'qbmm'                          : None,                     \
                    'dist_type'                     : None,                     \
                    'R0_type'                       : None,                     \
                    'nnode'                         : None,                     \
                    'sigR'                          : None,                     \
                    'sigV'                          : None,                     \
                    'rhoRV'                         : None,                     \
                    'patch_icpp(1)%r0'              : None,                     \
                    'patch_icpp(2)%r0'              : None,                     \
                    'patch_icpp(3)%r0'              : None,                     \
                    'patch_icpp(4)%r0'              : None,                     \
                    'patch_icpp(5)%r0'              : None,                     \
                    'patch_icpp(6)%r0'              : None,                     \
                    'patch_icpp(7)%r0'              : None,                     \
                    'patch_icpp(8)%r0'              : None,                     \
                    'patch_icpp(9)%r0'              : None,                     \
                    'patch_icpp(10)%r0'             : None,                     \
                    'patch_icpp(1)%v0'              : None,                     \
                    'patch_icpp(2)%v0'              : None,                     \
                    'patch_icpp(3)%v0'              : None,                     \
                    'patch_icpp(4)%v0'              : None,                     \
                    'patch_icpp(5)%v0'              : None,                     \
                    'patch_icpp(6)%v0'              : None,                     \
                    'patch_icpp(7)%v0'              : None,                     \
                    'patch_icpp(8)%v0'              : None,                     \
                    'patch_icpp(9)%v0'              : None,                     \
                    'patch_icpp(10)%v0'             : None,                      \
                    'patch_icpp(1)%p0'              : None,                     \
                    'patch_icpp(2)%p0'              : None,                     \
                    'patch_icpp(3)%p0'              : None,                     \
                    'patch_icpp(4)%p0'              : None,                     \
                    'patch_icpp(5)%p0'              : None,                     \
                    'patch_icpp(6)%p0'              : None,                     \
                    'patch_icpp(7)%p0'              : None,                     \
                    'patch_icpp(8)%p0'              : None,                     \
                    'patch_icpp(9)%p0'              : None,                     \
                    'patch_icpp(10)%p0'             : None,                      \
                    'patch_icpp(1)%m0'              : None,                     \
                    'patch_icpp(2)%m0'              : None,                     \
                    'patch_icpp(3)%m0'              : None,                     \
                    'patch_icpp(4)%m0'              : None,                     \
                    'patch_icpp(5)%m0'              : None,                     \
                    'patch_icpp(6)%m0'              : None,                     \
                    'patch_icpp(7)%m0'              : None,                     \
                    'patch_icpp(8)%m0'              : None,                     \
                    'patch_icpp(9)%m0'              : None,                     \
                    'patch_icpp(10)%m0'             : None                      \
}
# ==============================================================================


# Simulation Dictionary ========================================================
simulation_dict =                                                              \
    {                                                                          \
                    'case_dir'                      : None,                    \
                    'run_time_info'                 : None,                    \
                    't_step_old'                    : None,                    \
                    't_tol'                         : None,                    \
                    'debug'                         : None,                    \
                    'm'                             : None,                    \
                    'n'                             : None,                    \
                    'p'                             : None,                    \
                    'cyl_coord'                     : None,                    \
                    'dt'                            : None,                    \
                    't_step_start'                  : None,                    \
                    't_step_stop'                   : None,                    \
                    't_step_save'                   : None,                    \
                    'model_eqns'                    : None,                    \
                    'num_fluids'                    : None,                    \
                    'adv_alphan'                    : None,                    \
                    'mpp_lim'                       : None,                    \
                    'time_stepper'                  : None,                    \
                    'weno_vars'                     : None,                    \
                    'weno_order'                    : None,                    \
                    'weno_eps'                      : None,                    \
                    'char_decomp'                   : None,                    \
                    'mapped_weno'                   : None,                    \
                    'mp_weno'                       : None,                    \
                    'weno_avg'                      : None,                    \
                    'weno_Re_flux'                  : None,                    \
                    'riemann_solver'                : None,                    \
                    'wave_speeds'                   : None,                    \
                    'avg_state'                     : None,                    \
                    'commute_err'                   : None,                    \
                    'split_err'                     : None,                    \
                    'alt_crv'                       : None,                    \
                    'alt_soundspeed'                : None,                    \
                    'regularization'                : None,                    \
                    'reg_eps'                       : None,                    \
                    'null_weights'                  : None,                    \
                    'mixture_err'                   : None,                    \
                    'tvd_riemann_flux'              : None,                    \
                    'tvd_rhs_flux'                  : None,                    \
                    'tvd_wave_speeds'               : None,                    \
                    'flux_lim'                      : None,                    \
                    'We_riemann_flux'               : None,                    \
                    'We_rhs_flux'                   : None,                    \
                    'We_src'                        : None,                    \
                    'We_wave_speeds'                : None,                    \
                    'lsq_deriv'                     : None,                    \
                    'parallel_io'                   : None,                    \
                    'precision'                     : None,                    \
                    'bc_x%beg'                      : None,                    \
                    'bc_x%end'                      : None,                    \
                    'bc_y%beg'                      : None,                    \
                    'bc_y%end'                      : None,                    \
                    'bc_z%beg'                      : None,                    \
                    'bc_z%end'                      : None,                    \
                    'hypoelasticity'                : None,                    \
                    'fd_order'                      : None,                    \
                    'com_wrt'                       : None,                    \
                    'com_wrt(1)'                    : None,                    \
                    'com_wrt(2)'                    : None,                    \
                    'com_wrt(3)'                    : None,                    \
                    'com_wrt(4)'                    : None,                    \
                    'com_wrt(5)'                    : None,                    \
                    'com_wrt(6)'                    : None,                    \
                    'com_wrt(7)'                    : None,                    \
                    'com_wrt(8)'                    : None,                    \
                    'com_wrt(9)'                    : None,                    \
                    'com_wrt(10)'                   : None,                    \
                    'cb_wrt'                        : None,                    \
                    'cb_wrt(1)'                     : None,                    \
                    'cb_wrt(2)'                     : None,                    \
                    'cb_wrt(3)'                     : None,                    \
                    'cb_wrt(4)'                     : None,                    \
                    'cb_wrt(5)'                     : None,                    \
                    'cb_wrt(6)'                     : None,                    \
                    'cb_wrt(7)'                     : None,                    \
                    'cb_wrt(8)'                     : None,                    \
                    'cb_wrt(9)'                     : None,                    \
                    'cb_wrt(10)'                    : None,                    \
                    'num_probes'                    : None,                    \
                    'probe_wrt'                     : None,                    \
                    'probe(1)%x'                    : None,                    \
                    'probe(1)%y'                    : None,                    \
                    'probe(1)%z'                    : None,                    \
                    'probe(2)%x'                    : None,                    \
                    'probe(2)%y'                    : None,                    \
                    'probe(2)%z'                    : None,                    \
                    'probe(3)%x'                    : None,                    \
                    'probe(3)%y'                    : None,                    \
                    'probe(3)%z'                    : None,                    \
                    'probe_wrt(1)%x'                : None,                    \
                    'probe_wrt(1)%y'                : None,                    \
                    'probe_wrt(1)%z'                : None,                    \
                    'probe_wrt(2)%x'                : None,                    \
                    'probe_wrt(2)%y'                : None,                    \
                    'probe_wrt(2)%z'                : None,                    \
                    'probe_wrt(3)%x'                : None,                    \
                    'probe_wrt(3)%y'                : None,                    \
                    'probe_wrt(3)%z'                : None,                    \
                    'probe_wrt(4)%x'                : None,                    \
                    'probe_wrt(4)%y'                : None,                    \
                    'probe_wrt(4)%z'                : None,                    \
                    'probe_wrt(5)%x'                : None,                    \
                    'probe_wrt(5)%y'                : None,                    \
                    'probe_wrt(5)%z'                : None,                    \
                    'probe_wrt(6)%x'                : None,                    \
                    'probe_wrt(6)%y'                : None,                    \
                    'probe_wrt(6)%z'                : None,                    \
                    'probe_wrt(7)%x'                : None,                    \
                    'probe_wrt(7)%y'                : None,                    \
                    'probe_wrt(7)%z'                : None,                    \
                    'probe_wrt(8)%x'                : None,                    \
                    'probe_wrt(8)%y'                : None,                    \
                    'probe_wrt(8)%z'                : None,                    \
                    'probe_wrt(9)%x'                : None,                    \
                    'probe_wrt(9)%y'                : None,                    \
                    'probe_wrt(9)%z'                : None,                    \
                    'probe_wrt(10)%x'               : None,                    \
                    'probe_wrt(10)%y'               : None,                    \
                    'probe_wrt(10)%z'               : None,                    \
                    'threshold_mf'                  : None,                    \
                    'threshold_mf(1)'               : None,                    \
                    'threshold_mf(2)'               : None,                    \
                    'threshold_mf(3)'               : None,                    \
                    'threshold_mf(4)'               : None,                    \
                    'threshold_mf(5)'               : None,                    \
                    'moment_order'                  : None,                    \
                    'moment_order(1)'               : None,                    \
                    'moment_order(2)'               : None,                    \
                    'moment_order(3)'               : None,                    \
                    'moment_order(4)'               : None,                    \
                    'moment_order(5)'               : None,                    \
                    'fluid_pp(1)%gamma'             : None,                    \
                    'fluid_pp(1)%pi_inf'            : None,                    \
                    'fluid_pp(1)%Re(1)'             : None,                    \
                    'fluid_pp(1)%Re(2)'             : None,                    \
                    'fluid_pp(1)%We(2)'             : None,                    \
                    'fluid_pp(1)%We(3)'             : None,                    \
                    'fluid_pp(1)%We(4)'             : None,                    \
                    'fluid_pp(1)%We(5)'             : None,                    \
                    'fluid_pp(1)%We(6)'             : None,                    \
                    'fluid_pp(1)%We(7)'             : None,                    \
                    'fluid_pp(1)%We(8)'             : None,                    \
                    'fluid_pp(1)%We(9)'             : None,                    \
                    'fluid_pp(1)%We(10)'            : None,                    \
                    'fluid_pp(2)%gamma'             : None,                    \
                    'fluid_pp(2)%pi_inf'            : None,                    \
                    'fluid_pp(2)%Re(1)'             : None,                    \
                    'fluid_pp(2)%Re(2)'             : None,                    \
                    'fluid_pp(2)%We(1)'             : None,                    \
                    'fluid_pp(2)%We(3)'             : None,                    \
                    'fluid_pp(2)%We(4)'             : None,                    \
                    'fluid_pp(2)%We(5)'             : None,                    \
                    'fluid_pp(2)%We(6)'             : None,                    \
                    'fluid_pp(2)%We(7)'             : None,                    \
                    'fluid_pp(2)%We(8)'             : None,                    \
                    'fluid_pp(2)%We(9)'             : None,                    \
                    'fluid_pp(2)%We(10)'            : None,                    \
                    'fluid_pp(3)%gamma'             : None,                    \
                    'fluid_pp(3)%pi_inf'            : None,                    \
                    'fluid_pp(3)%Re(1)'             : None,                    \
                    'fluid_pp(3)%Re(2)'             : None,                    \
                    'fluid_pp(3)%We(1)'             : None,                    \
                    'fluid_pp(3)%We(2)'             : None,                    \
                    'fluid_pp(3)%We(4)'             : None,                    \
                    'fluid_pp(3)%We(5)'             : None,                    \
                    'fluid_pp(3)%We(6)'             : None,                    \
                    'fluid_pp(3)%We(7)'             : None,                    \
                    'fluid_pp(3)%We(8)'             : None,                    \
                    'fluid_pp(3)%We(9)'             : None,                    \
                    'fluid_pp(3)%We(10)'            : None,                    \
                    'fluid_pp(4)%gamma'             : None,                    \
                    'fluid_pp(4)%pi_inf'            : None,                    \
                    'fluid_pp(4)%Re(1)'             : None,                    \
                    'fluid_pp(4)%Re(2)'             : None,                    \
                    'fluid_pp(4)%We(1)'             : None,                    \
                    'fluid_pp(4)%We(2)'             : None,                    \
                    'fluid_pp(4)%We(3)'             : None,                    \
                    'fluid_pp(4)%We(5)'             : None,                    \
                    'fluid_pp(4)%We(6)'             : None,                    \
                    'fluid_pp(4)%We(7)'             : None,                    \
                    'fluid_pp(4)%We(8)'             : None,                    \
                    'fluid_pp(4)%We(9)'             : None,                    \
                    'fluid_pp(4)%We(10)'            : None,                    \
                    'fluid_pp(5)%gamma'             : None,                    \
                    'fluid_pp(5)%pi_inf'            : None,                    \
                    'fluid_pp(5)%Re(1)'             : None,                    \
                    'fluid_pp(5)%Re(2)'             : None,                    \
                    'fluid_pp(5)%We(1)'             : None,                    \
                    'fluid_pp(5)%We(2)'             : None,                    \
                    'fluid_pp(5)%We(3)'             : None,                    \
                    'fluid_pp(5)%We(4)'             : None,                    \
                    'fluid_pp(5)%We(6)'             : None,                    \
                    'fluid_pp(5)%We(7)'             : None,                    \
                    'fluid_pp(5)%We(8)'             : None,                    \
                    'fluid_pp(5)%We(9)'             : None,                    \
                    'fluid_pp(5)%We(10)'            : None,                    \
                    'fluid_pp(6)%gamma'             : None,                    \
                    'fluid_pp(6)%pi_inf'            : None,                    \
                    'fluid_pp(6)%Re(1)'             : None,                    \
                    'fluid_pp(6)%Re(2)'             : None,                    \
                    'fluid_pp(6)%We(1)'             : None,                    \
                    'fluid_pp(6)%We(2)'             : None,                    \
                    'fluid_pp(6)%We(3)'             : None,                    \
                    'fluid_pp(6)%We(4)'             : None,                    \
                    'fluid_pp(6)%We(5)'             : None,                    \
                    'fluid_pp(6)%We(7)'             : None,                    \
                    'fluid_pp(6)%We(8)'             : None,                    \
                    'fluid_pp(6)%We(9)'             : None,                    \
                    'fluid_pp(6)%We(10)'            : None,                    \
                    'fluid_pp(7)%gamma'             : None,                    \
                    'fluid_pp(7)%pi_inf'            : None,                    \
                    'fluid_pp(7)%Re(1)'             : None,                    \
                    'fluid_pp(7)%Re(2)'             : None,                    \
                    'fluid_pp(7)%We(1)'             : None,                    \
                    'fluid_pp(7)%We(2)'             : None,                    \
                    'fluid_pp(7)%We(3)'             : None,                    \
                    'fluid_pp(7)%We(4)'             : None,                    \
                    'fluid_pp(7)%We(5)'             : None,                    \
                    'fluid_pp(7)%We(6)'             : None,                    \
                    'fluid_pp(7)%We(8)'             : None,                    \
                    'fluid_pp(7)%We(9)'             : None,                    \
                    'fluid_pp(7)%We(10)'            : None,                    \
                    'fluid_pp(8)%gamma'             : None,                    \
                    'fluid_pp(8)%pi_inf'            : None,                    \
                    'fluid_pp(8)%Re(1)'             : None,                    \
                    'fluid_pp(8)%Re(2)'             : None,                    \
                    'fluid_pp(8)%We(1)'             : None,                    \
                    'fluid_pp(8)%We(2)'             : None,                    \
                    'fluid_pp(8)%We(3)'             : None,                    \
                    'fluid_pp(8)%We(4)'             : None,                    \
                    'fluid_pp(8)%We(5)'             : None,                    \
                    'fluid_pp(8)%We(6)'             : None,                    \
                    'fluid_pp(8)%We(7)'             : None,                    \
                    'fluid_pp(8)%We(9)'             : None,                    \
                    'fluid_pp(8)%We(10)'            : None,                    \
                    'fluid_pp(9)%gamma'             : None,                    \
                    'fluid_pp(9)%pi_inf'            : None,                    \
                    'fluid_pp(9)%Re(1)'             : None,                    \
                    'fluid_pp(9)%Re(2)'             : None,                    \
                    'fluid_pp(9)%We(1)'             : None,                    \
                    'fluid_pp(9)%We(2)'             : None,                    \
                    'fluid_pp(9)%We(3)'             : None,                    \
                    'fluid_pp(9)%We(4)'             : None,                    \
                    'fluid_pp(9)%We(5)'             : None,                    \
                    'fluid_pp(9)%We(6)'             : None,                    \
                    'fluid_pp(9)%We(7)'             : None,                    \
                    'fluid_pp(9)%We(8)'             : None,                    \
                    'fluid_pp(9)%We(10)'            : None,                    \
                    'fluid_pp(10)%gamma'            : None,                    \
                    'fluid_pp(10)%pi_inf'           : None,                    \
                    'fluid_pp(10)%Re(1)'            : None,                    \
                    'fluid_pp(10)%Re(2)'            : None,                    \
                    'fluid_pp(10)%We(1)'            : None,                     \
                    'fluid_pp(10)%We(2)'            : None,                     \
                    'fluid_pp(10)%We(3)'            : None,                     \
                    'fluid_pp(10)%We(4)'            : None,                     \
                    'fluid_pp(10)%We(5)'            : None,                     \
                    'fluid_pp(10)%We(6)'            : None,                     \
                    'fluid_pp(10)%We(7)'            : None,                     \
                    'fluid_pp(10)%We(8)'            : None,                     \
                    'fluid_pp(10)%We(9)'            : None,                     \
                    'fluid_pp(1)%mul0'              : None,                    \
                    'fluid_pp(1)%ss'                : None,                     \
                    'fluid_pp(1)%pv'                : None,                     \
                    'fluid_pp(1)%gamma_v'           : None,                     \
                    'fluid_pp(1)%M_v'               : None,                     \
                    'fluid_pp(1)%mu_v'              : None,                     \
                    'fluid_pp(1)%k_v'               : None,                     \
                    'fluid_pp(2)%mul0'              : None,                    \
                    'fluid_pp(2)%ss'                : None,                     \
                    'fluid_pp(2)%pv'                : None,                     \
                    'fluid_pp(2)%gamma_v'           : None,                     \
                    'fluid_pp(2)%M_v'               : None,                     \
                    'fluid_pp(2)%mu_v'              : None,                     \
                    'fluid_pp(2)%k_v'               : None,                     \
                    'fluid_pp(3)%mul0'              : None,                    \
                    'fluid_pp(3)%ss'                : None,                     \
                    'fluid_pp(3)%pv'                : None,                     \
                    'fluid_pp(3)%gamma_v'           : None,                     \
                    'fluid_pp(3)%M_v'               : None,                     \
                    'fluid_pp(3)%mu_v'              : None,                     \
                    'fluid_pp(3)%k_v'               : None,                     \
                    'fluid_pp(4)%mul0'              : None,                    \
                    'fluid_pp(4)%ss'                : None,                     \
                    'fluid_pp(4)%pv'                : None,                     \
                    'fluid_pp(4)%gamma_v'           : None,                     \
                    'fluid_pp(4)%M_v'               : None,                     \
                    'fluid_pp(4)%mu_v'              : None,                     \
                    'fluid_pp(4)%k_v'               : None,                     \
                    'fluid_pp(5)%mul0'              : None,                    \
                    'fluid_pp(5)%ss'                : None,                     \
                    'fluid_pp(5)%pv'                : None,                     \
                    'fluid_pp(5)%gamma_v'           : None,                     \
                    'fluid_pp(5)%M_v'               : None,                     \
                    'fluid_pp(5)%mu_v'              : None,                     \
                    'fluid_pp(5)%k_v'               : None,                     \
                    'fluid_pp(6)%mul0'              : None,                    \
                    'fluid_pp(6)%ss'                : None,                     \
                    'fluid_pp(6)%pv'                : None,                     \
                    'fluid_pp(6)%gamma_v'           : None,                     \
                    'fluid_pp(6)%M_v'               : None,                     \
                    'fluid_pp(6)%mu_v'              : None,                     \
                    'fluid_pp(6)%k_v'               : None,                     \
                    'fluid_pp(7)%mul0'              : None,                    \
                    'fluid_pp(7)%ss'                : None,                     \
                    'fluid_pp(7)%pv'                : None,                     \
                    'fluid_pp(7)%gamma_v'           : None,                     \
                    'fluid_pp(7)%M_v'               : None,                     \
                    'fluid_pp(7)%mu_v'              : None,                     \
                    'fluid_pp(7)%k_v'               : None,                     \
                    'fluid_pp(8)%mul0'              : None,                    \
                    'fluid_pp(8)%ss'                : None,                     \
                    'fluid_pp(8)%pv'                : None,                     \
                    'fluid_pp(8)%gamma_v'           : None,                     \
                    'fluid_pp(8)%M_v'               : None,                     \
                    'fluid_pp(8)%mu_v'              : None,                     \
                    'fluid_pp(8)%k_v'               : None,                     \
                    'fluid_pp(9)%mul0'              : None,                    \
                    'fluid_pp(9)%ss'                : None,                     \
                    'fluid_pp(9)%pv'                : None,                     \
                    'fluid_pp(9)%gamma_v'           : None,                     \
                    'fluid_pp(9)%M_v'               : None,                     \
                    'fluid_pp(9)%mu_v'              : None,                     \
                    'fluid_pp(9)%k_v'               : None,                     \
                    'fluid_pp(10)%mul0'              : None,                    \
                    'fluid_pp(10)%ss'                : None,                     \
                    'fluid_pp(10)%pv'                : None,                     \
                    'fluid_pp(10)%gamma_v'           : None,                     \
                    'fluid_pp(10)%M_v'               : None,                     \
                    'fluid_pp(10)%mu_v'              : None,                     \
                    'fluid_pp(10)%k_v'              : None,                     \
                    'fluid_pp(1)%G'                 : None,                     \
                    'fluid_pp(2)%G'                 : None,                     \
                    'fluid_pp(3)%G'                 : None,                     \
                    'fluid_pp(4)%G'                 : None,                     \
                    'fluid_pp(5)%G'                 : None,                     \
                    'fluid_pp(6)%G'                 : None,                     \
                    'fluid_pp(7)%G'                 : None,                     \
                    'fluid_pp(8)%G'                 : None,                     \
                    'fluid_pp(9)%G'                 : None,                     \
                    'fluid_pp(10)%G'                : None,                     \
                    'pref'                          : None,                     \
                    'rhoref'                        : None,                     \
                    'polydisperse'                  : None,                     \
                    'poly_sigma'                    : None,                     \
                    'bubbles'                       : None,                     \
                    'bubble_model'                  : None,                     \
                    'polytropic'                    : None,                     \
                    'thermal'                       : None,                     \
                    'R0ref'                         : None,                     \
                    'Ca'                            : None,                     \
                    'Web'                           : None,                     \
                    'Re_inv'                        : None,                     \
                    'nb'                            : None,                     \
                    'Monopole'                      : None,                     \
                    'num_mono'                      : None,                     \
                    'qbmm'                          : None,                     \
                    'R0_type'                       : None,                     \
                    'nnode'                         : None,                     \
                    'Mono(1)%loc(1)'                : None,                     \
                    'Mono(1)%loc(2)'                : None,                     \
                    'Mono(1)%loc(3)'                : None,                     \
                    'Mono(1)%mag'                   : None,                     \
                    'Mono(1)%length'                : None,                     \
                    'Mono(1)%dir'                   : None,                     \
                    'Mono(1)%npulse'                : None,                     \
                    'Mono(1)%pulse'                 : None,                     \
                    'Mono(1)%support'               : None,                     \
                    'Mono(2)%loc(1)'                : None,                     \
                    'Mono(2)%loc(2)'                : None,                     \
                    'Mono(2)%loc(3)'                : None,                     \
                    'Mono(2)%mag'                   : None,                     \
                    'Mono(2)%length'                : None,                     \
                    'Mono(2)%dir'                   : None,                     \
                    'Mono(2)%npulse'                : None,                     \
                    'Mono(2)%pulse'                 : None,                     \
                    'Mono(2)%support'               : None,                     \
                    'Mono(3)%loc(1)'                : None,                     \
                    'Mono(3)%loc(2)'                : None,                     \
                    'Mono(3)%loc(3)'                : None,                     \
                    'Mono(3)%mag'                   : None,                     \
                    'Mono(3)%length'                : None,                     \
                    'Mono(3)%dir'                   : None,                     \
                    'Mono(3)%npulse'                : None,                     \
                    'Mono(3)%pulse'                 : None,                     \
                    'Mono(3)%support'               : None,                     \
                    'Mono(4)%loc(1)'                : None,                     \
                    'Mono(4)%loc(2)'                : None,                     \
                    'Mono(4)%loc(3)'                : None,                     \
                    'Mono(4)%mag'                   : None,                     \
                    'Mono(4)%length'                : None,                     \
                    'Mono(4)%dir'                   : None,                     \
                    'Mono(4)%npulse'                : None,                     \
                    'Mono(4)%pulse'                 : None,                     \
                    'Mono(4)%support'               : None,                     \
                    'Mono(1)%delay'                 : None,                     \
                    'Mono(2)%delay'                 : None,                     \
                    'Mono(3)%delay'                 : None,                     \
                    'Mono(4)%delay'                 : None,                     \
                    'integral_wrt'                  : None,                     \
                    'num_integrals'                 : None,                     \
                    'integral(1)%xmin'              : None,                     \
                    'integral(1)%xmax'              : None,                     \
                    'integral(1)%ymin'              : None,                     \
                    'integral(1)%ymax'              : None,                     \
                    'integral(1)%zmin'              : None,                     \
                    'integral(1)%zmax'              : None,                     \
                    'integral(2)%xmin'              : None,                     \
                    'integral(2)%xmax'              : None,                     \
                    'integral(2)%ymin'              : None,                     \
                    'integral(2)%ymax'              : None,                     \
                    'integral(2)%zmin'              : None,                     \
                    'integral(2)%zmax'              : None,                     \
                    'integral(3)%xmin'              : None,                     \
                    'integral(3)%xmax'              : None,                     \
                    'integral(3)%ymin'              : None,                     \
                    'integral(3)%ymax'              : None,                     \
                    'integral(3)%zmin'              : None,                     \
                    'integral(3)%zmax'              : None,                     \
                    'integral(4)%xmin'              : None,                     \
                    'integral(4)%xmax'              : None,                     \
                    'integral(4)%ymin'              : None,                     \
                    'integral(4)%ymax'              : None,                     \
                    'integral(4)%zmin'              : None,                     \
                    'integral(4)%zmax'              : None,                     \
                    'integral(5)%xmin'              : None,                     \
                    'integral(5)%xmax'              : None,                     \
                    'integral(5)%ymin'              : None,                     \
                    'integral(5)%ymax'              : None,                     \
                    'integral(5)%zmin'              : None,                     \
                    'integral(5)%zmax'              : None                      \
}
# ==============================================================================


# Post-process Dictionary ======================================================
post_process_dict =                                                            \
    {                                                                          \
                    'case_dir'                      : None,                    \
                    'cyl_coord'                     : None,                    \
                    'm'                             : None,                    \
                    'n'                             : None,                    \
                    'p'                             : None,                    \
                    't_step_start'                  : None,                    \
                    't_step_stop'                   : None,                    \
                    't_step_save'                   : None,                    \
                    'model_eqns'                    : None,                    \
                    'num_fluids'                    : None,                    \
                    'adv_alphan'                    : None,                    \
                    'mpp_lim'                       : None,                    \
                    'weno_order'                    : None,                    \
                    'alt_soundspeed'                : None,                    \
                    'mixture_err'                   : None,                    \
                    'parallel_io'                   : None,                    \
                    'bc_x%beg'                      : None,                    \
                    'bc_x%end'                      : None,                    \
                    'bc_y%beg'                      : None,                    \
                    'bc_y%end'                      : None,                    \
                    'bc_z%beg'                      : None,                    \
                    'bc_z%end'                      : None,                    \
                    'hypoelasticity'                : None,                    \
                    'fluid_pp(1)%gamma'             : None,                    \
                    'fluid_pp(1)%pi_inf'            : None,                    \
                    'fluid_pp(2)%gamma'             : None,                    \
                    'fluid_pp(2)%pi_inf'            : None,                    \
                    'fluid_pp(3)%gamma'             : None,                    \
                    'fluid_pp(3)%pi_inf'            : None,                    \
                    'fluid_pp(4)%gamma'             : None,                    \
                    'fluid_pp(4)%pi_inf'            : None,                    \
                    'fluid_pp(5)%gamma'             : None,                    \
                    'fluid_pp(5)%pi_inf'            : None,                    \
                    'fluid_pp(6)%gamma'             : None,                    \
                    'fluid_pp(6)%pi_inf'            : None,                    \
                    'fluid_pp(7)%gamma'             : None,                    \
                    'fluid_pp(7)%pi_inf'            : None,                    \
                    'fluid_pp(8)%gamma'             : None,                    \
                    'fluid_pp(8)%pi_inf'            : None,                    \
                    'fluid_pp(9)%gamma'             : None,                    \
                    'fluid_pp(9)%pi_inf'            : None,                    \
                    'fluid_pp(10)%gamma'            : None,                    \
                    'fluid_pp(10)%pi_inf'           : None,                    \
                    'format'                        : None,                    \
                    'precision'                     : None,                    \
                    'coarsen_silo'                  : None,                    \
                    'fourier_decomp'                : None,                    \
                    'fourier_modes%beg'             : None,                    \
                    'fourier_modes%end'             : None,                    \
                    'alpha_rho_wrt'                 : None,                    \
                    'alpha_rho_wrt(1)'              : None,                    \
                    'alpha_rho_wrt(2)'              : None,                    \
                    'alpha_rho_wrt(3)'              : None,                    \
                    'alpha_rho_wrt(4)'              : None,                    \
                    'alpha_rho_wrt(5)'              : None,                    \
                    'alpha_rho_wrt(6)'              : None,                    \
                    'alpha_rho_wrt(7)'              : None,                    \
                    'alpha_rho_wrt(8)'              : None,                    \
                    'alpha_rho_wrt(9)'              : None,                    \
                    'alpha_rho_wrt(10)'             : None,                    \
                    'rho_wrt'                       : None,                    \
                    'mom_wrt'                       : None,                    \
                    'mom_wrt(1)'                    : None,                    \
                    'mom_wrt(2)'                    : None,                    \
                    'mom_wrt(3)'                    : None,                    \
                    'vel_wrt'                       : None,                    \
                    'vel_wrt(1)'                    : None,                    \
                    'vel_wrt(2)'                    : None,                    \
                    'vel_wrt(3)'                    : None,                    \
                    'flux_lim'                      : None,                    \
                    'flux_wrt'                      : None,                    \
                    'flux_wrt(1)'                   : None,                    \
                    'flux_wrt(2)'                   : None,                    \
                    'flux_wrt(3)'                   : None,                    \
                    'E_wrt'                         : None,                    \
                    'pres_wrt'                      : None,                    \
                    'alpha_wrt'                     : None,                    \
                    'alpha_wrt(1)'                  : None,                    \
                    'alpha_wrt(2)'                  : None,                    \
                    'alpha_wrt(3)'                  : None,                    \
                    'alpha_wrt(4)'                  : None,                    \
                    'alpha_wrt(5)'                  : None,                    \
                    'alpha_wrt(6)'                  : None,                    \
                    'alpha_wrt(7)'                  : None,                    \
                    'alpha_wrt(8)'                  : None,                    \
                    'alpha_wrt(9)'                  : None,                    \
                    'alpha_wrt(10)'                 : None,                    \
                    'kappa_wrt'                     : None,                    \
                    'kappa_wrt(1)'                  : None,                    \
                    'kappa_wrt(2)'                  : None,                    \
                    'kappa_wrt(3)'                  : None,                    \
                    'kappa_wrt(4)'                  : None,                    \
                    'kappa_wrt(5)'                  : None,                    \
                    'kappa_wrt(6)'                  : None,                    \
                    'kappa_wrt(7)'                  : None,                    \
                    'kappa_wrt(8)'                  : None,                    \
                    'kappa_wrt(9)'                  : None,                    \
                    'kappa_wrt(10)'                 : None,                    \
                    'gamma_wrt'                     : None,                    \
                    'heat_ratio_wrt'                : None,                    \
                    'pi_inf_wrt'                    : None,                    \
                    'pres_inf_wrt'                  : None,                    \
                    'cons_vars_wrt'                 : None,                    \
                    'prim_vars_wrt'                 : None,                    \
                    'c_wrt'                         : None,                    \
                    'omega_wrt'                     : None,                    \
                    'omega_wrt(1)'                  : None,                    \
                    'omega_wrt(2)'                  : None,                    \
                    'omega_wrt(3)'                  : None,                    \
                    'schlieren_wrt'                 : None,                    \
                    'schlieren_alpha'               : None,                    \
                    'schlieren_alpha(1)'            : None,                    \
                    'schlieren_alpha(2)'            : None,                    \
                    'schlieren_alpha(3)'            : None,                    \
                    'schlieren_alpha(4)'            : None,                    \
                    'schlieren_alpha(5)'            : None,                    \
                    'schlieren_alpha(6)'            : None,                    \
                    'schlieren_alpha(7)'            : None,                    \
                    'schlieren_alpha(8)'            : None,                    \
                    'schlieren_alpha(9)'            : None,                    \
                    'schlieren_alpha(10)'           : None,                    \
                    'fd_order'                      : None,                     \
                    'fluid_pp(1)%ss'                : None,                     \
                    'fluid_pp(1)%pv'                : None,                     \
                    'fluid_pp(1)%gamma_v'           : None,                     \
                    'fluid_pp(1)%M_v'               : None,                     \
                    'fluid_pp(1)%mu_v'              : None,                     \
                    'fluid_pp(1)%k_v'               : None,                     \
                    'fluid_pp(2)%mul0'              : None,                     \
                    'fluid_pp(2)%ss'                : None,                     \
                    'fluid_pp(2)%pv'                : None,                     \
                    'fluid_pp(2)%gamma_v'           : None,                     \
                    'fluid_pp(2)%M_v'               : None,                     \
                    'fluid_pp(2)%mu_v'              : None,                     \
                    'fluid_pp(2)%k_v'               : None,                     \
                    'fluid_pp(3)%mul0'              : None,                     \
                    'fluid_pp(3)%ss'                : None,                     \
                    'fluid_pp(3)%pv'                : None,                     \
                    'fluid_pp(3)%gamma_v'           : None,                     \
                    'fluid_pp(3)%M_v'               : None,                     \
                    'fluid_pp(4)%mu_v'              : None,                     \
                    'fluid_pp(4)%k_v'               : None,                     \
                    'fluid_pp(4)%mul0'              : None,                     \
                    'fluid_pp(4)%ss'                : None,                     \
                    'fluid_pp(4)%pv'                : None,                     \
                    'fluid_pp(4)%gamma_v'           : None,                     \
                    'fluid_pp(4)%M_v'               : None,                     \
                    'fluid_pp(4)%mu_v'              : None,                     \
                    'fluid_pp(4)%k_v'               : None,                     \
                    'fluid_pp(5)%mul0'              : None,                     \
                    'fluid_pp(5)%ss'                : None,                     \
                    'fluid_pp(5)%pv'                : None,                     \
                    'fluid_pp(5)%gamma_v'           : None,                     \
                    'fluid_pp(5)%M_v'               : None,                     \
                    'fluid_pp(5)%mu_v'              : None,                     \
                    'fluid_pp(5)%k_v'               : None,                     \
                    'fluid_pp(6)%mul0'              : None,                     \
                    'fluid_pp(6)%ss'                : None,                     \
                    'fluid_pp(6)%pv'                : None,                     \
                    'fluid_pp(6)%gamma_v'           : None,                     \
                    'fluid_pp(6)%M_v'               : None,                     \
                    'fluid_pp(6)%mu_v'              : None,                     \
                    'fluid_pp(6)%k_v'               : None,                     \
                    'fluid_pp(7)%mul0'              : None,                     \
                    'fluid_pp(7)%ss'                : None,                     \
                    'fluid_pp(7)%pv'                : None,                     \
                    'fluid_pp(7)%gamma_v'           : None,                     \
                    'fluid_pp(7)%M_v'               : None,                     \
                    'fluid_pp(7)%mu_v'              : None,                     \
                    'fluid_pp(7)%k_v'               : None,                     \
                    'fluid_pp(8)%mul0'              : None,                     \
                    'fluid_pp(8)%ss'                : None,                     \
                    'fluid_pp(8)%pv'                : None,                     \
                    'fluid_pp(8)%gamma_v'           : None,                     \
                    'fluid_pp(8)%M_v'               : None,                     \
                    'fluid_pp(8)%mu_v'              : None,                     \
                    'fluid_pp(8)%k_v'               : None,                     \
                    'fluid_pp(9)%mul0'              : None,                     \
                    'fluid_pp(9)%ss'                : None,                     \
                    'fluid_pp(9)%pv'                : None,                     \
                    'fluid_pp(9)%gamma_v'           : None,                     \
                    'fluid_pp(9)%M_v'               : None,                     \
                    'fluid_pp(9)%mu_v'              : None,                     \
                    'fluid_pp(9)%k_v'               : None,                     \
                    'fluid_pp(10)%mul0'             : None,                     \
                    'fluid_pp(10)%ss'               : None,                     \
                    'fluid_pp(10)%pv'               : None,                     \
                    'fluid_pp(10)%gamma_v'          : None,                     \
                    'fluid_pp(10)%M_v'              : None,                     \
                    'fluid_pp(10)%mu_v'             : None,                     \
                    'fluid_pp(10)%k_v'              : None,                     \
                    'fluid_pp(1)%G'                 : None,                     \
                    'fluid_pp(2)%G'                 : None,                     \
                    'fluid_pp(3)%G'                 : None,                     \
                    'fluid_pp(4)%G'                 : None,                     \
                    'fluid_pp(5)%G'                 : None,                     \
                    'fluid_pp(6)%G'                 : None,                     \
                    'fluid_pp(7)%G'                 : None,                     \
                    'fluid_pp(8)%G'                 : None,                     \
                    'fluid_pp(9)%G'                 : None,                     \
                    'fluid_pp(10)%G'                : None,                     \
                    'polydisperse'                  : None,                     \
                    'poly_sigma'                    : None,                     \
                    'polytropic'                    : None,                     \
                    'thermal'                       : None,                     \
                    'pref'                          : None,                     \
                    'Ca'                            : None,                     \
                    'Web'                           : None,                     \
                    'Re_inv'                        : None,                     \
                    'rhoref'                        : None,                     \
                    'bubbles'                       : None,                     \
                    'R0ref'                         : None,                     \
                    'nb'                            : None                      \
    }
# ==============================================================================

# PBS Dictionary ===============================================================
pbs_dict =                                                                     \
    {                                                                          \
                    'queue'                         : None,                    \
                    'nodes'                         : None,                    \
                    'ppn'                           : None,                    \
                    'walltime'                      : None,                    \
                    'mail_list'                     : None                     \
    }
# ==============================================================================



# CONTAINS =====================================================================

def f_execute_mfc_component_SHB(comp_name, case_dict, mfc_dir, engine, sub_name): # ----------
    # Description: The following function receives the name of the MFC component
    #              the user wishes to execute, the case dictionary, the location
    #              of the MFC folder and lastly, the configuration of the engine
    #              with which the component will be executed. Then, given this
    #              information, the function compiles the code for the selected
    #              component and then writes the component's input file. A batch
    #              file may also be generated, given the engine is configured in
    #              parallel. If this is the case, the function then runs the
    #              component's executable by submitting the batch file to PBS,
    #              otherwise, it runs the executable serially, directly from the
    #              command-line.


    # Enabling access to the MFC component and PBS dictionaries
    global pre_process_dict, simulation_dict, post_process_dict, pbs_dict


    # Checking the validity of the configuration of the engine
    if (engine != 'parallel') and (engine != 'serial'):
        print('\n' + comp_name + '>> Unsupported engine configuration. ' \
                                 'Exiting ...' + '\n')
        exit(0)


    # Checking whether the MFC component selected by the user exists
    if (comp_name != 'MFC_PreProcess' ) and \
       (comp_name != 'MFC_Simulation'  ) and \
       (comp_name != 'MFC_PostProcess'):
        print('\n' + 'Unsupported choice of MFC component to execute. ' \
                   + 'Exiting ...' + '\n')
        exit(0)


    # Checking the consistency of the case dictionary with respect to the MFC
    # component and PBS dictionaries
    for parameter in case_dict:
        if ( (parameter in pre_process_dict) == False) and \
           (  (parameter in simulation_dict) == False) and \
           ((parameter in post_process_dict) == False) and \
           (         (parameter in pbs_dict) == False):
               print('\n' + comp_name + '>> Unsupported parameter choice ' \
                          + parameter + '. Exiting ...' + '\n')
               exit(0)


    # Updating the values in the PBS dictionary using the values provided by the
    # user in the case dictionary
    for parameter in case_dict:
        if (parameter in pbs_dict) == True:
            pbs_dict[parameter] = case_dict[parameter]


    # Outputting the component's start-up message
    print('\n' + comp_name + '>> Preparing ' + engine + ' job ...' + '\n')


    # Setting the directory location for the MFC component
    comp_dir = mfc_dir + '/' + comp_name + '_code'

    makefile = 'makefile_richardson'
    if sys.platform == 'darwin':
        makefile = 'makefile'

    # Compiling the MFC component's code if necessary
    cmd_status = Popen('make -C ' + comp_dir + ' all', shell=True, stdout=PIPE)
    output, errors = cmd_status.communicate()


    # Generating input file to be read in by the MFC component's executable
    f_create_input_file(comp_name, case_dict)


     # If the engine is configured serially, the job, i.e. the executable of the
    # component, is run in the command-line, else, for a parallel configuration,
    # a bash script is generated and the job is submitted to a queue via PBS.
    if engine == 'serial':
        print('\n' + comp_name + '>> Serial job in progress ...' + '\n')
        #cmd_status = Popen('mpirun -n '+str(pbs_dict[ 'ppn' ])+' ./'+comp_dir+'/'+comp_name, shell=True, stdout=PIPE)

        cmd_status = Popen(f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{pathlib.Path(__file__).parent.resolve()}/../../.mfc/___current___/build/lib" mpirun -n {str(pbs_dict["ppn"])} "{mfc_dir}/../.mfc/___current___/build/bin/{comp_name}"', shell=True )
        output, errors = cmd_status.communicate()
        print('\n' + output)
        print(comp_name + '>> Serial job completed!' + '\n')
        cmd_status = Popen('rm -f '+ comp_name +'.inp', shell=True, stdout=PIPE)
        output, errors = cmd_status.communicate()
    #else if engine == 'interactive':
    #    print '\n' + comp_name + '>> Interactive job in progress ...' + '\n'
    #    cmd_status = Popen('./'+comp_dir+'/'+comp_name, shell=True, stdout=PIPE)
    #    output, errors = cmd_status.communicate()
    #    print '\n' + output
    #    print comp_name + '>> Serial job completed!' + '\n'
    #    cmd_status = Popen('rm -f '+ comp_name +'.inp', shell=True, stdout=PIPE)
    #    output, errors = cmd_status.communicate()
    else:
        f_create_batch_file_SHB(comp_name, case_dict, mfc_dir, sub_name)
        # Submit job to queue (Hooke/Thomson/Darter)
        #cmd_status = Popen('qsub ' + comp_name + '.sh', shell=True, stdout=PIPE)
        # submit job to queue (Stampede)
        cmd_status = Popen('sbatch ' + comp_name + '.sh', shell=True, stdout=PIPE)
        output, errors = cmd_status.communicate()
        print('\n' + output)
        print(comp_name + '>> Parallel job submitted to queue!' + '\n')
# END: def f_execute_mfc_component ---------------------------------------------

def f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine): # ----------
    # Description: The following function receives the name of the MFC component
    #              the user wishes to execute, the case dictionary, the location
    #              of the MFC folder and lastly, the configuration of the engine
    #              with which the component will be executed. Then, given this
    #              information, the function compiles the code for the selected
    #              component and then writes the component's input file. A batch
    #              file may also be generated, given the engine is configured in
    #              parallel. If this is the case, the function then runs the
    #              component's executable by submitting the batch file to PBS,
    #              otherwise, it runs the executable serially, directly from the
    #              command-line.


    # Enabling access to the MFC component and PBS dictionaries
    global pre_process_dict, simulation_dict, post_process_dict, pbs_dict


    # Checking the validity of the configuration of the engine
    if (engine != 'parallel') and (engine != 'serial'):
        print('\n' + comp_name + '>> Unsupported engine configuration. Exiting ...' + '\n')
        exit(0)


    # Checking whether the MFC component selected by the user exists
    if (comp_name != 'MFC_PreProcess' ) and \
       (comp_name != 'MFC_Simulation'  ) and \
       (comp_name != 'MFC_PostProcess'):
        print( '\n' + 'Unsupported choice of MFC component to execute. ' \
                   + 'Exiting ...' + '\n')
        exit(0)


    # Checking the consistency of the case dictionary with respect to the MFC
    # component and PBS dictionaries
    for parameter in case_dict:
        if  (  parameter not in  simulation_dict ) and \
            (  parameter not in  post_process_dict ) and \
            (  parameter not in  pre_process_dict ) and \
            (  parameter not in  pbs_dict ):
               print( '\n' + comp_name + '>> Unsupported parameter choice ' \
                          + parameter + '. Exiting ...' + '\n')
               exit(0)


    # Updating the values in the PBS dictionary using the values provided by the
    # user in the case dictionary
    for parameter in case_dict:
        if (parameter in pbs_dict):
            pbs_dict[parameter] = case_dict[parameter]


    # Checking whether the engine configuration is compatible with the values of
    # the parameters located in the PBS dictionary
#    if engine == 'serial':
#        for parameter in pbs_dict:
#            if pbs_dict[parameter] is not None:
#                print '\n' + comp_name + '>> Serial engine configuration '  \
#                                       + 'incompatible with value(s) of '   \
#                                       + 'parameter(s) in PBS dictionary. ' \
#                                       + 'Exiting ...' + '\n'
#                exit(0)
#    else:
#        for parameter in pbs_dict:
#            if pbs_dict[parameter] is None:
#                print '\n' + comp_name + '>> Parallel engine configuration ' \
#                                       + 'incompatible with value(s) of '    \
#                                       + 'parameter(s) in PBS dictionary. '  \
#                                       + 'Exiting ...' + '\n'
#                exit(0)


    # Outputting the component's start-up message
    print( '\n' + comp_name + '>> Preparing ' + engine + ' job ...' + '\n')


    # Setting the directory location for the MFC component
    comp_dir = mfc_dir + '/' + comp_name + '_code'

    # Compiling the MFC component's code if necessary
    #cmd_status = Popen(f'./mfc.py --build {comp_dir} all', shell=True, stdout=PIPE)
    #output, errors = cmd_status.communicate()


    # Generating input file to be read in by the MFC component's executable
    f_create_input_file(comp_name, case_dict)


    # If the engine is configured serially, the job, i.e. the executable of the
    # component, is run in the command-line, else, for a parallel configuration,
    # a bash script is generated and the job is submitted to a queue via PBS.
    if engine == 'serial':
        print( '\n' + comp_name + '>> Serial job in progress ...' + '\n')
        #cmd_status = Popen('./'+comp_dir+'/'+comp_name, shell=True, stdout=PIPE)

        cmd_status = Popen(f'LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{pathlib.Path(__file__).parent.resolve()}/../../.mfc/___current___/build/lib" mpirun -n {str(pbs_dict["ppn"])} "{mfc_dir}/../.mfc/___current___/build/bin/{comp_name}"', shell=True )
        output, errors = cmd_status.communicate()
        #print '\n' + output
        print( comp_name + '>> Serial job completed!' + '\n')
        #cmd_status = Popen('rm -f '+ comp_name +'.inp', shell=True, stdout=PIPE)
        #output, errors = cmd_status.communicate()
    else:
        f_create_batch_file(comp_name, case_dict, mfc_dir)
        # Submit job to queue (qsub)
        # cmd_status = Popen('qsub ' + comp_name + '.sh', shell=True, stdout=PIPE)
        # submit job to queue (sbatch)
        cmd_status = Popen('sbatch ' + comp_name + '.sh', shell=True, stdout=PIPE)
        output, errors = cmd_status.communicate()
        print( '\n' + output)
        print( comp_name + '>> Parallel job submitted to queue!' + '\n')
# END: def f_execute_mfc_component ---------------------------------------------


def f_create_input_file(comp_name, case_dict): # -------------------------------
    # Description: The following function generates an input file given the name
    #              of the MFC component for which the file is to be created and
    #              the case dictionary from which the parameters will be used to
    #              populate the file.


    # Enabling access to the MFC component dictionaries
    global pre_process_dict, simulation_dict, post_process_dict


    # Updating the values in the relevant MFC component dictionary using the
    # values provided by the user in the case dictionary
    if comp_name == 'MFC_PreProcess':
        for parameter in case_dict:
            if parameter in pre_process_dict:
                pre_process_dict[parameter] = case_dict[parameter]
        comp_dict = pre_process_dict
    elif comp_name == 'MFC_Simulation':
        for parameter in case_dict:
            if parameter in simulation_dict:
                simulation_dict[parameter] = case_dict[parameter]
        comp_dict = simulation_dict
    else:
        for parameter in case_dict:
            if parameter in post_process_dict:
                post_process_dict[parameter] = case_dict[parameter]
        comp_dict = post_process_dict


    # Setting the location of the input file
    file_loc = comp_name + '.inp'


    # Opening and obtaining a handle for it
    file_id = open(file_loc, 'w')


    # Populating the file's header information
    file_id.write('&user_inputs\n')


    # Populating the body of the input file
    for parameter in comp_dict:
        if comp_dict[parameter] is not None:
            file_id.write(parameter + ' = ' + str(comp_dict[parameter]) + '\n')

    # OTS: end statement for namelist [gfortran compatibility]
    file_id.write('&end\n')

    # Populating the file's footer information
    file_id.write('/')


    # Closing the input file
    file_id.close()
# END: def f_create_input_file -------------------------------------------------


def f_create_batch_file_SHB(comp_name, case_dict, mfc_dir,sub_name): # ----------------------
    # Description: The following function generates a batch file given the name
    #              of the MFC component for which the file is to be created, the
    #              case dictionary from which the parameters will be used to
    #              populate the file, and the location of the MFC folder.


    # Enabling access to the PBS dictionary
    global pbs_dict


    # Setting the location of the batch file
    file_loc = comp_name + '.sh'


    # Opening and obtaining a handle for it
    file_id = open(file_loc, 'w')


    # Populating Batch File  ===================================================
    file_id.write(                                                             \
                                                                               \
        # Script interpreter
        '#!/bin/sh'                                                     + '\n' \
                                                                               \
        # Account to be charged for the job:
        # (Darter/Gordon)
        # '#PBS -A TG-CTS120005'                                          + '\n' \
        # (Stampede)
        #'#SBATCH -A TG-CTS120005'                                       + '\n' \
                                                                               \
        # Name of the queue to which the job should be submitted:
        # (Hooke/Thomson/Darter/Gordon)
        # '#PBS -q ' + str(pbs_dict['queue'])                             + '\n' \
        # (Stampede)
        # '#SBATCH -p ' + str(pbs_dict['queue'])                          + '\n' \
        # (Comet)
        '#SBATCH --partition=' + str(pbs_dict['queue'])                 + '\n' \
                                                                               \
        # Name of the job to be submitted to the scheduler:
        # (Hooke/Thomson/Darter/Gordon)
        # '#PBS -N ' + comp_name                                          + '\n' \
        # (Stampede)
        # '#SBATCH -J ' + comp_name                                       + '\n' \
        # (Comet/Richardson)
        '#SBATCH --job-name=' + sub_name                               + '\n' \
                                                                               \
        # Node(s) and processor(s) per node (ppn) for job:
        # (Thomson)
        # '#PBS -l nodes=' + str(pbs_dict['nodes'])                              \
        #        + ':ppn=' + str(pbs_dict[ 'ppn' ])                       + '\n' \
        # (Hooke)
        # '#PBS -l nodes=0' + str(pbs_dict['nodes'])                             \
        #        + ':ppn=' + str(pbs_dict[ 'ppn' ])                       + '\n' \
        # (Darter)
        # '#PBS -l size=' + str( pbs_dict['nodes']*pbs_dict['ppn']               \
        #                      + min(1,( pbs_dict['nodes']                       \
        #                              * pbs_dict[ 'ppn' ] )%16)                 \
        #                      * (16 - ( pbs_dict['nodes']                       \
        #                              * pbs_dict[ 'ppn' ] )%16) )        + '\n' \
        # (Stampede)
        # '#SBATCH -n ' + str( pbs_dict['nodes']*pbs_dict['ppn']                 \
        #                      + min(1,( pbs_dict['nodes']                       \
        #                              * pbs_dict[ 'ppn' ] )%16)                 \
        #                      * (16 - ( pbs_dict['nodes']                       \
        #                              * pbs_dict[ 'ppn' ] )%16) )        + '\n' \
        # (Gordon)
        # '#PBS -l nodes=' + str(pbs_dict['nodes'])                              \
        #        + ':ppn=' + str(pbs_dict[ 'ppn' ]) + ':native'           + '\n' \
        # (Comet)
        #'#SBATCH --nodes=' + str(pbs_dict['nodes'])                      + '\n' \
        #'#SBATCH --ntasks-per-node=' + str(pbs_dict['ppn'])              + '\n' \
        #'#SBATCH --switches=' + '1'                                      + '\n' \
        #                                                                       \
        # (Richardson)
        #'#SBATCH -n ' + str( pbs_dict['nodes']*pbs_dict['ppn'])                \
                #              + ' --ntasks-per-node=' + str(pbs_dict['ppn'])    + '\n' \
        '#SBATCH -n ' + str( pbs_dict['nodes']*pbs_dict['ppn'])         + '\n' \
                                                                               \
        # (Richardson)
        #'#SBATCH --mem-per-cpu ' + str(5) + 'G'                         + '\n' \
                                                                               \
        # Maximum amount of time to commit to the execution of the job:
        # (Hooke/Thomson/Gordon)
        # '#PBS -l walltime=' + str(pbs_dict['walltime'])                 + '\n' \
        # (Stampede/Comet)
        '#SBATCH -t ' + str(pbs_dict['walltime'])                       + '\n' \
                                                                               \
        # Declare the job rerunable (y) or non-rerunable (n)
        # (Hooke/Thomson)
        # '#PBS -r n'                                                     + '\n' \
                                                                               \
        # Output standard output and error in a single file
        # (Hooke/Thomson/Darter/Gordon?)
        # '#PBS -j oe'                                                    + '\n' \
        # (Stampede/Comet/Richardson)
        '#SBATCH -o ' + comp_name + '.o%j'                              + '\n' \
        '#SBATCH -e ' + comp_name + '.o%j'                              + '\n' \
                                                                               \
        # Notify by email when job begins (b), aborts (a), and/or ends (e):
        # (Hooke/Thomson/Darter/Gordon)
        # '#PBS -m bae'                                                   + '\n' \
        # '#PBS -M ' + str(pbs_dict['mail_list'])                         + '\n' \
        # (Stampede/Comet/Richardson)
        '#SBATCH --mail-type=all'                                       + '\n' \
        '#SBATCH --mail-user=' + str(pbs_dict['mail_list'])             + '\n' \
                                                                               \
        # Total number of processor(s) allocated for job execution
        # (Hooke/Thomson/Darter/Gordon?)
        # 'num_procs=$(cat $PBS_NODEFILE | wc -l)'                        + '\n' \

        # Moving to the case directory
        # (Hooke/Thomson/Darter/Gordon?)
        # 'cd $PBS_O_WORKDIR'                                             + '\n' \
                                                                               \
        # Setting up environment variables for MPI I/O on Cray systems (Darter)
        # 'export MPICH_PTL_UNEX_EVENTS=400000'                           + '\n' \
                                                                               \
        # 'export MPICH_PTL_OTHER_EVENTS=100000'                          + '\n' \
                                                                               \
        # 'export MPICH_UNEX_BUFFER_SIZE=536870912'                       + '\n' \
                                                                               \
        # 'export MPICH_MPIIO_HINTS=*:romio_ds_write=disable'             + '\n' \
                                                                               \
        # Setting up the output file's header information:
        # (Hooke/Thomson/Darter/Gordon?)
        # 'echo MFC v3.0 - Cases - ' + basename(getcwd())                        \
        #                           + ': $PBS_JOBNAME.o${PBS_JOBID:0:7}' + '\n' \
        # 'echo Description: $PBS_JOBID executed on $num_procs '                 \
        # (Stampede/Comet/Richardson)
        'echo MFC v3.0 - Cases - ' + basename(getcwd())                        \
                               + ': $SLURM_JOB_NAME.o%j'                 + '\n' \
        'echo Description: %j executed on $num_procs '                         \

                         + 'processor\'(s)\'. The' + '\n' + 'echo '            \
                         + '\'            \' command-line output '             \
                         + 'information may be found below.'            + '\n' \
        'echo Author: Vedran Coralic'                                   + '\n' \
        'echo Start-date: `date +%D`'                                   + '\n' \
        'echo Start-time: `date +%T`'                                   + '\n' \
        'echo' + '\n' + 'echo'                                          + '\n' \
        'echo \'================================ Terminal Output '             \
             + '===============================\'' + '\n' + 'echo'      + '\n' \
                                                                               \
        # Starting the timer for the job execution
        't_start=$(date +%s)'                                           + '\n' \
        #'t_start=$(date +"%T.%3N")'                                           + '\n' \
        # Executing job:
        # (Hooke)
        # '/opt/mvapich2/ch3_mrail_gen2-intel12/bin/mpirun '                     \
        #                                + mfc_dir + '/' + comp_name             \
        #                                + '_code' + '/' + comp_name      + '\n' \
        # (Darter)
        # 'aprun -n ' + str(pbs_dict['nodes']*pbs_dict['ppn']) + ' '             \
        #                                + mfc_dir + '/' + comp_name             \
        #                                + '_code' + '/' + comp_name      + '\n' \
        # (Thomson)
        # '/share/apps/openmpi-1.4.3/nag_fort/bin/mpirun '                       \
        #                                + mfc_dir + '/' + comp_name             \
        #                                + '_code' + '/' + comp_name      + '\n' \
        # (Stampede/Comet)
        #'ibrun '                                                               \
        #                               + mfc_dir + '/' + comp_name             \
        #                               + '_code' + '/' + comp_name      + '\n' \
        # (Gordon)
        # 'mpirun_rsh -np ' + str(pbs_dict['nodes']*pbs_dict['ppn']) + ' '       \
        #                               + '-hostfile $PBS_NODEFILE '            \
        #                               + mfc_dir + '/' + comp_name             \
        #                               + '_code' + '/' + comp_name      + '\n' \
        # (Richardson)
        'mpirun '                                                               \
                                       + mfc_dir + '/' + comp_name             \
                                       + '_code' + '/' + comp_name      + '\n' \
        # Stopping the timer for the job
        't_stop=$(date +%s)' + '\n' + 'echo'                            + '\n' \
        #'t_stop=$(date +"%T.%3N")' + '\n' + 'echo'                            + '\n' \
                                                                               \
        # Setting up the PBS output file's footer information
        'echo \'================================================='             \
             + '===============================\''                      + '\n' \
        'echo' + '\n' + 'echo'                                          + '\n' \
        'echo End-date: `date +%D`'                                     + '\n' \
        'echo End-time: `date +%T`' + '\n' + 'echo'                     + '\n' \
        'echo Total-time: $(expr $t_stop - $t_start)s'                  + '\n' \
                                                                               \
        # Removing the input file
        #'rm -f ' + comp_name + '.inp'                                   + '\n' \
                                                                               \
        # Removing the batch file
        'rm -f ' + comp_name + '.sh'                                           )
    # END: Populating Batch File ===============================================


    # Closing the batch file
    file_id.close()


    # Giving the batch file the permission to be executed
    cmd_status = Popen('chmod +x ' + comp_name + '.sh', shell=True, stdout=PIPE)
    output, errors = cmd_status.communicate()


def f_create_batch_file(comp_name, case_dict, mfc_dir): # ----------------------
    # Description: The following function generates a batch file given the name
    #              of the MFC component for which the file is to be created, the
    #              case dictionary from which the parameters will be used to
    #              populate the file, and the location of the MFC folder.


    # Enabling access to the PBS dictionary
    global pbs_dict


    # Setting the location of the batch file
    file_loc = comp_name + '.sh'


    # Opening and obtaining a handle for it
    file_id = open(file_loc, 'w')


    # Populating Batch File  ===================================================
    file_id.write(                                                             \
                                                                               \
        # Script interpreter
        '#!/bin/sh'                                                     + '\n' \
                                                                               \
        # Account to be charged for the job:
        # (PBS)
        # '#PBS -A xxx'                                          + '\n' \
        # (Slurm)
        # '#SBATCH -A xxx'                                       + '\n' \
                                                                               \
        # Name of the queue to which the job should be submitted:
        # (PBS)
        # '#PBS -q ' + str(pbs_dict['queue'])                             + '\n' \
        # (Slurm)
        '#SBATCH -p ' + str(pbs_dict['queue'])                          + '\n' \
                                                                               \
        # Name of the job to be submitted to the scheduler:
        # (PBS)
        # '#PBS -N ' + comp_name                                          + '\n' \
        # (Slurm)
        '#SBATCH -J ' + comp_name                                       + '\n' \
                                                                               \
        # Node(s) and processor(s) per node (ppn) for job:
        # (PBS)
        # '#PBS -l nodes=0' + str(pbs_dict['nodes'])                             \
        #        + ':ppn=' + str(pbs_dict[ 'ppn' ])                       + '\n' \
        # (Slurm)
        '#SBATCH --nodes=' + str(pbs_dict['nodes'])                     + '\n' \
        '#SBATCH --ntasks-per-node=' + str(pbs_dict['ppn'])             + '\n' \
                                                                               \
        # Maximum amount of time to commit to the execution of the job:
        # (PBS)
        # '#PBS -l walltime=' + str(pbs_dict['walltime'])                 + '\n' \
        # (Slurm)
        '#SBATCH -t ' + str(pbs_dict['walltime'])                       + '\n' \
                                                                               \
        # Declare the job rerunable (y) or non-rerunable (n)
        # (PBS)
        # '#PBS -r n'                                                     + '\n' \
        #                                                                        \
        # Output standard output and error in a single file
        # (PBS)
        # '#PBS -j oe'                                                    + '\n' \
        # (Slurm)
        '#SBATCH -o ' + comp_name + '.o%j'                              + '\n' \
        '#SBATCH -e ' + comp_name + '.o%j'                              + '\n' \
                                                                               \
        # Notify by email when job begins (b), aborts (a), and/or ends (e):
        # (PBS)
        # '#PBS -m bae'                                                   + '\n' \
        # '#PBS -M ' + str(pbs_dict['mail_list'])                         + '\n' \
        # (Slurm)
        '#SBATCH --mail-type=all'                                       + '\n' \
        '#SBATCH --mail-user=' + str(pbs_dict['mail_list'])             + '\n' \
                                                                               \
        #'sleep 30s'                                                     + '\n' \
        # Total number of processor(s) allocated for job execution
        # (PBS)
        # 'num_procs=$(cat $PBS_NODEFILE | wc -l)'                        + '\n' \
                                                                               \
        # Moving to the case directory
        # (PBS)
        # 'cd $PBS_O_WORKDIR'                                             + '\n' \
                                                                               \
        # Setting up environment variables for MPI I/O on Cray systems
        # 'export MPICH_PTL_UNEX_EVENTS=400000'                           + '\n' \
        #                                                                        \
        # 'export MPICH_PTL_OTHER_EVENTS=100000'                          + '\n' \
        #                                                                        \
        # 'export MPICH_UNEX_BUFFER_SIZE=536870912'                       + '\n' \
        #                                                                        \
        # 'export MPICH_MPIIO_HINTS=*:romio_ds_write=disable'             + '\n' \
        #                                                                        \
        # Setting up the output file's header information:
        # (PBS)
        # 'echo MFC - Cases - ' + basename(getcwd())                        \
        #                            + ': $PBS_JOBNAME.o${PBS_JOBID:0:7}' + '\n' \
        # 'echo Description: $PBS_JOBID executed on $num_procs '                 \
        # (Slurm)
        'echo MFC - Cases - ' + basename(getcwd())                        \
                              + ': $SLURM_JOB_NAME.o$SLURM_JOB_ID'      + '\n' \
        'echo Description: $SLURM_JOB_ID executed on $SLURM_NTASKS '           \

                         + 'processor\'(s)\'. The' + '\n' + 'echo '            \
                         + '\'            \' command-line output '             \
                         + 'information may be found below.'            + '\n' \
        'echo Start-date: `date +%D`'                                   + '\n' \
        'echo Start-time: `date +%T`'                                   + '\n' \
        'echo' + '\n' + 'echo'                                          + '\n' \
        'echo \'================================ Terminal Output '             \
             + '===============================\'' + '\n' + 'echo'      + '\n' \
                                                                               \
        # Starting the timer for the job execution
         't_start=$(date +%s)'                                          + '\n' \
                                                                               \
        # Executing job:
        'mpirun '                                                               \
                                       + mfc_dir + '/' + comp_name             \
                                       + '_code' + '/' + comp_name      + '\n' \
        # Stopping the timer for the job
        't_stop=$(date +%s)' + '\n' + 'echo'                            + '\n' \
                                                                               \
        # Setting up the PBS output file's footer information
        'echo \'================================================='             \
             + '===============================\''                      + '\n' \
        'echo' + '\n' + 'echo'                                          + '\n' \
        'echo End-date: `date +%D`'                                     + '\n' \
        'echo End-time: `date +%T`' + '\n' + 'echo'                     + '\n' \
        'echo Total-time: $(expr $t_stop - $t_start)s'                  + '\n' \
                                                                               \
        # Removing the input file
        'rm -f ' + comp_name + '.inp'                                   + '\n' \
                                                                               \
        # Removing the batch file
        'rm -f ' + comp_name + '.sh'                                           )
    # END: Populating Batch File ===============================================


    # Closing the batch file
    file_id.close()


    # Giving the batch file the permission to be executed
    cmd_status = Popen('chmod +x ' + comp_name + '.sh', shell=True, stdout=PIPE)
    output, errors = cmd_status.communicate()
# END: def f_create_batch_file -------------------------------------------------


# END: CONTAINS ================================================================
