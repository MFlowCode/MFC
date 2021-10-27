#!/usr/bin/python

## @brief The following script describes the steps necessary to configure
##		and interconnect the pre-process, simulation and post-process
##		components of the MFC for the analysis of a sample case. Note
##		that this file is meant to serve as a template and so should
##		not be edited in this location. Suggestions on how to properly
##		customize this script for practical use are described at the end
##		of the file.

# Dependencies and Logistics ===================================================

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path

# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../..'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
engine = 'parallel'
# ==============================================================================


# Case Analysis Configuration ==================================================

# MFC component to be executed, inputted by the user at the command-line. The
# three options include the 'pre_process', 'simulation' and 'post_process'. 
comp_name = argv[1].strip()

# Case-sensitive dictionary needing only the specification of parameters that
# are necessary to execute the user selected MFC component. Note that providing
# parameters that are not used by the chosen component is not a problem as long
# as those parameters are needed in the execution of one of the two remaining
# MFC components or the use of the portable batch system (PBS).
case_dict = {}

# Executing the MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================


# Modification Notes ===========================================================

# In order to setup a case in the cases folder, the user should first create the
# appropriate directory in cases and then proceed to make a copy of the template
# to this new location. The module m_python_proxy.py does NOT need to be copied.
# The new script file may then be edited by the user at will. However, it is
# recommended that the user only modify the case dictionary variable. The latter
# should contain all of the parameters necessary to execute the selected MFC
# component for the new case.

# ==============================================================================
