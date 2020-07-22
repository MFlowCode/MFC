# How to add variable to MFC

* Add variable to input.py
* Add variable to src/master_scripts/m_python_proxy.py dictionary as appropriate (pre_process, simulation, and/or post_process)
    * This will then write this variable to the \*.in input file that gets created by input.py
* Go into the code that you added the variable for (pre_process / simulation / post_process)
    * m_startup.f90
      * Add variable to NAMELIST user_inputs in subroutine s_read_input_file
      * Add appropriate checks to s_check_input_file 
    * m_global_parameters.f90
      * Make your new variable/parameter global so that the whole code can see it
      * Add it to the header of this module (document as appropriate)
      * In s_assign_default_values_to_user_inputs set the default value for this variable
    * m_mpi_proxy.f90
      * s_mpi_bcast_user_inputs
      * This is where the variable is broadcasted to all processors 
      * Add to list of broadcasted variables 
      * e.g. CALL MPI_BCAST(weno_nn, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr  )
* You can now use your variable (e.g. weno_nn) in the code. At least anywhere that m_global_parameters is used (see header of module)

