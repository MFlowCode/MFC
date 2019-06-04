include ../../Makefile.in

# List of the object files that are members of the pre-process code
object_files = m_derived_types.o \
                           m_compile_specific.o \
                           m_global_parameters.o \
                           m_mpi_proxy.o \
                           m_variables_conversion.o \
                           m_data_output.o \
                           m_start_up.o \
                           m_grid.o \
                           m_initial_condition.o \
                           p_main.o

# ==============================================================================


# Rules ========================================================================

# Creation of the executable
all:	 $(object_files)
		$(compiler) $(object_files) -o ./pre_process

# Clean-up of the working directory
clean:
		rm -f *.o *.mod *.f90~
		rm -f pre_process

# Treatment of the object files
%.o:	%.f90
		$(compiler) $(f90_flags) $<

# ==============================================================================
