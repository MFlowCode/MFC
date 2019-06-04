include ../../Makefile.in

# List of the object files that constitute the simulation code
object_files = m_derived_types.o \
	       			m_compile_specific.o \
			   m_global_parameters.o \
			   m_mpi_proxy.o \
			   m_variables_conversion.o \
			   m_start_up.o \
			   m_fftw.o \
			   m_weno.o \
			   m_riemann_solvers.o \
			   m_cbc.o \
			   m_bubbles.o \
			   m_rhs.o \
			   m_data_output.o \
			   m_time_steppers.o \
			   m_derived_variables.o \
			   p_main.o

# ==============================================================================


# Rules ========================================================================

# Creation of the executable
all:	$(object_files)
		$(compiler) $(object_files) $(ld_sim_flags) -o ./simulation

# Clean-up of the working directory
clean:
		rm -f *.o *.mod *.f90~
		rm -f simulation

# Treatment of the object files
%.o:	%.f90
		$(compiler) $(f90_flags) $<

# ===============================================================================
