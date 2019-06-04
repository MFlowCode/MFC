include ../../Makefile.in


# List of the object files that are members of the post-process code
object_files = m_derived_types.o \
	       			m_compile_specific.o \
			   m_global_parameters.o \
			   m_mpi_proxy.o \
               m_fftw.o \
			   m_variables_conversion.o \
			   m_start_up.o \
			   m_data_input.o \
			   m_derived_variables.o \
			   m_data_output.o\
			   p_main.o

# ==============================================================================


# Rules ========================================================================

# Creation of the executable
all:	$(object_files)
		$(compiler) $(object_files) $(ld_post_flags) -o ./post_process

# Clean-up of the working directory
clean:
		rm -f *.o *.mod *.f90~
		rm -f post_process

# Treatment of the object files
%.o:	%.f90
		$(compiler) $(f90_post_flags) $<

# ==============================================================================
