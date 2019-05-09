# MFC v3.0 - Simulation Code: makefile
# Description: Definitions and rules for compilation of simulation code
# Author: Vedran Coralic
# Date: 06/14/12


# Definitions ==================================================================

#For compiler-specific modules, options (gcc,intel)
mycompiler="gcc"

ifeq ($(mycompiler),"gcc")
myfile = m_compile_specific_gcc.f90
else
myfile = m_compile_specific_intel.f90
endif

# MKLROOT
# MKLROOT = /opt/intel/mkl/10.3

# MKLROOT (Heaviside)
MKLROOT = /opt/intel/mkl

# Linker flags for FFTW library (Comet)
#ld_flags = -L/usr/lib64 -lfftw3 -lm

# Linker flags for FFTW library (Richardson)
ld_flags = -L/usr/local/lib -lfftw3 -lm

# Header flags for FFTW library (Comet)
#hd_flags = -I/usr/include

# Header flags for FFTW library (Richardson)
#hd_flags = -I/usr/local/include
hd_flags = -I/usr/local/Cellar/fftw/3.3.8/include/fftw3.f03

# Linker flags (Hooke)
#ld_flags = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

# Header flags (Hooke/Heaviside)
#hd_flags =  -I$(MKLROOT)/include

# Intel compiler flags (development)
# f90_flags = -O0 -g -CA -CB -CS -CV -c -C -traceback -debug all -ftrapuv -check all,noarg_temp_created -WB -warn all -nogen-interfaces $(hd_flags)
# f90_flags = -c -C -fpe-all=0 -O0 -traceback -warn $(hd_flags)
# f90_flags = -c -ip -O0 -w -traceback $(hd_flags)


# Intel compiler flags (production)
# f90_flags = -c -ip -O3 -w -heap-arrays $(hd_flags)
#f90_flags = -c -ip -O3 -w $(hd_flags)
# f90_flags = -c -ip -O3 -w $(hd_flags)
# f90_flags = -c -ip -O1 -w $(hd_flags)

# NAG compiler flags (development)
# f90_flags = -c -C=all -colour -f2003 -gline -mismatch -mtrace=all -O0 -v

# NAG compiler flags (production)
# f90_flags = -c -colour -f2003 -mismatch -O4 -w
# f90_flags = -c -colour -f2003 -mismatch -O0 -w

# PGI compiler flags (development)
# f90_flags = -c -C -Mchkptr -Minform=inform -O0 -traceback -v

# PGI compiler flags (production)
# f90_flags = -c -Minform=severe -O4 -w

# GCC compiler flags (development)
#f90_flags = -c -C -O0 -Wall -g3 -fbacktrace -v

# GCC compiler flags (production)
f90_flags = -ffree-line-length-none -c -O3 -w 

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
			   m_anti_diffusion.o \
			   p_main.o

# Compiler choice (Hooke/Intel compiler)
# compiler = /opt/mvapich2/ch3_mrail_gen2-intel12/bin/mpif90

# Compiler choice (Darter/Intel compiler)
# export XTPE_INFO_MESSAGE_OFF = yes
# compiler = ftn

# Compiler choice (Stampede/Intel compiler)
# compiler = /opt/apps/intel13/mvapich2/1.9/bin/mpif90

# Compiler choice (Thomson/NAG compiler)
# compiler = /share/apps/openmpi-1.4.3/nag_fort/bin/mpif90

# Compiler choice (Laptop compiler)
# compiler = /opt/local/bin/openmpif90
# compiler = mpif90

# Compiler choice (Comet)
# compiler = mpif90

# Compiler choice (Heaviside)
# compiler = mpiifort

# Compiler choice (Richardson)
compiler = mpif90

# ==============================================================================


# Rules ========================================================================

# Creation of the executable
all:	$(object_files)
		$(compiler) $(object_files) $(ld_flags) -o ./simulation

# Clean-up of the working directory
clean:
		rm -f *.o *.mod *.f90~
		rm -f simulation

# Treatment of the object files
m_compile_specific.o:	$(myfile)
		$(compiler) $(f90_flags) -o m_compile_specific.o $<

%.o:	%.f90
		$(compiler) $(f90_flags) $<

# ===============================================================================
