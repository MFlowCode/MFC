# MFC v3.0 - Post-process Code: makefile
# Description: Definitions and rules for compilation of post-process code
# Author: Vedran Coralic
# Date: 06/14/12


# Definitions ==================================================================

mycompiler="intel"

ifeq ($(mycompiler),"gcc")
myfile = m_compile_specific_gcc.f90
else
myfile = m_compile_specific_intel.f90
endif

# Location of MKL (Heaviside)
MKLROOT = /opt/intel/mkl

# Location of the Silo directory (Hooke)
# silo_dir = /home/maeda/silo-4.7.2

# Location of the Silo directory (Darter)
# silo_dir = $(SILO_DIR)

# Location of the Silo directory (Stampede)
# silo_dir = /home1/02670/jcmeng/silo-4.9.1/

# Location of the Silo directory (Thomson)
# silo_dir = /share/apps/silo-4.7.2

# Location of the Silo directory (Heaviside)
# silo_dir = /home/maeda/Softwares/silo
# silo_dir = /usr/lib/silo-4.10.2

# Location of the Silo directory (Richardson)
#silo_dir = /home/kevinsch/software/silo-4.10.2
silo_dir = /home/spencer/packages/silo-4.10.2
#silo_dir = /home/spencer/packages/visit2.13.2/2.13.2/linux-x86_64/lib

# Location of the Silo directory (Laptop)
# silo_dir = /Users/kevinschmidmayer/Downloads/software/silo-4.10.2

# Linker flags (Hooke)
# ld_flags = -L$(silo_dir)/lib -lsiloh5 -lm -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

# Linker flags (Thomson/Stampede)
ld_flags = -L$(silo_dir)/lib -lsiloh5 -lm -L/usr/lib64 -lfftw3 -lm
# ld_flags = -L$(silo_dir)/lib -lsiloh5 -lm -L/usr/local/lib64 -lfftw3 -lm
# ld_flags = -L$(silo_dir)/lib -lsiloh5 -lz -lm -L/usr/local/lib64 -lfftw3 -lm
#ld_flags = -L/home/spencer/packages/visit2.13.2/2.13.2/linux-x86_64/lib -lsiloh5 -lz -lm -L/usr/local/lib64 -lfftw3 -lm

# Linker flags (Darter)
# ld_flags = $(SILO_LIB)

# Header flags (Hooke)
# hd_flags = -I$(silo_dir)/include -I$(MKLROOT)/include -I/usr/include

# Header flags (Richardson)
hd_flags = -I$(silo_dir)/include -I$(MKLROOT)/include -I/usr/include
#hd_flags = -I/home/spencer/packages/visit2.13.2/2.13.2/linux-x86_64/include/silo/include -I$(MKLROOT)/include -I/usr/include

# Header flags (Laptop)
# hd_flags = -I$(silo_dir)/include -I/usr/include

# Header flags (Thomson/Darter/Stampede)
# hd_flags = -I$(silo_dir)/include

# Intel compiler flags (development)
# f90_flags = -c -C -fpe-all=0 -O0 -traceback -warn $(hd_flags)
# f90_flags = -O0 -g -CA -CB -CS -CV -c -C -traceback -debug all -ftrapuv -check all,noarg_temp_created -WB -warn all -nogen-interfaces $(hd_flags)
# f90_flags = -c -C -fpe-all=0 -O0 -traceback $(hd_flags)

# Intel compiler flags (production)
f90_flags = -c -ip -O3 -w $(hd_flags)

# NAG compiler flags (development)
# f90_flags = -c -C=all -colour -f2003 -gline -mismatch -mtrace=all -O0 -v \
#             $(hd_flags)

# NAG compiler flags (production)
# f90_flags = -c -colour -f2003 -mismatch -O4 -w $(hd_flags)

# PGI compiler flags (development)
# f90_flags = -c -C -Mchkptr -Minform=inform -O0 -traceback -v $(hd_flags)

# PGI compiler flags (production)
# f90_flags = -c -Minform=severe -O4 -w $(hd_flags)

# GCC compiler flags (development)
# f90_flags = -c -C -O0 -Wall -g3 -fbacktrace -v

# GCC compiler flags (production)
# f90_flags = -ffree-line-length-none -c -O3 -w $(hd_flags)

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

# Compiler choice (Hooke/Intel compiler)
# compiler = /opt/mvapich2/ch3_mrail_gen2-intel12/bin/mpif90

# Compiler choice (Darter/Intel compiler)
# export XTPE_INFO_MESSAGE_OFF = yes
# compiler = ftn

# Compiler choice (Stampede/Intel compiler)
# compiler = /opt/apps/intel13/mvapich2/1.9/bin/mpif90

# Compiler choice (Thomson/NAG compiler)
# compiler = /share/apps/openmpi-1.4.3/nag_fort/bin/mpif90

# Compiler choice (Richardson/Intel compiler)
#compiler = mpiifort
compiler = mpifort

# Compiler choice (Laptop compiler)
# compiler = /opt/local/bin/openmpif90
# compiler = mpif90

# Compiler choice (Heaviside/Intel compiler)
# compiler = mpiifort

# ==============================================================================


# Rules ========================================================================

# Creation of the executable
all:	$(object_files)
		$(compiler) $(object_files) $(ld_flags) -o ./post_process

# Clean-up of the working directory
clean:
		rm -f *.o *.mod *.f90~
		rm -f post_process

# Treatment of the object files
m_compile_specific.o:	$(myfile)
		$(compiler) $(f90_flags) -o m_compile_specific.o $<
%.o:	%.f90
		$(compiler) $(f90_flags) $<

# ==============================================================================
