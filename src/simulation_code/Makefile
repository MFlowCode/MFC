FYPPS=$(wildcard *.fpp)
F90S=$(patsubst %.fpp,%.f90,$(wildcard *.f90 *.fpp))
OBJS=$(patsubst %.f90,%.o,$(F90S))

MOD_OBJS=$(patsubst %.fpp,%.o,$(patsubst %.f90,%.o,$(wildcard m*.f90 m*.fpp)))

EXEC = $(INSTALL_PATH)/bin/simulation
PROGRAM = p_main
PRG_OBJ = $(PROGRAM).o

all : $(PROGRAM)

define CREATE_FYPP_RULE
$(2) : $(1)
	fypp $(FYPPFLAGS) $(1) $(2)
endef

define CREATE_FC_RULE
$(2) : $(1)
	$(FC) $(FFLAGS) -c -o $(2) $< $(1) $(LNKFLAGS)
endef

$(foreach l,$(FYPPS),$(eval $(call CREATE_FYPP_RULE,$(l),$(patsubst %.fpp,%.f90,$(l)))))
$(foreach l,$(F90S), $(eval $(call CREATE_FC_RULE,  $(l),$(patsubst %.f90,%.o,  $(l)))))

$(PROGRAM) : $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) -o $(EXEC) $(LNKFLAGS)

$(PRG_OBJ) : $(MOD_OBJS)

.PHONY: all clean

clean:
	rm -rf $(OBJS) $(PROGRAM) $(EXEC) $(patsubst %.o,%.mod,$(MOD_OBJS))

# Modules depend on each other
m_rhs.o :                  m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o m_weno.o m_riemann_solvers.o m_cbc.o m_bubbles.o m_qbmm.o nvtx.o
m_start_up.o :             m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o m_compile_specific.o
m_data_output.o :          m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o m_compile_specific.o
m_riemann_solvers.o :      m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o m_bubbles.o
m_bubbles.o :              m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o 
m_qbmm.o :                 m_derived_types.o m_global_parameters.o m_variables_conversion.o m_mpi_proxy.o 

m_cbc.o :                  m_derived_types.o m_global_parameters.o m_variables_conversion.o
m_weno.o :                 m_derived_types.o m_global_parameters.o m_variables_conversion.o

m_time_steppers.o :        m_derived_types.o m_global_parameters.o m_mpi_proxy.o m_rhs.o m_data_output.o m_bubbles.o
m_derived_variables.o :    m_derived_types.o m_global_parameters.o m_mpi_proxy.o m_data_output.o m_time_steppers.o m_rhs.o
m_variables_conversion.o : m_derived_types.o m_global_parameters.o m_mpi_proxy.o nvtx.o

m_mpi_proxy.o :            m_derived_types.o m_global_parameters.o
m_global_parameters.o :    m_derived_types.o
