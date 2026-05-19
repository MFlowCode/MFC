FC      = ftn
# Flags matching MFC's AMD GPU build (CMakeLists.txt, gfx90a / Frontier MI250X)
FFLAGS  = -fopenmp --offload-arch=gfx90a -fopenmp-target-fast \
          -fopenmp-assume-threads-oversubscription \
          -fopenmp-assume-teams-oversubscription

TESTS = test01_dt_elem_access \
        test02_dt_array_constructor \
        test03_dt_whole_array_ops \
        test04_dt_zero_init \
        test05_device_routine_int_arg \
        test06_dt_member_as_array_arg \
        test07_module_param_array \
        test08_dt_private_array_elems \
        test09_dt_arr7_arithmetic \
        test10_hlld_flux

.PHONY: all clean run

all: $(TESTS)

%: %.f90
	$(FC) $(FFLAGS) -o $@ $<

run: $(TESTS)
	@echo "=== Running AMD OMP target miniapps ==="
	@for t in $(TESTS); do \
	    echo "--- $$t ---"; \
	    ./$$t; \
	done
	@echo "=== Done ==="

clean:
	rm -f $(TESTS) *.mod *.o
