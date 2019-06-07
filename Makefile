include Makefile.user

all: pre_process simulation post_process

pre_process:
	@$(MAKE) -C src/pre_process_code -f makefile

simulation:

ifneq ("$(wildcard $(fftw_include_dir)/fftw3.*)","")
ifneq ("$(wildcard $(fftw_lib_dir)/libfftw*.la)","")
	@$(MAKE) -C src/simulation_code -f makefile
else
	@echo "==================================="
	@echo "Error: FFTW library files not found"
	@echo "==================================="
endif
else
	@echo "==================================="
	@echo "Error: FFTW include files not found"
	@echo "==================================="
endif

post_process:

ifneq ("$(wildcard $(silo_include_dir)/silo_*.inc)","")
ifneq ("$(wildcard $(silo_lib_dir)/libsilo*.*a)","")
	@$(MAKE) -C src/post_process_code -f makefile
else
	@echo "==================================="
	@echo "Error: Silo library files not found"
	@echo "==================================="
endif
else
	@echo "==================================="
	@echo "Error: Silo include files not found"
	@echo "==================================="
endif

check:
	./tests/checks.sh

.PHONY: clean
clean:
	@$(MAKE) -C src/pre_process_code -f makefile clean
	@$(MAKE) -C src/simulation_code -f makefile clean
	@$(MAKE) -C src/post_process_code -f makefile clean
