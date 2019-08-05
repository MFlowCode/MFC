include Makefile.user 

default: all

include installers/Makefile.messages

all: pre_process simulation post_process

pre_process:
	@$(MAKE) -C src/pre_process_code -f makefile

simulation:

ifneq ("$(wildcard $(fftw_include_dir)/fftw3.*)","")
ifneq ("$(wildcard $(fftw_lib_dir)/libfftw*.la)","")
	@$(MAKE) -C src/simulation_code -f makefile
else
	@echo "$$FFTW_LIB_ERR"
endif
else
	@echo "$$FFTW_INC_ERR"
endif

post_process:

ifneq ("$(wildcard $(silo_include_dir)/silo*.inc)","")
ifneq ("$(wildcard $(silo_lib_dir)/libsilo*.*a)","")
	@$(MAKE) -C src/post_process_code -f makefile
else
	@echo "$$SILO_LIB_ERR"
endif
else
	@echo "$$SILO_INC_ERR"
endif

test:
	./tests/checks.sh

.PHONY: clean
clean:
	@$(MAKE) -C src/pre_process_code -f makefile clean
	@$(MAKE) -C src/simulation_code -f makefile clean
	@$(MAKE) -C src/post_process_code -f makefile clean
