all:
	@$(MAKE) -C src/pre_process_code -f makefile
	@$(MAKE) -C src/simulation_code -f makefile
	@$(MAKE) -C src/post_process_code -f makefile

pre_process:
	@$(MAKE) -C src/pre_process_code -f makefile

simulation:
	@$(MAKE) -C src/simulation_code -f makefile

post_process:
	@$(MAKE) -C src/post_process_code -f makefile

.PHONY: clean
clean:
	@$(MAKE) -C src/pre_process_code -f makefile clean
	@$(MAKE) -C src/simulation_code -f makefile clean
	@$(MAKE) -C src/post_process_code -f makefile clean
