all:
	@$(MAKE) -C pre_process_code -f makefile
	@$(MAKE) -C simulation_code -f makefile
	@$(MAKE) -C post_process_code -f makefile

pre_process:
	@$(MAKE) -C pre_process_code -f makefile

simulation:
	@$(MAKE) -C simulation_code -f makefile

post_process:
	@$(MAKE) -C post_process_code -f makefile

.PHONY: clean
clean:
	@$(MAKE) -C pre_process_code -f makefile clean
	@$(MAKE) -C simulation_code -f makefile clean
	@$(MAKE) -C post_process_code -f makefile clean
