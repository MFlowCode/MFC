
IF (MFC_BUILD_PRE)
    ADD_SUBDIRECTORY(pre_process_code)
ENDIF()

IF (MFC_BUILD_SIM)
    ADD_SUBDIRECTORY(simulation_code)
ENDIF()

IF (MFC_BUILD_POST)
    ADD_SUBDIRECTORY(post_process_code)
ENDIF()
