# Initial configuration
CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

PROJECT(MFC LANGUAGES Fortran)


# Imports
INCLUDE(GNUInstallDirs)
INCLUDE(ExternalProject)


# MPI
FIND_PACKAGE(MPI REQUIRED COMPONENTS Fortran)

ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_DIRS})
LINK_DIRECTORIES(${MPI_Fortran_LIBRARIES})
LINK_LIBRARIES(${MPI_Fortran_LIBRARIES})
LIST(APPEND CMAKE_EXE_LINKER_FLAGS ${MPI_Fortran_LINK_FLAGS})


# Options
OPTION(BUILD_PRE  "Build pre_process"  OFF)
OPTION(BUILD_SIM  "Build simulation"   OFF)
OPTION(BUILD_POST "Build post_process" OFF)
OPTION(BUILD_ALL  "Build pre_process, simulation, and post_process" OFF)

IF(BUILD_ALL)
    SET(BUILD_PRE  ON FORCE)
    SET(BUILD_SIM  ON FORCE)
    SET(BUILD_POST ON FORCE)
ENDIF()


# Pre Process
IF(BUILD_PRE)
    FILE(GLOB pre_process_f90s "src/pre_process_code/*.f90")
    ADD_EXECUTABLE(pre_process "${pre_process_f90s}")

    INSTALL(TARGETS pre_process RUNTIME DESTINATION bin)
ENDIF()


if(BUILD_PRE OR BUILD_SIM)
    # FFTW3
    ExternalProject_Add(fftw3
        URL        "http://www.fftw.org/fftw-3.3.10.tar.gz"
        CMAKE_ARGS -DBUILD_TESTS=OFF
                   -DBUILD_SHARED_LIBS=OFF
                   -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    )
endif()


# Simulation
IF(BUILD_SIM)
    # Run FYPP
    FILE(GLOB simulation_f90s "src/simulation_code/*.f90")
    FILE(GLOB simulation_fpps "src/simulation_code/*.fpp")

    FOREACH(in_fpp_file IN LISTS simulation_fpps)
        STRING(REGEX REPLACE ".fpp\$" ".f90" out_fpp_file "${in_fpp_file}")
        
        ADD_CUSTOM_COMMAND(
            OUTPUT          "${out_fpp_file}"
            COMMAND          fypp "${in_fpp_file}" "${out_fpp_file}"
            MAIN_DEPENDENCY "${in_fpp_file}"
            VERBATIM)
        
        SET(simulation_f90s ${simulation_f90s} "${out_fpp_file}")
    ENDFOREACH(in_fpp_file)

    ADD_EXECUTABLE(simulation "${simulation_f90s}")

    # ZFP
    #ExternalProject_Add(zfp
    #    GIT_REPOSITORY "https://github.com/LLNL/zfp"
    #    GIT_TAG        980155a8e83bd6dabbcc600a4059df4d54521bfe
    #    CMAKE_ARGS     -DZFP_WITH_CUDA=ON
    #                   -DZFP_WITH_OPENMP=OFF
    #                   -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    #)
    
    TARGET_INCLUDE_DIRECTORIES(simulation PRIVATE "${CMAKE_INSTALL_PREFIX}/include")
    TARGET_LINK_LIBRARIES     (simulation PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/libfftw3.a")
    
    ADD_DEPENDENCIES(simulation fftw3)

    INSTALL(TARGETS simulation RUNTIME DESTINATION bin)
ENDIF()

# Post Process
IF(BUILD_POST)
    # HDF5
    ExternalProject_Add(hdf5
        GIT_REPOSITORY "https://github.com/HDFGroup/hdf5"
        GIT_TAG        hdf5-1_12_1
        CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
                       -DBUILD_SHARED_LIBS=OFF
                       -DFORTRAN_LIBRARIES=ON
                       -DBUILD_TESTING=OFF
                       -DHDF5_BUILD_UTILS=OFF
                       -DHDF5_BUILD_TOOLS=OFF
                       -DHDF5_BUILD_EXAMPLES=OFF
    )

    # SILO
    ExternalProject_Add(silo
        GIT_REPOSITORY "https://github.com/henryleberre/Silo"
        GIT_TAG        5c69eacd526002ea4cab063f8231e3c7ce0a9fb9
        CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                       -DSILO_HDF5_DIR=${CMAKE_INSTALL_PREFIX}
                       -DSILO_ENABLE_BROWSER=OFF
    )

    ADD_DEPENDENCIES(silo hdf5)

    FILE(GLOB post_process_f90s "src/post_process_code/*.f90")
    ADD_EXECUTABLE(post_process "${post_process_f90s}")

    TARGET_INCLUDE_DIRECTORIES(post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/include")
    TARGET_LINK_LIBRARIES     (post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/libfftw3.a")
    TARGET_LINK_LIBRARIES     (post_process PRIVATE "${CMAKE_INSTALL_PREFIX}/lib/libsiloh5.so")
    
    ADD_DEPENDENCIES(post_process silo hdf5 fftw3)

    INSTALL(TARGETS post_process RUNTIME DESTINATION bin)
ENDIF()
