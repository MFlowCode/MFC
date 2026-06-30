!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic
    use m_constants, only: recon_type_weno, recon_type_muscl, muscl_order_first_order, time_stepper_rk3

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        call s_check_inputs_time_stepping

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")
        @:PROHIBIT(many_ib_patch_parallelism .and. .not. ib, "many_ib_patch_parallelism requires ib to be enabled")

        if (active_box) then
            @:PROHIBIT(recon_type /= recon_type_weno, "active_box requires WENO reconstruction")
            @:PROHIBIT(ib, "active_box is incompatible with immersed boundaries")
            @:PROHIBIT(acoustic_source, "active_box is incompatible with acoustic sources")
            @:PROHIBIT(bodyForces, "active_box is incompatible with body forces")
            @:PROHIBIT(bubbles_lagrange, "active_box is incompatible with Lagrangian bubbles")
            @:PROHIBIT(relax, "active_box is incompatible with phase change")
            @:PROHIBIT(igr, "active_box is incompatible with the IGR solver")
            @:PROHIBIT(time_stepper /= time_stepper_rk3, "active_box requires time_stepper = 3 (SSP-RK3)")
            @:PROHIBIT(viscous, "active_box is incompatible with viscous (no finite domain of dependence for the frozen exterior)")
            @:PROHIBIT(surface_tension, &
                       & "active_box is incompatible with surface_tension (nonlocal curvature coupling violates the static-uniform-exterior assumption)")
            @:PROHIBIT(cyl_coord, &
                       & "active_box is incompatible with cyl_coord (geometric source terms are nonzero for uniform flow; exterior is not static)")
            @:PROHIBIT(hypoelasticity, &
                       & "active_box is incompatible with hypoelasticity (stress source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(hyperelasticity, &
                       & "active_box is incompatible with hyperelasticity (stress source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(mhd, &
                       & "active_box is incompatible with mhd (magnetic field source terms violate the static-uniform-exterior assumption)")
            @:PROHIBIT(chemistry, &
                       & "active_box is incompatible with chemistry (reactive source terms violate the static-uniform-exterior assumption)")
        end if

        @:PROHIBIT(sfc_partition_wrt .and. partition_tile_size < 1, "partition_tile_size must be >= 1")
        @:PROHIBIT(load_balance .and. .not. parallel_io, "load_balance requires parallel_io = T")
        @:PROHIBIT(load_balance .and. num_procs == 1, "load_balance requires more than one MPI rank")
        @:PROHIBIT(hybrid_weno .and. recon_type /= recon_type_weno, "hybrid_weno requires WENO reconstruction")
        @:PROHIBIT(hybrid_weno .and. hybrid_weno_eps <= 0._wp, "hybrid_weno_eps must be > 0")

        if (num_particle_clouds > 0) then
            call s_check_inputs_particle_clouds
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(muscl_order == muscl_order_first_order .and. int_comp > 0, &
                   & "int_comp requires muscl_order >= 2 (muscl_order=1 leaves the reconstruction workspace uninitialised)")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    !> Checks that each active particle cloud has a valid packing_method specified
    impure subroutine s_check_inputs_particle_clouds

        integer          :: i
        character(len=5) :: idxStr

        do i = 1, num_particle_clouds
            call s_int_to_str(i, idxStr)
            @:PROHIBIT(particle_cloud(i)%packing_method == dflt_int, &
                       & "particle_cloud("//trim(idxStr) &
                       & //")%packing_method must be specified (1 = rejection sampling, 2 = lattice)")
            @:PROHIBIT(particle_cloud(i)%packing_method /= 1 .and. particle_cloud(i)%packing_method /= 2, &
                       & "particle_cloud("//trim(idxStr) //")%packing_method must be 1 (rejection sampling) or 2 (lattice)")
        end do

    end subroutine s_check_inputs_particle_clouds

end module m_checker
