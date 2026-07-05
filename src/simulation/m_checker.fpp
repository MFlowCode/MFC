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
    use m_constants, only: recon_type_weno, recon_type_muscl, muscl_order_first_order, time_stepper_rk3, riemann_solver_hllc

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
        @:PROHIBIT(hybrid_weno .and. weno_order == 1, "hybrid_weno requires weno_order > 1")
        @:PROHIBIT(hybrid_weno .and. hybrid_weno_eps <= 0._wp, "hybrid_weno_eps must be > 0")
        @:PROHIBIT(hybrid_weno .and. igr, "hybrid_weno is incompatible with the IGR solver")
        @:PROHIBIT(hybrid_riemann .and. riemann_solver /= riemann_solver_hllc, "hybrid_riemann requires riemann_solver = 2 (HLLC)")
        @:PROHIBIT(hybrid_riemann .and. recon_type /= recon_type_weno, &
                   & "hybrid_riemann requires WENO reconstruction (the shared sensor lives in the WENO module)")
        @:PROHIBIT(hybrid_riemann .and. weno_order == 1, "hybrid_riemann requires weno_order > 1")
        @:PROHIBIT(hybrid_riemann .and. .not. (model_eqns == model_eqns_5eq .or. model_eqns == model_eqns_6eq), &
                   & "hybrid_riemann supports only the 5- and 6-equation models")
        @:PROHIBIT(hybrid_riemann .and. (hybrid_smooth_flux < 1 .or. hybrid_smooth_flux > 2), &
                   & "hybrid_smooth_flux must be 1 (central) or 2 (Rusanov)")
        @:PROHIBIT(hybrid_riemann .and. igr, "hybrid_riemann is incompatible with the IGR solver")
        @:PROHIBIT(hybrid_riemann .and. (viscous .or. surface_tension .or. hypoelasticity .or. hyperelasticity .or. elasticity), &
                   & "hybrid_riemann does not support viscous/elastic/surface-tension physics")
        @:PROHIBIT(hybrid_riemann .and. (bubbles_euler .or. bubbles_lagrange .or. qbmm), &
                   & "hybrid_riemann does not support bubble models")
        @:PROHIBIT(hybrid_riemann .and. chemistry, "hybrid_riemann does not support chemistry")
        @:PROHIBIT(hybrid_riemann .and. (cyl_coord .or. mhd), "hybrid_riemann does not support cylindrical/axisymmetric or MHD")
        @:PROHIBIT(hybrid_riemann .and. low_Mach /= 0, &
                   & "hybrid_riemann (cheap central/Rusanov flux) is incompatible with the low_Mach correction")

        if (amr) then
            @:PROHIBIT(recon_type /= recon_type_weno, "amr requires WENO reconstruction")
            @:PROHIBIT(time_stepper /= time_stepper_rk3, "amr requires time_stepper = 3 (SSP-RK3)")
            @:PROHIBIT(model_eqns /= 2, "amr requires model_eqns = 2 (5-equation)")
            @:PROHIBIT(num_fluids > 1 .and. .not. mpp_lim, &
                       & "amr with num_fluids > 1 requires mpp_lim (its volume-fraction clamp+renormalize maintains coarse/fine alpha consistency)")
            @:PROHIBIT(surface_tension, &
                       & "amr does not support surface_tension: the capillary force depends on the interface normal (grad-c direction), which the conservative-linearly-prolonged fine ghost color cannot reproduce consistently with the coarse solver across a 2:1 coarse/fine boundary - a growing spurious seam current results")
            @:PROHIBIT(hypoelasticity .or. hyperelasticity .or. mhd, "amr does not support elastic/MHD")
            @:PROHIBIT(bubbles_lagrange .or. igr .or. cyl_coord, "amr does not support Lagrangian bubbles/IGR/cylindrical")
            @:PROHIBIT(qbmm .and. .not. polytropic, &
                       & "amr does not support non-polytropic QBMM: its pb/mv quadrature side-state evolves as a global array that the fine advance would corrupt through the swap. Polytropic QBMM (pb/mv inert) is supported: its bubble moments live entirely in q_cons and are injected piecewise-constant at prolongation to preserve CHyQMOM realizability")
            ! static-body IB AMR (SP20) + prescribed-motion moving bodies (SP21): fixed or analytically-moving
            ! (moving_ibm==1) bodies resolved on a static fine block. Multi-body (num_ibs>1) is supported - every body
            ! shares the one static block and reuses the multi-body-capable core IB setup. Force/torque-driven motion
            ! (moving_ibm==2), STL geometry, and dynamic regrid with IB remain gated (unvalidated / recompute-on-move).
            @:PROHIBIT(ib .and. any(patch_ib(1:num_ibs)%moving_ibm == 2), &
                       & "amr with ib supports static or prescribed-motion (moving_ibm=1) bodies only; force-driven moving IB (moving_ibm=2) under amr is not yet validated")
            @:PROHIBIT(ib .and. any(patch_ib(1:num_ibs)%geometry == 12), &
                       & "amr with ib does not support STL-model geometry (not yet validated)")
            @:PROHIBIT(ib .and. amr_regrid_int > 0, &
                       & "amr with ib requires a static block (amr_regrid_int = 0); dynamic regrid with IB is not yet validated")
            @:PROHIBIT(active_box, "amr is incompatible with active_box (unvalidated combination)")
            @:PROHIBIT(hybrid_weno, "amr is incompatible with hybrid_weno (unvalidated combination)")
            @:PROHIBIT(hybrid_riemann, "amr is incompatible with hybrid_riemann (unvalidated combination)")
            @:PROHIBIT(acoustic_source, &
                       & "amr is incompatible with acoustic_source (dt-dependent RHS source; unvalidated with the fine-level advance)")
            @:PROHIBIT(any(amr_block_beg(1:num_dims) < 0), "amr_block_beg must be >= 0")
            @:PROHIBIT(amr_block_end(1) > m_glb .or. (num_dims >= 2 .and. amr_block_end(2) > n_glb) .or. (num_dims >= 3 &
                       & .and. amr_block_end(3) > p_glb), "amr_block_end must be <= global cell max per axis")
            @:PROHIBIT(any(amr_block_end(1:num_dims) <= amr_block_beg(1:num_dims)), &
                       & "amr_block_end must exceed amr_block_beg on each active axis")
            @:PROHIBIT(2*(amr_block_end(1) - amr_block_beg(1) + 1) - 1 > m_glb, &
                       & "amr fine x-extent exceeds the base grid (module scratch is sized to the base)")
            @:PROHIBIT(num_dims >= 2 .and. 2*(amr_block_end(2) - amr_block_beg(2) + 1) - 1 > n_glb, &
                       & "amr fine y-extent exceeds the base grid")
            @:PROHIBIT(num_dims >= 3 .and. 2*(amr_block_end(3) - amr_block_beg(3) + 1) - 1 > p_glb, &
                       & "amr fine z-extent exceeds the base grid")
            @:PROHIBIT(amr_regrid_int < 0, "amr_regrid_int must be >= 0")
            @:PROHIBIT(amr_regrid_int > 0 .and. amr_tag_eps <= 0._wp, "amr_tag_eps must be > 0 when regridding")
            @:PROHIBIT(amr_regrid_int > 0 .and. amr_buf < 1, "amr_buf must be >= 1 when regridding")
            @:PROHIBIT(amr_max_blocks < 1, "amr_max_blocks must be >= 1")
            @:PROHIBIT(amr_cluster_eff <= 0._wp .or. amr_cluster_eff > 1._wp, &
                       & "amr_cluster_eff must satisfy 0 < amr_cluster_eff <= 1")
        end if
        @:PROHIBIT(.not. amr .and. amr_regrid_int > 0, "amr_regrid_int requires amr")
        @:PROHIBIT(amr_subcycle .and. .not. amr, "amr_subcycle requires amr")
        @:PROHIBIT(amr_subcycle .and. cfl_dt, "amr_subcycle requires a fixed dt (cfl_dt not supported)")

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
