!>
!! @file
!! @brief Contains module m_rhs

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'graph_macros.fpp'

!> @brief Assemble the right-hand side of the governing equations using finite-volume flux differencing, Riemann solvers, and
!! physical source terms
module m_rhs

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_eos
    use m_weno
    use m_constants, only: riemann_solver_hll, riemann_solver_hlld, model_eqns_6eq, int_comp_mthinc, recon_type_weno, &
        & recon_type_muscl
#if defined(MFC_FUSED_RHS) || defined(MFC_DIRECT_RIEMANN_FLUX) || defined(MFC_CUDA_GRAPHS)
    use m_constants, only: riemann_solver_hllc, model_eqns_5eq, int_comp_thinc, avg_state_arithmetic, wave_speeds_direct, &
        & riemann_solver_lax_friedrichs, BC_CHAR_SLIP_WALL, BC_CHAR_SUP_OUTFLOW
#endif
#ifdef MFC_FUSED_RHS
    use m_constants, only: wave_speeds_pressure
    #:if MFC_CASE_OPTIMIZATION and hllc_alt_soundspeed
        use m_fused_rhs_openacc, only: f_hllc_kapila
    #:endif
#endif
    use m_muscl
    use m_riemann_solvers
#if defined(MFC_DIRECT_RIEMANN_FLUX) || defined(MFC_FUSED_RHS)
    use m_riemann_state, only: flux_rsx_vf, flux_src_rsx_vf
#endif
#if defined(MFC_FUSED_RHS)
    use m_riemann_state, only: s_fused_viscous_source_flux, s_fused_source_flux_state, isx, isy, isz, vel_src_rsx_vf
#endif
    use m_cbc
    use m_bubbles_EE
    use m_bubbles_EL
    use m_qbmm
    use m_hypoelastic
    use m_acoustic_src
    use m_viscous
    use m_ibm
    use m_nvtx
    use m_boundary_common
    use m_helper
    use m_surface_tension
    use m_body_forces
    use m_chemistry
    use m_conduction
    use m_reactive_burn
    #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous))
#ifdef MFC_FUSED_RHS
        use m_fused_rhs_openacc, only: s_fused_rhs_openacc, s_fused_rhs_openacc_overlap_bounds
        #:if model_eqns == 2
            use m_fused_rhs_openacc, only: s_fused_rhs_openacc_boundary_rk, s_fused_rhs_openacc_update_rk
        #:endif
#endif
    #:endif
    use m_igr
    use m_thinc
    use m_pressure_relaxation
#ifdef MFC_CUDA_GRAPHS
    use m_cuda_graph
#endif
#ifdef MFC_ACTIVE_BOX
    use m_active_box, only: ab_x, ab_y, ab_z, ab_active
#endif

    implicit none

    private; public :: s_initialize_rhs_module, s_associate_rhs_flux_views, s_compute_rhs, s_finalize_rhs_module
#ifdef MFC_FUSED_RHS
    public :: f_use_fused_rhs_openacc_final_rk, f_use_fused_rhs_openacc_rk, s_update_fused_rhs_openacc_state
#endif

    type(vector_field) :: q_cons_qp  !< WENO-reconstructed cell-average conservative variables at quadrature points
    $:GPU_DECLARE(create='[q_cons_qp]')

    type(vector_field) :: q_prim_qp  !< Primitive variables at cell-interior quadrature points
    $:GPU_DECLARE(create='[q_prim_qp]')
#if defined(MFC_FUSED_RHS)
    type(vector_field) :: q_prim_next  !< Next-stage primitive variables for the fused RK path
    $:GPU_DECLARE(create='[q_prim_next]')

#ifdef MFC_ACTIVE_BOX
    logical :: ab_prim_refreshed = .false.  !< Both primitive arrays fully converted since the active box engaged
#endif
#endif

    !> @name The first-order spatial derivatives of the primitive variables at cell- interior Gaussian quadrature points. These are
    !! WENO-reconstructed from their respective cell-average values, obtained through the application of the divergence theorem on
    !! the integral-average cell-boundary values of the primitive variables, located in qK_prim_n, where K = L or R.
    !> @{
    type(vector_field), allocatable, dimension(:) :: dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp
    $:GPU_DECLARE(create='[dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp]')
    !> @}

    !> @name The left and right WENO-reconstructed cell-boundary values of the cell- average first-order spatial derivatives of the
    !! primitive variables. The cell-average of the first-order spatial derivatives may be found in the variables dq_prim_ds_qp,
    !! where s = x, y or z.
    !> @{
    type(vector_field), allocatable, dimension(:) :: dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n
    type(vector_field), allocatable, dimension(:) :: dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n]')
    $:GPU_DECLARE(create='[dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n]')
#endif
    !> @}

    type(scalar_field), allocatable, dimension(:) :: tau_Re_vf
    $:GPU_DECLARE(create='[tau_Re_vf]')

    !> @name The cell-boundary values of the fluxes (src - source, gsrc - geometrical source). These are computed by applying the
    !! chosen Riemann problem solver on the left and right cell-boundary values of the primitive variables
    !> @{
    type(vector_field), allocatable, dimension(:) :: flux_n
    type(vector_field), allocatable, dimension(:) :: flux_src_n
    type(vector_field), allocatable, dimension(:) :: flux_gsrc_n

#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[flux_n, flux_src_n, flux_gsrc_n]')
#endif
    !> @}

    type(vector_field), allocatable, dimension(:) :: nc_iface_vel_n
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[nc_iface_vel_n]')
#endif

    !> hat_R-pass interface velocities for fused dual-pass HLLD: the axisymmetric geometric source consumes both passes' face
    !! velocities after the direction loop, so the hat_R values need their own field set (hat_L uses nc_iface_vel_n).
    type(vector_field), allocatable, dimension(:) :: nc_iface_vel_hatR_n
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[nc_iface_vel_hatR_n]')
#endif

    type(vector_field), allocatable, dimension(:) :: qL_prim, qR_prim
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[qL_prim, qR_prim]')
#endif

    !> Reconstruction and loop bounds, host-only; kernels receive them as arguments
    type(int_bounds_info) :: iv  !< Vector field indical bounds

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(int_bounds_info) :: irx, iry, irz
    $:GPU_DECLARE(create='[irx, iry, irz]')

    type(int_bounds_info) :: is1, is2, is3
    !> @}

    !> @name Saved fluxes for testing
    !> @{
    type(scalar_field) :: alf_sum
    !> @}
    $:GPU_DECLARE(create='[alf_sum]')

    real(wp), allocatable, dimension(:,:,:)   :: blkmod1, blkmod2, alpha1, alpha2, Kterm
    real(wp), allocatable, dimension(:,:,:,:) :: qL_rsx_vf, qR_rsx_vf
    real(wp), allocatable, dimension(:,:,:,:) :: dqL_rsx_vf, dqR_rsx_vf
    $:GPU_DECLARE(create='[blkmod1, blkmod2, alpha1, alpha2, Kterm]')
    $:GPU_DECLARE(create='[qL_rsx_vf, qR_rsx_vf]')
    $:GPU_DECLARE(create='[dqL_rsx_vf, dqR_rsx_vf]')

    !> @name Partial right-hand sides of the two anchored passes of dual-pass HLLD
    !> @{
    type(scalar_field), allocatable, dimension(:) :: rhs_hatL_vf, rhs_hatR_vf
    !> @}
    $:GPU_DECLARE(create='[rhs_hatL_vf, rhs_hatR_vf]')

#if defined(MFC_FUSED_RHS) || defined(MFC_NATIVE_HALO_OVERLAP)
    !> Overlap the halo exchange with the interior x columns, then finish the halo and boundary columns
    logical :: halo_overlap = .false.
    logical :: native_halo_overlap_active = .false.
#endif

contains

#if defined(MFC_FUSED_RHS) || defined(MFC_CUDA_GRAPHS)
    logical function f_has_characteristic_bc()

        integer, dimension(6) :: gate_bcs

        gate_bcs = (/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/)
        f_has_characteristic_bc = any(gate_bcs <= BC_CHAR_SLIP_WALL .and. gate_bcs >= BC_CHAR_SUP_OUTFLOW)

    end function f_has_characteristic_bc
#endif

#ifdef MFC_FUSED_RHS
    logical function f_use_fused_rhs_openacc()

        #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and not hllc_hll_u_interface and (not surface_tension or (hllc_riemann == 2 and not viscous)) and (not heat_conduction or (not viscous and not chemistry and not surface_tension and not bubbles_euler and not hypoelasticity and not mhd and not relativity))
            f_use_fused_rhs_openacc = .not. igr .and. sys_size == ${3 if model_eqns == 3 else 2}$*num_fluids &
                & + ${'3' if mhd else 'num_dims'}$ &
                & + 1${' + 1' if surface_tension else ''}$${' + ' + str(nb*(6 if qbmm else 2 if polytropic else 4) + int(adv_n)) if bubbles_euler else ''}$${' + ' + str(num_dims*(num_dims + 1)//2 + int(cyl_coord)) if hypoelasticity else ''}$${' + 1' if cont_damage else ''}$${' + ' + str(hllc_num_species) if chemistry else ''}$${' + ' + str((2 if num_dims == 1 else 3) + int(hyper_cleaning)) if mhd else ''}$ .and. (riemann_solver == riemann_solver_hllc .or. riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_lax_friedrichs${' .or. (riemann_solver == riemann_solver_hlld .and. .not. f_has_characteristic_bc())' if mhd else ''}$) .and. (.not. alt_soundspeed .or. ((num_fluids == 2 .or. num_fluids == 3) .and. riemann_solver /= riemann_solver_hll)) .and. (.not. heat_conduction .or. (.not. alt_soundspeed .and. grid_geometry /= 3)) .and. (wave_speeds == wave_speeds_direct .or. (wave_speeds == wave_speeds_pressure .and. avg_state == avg_state_arithmetic)) .and. ${'.not. chemistry .and. ' if not chemistry else ''}$ (.not. cont_damage .or. ${'.true.' if hllc_hypoelastic else '.false.'}$) .and. (.not. f_has_characteristic_bc() .or. ${'min(m, n, p)' if num_dims == 3 else 'min(m, n)' if num_dims == 2 else 'm'}$ > 2*buff_size + 2)
        #:else
            f_use_fused_rhs_openacc = .false.
        #:endif

    end function f_use_fused_rhs_openacc

    logical function f_use_hllc_halo_overlap()

        integer :: first, last

        f_use_hllc_halo_overlap = .false.
#ifdef MFC_NATIVE_HALO_OVERLAP
        #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous)) and not viscous and not qbmm and not surface_tension and hllc_int_comp != 2
#if !defined(MFC_NVSHMEM) || !defined(__NVCOMPILER)
            if (.not. rdma_mpi) return
#endif
            call s_fused_rhs_openacc_overlap_bounds(first, last)
            f_use_hllc_halo_overlap = halo_overlap .and. f_use_fused_rhs_openacc() .and. num_procs > 1 .and. num_dims > 1 &
                & .and. first <= last .and. ((recon_type == recon_type_weno .and. (weno_order == 1 .or. weno_order == 3 &
                & .or. weno_order == 5 .or. weno_order == 7) .and. (wenojs .or. mapped_weno .or. wenoz .or. teno)) &
                & .or. (recon_type == recon_type_muscl .and. (muscl_order == 1 .or. muscl_order == 2))) .and. .not. weno_Re_flux &
                & .and. (int_comp == 0 .or. (num_fluids == 2 .and. int_comp == int_comp_thinc))
        #:endif
#endif

    end function f_use_hllc_halo_overlap

    logical function f_use_fused_rhs_openacc_final_rk()

        #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and model_eqns == 2 and not bubbles_euler and (not mhd or (hllc_mhd and not hyper_cleaning)) and not relativity and not hypoelasticity and (not viscous or hllc_viscous_stress) and not surface_tension and not chemistry and not heat_conduction
            ! Heat conduction keeps the stepper-side RK update because its source is added to rhs_vf after the sweep
            f_use_fused_rhs_openacc_final_rk = f_use_fused_rhs_openacc() .and. ((recon_type == recon_type_weno &
                & .and. (weno_order == 1 .or. weno_order == 3 .or. weno_order == 5 .or. weno_order == 7) .and. (wenojs &
                & .or. mapped_weno .or. wenoz .or. teno)) .or. (recon_type == recon_type_muscl .and. (muscl_order == 1 &
                & .or. muscl_order == 2))) .and. (.not. weno_Re_flux .or. viscous) .and. (int_comp == 0 .or. (num_fluids == 2 &
                & .and. (int_comp == int_comp_thinc .or. int_comp == int_comp_mthinc))) .and. .not. bodyForces &
                & .and. .not. synthetic_turbulence .and. .not. ib .and. .not. relax .and. .not. bubbles_lagrange &
                & .and. .not. acoustic_source .and. .not. probe_wrt .and. .not. igr .and. grid_geometry /= 3 .and. .not. (viscous &
                & .and. cyl_coord .and. (bc_y%beg == -2 .or. bc_y%beg == -14))
        #:else
            f_use_fused_rhs_openacc_final_rk = .false.
        #:endif

    end function f_use_fused_rhs_openacc_final_rk

    logical function f_use_fused_rhs_openacc_rk()

        f_use_fused_rhs_openacc_rk = f_use_fused_rhs_openacc_final_rk()
        #:if MFC_CASE_OPTIMIZATION and model_eqns == 2 and not bubbles_euler and not hllc_reactive_burn and not chemistry and not relativity
            ! Sources and probes consume the stage state before this separate RK update.
            f_use_fused_rhs_openacc_rk = f_use_fused_rhs_openacc_rk .or. (f_use_fused_rhs_openacc() &
            #:if not ((viscous and not hllc_viscous_stress) or surface_tension or hllc_hypoelastic or hllc_mhd)
                & .and. (acoustic_source .or. probe_wrt) &
            #:endif
            & .and. .not. bodyForces .and. .not. synthetic_turbulence .and. .not. ib .and. .not. relax .and. .not. bubbles_lagrange)
        #:endif

    end function f_use_fused_rhs_openacc_rk

    subroutine s_update_fused_rhs_openacc_state(q_cons, q_stor, q_output, rhs, save_state, q_coef, stor_coef, rhs_coef, den_coef)

        type(scalar_field), intent(inout) :: q_cons(:), q_stor(:), q_output(:)
        type(scalar_field), intent(in)    :: rhs(:)
        logical, intent(in)               :: save_state
        real(wp), intent(in)              :: q_coef, stor_coef, rhs_coef, den_coef

        #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and model_eqns == 2 and not bubbles_euler and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous))
            call s_fused_rhs_openacc_update_rk(q_cons, q_stor, q_output, q_prim_qp%vf, q_prim_qp%vf, rhs, save_state, q_coef, &
                                               & stor_coef, rhs_coef, den_coef)
        #:endif

    end subroutine s_update_fused_rhs_openacc_state
#endif

    !> Initialize the RHS module
    impure subroutine s_initialize_rhs_module

        integer :: i, j, k, l, id  !< Generic loop iterators

#ifdef MFC_NATIVE_HALO_OVERLAP
        character(len=8) :: halo_overlap_env
        integer          :: halo_overlap_env_status

        halo_overlap = .false.
        call get_environment_variable('MFC_HALO_OVERLAP', halo_overlap_env, status=halo_overlap_env_status)
        if (halo_overlap_env_status == 0) halo_overlap = trim(halo_overlap_env) == '1'
#endif

        $:GPU_ENTER_DATA(copyin='[idwbuff]')
        $:GPU_UPDATE(device='[idwbuff]')

        @:ALLOCATE(q_cons_qp%vf(1:sys_size))
        @:ALLOCATE(q_prim_qp%vf(1:sys_size))
#if defined(MFC_FUSED_RHS)
        @:ALLOCATE(q_prim_next%vf(1:sys_size))
#endif

        if (.not. igr) then
            do l = 1, sys_size
                @:ALLOCATE(q_cons_qp%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            do l = eqn_idx%mom%beg, eqn_idx%E
                @:ALLOCATE(q_prim_qp%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
#if defined(MFC_FUSED_RHS)
            do l = 1, sys_size
                @:ALLOCATE(q_prim_next%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
#endif
        end if

        if (surface_tension) then
            do l = eqn_idx%adv%end + 1, eqn_idx%c - 1
                @:ALLOCATE(q_prim_qp%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        else
            do l = eqn_idx%adv%end + 1, sys_size
                @:ALLOCATE(q_prim_qp%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        end if

        if (.not. igr) then
            @:ACC_SETUP_VFs(q_cons_qp, q_prim_qp)
#if defined(MFC_FUSED_RHS)
            @:ACC_SETUP_VFs(q_prim_next)
#endif

            do l = 1, eqn_idx%cont%end
                if (relativity) then
                    ! Cons and Prim densities are different for relativity
                    @:ALLOCATE(q_prim_qp%vf(l)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                else
                    q_prim_qp%vf(l)%sf => q_cons_qp%vf(l)%sf
                    $:GPU_ENTER_DATA(copyin='[q_prim_qp%vf(l)%sf]')
                    $:GPU_ENTER_DATA(attach='[q_prim_qp%vf(l)%sf]')
                end if
            end do

            do l = eqn_idx%adv%beg, eqn_idx%adv%end
                q_prim_qp%vf(l)%sf => q_cons_qp%vf(l)%sf
                $:GPU_ENTER_DATA(copyin='[q_prim_qp%vf(l)%sf]')
                $:GPU_ENTER_DATA(attach='[q_prim_qp%vf(l)%sf]')
            end do
        end if

        if (surface_tension) then
            q_prim_qp%vf(eqn_idx%c)%sf => q_cons_qp%vf(eqn_idx%c)%sf
            $:GPU_ENTER_DATA(copyin='[q_prim_qp%vf(eqn_idx%c)%sf]')
            $:GPU_ENTER_DATA(attach='[q_prim_qp%vf(eqn_idx%c)%sf]')
        end if

        if (hyper_cleaning) then
            q_prim_qp%vf(eqn_idx%psi)%sf => q_cons_qp%vf(eqn_idx%psi)%sf
            $:GPU_ENTER_DATA(copyin='[q_prim_qp%vf(eqn_idx%psi)%sf]')
            $:GPU_ENTER_DATA(attach='[q_prim_qp%vf(eqn_idx%psi)%sf]')
        end if

        if (.not. igr) then
            @:ALLOCATE(flux_n(1:num_dims))
            @:ALLOCATE(flux_src_n(1:num_dims))
            @:ALLOCATE(flux_gsrc_n(1:num_dims))

            do i = 1, num_dims
                @:ALLOCATE(flux_n(i)%vf(1:sys_size))
                @:ALLOCATE(flux_src_n(i)%vf(1:sys_size))
                @:ALLOCATE(flux_gsrc_n(i)%vf(1:sys_size))

                if (i == 1) then
#ifndef MFC_DIRECT_RIEMANN_FLUX
                    do l = 1, sys_size
                        @:ALLOCATE(flux_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do
#endif
                    do l = 1, sys_size
                        @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do

                    if (viscous .or. surface_tension .or. heat_conduction) then
                        do l = eqn_idx%mom%beg, eqn_idx%E
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if

#ifndef MFC_DIRECT_RIEMANN_FLUX
                    ! flux_src at adv%beg:adv%end is overloaded by design:
                    !
                    ! adv_src_mode_alpha_iface: flux_src(adv%beg:adv%end) = per-fluid interface alpha values. Each fluid
                    ! gets a separate 3D array.
                    !
                    ! adv_src_mode_vel_iface: flux_src(adv%beg) = one shared face-normal interface velocity.
                    ! adv%beg+1:adv%end are pointer-aliased to adv%beg so loops over adv%beg:adv%end can keep fluid indexing
                    ! while still reading one value. Saves (num_fluids - 1) fields.
                    !
                    ! adv_src_mode_none: flux_src(adv%beg:adv%end) = zeros. No separate NC advection source term. Allocated
                    ! for structural consistency with s_finalize_riemann_solver.
                    @:ALLOCATE(flux_src_n(i)%vf(eqn_idx%adv%beg)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))

                    if (adv_src_mode == adv_src_mode_alpha_iface .or. adv_src_mode == adv_src_mode_none) then
                        ! Alpha-interface needs separate per-fluid arrays. HLLD (adv_src_mode_none) allocates for structural
                        ! consistency with s_finalize_riemann_solver.
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if
#endif

                    if (chemistry) then
                        do l = eqn_idx%species%beg, eqn_idx%species%end
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                        if (chem_params%diffusion .and. .not. (viscous .or. surface_tension .or. heat_conduction)) then
                            @:ALLOCATE(flux_src_n(i)%vf(eqn_idx%E)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                       & idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
                        end if
                    end if
                else
                    do l = 1, sys_size
                        @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do
                end if

                @:ACC_SETUP_VFs(flux_n(i))
                @:ACC_SETUP_VFs(flux_src_n(i), flux_gsrc_n(i))

                if (i == 1) then
#ifndef MFC_DIRECT_RIEMANN_FLUX
                    if (adv_src_mode == adv_src_mode_vel_iface) then
                        ! u-interface: flux_src(adv%beg) holds one shared face-normal velocity. Pointer-alias adv%beg+1:adv%end to
                        ! the same memory so loops over adv%beg:adv%end can keep fluid indexing while still reading one value. This
                        ! saves (num_fluids - 1) 3D field allocations
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            flux_src_n(i)%vf(l)%sf => flux_src_n(i)%vf(eqn_idx%adv%beg)%sf
                            $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
                        end do
                    end if
#endif
                else
                    do l = 1, sys_size
#ifndef MFC_DIRECT_RIEMANN_FLUX
                        flux_n(i)%vf(l)%sf => flux_n(1)%vf(l)%sf
                        $:GPU_ENTER_DATA(attach='[flux_n(i)%vf(l)%sf]')
#endif
#ifdef MFC_DIRECT_RIEMANN_FLUX
                        if (l < eqn_idx%adv%beg .or. l > eqn_idx%adv%end) then
#endif
                            flux_src_n(i)%vf(l)%sf => flux_src_n(1)%vf(l)%sf
                            $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
#ifdef MFC_DIRECT_RIEMANN_FLUX
                        end if
#endif
                    end do
                end if
            end do
        end if

        if ((.not. igr)) then
            @:ALLOCATE(dq_prim_dx_qp(1:1))
            @:ALLOCATE(dq_prim_dy_qp(1:1))
            @:ALLOCATE(dq_prim_dz_qp(1:1))

            @:ALLOCATE(qL_prim(1:num_dims))
            @:ALLOCATE(qR_prim(1:num_dims))

            @:ALLOCATE(dqL_prim_dx_n(1:num_dims))
            @:ALLOCATE(dqL_prim_dy_n(1:num_dims))
            @:ALLOCATE(dqL_prim_dz_n(1:num_dims))
            @:ALLOCATE(dqR_prim_dx_n(1:num_dims))
            @:ALLOCATE(dqR_prim_dy_n(1:num_dims))
            @:ALLOCATE(dqR_prim_dz_n(1:num_dims))

            do i = 1, num_dims
                @:ALLOCATE(qL_prim(i)%vf(1:sys_size))
                @:ALLOCATE(qR_prim(i)%vf(1:sys_size))
                do l = eqn_idx%mom%beg, eqn_idx%mom%end
                    @:ALLOCATE(qL_prim(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ALLOCATE(qR_prim(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                end do
                @:ACC_SETUP_VFs(qL_prim(i), qR_prim(i))
            end do

            @:ALLOCATE(qL_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:sys_size))
            @:ALLOCATE(qR_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                       & 1:sys_size))

            if (.not. viscous) then
                do i = 1, num_dims
                    @:ALLOCATE(dqL_prim_dx_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqL_prim_dy_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqL_prim_dz_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dx_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dy_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dz_n(i)%vf(1:sys_size))

                    do l = eqn_idx%mom%beg, eqn_idx%mom%end
                        @:ALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                        @:ALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                        @:ALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                        @:ALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                        @:ALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                        @:ALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf(1:1, 1:1, 1:1))
                    end do
                    @:ACC_SETUP_VFs(dqL_prim_dx_n(i), dqL_prim_dy_n(i), dqL_prim_dz_n(i))
                    @:ACC_SETUP_VFs(dqR_prim_dx_n(i), dqR_prim_dy_n(i), dqR_prim_dz_n(i))
                end do
            end if

            ! The cylindrical axis source reads tau_Re_vf(mom%beg:E), so conduction needs it too,
            ! but without the viscous gradient fields below.
            if (viscous .or. heat_conduction) then
                @:ALLOCATE(tau_Re_vf(1:sys_size))
                do i = 1, num_dims
                    @:ALLOCATE(tau_Re_vf(eqn_idx%cont%end + i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(tau_Re_vf(eqn_idx%cont%end + i))
                end do
                @:ALLOCATE(tau_Re_vf(eqn_idx%E)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(tau_Re_vf(eqn_idx%E))
            end if

            if (viscous) then
                @:ALLOCATE(dq_prim_dx_qp(1)%vf(1:sys_size))
                @:ALLOCATE(dq_prim_dy_qp(1)%vf(1:sys_size))
                @:ALLOCATE(dq_prim_dz_qp(1)%vf(1:sys_size))

                do l = eqn_idx%mom%beg, eqn_idx%mom%end
                    @:ALLOCATE(dq_prim_dx_qp(1)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                end do

                @:ACC_SETUP_VFs(dq_prim_dx_qp(1))

                if (n > 0) then
                    do l = eqn_idx%mom%beg, eqn_idx%mom%end
                        @:ALLOCATE(dq_prim_dy_qp(1)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do

                    @:ACC_SETUP_VFs(dq_prim_dy_qp(1))

                    if (p > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:ALLOCATE(dq_prim_dz_qp(1)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                        @:ACC_SETUP_VFs(dq_prim_dz_qp(1))
                    end if
                end if

                do i = 1, num_dims
                    @:ALLOCATE(dqL_prim_dx_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqL_prim_dy_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqL_prim_dz_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dx_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dy_n(i)%vf(1:sys_size))
                    @:ALLOCATE(dqR_prim_dz_n(i)%vf(1:sys_size))
                end do

                do i = 1, num_dims
                    do l = eqn_idx%mom%beg, eqn_idx%mom%end
                        @:ALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                        @:ALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do

                    if (n > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:ALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                            @:ALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if

                    if (p > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:ALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                            @:ALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if

                    @:ACC_SETUP_VFs(dqL_prim_dx_n(i), dqL_prim_dy_n(i), dqL_prim_dz_n(i))
                    @:ACC_SETUP_VFs(dqR_prim_dx_n(i), dqR_prim_dy_n(i), dqR_prim_dz_n(i))
                end do

                if (weno_Re_flux) then
                    @:ALLOCATE(dqL_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end, eqn_idx%mom%beg:eqn_idx%mom%end))
                    @:ALLOCATE(dqR_rsx_vf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end, eqn_idx%mom%beg:eqn_idx%mom%end))
                end if
            else
                @:ALLOCATE(dq_prim_dx_qp(1)%vf(1:sys_size))
                @:ALLOCATE(dq_prim_dy_qp(1)%vf(1:sys_size))
                @:ALLOCATE(dq_prim_dz_qp(1)%vf(1:sys_size))

                do l = eqn_idx%mom%beg, eqn_idx%mom%end
                    @:ALLOCATE(dq_prim_dx_qp(1)%vf(l)%sf(0, 0, 0))
                    @:ACC_SETUP_VFs(dq_prim_dx_qp(1))
                    if (n > 0) then
                        @:ALLOCATE(dq_prim_dy_qp(1)%vf(l)%sf(0, 0, 0))
                        @:ACC_SETUP_VFs(dq_prim_dy_qp(1))
                        if (p > 0) then
                            @:ALLOCATE(dq_prim_dz_qp(1)%vf(l)%sf(0, 0, 0))
                            @:ACC_SETUP_VFs(dq_prim_dz_qp(1))
                        end if
                    end if
                end do
            end if

            $:GPU_PARALLEL_LOOP(private='[i, j, k, l, id]', collapse=4)
            do id = 1, num_dims
                do i = 1, sys_size
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                flux_gsrc_n(id)%vf(i)%sf(j, k, l) = 0._wp
                            end do
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (qbmm) then
            @:ALLOCATE(mom_sp(1:nmomsp), mom_3d(0:2, 0:2, nb))

            do i = 0, 2
                do j = 0, 2
                    do k = 1, nb
                        @:ALLOCATE(mom_3d(i, j, k)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                        @:ACC_SETUP_SFs(mom_3d(i, j, k))
                    end do
                end do
            end do

            do i = 1, nmomsp
                @:ALLOCATE(mom_sp(i)%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(mom_sp(i))
            end do
        end if

        if (mpp_lim .and. bubbles_euler) then
            @:ALLOCATE(alf_sum%sf(state_x_storage_beg:state_x_storage_end, idwbuff(2)%beg:idwbuff(2)%end, &
                       & idwbuff(3)%beg:idwbuff(3)%end))
        end if
        if (use_nc_iface_vel) then
            @:ALLOCATE(nc_iface_vel_n(1:num_dims))
            do i = 1, num_dims
                @:ALLOCATE(nc_iface_vel_n(i)%vf(1:num_dims))
                do l = 1, num_dims
                    @:ALLOCATE(nc_iface_vel_n(i)%vf(l)%sf( idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                end do
            end do
            if (num_dims == 2) then
                @:ACC_SETUP_VFs(nc_iface_vel_n(1), nc_iface_vel_n(2))
            else if (num_dims == 3) then
                @:ACC_SETUP_VFs(nc_iface_vel_n(1), nc_iface_vel_n(2), nc_iface_vel_n(3))
            else
                @:ACC_SETUP_VFs(nc_iface_vel_n(1))
            end if
            if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
                @:ALLOCATE(nc_iface_vel_hatR_n(1:num_dims))
                do i = 1, num_dims
                    @:ALLOCATE(nc_iface_vel_hatR_n(i)%vf(1:num_dims))
                    do l = 1, num_dims
                        @:ALLOCATE(nc_iface_vel_hatR_n(i)%vf(l)%sf( idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do
                end do
                if (num_dims == 2) then
                    @:ACC_SETUP_VFs(nc_iface_vel_hatR_n(1), nc_iface_vel_hatR_n(2))
                else if (num_dims == 3) then
                    @:ACC_SETUP_VFs(nc_iface_vel_hatR_n(1), nc_iface_vel_hatR_n(2), nc_iface_vel_hatR_n(3))
                else
                    @:ACC_SETUP_VFs(nc_iface_vel_hatR_n(1))
                end if
            end if
        end if

        if (alt_soundspeed) then
            @:ALLOCATE(blkmod1(0:m, 0:n, 0:p), blkmod2(0:m, 0:n, 0:p), alpha1(0:m, 0:n, 0:p), alpha2(0:m, 0:n, 0:p), Kterm(0:m, &
                       & 0:n, 0:p))
        end if

        if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
            @:ALLOCATE(rhs_hatL_vf(1:sys_size))
            @:ALLOCATE(rhs_hatR_vf(1:sys_size))
            @:PREFER_GPU(rhs_hatL_vf)
            @:PREFER_GPU(rhs_hatR_vf)

            do i = 1, sys_size
                @:ALLOCATE(rhs_hatL_vf(i)%sf(0:m, 0:n, 0:p))
                @:ALLOCATE(rhs_hatR_vf(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(rhs_hatL_vf(i))
                @:ACC_SETUP_SFs(rhs_hatR_vf(i))
            end do
        end if

        ! Spatial body force source arrays - sized to include ghost cells so the
        ! same indexing as q_*_vf is valid; iteration is restricted to interior
        ! cells via the `bounds` argument in m_body_forces. Only allocated when the
        ! feature is active (the write/read sites are likewise guarded).
        if (bf_spatial_support) then
            if (n > 0) then
                if (p > 0) then
                    @:ALLOCATE(spbf_source_x(-buff_size:buff_size + m, -buff_size:buff_size + n, -buff_size:buff_size + p))
                    @:ALLOCATE(spbf_source_y(-buff_size:buff_size + m, -buff_size:buff_size + n, -buff_size:buff_size + p))
                else
                    @:ALLOCATE(spbf_source_x(-buff_size:buff_size + m, -buff_size:buff_size + n, 0:0))
                    @:ALLOCATE(spbf_source_y(-buff_size:buff_size + m, -buff_size:buff_size + n, 0:0))
                end if
            else
                @:ALLOCATE(spbf_source_x(-buff_size:buff_size + m, 0:0, 0:0))
                @:ALLOCATE(spbf_source_y(-buff_size:buff_size + m, 0:0, 0:0))
            end if
            @:PREFER_GPU(spbf_source_x)
            @:PREFER_GPU(spbf_source_y)

            spbf_source_x = 0._wp
            spbf_source_y = 0._wp
            $:GPU_UPDATE(device='[spbf_source_x, spbf_source_y]')
        end if

    end subroutine s_initialize_rhs_module

    !> Associate per-direction flux fields with the contiguous Riemann output
    impure subroutine s_associate_rhs_flux_views

        ! The guard sits inside the routine: only m_rhs and m_riemann_state compile with
        ! MFC_DIRECT_RIEMANN_FLUX, so a guard at the call site would compile the call out.
#ifdef MFC_DIRECT_RIEMANN_FLUX
        integer :: i, l

        do i = 1, num_dims
            do l = 1, sys_size
                flux_n(i)%vf(l)%sf(-1:m,-1:n,-1:p) => flux_rsx_vf(:,:,:,l)
                $:GPU_ENTER_DATA(attach='[flux_n(i)%vf(l)%sf]')
            end do

            ! HLL Method 1 and HLLD map per-fluid rows; HLLC and HLL Method 2 (u-interface) share the adv%beg row
            if ((riemann_solver == riemann_solver_hll .and. .not. hll_u_interface) .or. riemann_solver == riemann_solver_hlld) then
                do l = eqn_idx%adv%beg, eqn_idx%adv%end
                    flux_src_n(i)%vf(l)%sf(-1:m,-1:n,-1:p) => flux_src_rsx_vf(:,:,:,l)
                    $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
                end do
            else
                do l = eqn_idx%adv%beg, eqn_idx%adv%end
                    flux_src_n(i)%vf(l)%sf(-1:m,-1:n,-1:p) => flux_src_rsx_vf(:,:,:,eqn_idx%adv%beg)
                    $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
                end do
            end if
        end do
#endif

    end subroutine s_associate_rhs_flux_views

    !> Compute the right-hand side of the semi-discrete governing equations for a single time stage
    #:set FUSED_RK_DUMMIES = ', primitive_ready, primitive_next_active, fuse_final_rk, q_stor_vf, q_output_vf, save_rk_state, ' &
        & + 'rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef' if MFC_FUSED_RHS else ''
    impure subroutine s_compute_rhs(q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                                    & stage${FUSED_RK_DUMMIES}$)

        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_cons_vf
        type(scalar_field), intent(inout)                                                          :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                 :: bc_type
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: rhs_vf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in

        ! TODO :: I think these other two variables need to be stp as well, but it doesn't compile like that right now
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout)  :: rhs_pb
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: mv_in
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout)  :: rhs_mv
        integer, intent(in)                                                                        :: t_step
        integer, intent(in)                                                                        :: stage

#if defined(MFC_FUSED_RHS)
        logical, intent(in)                                    :: primitive_ready, primitive_next_active, fuse_final_rk
        type(scalar_field), dimension(sys_size), intent(inout) :: q_stor_vf, q_output_vf
        logical, intent(in)                                    :: save_rk_state
        real(wp), intent(in)                                   :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef
        logical                                                :: copy_q_cons_qp
#endif
#ifdef MFC_CUDA_GRAPHS
        integer         :: rhs_graph_action, bubble_graph_action
        integer(kind=8) :: rhs_graph_signature, bubble_graph_signature
        logical         :: rhs_graph_capturable, transport_graph_capturable
        logical         :: rhs_graph_comm, dynamic_sources
#endif
        integer         :: id
        integer(kind=8) :: i, j, k, l, q             !< Generic loop iterators
        logical         :: fused_viscous_stress_ran  !< The fused kernel assembles and differences the viscous stress in this sweep

        ! RHS: halo exchange -> reconstruct -> Riemann solve -> flux difference -> source terms

        call nvtxStartRange("COMPUTE-RHS")

        fused_viscous_stress_ran = .false.
#ifdef MFC_FUSED_RHS
        native_halo_overlap_active = f_use_hllc_halo_overlap()
        if (cfl_dt) then
            if (mytime >= t_stop) native_halo_overlap_active = .false.
        else
            if (t_step == t_step_stop) native_halo_overlap_active = .false.
        end if
        #:if MFC_CASE_OPTIMIZATION and hllc_viscous_stress
            fused_viscous_stress_ran = viscous .and. f_use_fused_rhs_openacc()
        #:endif
        copy_q_cons_qp = .not. f_use_fused_rhs_openacc() .or. bubbles_euler .or. bubbles_lagrange &
            & .or. acoustic_source${' .or. relativity' if MFC_CASE_OPTIMIZATION and relativity else ''}$
#endif

        if (.not. igr) then
            ! Association/Population of Working Variables
#if defined(MFC_FUSED_RHS)
            if (.not. primitive_ready .and. copy_q_cons_qp) then
#endif
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do i = 1, sys_size
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                q_cons_qp%vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
#if defined(MFC_FUSED_RHS)
            end if
#endif

            ! Converting Conservative to Primitive Variables

            if (mpp_lim .and. bubbles_euler) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            alf_sum%sf(j, k, l) = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end - 1
                                alf_sum%sf(j, k, l) = alf_sum%sf(j, k, l) + q_cons_qp%vf(i)%sf(j, k, l)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%adv%beg, eqn_idx%adv%end - 1
                                q_cons_qp%vf(i)%sf(j, k, l) = q_cons_qp%vf(i)%sf(j, k, &
                                             & l)*(1._wp - q_cons_qp%vf(eqn_idx%alf)%sf(j, k, l))/alf_sum%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

#ifdef MFC_CUDA_GRAPHS
        ! Captured kernels read the permuted bounds as launch bounds, fixed per node at capture. A kernel that reads
        ! them in its body sets rhs_graph_capturable = .false.
        rhs_graph_action = graph_plain
        rhs_graph_capturable = .not. igr .or. num_procs == 1
        if (cfl_adap_dt) rhs_graph_capturable = .false.
        ! Lagrangian transport stays outside graph capture.
        if (bubbles_lagrange) rhs_graph_capturable = .false.
#if defined(MFC_FUSED_RHS) && defined(MFC_CBC_NO_HOST_SYNC)
        if (f_has_characteristic_bc() .and. .not. f_use_fused_rhs_openacc()) rhs_graph_capturable = .false.
#else
        if (f_has_characteristic_bc()) rhs_graph_capturable = .false.
#endif
#ifdef MFC_FUSED_RHS
        if (viscous .and. .not. f_use_fused_rhs_openacc()) rhs_graph_capturable = .false.
#else
        if (viscous) rhs_graph_capturable = .false.
#endif
        ! Dynamic source and diagnostic work stays outside the transport graph
        dynamic_sources = acoustic_source .or. probe_wrt .or. bubbles_lagrange .or. cont_damage
        transport_graph_capturable = rhs_graph_capturable .and. dynamic_sources .and. .not. igr
        if (dynamic_sources) rhs_graph_capturable = .false.
        ! With the device-resident NVSHMEM epoch exchange every halo launch is parameter-free, so the capture starts
        ! before the pack and the whole stage replays as one graph. Regions that run uncaptured keep the capture after
        ! the communication.
#if defined(MFC_HALO_SINGLE_ROUND) && defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
        rhs_graph_comm = rhs_graph_capturable .and. halo_26_epoch_device .and. num_procs > 1
#else
        rhs_graph_comm = .false.
#endif
        if (qbmm .or. chemistry .or. viscous .or. surface_tension .or. bubbles_euler) rhs_graph_comm = .false.
        if (int_comp == int_comp_mthinc .and. n > 0) rhs_graph_comm = .false.
        if (rhs_graph_capturable .or. transport_graph_capturable) then
            rhs_graph_signature = 0_8
            call s_graph_signature(rhs_graph_signature, int(stage, 8))
            if (igr) call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, t_step == t_step_start))
#if defined(MFC_FUSED_RHS)
            call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, native_halo_overlap_active))
            call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, primitive_next_active))
            call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, fuse_final_rk))
            call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, save_rk_state))
            call s_graph_signature(rhs_graph_signature, transfer(real(rk_q_coef, 8), 0_8))
            call s_graph_signature(rhs_graph_signature, transfer(real(rk_stor_coef, 8), 0_8))
            call s_graph_signature(rhs_graph_signature, transfer(real(rk_rhs_coef, 8), 0_8))
            call s_graph_signature(rhs_graph_signature, transfer(real(rk_den_coef, 8), 0_8))
#ifdef MFC_ACTIVE_BOX
            ! Engaging the window changes the captured window guards on every rank, so it forces a recapture
            call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, ab_active))
            if (ab_active) then
                call s_graph_signature(rhs_graph_signature, int(ab_x%beg, 8))
                call s_graph_signature(rhs_graph_signature, int(ab_x%end, 8))
                call s_graph_signature(rhs_graph_signature, int(ab_y%beg, 8))
                call s_graph_signature(rhs_graph_signature, int(ab_y%end, 8))
                call s_graph_signature(rhs_graph_signature, int(ab_z%beg, 8))
                call s_graph_signature(rhs_graph_signature, int(ab_z%end, 8))
            end if
#endif
#endif
            if (rhs_graph_comm) then
                ! Launch bounds and branch shape of the exchange; the trim tables are device arrays read at every replay
#if defined(MFC_HALO_SINGLE_ROUND) && defined(MFC_NVSHMEM) && defined(__NVCOMPILER)
                call s_graph_signature(rhs_graph_signature, int(halo_26_pack_points, 8))
                call s_graph_signature(rhs_graph_signature, int(halo_26_recv_points, 8))
                call s_graph_signature(rhs_graph_signature, int(halo_26_neighbor_count, 8))
                call s_graph_signature(rhs_graph_signature, merge(1_8, 0_8, halo_26_trim_active))
#endif
            end if
            if (hyper_cleaning) then
                call s_graph_signature(rhs_graph_signature, transfer(real(hyper_cleaning_speed, 8), 0_8))
                call s_graph_signature(rhs_graph_signature, transfer(real(hyper_cleaning_tau, 8), 0_8))
            end if
            if (rhs_graph_comm) rhs_graph_action = f_graph_action(graph_slot_rhs, rhs_graph_signature)
        end if
#endif
        $:GPU_GRAPH_SKIP_BEGIN('rhs_graph_action')
        if (igr) then
            call nvtxStartRange("RHS-COMMUNICATION")
            call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)
            call nvtxEndRange
        end if
        if (.not. igr) then
#if defined(MFC_FUSED_RHS)
            if (.not. primitive_ready) then
                call nvtxStartRange("RHS-CONVERT")
                if (copy_q_cons_qp) then
                    call s_convert_conservative_to_primitive_variables(q_cons_qp%vf, q_T_sf, q_prim_qp%vf, idwint)
                else
                    call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_qp%vf, idwint)
                end if
                call nvtxEndRange
            end if
#ifdef MFC_ACTIVE_BOX
            ! The fused kernel updates the primitive ping-pong arrays inside the active-box window, while column
            ! staging and the halo exchange read beyond it. The state outside the window is frozen, so one full
            ! conversion of both arrays on engagement keeps those reads valid.
            if (ab_active .and. .not. ab_prim_refreshed) then
                ab_prim_refreshed = .true.
                call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_qp%vf, idwint)
                call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_next%vf, idwint)
            end if
#endif
#else
            call nvtxStartRange("RHS-CONVERT")
            call s_convert_conservative_to_primitive_variables(q_cons_qp%vf, q_T_sf, q_prim_qp%vf, idwint)
            call nvtxEndRange
#endif

            call nvtxStartRange("RHS-COMMUNICATION")
#if defined(MFC_FUSED_RHS)
            if (primitive_next_active) then
                call s_populate_variables_buffers(bc_type, q_prim_next%vf, pb_in, mv_in, q_T_sf, merge(1, 0, &
                                                  & native_halo_overlap_active))
            else
                call s_populate_variables_buffers(bc_type, q_prim_qp%vf, pb_in, mv_in, q_T_sf, merge(1, 0, &
                                                  & native_halo_overlap_active))
            end if
#else
            call s_populate_variables_buffers(bc_type, q_prim_qp%vf, pb_in, mv_in, q_T_sf)
#endif
            call nvtxEndRange
#ifdef MFC_FUSED_RHS
            #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous)) and not viscous and not qbmm and not surface_tension and hllc_int_comp != 2
                if (native_halo_overlap_active) then
                    call nvtxStartRange("RHS-INTERIOR-X")
                    if (primitive_next_active) then
                        call s_weno_fused_rhs_openacc(q_prim_next%vf, rhs_vf, 1, q_cons_vf, q_stor_vf, q_output_vf, q_prim_qp%vf, &
                                                      & fuse_final_rk .and. num_dims == 1, save_rk_state, rk_q_coef, &
                                                      & rk_stor_coef, rk_rhs_coef, rk_den_coef, overlap_region=1)
                    else
                        call s_weno_fused_rhs_openacc(q_prim_qp%vf, rhs_vf, 1, q_cons_vf, q_stor_vf, q_output_vf, q_prim_next%vf, &
                                                      & fuse_final_rk .and. num_dims == 1, save_rk_state, rk_q_coef, &
                                                      & rk_stor_coef, rk_rhs_coef, rk_den_coef, overlap_region=1)
                    end if
                    call nvtxEndRange
                    call nvtxStartRange("RHS-COMMUNICATION-COMPLETE")
                    if (primitive_next_active) then
                        call s_populate_variables_buffers(bc_type, q_prim_next%vf, pb_in, mv_in, q_T_sf, 2)
                    else
                        call s_populate_variables_buffers(bc_type, q_prim_qp%vf, pb_in, mv_in, q_T_sf, 2)
                    end if
                    call nvtxEndRange
                end if
            #:endif
#endif
        end if

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (qbmm) call s_mom_inv(q_cons_qp%vf, q_prim_qp%vf, mom_sp, mom_3d, pb_in, rhs_pb, mv_in, rhs_mv, idwbuff(1), &
            & idwbuff(2), idwbuff(3))

        if ((viscous .and. .not. igr)) then
            call nvtxStartRange("RHS-VISCOUS")
#if defined(MFC_FUSED_RHS)
            ! With the primitive ping-pong active the gradients read q_prim_next, as the halo exchange above does
            #:if MFC_CASE_OPTIMIZATION and viscous and hllc_vel_gradients
                ! The fused kernels build dq_prim_d?_qp from the primitive state; the host gradient pass runs for
                ! the configurations the fused WENO path does not cover.
                if (f_use_fused_rhs_openacc() .and. (wenojs .or. mapped_weno .or. wenoz .or. teno) .and. (int_comp == 0 &
                    & .or. (num_fluids == 2 .and. (int_comp == int_comp_thinc .or. int_comp == int_comp_mthinc)))) then
                    if (primitive_next_active) then
                        call s_compute_fused_velocity_gradients(q_prim_next%vf, &
                                                                & ${', '.join(f'dq_prim_d{a}_qp(1)%vf' for a in 'xyz'[:num_dims])}$)
                    else
                        call s_compute_fused_velocity_gradients(q_prim_qp%vf, &
                                                                & ${', '.join(f'dq_prim_d{a}_qp(1)%vf' for a in 'xyz'[:num_dims])}$)
                    end if
                else if (primitive_next_active) then
                    call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                                       & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_next, dq_prim_dx_qp, dq_prim_dy_qp, &
                                       & dq_prim_dz_qp, idwbuff(1), idwbuff(2), idwbuff(3), fused_viscous_stress_ran)
                else
                    call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                                       & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_qp, dq_prim_dx_qp, dq_prim_dy_qp, &
                                       & dq_prim_dz_qp, idwbuff(1), idwbuff(2), idwbuff(3), fused_viscous_stress_ran)
                end if
            #:else
                if (primitive_next_active) then
                    call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                                       & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_next, dq_prim_dx_qp, dq_prim_dy_qp, &
                                       & dq_prim_dz_qp, idwbuff(1), idwbuff(2), idwbuff(3), fused_viscous_stress_ran)
                else
                    call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                                       & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_qp, dq_prim_dx_qp, dq_prim_dy_qp, &
                                       & dq_prim_dz_qp, idwbuff(1), idwbuff(2), idwbuff(3), fused_viscous_stress_ran)
                end if
            #:endif
#else
            call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                               & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_qp, dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                               & idwbuff(1), idwbuff(2), idwbuff(3), .false.)
#endif
            call nvtxEndRange
        end if

        if (surface_tension) then
            call nvtxStartRange("RHS-SURFACE-TENSION")
#if defined(MFC_FUSED_RHS)
            ! Read the same primitive buffer as s_get_viscous above
            if (primitive_next_active) then
                call s_get_capillary(q_prim_next%vf, bc_type)
            else
                call s_get_capillary(q_prim_qp%vf, bc_type)
            end if
#else
            call s_get_capillary(q_prim_qp%vf, bc_type)
#endif
            call nvtxEndRange
        end if

        if (int_comp == int_comp_mthinc .and. n > 0) then
            call nvtxStartRange("RHS-COMPRESSION-NORMALS")
#if defined(MFC_FUSED_RHS)
            if (primitive_next_active) then
                call s_compute_mthinc_normals(q_prim_next%vf)
            else
                call s_compute_mthinc_normals(q_prim_qp%vf)
            end if
#else
            call s_compute_mthinc_normals(q_prim_qp%vf)
#endif
            call nvtxEndRange
        end if

#ifdef MFC_CUDA_GRAPHS
        ! With rhs_graph_comm the capture began before the pack; otherwise it begins here, after the communication
        if (.not. rhs_graph_comm .and. rhs_graph_capturable) then
            rhs_graph_action = f_graph_action(graph_slot_rhs, rhs_graph_signature)
        end if
#endif
        $:GPU_GRAPH_SKIP_BEGIN('rhs_graph_action')
        $:GPU_CAPTURE_BEGIN(slot='graph_slot_rhs', signature='rhs_graph_signature', enabled='transport_graph_capturable')
        if (hypo_nc_mode /= hypo_nc_mode_dual_pass) then
            do id = 1, num_dims
                if (igr) then
                    if (id == 1) then
                        ! `i` outermost, matching the RK update that reads it
                        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                        do i = 1, sys_size
                            do l = -min(1, p), p + min(1, p)
                                do k = -1, n + 1
                                    do j = -1, m + 1
                                        rhs_vf(i)%sf(j, k, l) = 0._stp
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        #:if MFC_CASE_OPTIMIZATION and num_dims > 2 and defined('MFC_IGR_TUNED_SWEEPS') and MFC_IGR_TUNED_SWEEPS
                            ! Gradient divided differences for the IGR sweeps, computed once per RHS evaluation
                            call s_igr_precompute_dd(q_cons_vf)
                        #:endif
                    end if

                    call nvtxStartRange("IGR_RIEMANN")
                    call s_igr_riemann_solver(q_cons_vf, rhs_vf, id)
                    call nvtxEndRange

                    if (id == 1) then
                        call nvtxStartRange("IGR_Jacobi")
                        call s_igr_iterative_solve(q_cons_vf, bc_type, t_step)
                        call nvtxEndRange

                        call nvtxStartRange("IGR_SIGMA")
                        call s_igr_sigma_x(q_cons_vf, rhs_vf)
                        call nvtxEndRange
                    end if
                end if
                if (.not. igr) then
#ifdef MFC_FUSED_RHS
                    if (f_use_fused_rhs_openacc()) then
                        call s_fused_directional_rhs(id, rhs_vf, q_cons_vf, primitive_next_active, fuse_final_rk, q_stor_vf, &
                                                     & q_output_vf, save_rk_state, rk_q_coef, rk_stor_coef, rk_rhs_coef, &
                                                     & rk_den_coef)
                    else
#endif
                        call s_reconstruct_riemann_states(id)

                        call s_compute_directional_rhs(id, rhs_vf, .false.)
#ifdef MFC_FUSED_RHS
                    end if
#endif

                    ! RHS additions for hypoelasticity
                    if (hypo_nc_mode == hypo_nc_mode_finite_diff) then
                        call nvtxStartRange("RHS-HYPOELASTICITY-FD-PER-SWEEP")
                        call s_compute_hypoelastic_rhs_finite_diff_per_sweep(id, q_prim_qp%vf, rhs_vf)
                        call nvtxEndRange
                    end if

                    ! RHS for diffusion
                    if (chemistry .and. chem_params%diffusion) then
                        call nvtxStartRange("RHS-CHEM-DIFFUSION")
#ifdef MFC_FUSED_RHS
                        #:if MFC_CASE_OPTIMIZATION and chemistry and hllc_chemistry
                            ! The fused path leaves the energy and species rows of flux_src unset, so zero them first
                            if (f_use_fused_rhs_openacc() .and. (.not. viscous .or. fused_viscous_stress_ran)) then
                                call s_zero_fused_diffusion_rows(id)
                            end if
                        #:endif
#endif
                        call s_compute_chemistry_diffusion_flux(id, q_prim_qp%vf, flux_src_n(id)%vf, irx, iry, irz, q_T_sf)
                        call nvtxEndRange
                    end if

                    ! RHS for Fourier heat conduction
                    if (heat_conduction) then
                        call nvtxStartRange("RHS-CONDUCTION")
#ifdef MFC_FUSED_RHS
                        call s_compute_conduction_source_flux(id, q_prim_qp%vf, q_T_sf, flux_src_n(id)%vf, irx, iry, irz, &
                                                              & replace_e_row=f_use_fused_rhs_openacc())
#else
                        call s_compute_conduction_source_flux(id, q_prim_qp%vf, q_T_sf, flux_src_n(id)%vf, irx, iry, irz)
#endif
                        call nvtxEndRange
                    end if

                    ! Viscous stress contribution to RHS, unless the fused kernel already differenced it
                    if ((viscous .and. .not. fused_viscous_stress_ran) .or. surface_tension .or. chem_params%diffusion &
                        & .or. heat_conduction .or. (viscous .and. fused_viscous_stress_ran .and. id == 2 .and. cyl_coord &
                        & .and. ((bc_y%beg == -2) .or. (bc_y%beg == -14)))) then
                        call nvtxStartRange("RHS-ADD-PHYSICS")
                        call s_compute_additional_physics_rhs(id, q_prim_qp%vf, q_T_sf, rhs_vf, flux_src_n(id)%vf, &
                                                              & dq_prim_dx_qp(1)%vf, dq_prim_dy_qp(1)%vf, dq_prim_dz_qp(1)%vf, &
                                                              & fused_viscous_stress_ran)
                        call nvtxEndRange
                    end if

                    ! Bubble dynamics source terms
                    if (bubbles_euler .and. .not. qbmm .and. .not. f_use_fused_bubble_source()) then
                        call nvtxStartRange("RHS-BUBBLES-COMPUTE")
                        call s_compute_bubbles_EE_rhs(id, q_prim_qp%vf, divu)
                        call nvtxEndRange
                    end if

                    ! RHS additions for qbmm bubbles
                    if (qbmm) then
                        call nvtxStartRange("RHS-QBMM")
                        call s_compute_qbmm_rhs(id, q_cons_qp%vf, q_prim_qp%vf, rhs_vf, flux_n(id)%vf, pb_in, rhs_pb)
                        call nvtxEndRange
                    end if
                    ! END: Additional physics and source terms

                    if (hyper_cleaning) then
                        $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                        do l = 0, p
                            do k = 0, n
                                do j = 0, m
                                    rhs_vf(eqn_idx%psi)%sf(j, k, l) = rhs_vf(eqn_idx%psi)%sf(j, k, &
                                           & l) - q_prim_vf(eqn_idx%psi)%sf(j, k, l)/hyper_cleaning_tau
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if

                    ! END: Additional physics and source terms
                end if
            end do
        else
            ! Fused dual-pass HLLD: per direction, the Riemann states are reconstructed once and ONE fused solve computes both
            ! anchored flux sets (hat_L into flux_n/nc_iface_vel_n through the solver's own finalize; hat_R into the flux_hatR_rs*
            ! buffers). The hat_L partial RHS is assembled first; the hat_R set is then finalized into the same (already consumed)
            ! flux field arrays and assembled. The full RHS is the sum of the two anchored partial RHS's; the axisymmetric
            ! cylindrical completion for the volume fractions and stresses is applied once on the summed RHS, using the cell's own
            ! two anchored radial face traces (hat_L outer face, hat_R inner face).
            do id = 1, num_dims
                call s_reconstruct_riemann_states(id)
                call s_compute_directional_rhs(id, rhs_hatL_vf, .true.)
                call nvtxStartRange("RHS-RIEMANN-FINALIZE-HATR")
                call s_finalize_riemann_solver_hatR(flux_n(id)%vf, flux_gsrc_n(id)%vf, id)
                if (use_nc_iface_vel) then
                    call s_finalize_nc_iface_vel_hatR(nc_iface_vel_hatR_n(id)%vf, id)
                end if
                call nvtxEndRange
                call nvtxStartRange("RHS-ADVECTION-SRC")
                call s_compute_advection_source_term(id, rhs_hatR_vf, q_cons_qp, q_prim_qp, flux_src_n(id), .false.)
                call nvtxEndRange
            end do

            $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(i)%sf(j, k, l) = rhs_hatL_vf(i)%sf(j, k, l) + rhs_hatR_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

            if (grid_geometry == 2) then
                call nvtxStartRange("RHS-HYPOELASTICITY-AXISYM-HLLD")
                call s_compute_hypoelastic_rhs_axisym_geom_dual_pass(q_prim_qp%vf, rhs_vf, nc_iface_vel_n(2)%vf, &
                    & nc_iface_vel_hatR_n(2)%vf)
                call nvtxEndRange
            end if
        end if
        ! END: Dimensional Splitting Loop

        ! RHS additions for hypoelasticity (interface-consistent path, after all sweeps)
        if (hypo_nc_mode == hypo_nc_mode_interface) then
            call nvtxStartRange("RHS-HYPOELASTICITY-IFACE")
            call s_compute_hypoelastic_rhs_iface(q_prim_qp%vf, rhs_vf, nc_iface_vel_n)
            call nvtxEndRange
        end if
        if (ib) then
            $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ib_markers%sf(j, k, l) /= 0) then
                            do i = 1, sys_size
                                rhs_vf(i)%sf(j, k, l) = 0._wp
                            end do
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_CAPTURE_END(slot='graph_slot_rhs')

        ! Additional Physics and Source Terms Additions for acoustic_source
        if (acoustic_source) then
            call nvtxStartRange("RHS-ACOUSTIC-SRC")
            call s_acoustic_src_calculations(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf)
            call nvtxEndRange
        end if

        ! Add bubbles source term
        if (bubbles_euler .and. (.not. adap_dt) .and. (.not. qbmm)) then
            call nvtxStartRange("RHS-BUBBLES-SRC")
#ifdef MFC_CUDA_GRAPHS
            bubble_graph_action = graph_plain
            if (.not. rhs_graph_capturable .and. f_use_fused_bubble_source()) then
                bubble_graph_signature = 0_8
                call s_graph_signature(bubble_graph_signature, int(bubble_model, 8))
                call s_graph_signature(bubble_graph_signature, int(nb, 8))
                call s_graph_signature(bubble_graph_signature, merge(1_8, 0_8, polytropic))
                bubble_graph_action = f_graph_action(graph_slot_bubble_ee, bubble_graph_signature)
            end if
#endif
            $:GPU_GRAPH_SKIP_BEGIN('bubble_graph_action')
            call s_compute_bubble_EE_source(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf, divu)
#ifdef MFC_CUDA_GRAPHS
            if (bubble_graph_action == graph_capture) then
                if (f_graph_commit(graph_slot_bubble_ee) /= 0) &
                    & call s_mpi_abort('CUDA-graph capture of the bubble source failed.  Rebuild without MFC_CUDA_GRAPHS.')
            end if
#endif
            $:GPU_GRAPH_SKIP_END()
            call nvtxEndRange
        end if

        ! ptil is a probe-output diagnostic only, so evaluate it only when probes are on.
        ! This sits outside the source-term branch above so that the qbmm and adap_dt
        ! paths, which skip that branch, still get a filled ptil.
        if (bubbles_euler .and. probe_wrt) then
            call s_compute_ptilde(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size))
        end if

        if (bubbles_lagrange) then
            if (.not. adap_dt) then
                call nvtxStartRange("RHS-EL-BUBBLES-DYN")
                call s_compute_bubble_EL_dynamics(q_prim_qp%vf(1:sys_size), bc_type, stage)
                call nvtxEndRange
            end if

            if (lag_params%solver_approach == 2) then
                call nvtxStartRange("RHS-EL-BUBBLES-SRC")
                call s_compute_bubbles_EL_source(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf)
                call nvtxEndRange
            end if
        end if

        ! When reaction_substeps > 0 the reaction source is integrated by operator splitting
        ! after the flow update (s_chemistry_reaction_substep), not added to the flow RHS here.
        if (chemistry .and. chem_params%reactions .and. chem_params%reaction_substeps == 0) then
            call nvtxStartRange("RHS-CHEM-REACTIONS")
            call s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp%vf, q_T_sf, q_prim_qp%vf, idwint)
            call nvtxEndRange
        end if

        ! With rburn%substeps > 0 the burn is integrated by operator splitting after the flow
        ! update (s_reactive_burn_substep), not added to the flow RHS here.
        if (reactive_burn .and. rburn%substeps == 0) then
#if defined(MFC_FUSED_RHS)
            ! The fused final RK update applies the burn in-kernel
            if (.not. fuse_final_rk) then
#endif
                call nvtxStartRange("RHS-REACTIVE-BURN")
                call s_compute_reactive_burn(rhs_vf, q_cons_qp%vf, q_prim_qp%vf, idwint)
                call nvtxEndRange
#if defined(MFC_FUSED_RHS)
            end if
#endif
        end if

#if defined(MFC_FUSED_RHS)
        ! The fused path skips the q_cons_qp refresh, so the damage source reads the stage state and the
        ! primitive buffer this stage's sweep read
        if (primitive_next_active) then
            if (cont_damage) call s_compute_damage_state(q_cons_vf, q_prim_next%vf, rhs_vf)
        else
            if (cont_damage) call s_compute_damage_state(q_cons_vf, q_prim_qp%vf, rhs_vf)
        end if
#else
        if (cont_damage) call s_compute_damage_state(q_cons_qp%vf, q_prim_qp%vf, rhs_vf)
#endif

        ! END: Additional physics and source terms

        if (run_time_info .or. probe_wrt .or. ib .or. bubbles_lagrange) then
            if (.not. igr) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do i = 1, sys_size
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
#if defined(MFC_FUSED_RHS)
                                ! Diagnostics read the primitive buffer this stage's sweep read
                                if (primitive_next_active) then
                                    q_prim_vf(i)%sf(j, k, l) = q_prim_next%vf(i)%sf(j, k, l)
                                else
                                    q_prim_vf(i)%sf(j, k, l) = q_prim_qp%vf(i)%sf(j, k, l)
                                end if
#else
                                q_prim_vf(i)%sf(j, k, l) = q_prim_qp%vf(i)%sf(j, k, l)
#endif
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

#ifdef MFC_CUDA_GRAPHS
        if (rhs_graph_action == graph_capture) then
            if (f_graph_commit(graph_slot_rhs) /= 0) &
                & call s_mpi_abort('CUDA-graph capture of the RHS device region failed.  Rebuild without MFC_CUDA_GRAPHS.')
        end if
#endif
        $:GPU_GRAPH_SKIP_END()
        $:GPU_GRAPH_SKIP_END()

        call nvtxEndRange

    end subroutine s_compute_rhs

    !> Reconstructs the left/right Riemann states (cell-boundary values) for one sweep direction into the module face buffers
    !! (qL/qR_rs*_vf and, for WENO-reconstructed viscous fluxes, dqL/dqR_rs*_vf). Depends only on q_prim_qp, so it runs once per
    !! sweep direction; both anchored passes of the fused dual-pass HLLD solve in that direction reuse these buffers.
    subroutine s_reconstruct_riemann_states(id)

        integer, intent(in) :: id
        integer             :: i, j, k, l  !< Generic loop iterators

        call nvtxStartRange("RHS-RECONSTRUCTION")

        if (.not. surface_tension) then
            if ((.not. weno_Re_flux) .or. int_comp > 0) then
                ! Reconstruct densitiess
                iv%beg = 1; iv%end = sys_size
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(1:sys_size), qL_rsx_vf, qR_rsx_vf, id)
            else
                iv%beg = 1; iv%end = eqn_idx%cont%end
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)

                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = eqn_idx%mom%beg, eqn_idx%mom%end
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                qL_rsx_vf(j, k, l, i) = qL_prim(id)%vf(i)%sf(j, k, l)
                                qR_rsx_vf(j, k, l, i) = qR_prim(id)%vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                iv%beg = eqn_idx%E; iv%end = sys_size
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)
            end if
        else
            if (int_comp > 0) then
                ! THINC reads cont and adv from v_rs_ws; must reconstruct full sys_size range to populate both
                iv%beg = 1; iv%end = sys_size
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(1:sys_size), qL_rsx_vf, qR_rsx_vf, id)
                ! Surface tension requires first-order energy; overwrite the higher-order result from the full pass above
                iv%beg = eqn_idx%E; iv%end = eqn_idx%E
                call s_reconstruct_cell_boundary_values_first_order(q_prim_qp%vf(eqn_idx%E), qL_rsx_vf, qR_rsx_vf, id)
            else if ((.not. weno_Re_flux)) then
                iv%beg = 1; iv%end = eqn_idx%E - 1
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)

                iv%beg = eqn_idx%E; iv%end = eqn_idx%E
                call s_reconstruct_cell_boundary_values_first_order(q_prim_qp%vf(eqn_idx%E), qL_rsx_vf, qR_rsx_vf, id)

                iv%beg = eqn_idx%E + 1; iv%end = sys_size
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)
            else
                iv%beg = 1; iv%end = eqn_idx%cont%end
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)

                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = eqn_idx%mom%beg, eqn_idx%mom%end
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                qL_rsx_vf(j, k, l, i) = qL_prim(id)%vf(i)%sf(j, k, l)
                                qR_rsx_vf(j, k, l, i) = qR_prim(id)%vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()

                iv%beg = eqn_idx%E; iv%end = eqn_idx%E
                call s_reconstruct_cell_boundary_values_first_order(q_prim_qp%vf(eqn_idx%E), qL_rsx_vf, qR_rsx_vf, id)

                iv%beg = eqn_idx%E + 1; iv%end = sys_size
                call s_reconstruct_cell_boundary_values(q_prim_qp%vf(iv%beg:iv%end), qL_rsx_vf, qR_rsx_vf, id)
            end if
        end if

        ! Reconstruct viscous derivatives for viscosity
        if (weno_Re_flux) then
            iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end
            call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dx_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, dqR_rsx_vf, id, &
                & dqL_prim_dx_n(id)%vf(iv%beg:iv%end), dqR_prim_dx_n(id)%vf(iv%beg:iv%end), idwbuff(1), idwbuff(2), idwbuff(3))
            if (n > 0) then
                call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dy_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, dqR_rsx_vf, &
                    & id, dqL_prim_dy_n(id)%vf(iv%beg:iv%end), dqR_prim_dy_n(id)%vf(iv%beg:iv%end), idwbuff(1), idwbuff(2), &
                    & idwbuff(3))
                if (p > 0) then
                    call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dz_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, &
                        & dqR_rsx_vf, id, dqL_prim_dz_n(id)%vf(iv%beg:iv%end), dqR_prim_dz_n(id)%vf(iv%beg:iv%end), idwbuff(1), &
                        & idwbuff(2), idwbuff(3))
                end if
            end if
        end if

        call nvtxEndRange

    end subroutine s_reconstruct_riemann_states

    !> Computes one sweep direction's contribution to the RHS: the Riemann solve on the reconstructed cell-boundary states followed
    !! by the advection source term. For dual-pass HLLD the fused solve computes BOTH anchored flux sets in one call; this routine
    !! assembles the hat_L (is_hat_L=.true.) partial RHS, and the caller finalizes + assembles the hat_R set afterwards. is_hat_L
    !! selects which one-sided flux difference the advection source applies into rhs_vf.
    subroutine s_compute_directional_rhs(id, rhs_vf, is_hat_L)

        integer, intent(in)                                    :: id
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        logical, intent(in)                                    :: is_hat_L

        ! Configuring Coordinate Direction Indexes

        if (id == 1) then
            irx%beg = -1; iry%beg = 0; irz%beg = 0
        else if (id == 2) then
            irx%beg = 0; iry%beg = -1; irz%beg = 0
        else
            irx%beg = 0; iry%beg = 0; irz%beg = -1
        end if
        irx%end = m; iry%end = n; irz%end = p

        ! Computing Riemann Solver Flux and Source Flux
        call nvtxStartRange("RHS-RIEMANN-SOLVER")
        call s_riemann_solver(qR_rsx_vf, dqR_prim_dx_n(id)%vf, dqR_prim_dy_n(id)%vf, dqR_prim_dz_n(id)%vf, qR_prim(id)%vf, &
                              & qL_rsx_vf, dqL_prim_dx_n(id)%vf, dqL_prim_dy_n(id)%vf, dqL_prim_dz_n(id)%vf, qL_prim(id)%vf, &
                              & q_prim_qp%vf, flux_n(id)%vf, flux_src_n(id)%vf, flux_gsrc_n(id)%vf, id, irx, iry, irz)
        call nvtxEndRange

        if (use_nc_iface_vel) then
            call s_finalize_nc_iface_vel(nc_iface_vel_n(id)%vf, id)
        end if

        ! Additional physics and source terms RHS addition for advection source
        call nvtxStartRange("RHS-ADVECTION-SRC")
        call s_compute_advection_source_term(id, rhs_vf, q_cons_qp, q_prim_qp, flux_src_n(id), is_hat_L)
        call nvtxEndRange

    end subroutine s_compute_directional_rhs

#ifdef MFC_FUSED_RHS
    !> One fused sweep direction: the reconstruction the kernel does not do itself, the fused Riemann solve and flux difference, and
    !! the characteristic-boundary band
    subroutine s_fused_directional_rhs(id, rhs_vf, q_cons_vf, primitive_next_active, fuse_final_rk, q_stor_vf, q_output_vf, &
                                       & save_rk_state, rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef)

        integer, intent(in)                                    :: id
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf, q_cons_vf, q_stor_vf, q_output_vf
        logical, intent(in)                                    :: primitive_next_active, fuse_final_rk, save_rk_state
        real(wp), intent(in)                                   :: rk_q_coef, rk_stor_coef, rk_rhs_coef, rk_den_coef

        #:if MFC_CASE_OPTIMIZATION and 1 <= num_dims <= 3 and 1 <= num_fluids <= 4 and (model_eqns == 2 or (model_eqns == 3 and hllc_riemann == 2 and hllc_low_mach == 0 and not hllc_alt_soundspeed)) and (not bubbles_euler or hllc_bubbles) and (not chemistry or hllc_chemistry) and (not mhd or hllc_mhd) and (not relativity or hllc_mhd) and (not hypoelasticity or hllc_hypoelastic) and (not surface_tension or (hllc_riemann == 2 and not viscous))
            logical :: use_fused_rhs_openacc, use_weno_fused_rhs_openacc

            use_fused_rhs_openacc = f_use_fused_rhs_openacc()
            use_weno_fused_rhs_openacc = use_fused_rhs_openacc .and. ((recon_type == recon_type_weno .and. (weno_order == 1 &
                & .or. weno_order == 3 .or. weno_order == 5 .or. weno_order == 7) .and. (wenojs .or. mapped_weno .or. wenoz &
                & .or. teno)) .or. (recon_type == recon_type_muscl .and. (muscl_order == 1 .or. muscl_order == 2))) &
                & .and. (.not. weno_Re_flux .or. viscous) .and. (int_comp == 0 .or. (num_fluids == 2 &
                & .and. (int_comp == int_comp_thinc .or. int_comp == int_comp_mthinc)))

            if (.not. use_weno_fused_rhs_openacc) call s_reconstruct_riemann_states(id)

            ! Configuring Coordinate Direction Indexes
            if (id == 1) then
                irx%beg = -1; iry%beg = 0; irz%beg = 0
            else if (id == 2) then
                irx%beg = 0; iry%beg = -1; irz%beg = 0
            else
                irx%beg = 0; iry%beg = 0; irz%beg = -1
            end if
            irx%end = m; iry%end = n; irz%end = p

            call nvtxStartRange("RHS-RIEMANN-SOLVER")
            #:if viscous
                if (use_fused_rhs_openacc .and. weno_Re_flux) then
                    #:if hllc_grad_weno
                        ! The kernel WENO-reconstructs the face gradients in-column; the non-WENO fused path reads
                        ! the averaged host face gradients
                        #:if hllc_viscous_stress
                            if (.not. use_weno_fused_rhs_openacc) call s_average_fused_face_gradients(id)
                        #:endif
                    #:else
                        if (use_weno_fused_rhs_openacc) call s_reconstruct_viscous_derivatives(id)
                        #:if hllc_viscous_stress
                            call s_average_fused_face_gradients(id)
                        #:endif
                    #:endif
                end if
            #:endif
            if (use_weno_fused_rhs_openacc) then
                if (primitive_next_active) then
                    call s_weno_fused_rhs_openacc(q_prim_next%vf, rhs_vf, id, q_cons_vf, q_stor_vf, q_output_vf, q_prim_qp%vf, &
                                                  & fuse_final_rk .and. id == num_dims, save_rk_state, rk_q_coef, rk_stor_coef, &
                                                  & rk_rhs_coef, &
                                                  & rk_den_coef${', ' + ', '.join((f'dq_prim_d{a}_qp(1)%vf' if hllc_grad_weno else f'dqR_prim_d{a}_n(id)%vf') for a in 'xyz'[:num_dims]) if viscous and hllc_viscous_stress else ''}$, &
                                                  & overlap_region=merge(2, 0, native_halo_overlap_active .and. id == 1))
                else
                    call s_weno_fused_rhs_openacc(q_prim_qp%vf, rhs_vf, id, q_cons_vf, q_stor_vf, q_output_vf, q_prim_next%vf, &
                                                  & fuse_final_rk .and. id == num_dims, save_rk_state, rk_q_coef, rk_stor_coef, &
                                                  & rk_rhs_coef, &
                                                  & rk_den_coef${', ' + ', '.join((f'dq_prim_d{a}_qp(1)%vf' if hllc_grad_weno else f'dqR_prim_d{a}_n(id)%vf') for a in 'xyz'[:num_dims]) if viscous and hllc_viscous_stress else ''}$, &
                                                  & overlap_region=merge(2, 0, native_halo_overlap_active .and. id == 1))
                end if
            else
                call s_fused_rhs_openacc(q_prim_qp%vf, rhs_vf, id, &
                                         & ${(', '.join(f'dqR_prim_d{a}_n(id)%vf' for a in 'xyz'[:num_dims]) + ', ') if viscous and hllc_viscous_stress else ''}$qL_rsx_vf, &
                                         & qR_rsx_vf)
            end if
            #:if surface_tension
                if (use_fused_rhs_openacc) then
                    call s_fused_source_flux_state(flux_src_n(id)%vf, id, irx, iry, irz)
                    call s_compute_capillary_source_flux(vel_src_rsx_vf, flux_src_n(id)%vf, id, isx, isy, isz)
                end if
            #:endif
            #:if viscous and not hllc_viscous_stress
                if (use_fused_rhs_openacc) then
                    if (primitive_next_active) then
                        call s_fused_viscous_source_flux(flux_src_n(id)%vf, q_prim_next%vf, qR_prim(id)%vf, qL_prim(id)%vf, &
                                                         & dqR_prim_dx_n(id)%vf, dqR_prim_dy_n(id)%vf, dqR_prim_dz_n(id)%vf, &
                                                         & dqL_prim_dx_n(id)%vf, dqL_prim_dy_n(id)%vf, dqL_prim_dz_n(id)%vf, id, &
                                                         & irx, iry, irz)
                    else
                        call s_fused_viscous_source_flux(flux_src_n(id)%vf, q_prim_qp%vf, qR_prim(id)%vf, qL_prim(id)%vf, &
                                                         & dqR_prim_dx_n(id)%vf, dqR_prim_dy_n(id)%vf, dqR_prim_dz_n(id)%vf, &
                                                         & dqL_prim_dx_n(id)%vf, dqL_prim_dy_n(id)%vf, dqL_prim_dz_n(id)%vf, id, &
                                                         & irx, iry, irz)
                    end if
                end if
            #:endif
            call nvtxEndRange

            call nvtxStartRange("RHS-ADVECTION-SRC")
            if (primitive_next_active) then
                call s_compute_fused_advection_boundary(id, rhs_vf, q_prim_next%vf, flux_src_n(id))
            else
                call s_compute_fused_advection_boundary(id, rhs_vf, q_prim_qp%vf, flux_src_n(id))
            end if
            #:if hllc_cbc and model_eqns == 2
                if (fuse_final_rk .and. id == num_dims) then
                    if (primitive_next_active) then
                        call s_fused_rhs_openacc_boundary_rk(q_cons_vf, q_stor_vf, q_output_vf, q_prim_qp%vf, q_prim_next%vf, &
                                                             & rhs_vf, save_rk_state, rk_q_coef, rk_stor_coef, rk_rhs_coef, &
                                                             & rk_den_coef)
                    else
                        call s_fused_rhs_openacc_boundary_rk(q_cons_vf, q_stor_vf, q_output_vf, q_prim_next%vf, q_prim_qp%vf, &
                                                             & rhs_vf, save_rk_state, rk_q_coef, rk_stor_coef, rk_rhs_coef, &
                                                             & rk_den_coef)
                    end if
                end if
            #:endif
            call nvtxEndRange
        #:endif

    end subroutine s_fused_directional_rhs
#endif

#if defined(MFC_FUSED_RHS) && ${int(MFC_CASE_OPTIMIZATION and chemistry and hllc_chemistry)}$
    !> Zero the energy and species rows of flux_src that the fused path leaves unset ahead of the diffusion flux
    subroutine s_zero_fused_diffusion_rows(id)

        integer, intent(in) :: id
        integer             :: i, j, k, l, ixb, ixe, iyb, iye, izb, ize

        ixb = irx%beg; ixe = irx%end
        iyb = iry%beg; iye = iry%end
        izb = irz%beg; ize = irz%end
        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, j, k, l]')
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe
                    flux_src_n(id)%vf(eqn_idx%E)%sf(j, k, l) = 0._stp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        flux_src_n(id)%vf(i)%sf(j, k, l) = 0._stp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_zero_fused_diffusion_rows
#endif

#if defined(MFC_FUSED_RHS)
    !> Apply CBCs and finish only the normal cells intentionally omitted by the fused kernel.
    subroutine s_compute_fused_advection_boundary(idir, rhs_vf, q_prim_vf, flux_src_n_vf)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf, q_prim_vf
        type(vector_field), intent(inout)                      :: flux_src_n_vf
        logical                                                :: cbc_beg, cbc_end
        integer                                                :: extent

        select case (idir)
        case (1)
            cbc_beg = bc_x%beg <= BC_CHAR_SLIP_WALL .and. bc_x%beg >= BC_CHAR_SUP_OUTFLOW
            cbc_end = bc_x%end <= BC_CHAR_SLIP_WALL .and. bc_x%end >= BC_CHAR_SUP_OUTFLOW
            extent = m
        case (2)
            cbc_beg = bc_y%beg <= BC_CHAR_SLIP_WALL .and. bc_y%beg >= BC_CHAR_SUP_OUTFLOW
            cbc_end = bc_y%end <= BC_CHAR_SLIP_WALL .and. bc_y%end >= BC_CHAR_SUP_OUTFLOW
            extent = n
        case default
            cbc_beg = bc_z%beg <= BC_CHAR_SLIP_WALL .and. bc_z%beg >= BC_CHAR_SUP_OUTFLOW
            cbc_end = bc_z%end <= BC_CHAR_SLIP_WALL .and. bc_z%end >= BC_CHAR_SUP_OUTFLOW
            extent = p
        end select

        if (.not. cbc_beg .and. .not. cbc_end) return

#ifndef MFC_DIRECT_RIEMANN_FLUX
        if (cbc_beg) call s_copy_fused_boundary_fluxes(idir, -1, buff_size + 1, flux_src_n_vf)
        if (cbc_end) call s_copy_fused_boundary_fluxes(idir, extent - buff_size - 2, extent, flux_src_n_vf)
#endif
        if (cbc_beg) call s_cbc(q_prim_vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
        if (cbc_end) call s_cbc(q_prim_vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
        if (cbc_beg) call s_update_fused_boundary_band(idir, 0, buff_size + 1, rhs_vf, q_prim_vf, flux_src_n_vf)
        if (cbc_end) call s_update_fused_boundary_band(idir, extent - buff_size - 1, extent, rhs_vf, q_prim_vf, flux_src_n_vf)

    end subroutine s_compute_fused_advection_boundary

#ifndef MFC_DIRECT_RIEMANN_FLUX
    !> Publish the native column's face buffers to the separate CBC fields.
    subroutine s_copy_fused_boundary_fluxes(idir, face_beg, face_end, flux_src_n_vf)

        integer, intent(in)               :: idir, face_beg, face_end
        type(vector_field), intent(inout) :: flux_src_n_vf
        integer                           :: i, j, k, l, jlo, jhi, klo, khi, llo, lhi, source_end

        jlo = 0; jhi = m
        klo = 0; khi = n
        llo = 0; lhi = p
        select case (idir)
        case (1)
            jlo = face_beg; jhi = face_end
        case (2)
            klo = face_beg; khi = face_end
        case (3)
            llo = face_beg; lhi = face_end
        end select
        source_end = eqn_idx%adv%beg
        if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) source_end = eqn_idx%adv%end

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do l = llo, lhi
                do k = klo, khi
                    do j = jlo, jhi
                        flux_n(idir)%vf(i)%sf(j, k, l) = flux_rsx_vf(j, k, l, i)
                        if (i >= eqn_idx%adv%beg .and. i <= source_end) flux_src_n_vf%vf(i)%sf(j, k, l) = flux_src_rsx_vf(j, k, &
                            & l, i)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_copy_fused_boundary_fluxes
#endif

    !> Accumulate a compact normal band from CBC-adjusted face fluxes.
    subroutine s_update_fused_boundary_band(idir, normal_beg, normal_end, rhs_vf, q_source_vf, flux_src_n_vf)

        integer, intent(in) :: idir, normal_beg, normal_end
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_source_vf
        type(vector_field), intent(in) :: flux_src_n_vf
        integer :: i, j, k, l
        integer :: jlo, jhi, klo, khi, llo, lhi
        real(wp) :: inv_ds, flux_left, flux_right, source_left, source_right, directional_rhs

        #:if MFC_CASE_OPTIMIZATION and hllc_alt_soundspeed
            real(wp) :: source_alpha, kapila_term
        #:endif
        #:set KAPILA = MFC_CASE_OPTIMIZATION and hllc_alt_soundspeed
        #:set SOURCE_ALPHA = 'source_alpha' if KAPILA else 'q_source_vf(i)%sf(j, k, l)'
        ! Inside the active-box window only, matching the boundary RK update that reads it
        jlo = 0; jhi = m
        klo = 0; khi = n
        llo = 0; lhi = p
#ifdef MFC_ACTIVE_BOX
        if (ab_active) then
            jlo = max(jlo, ab_x%beg); jhi = min(jhi, ab_x%end)
            klo = max(klo, ab_y%beg); khi = min(khi, ab_y%end)
            llo = max(llo, ab_z%beg); lhi = min(lhi, ab_z%end)
            ! Skip the launch for an empty transverse window
            if (jlo > jhi .or. klo > khi .or. llo > lhi) return
        end if
#endif
        #:def fused_boundary_alpha()
            #:if KAPILA
                source_alpha = 0._wp
                if (i == eqn_idx%adv%beg .or. i == eqn_idx%adv%end) then
                    kapila_term = f_hllc_kapila(real(q_source_vf(eqn_idx%E)%sf(j, k, l), wp), &
                                                & real(q_source_vf(eqn_idx%adv%beg)%sf(j, k, l), wp), &
                                                & real(q_source_vf(eqn_idx%adv%end)%sf(j, k, l), wp), &
                                                & real(q_source_vf(eqn_idx%cont%beg)%sf(j, k, l), wp), &
                                                & real(q_source_vf(eqn_idx%cont%end)%sf(j, k, l), wp), gammas(1), gammas(2), &
                                                & pi_infs(1), pi_infs(2))
                    source_alpha = real(q_source_vf(i)%sf(j, k, l), wp) + merge(kapila_term, -kapila_term, i == eqn_idx%adv%beg)
                end if
            #:endif
        #:enddef

        select case (idir)
        case (1)
            $:GPU_PARALLEL_LOOP(collapse=4, &
                                & private='[i, j, k, l, inv_ds, flux_left, flux_right, source_left, source_right, directional_rhs' + (', source_alpha, kapila_term' if KAPILA else '') + ']')
            do i = 1, sys_size
                do l = llo, lhi
                    do k = klo, khi
                        do j = normal_beg, normal_end
                            inv_ds = 1._wp/dx(j)
                            flux_left = flux_n(1)%vf(i)%sf(j - 1, k, l)
                            flux_right = flux_n(1)%vf(i)%sf(j, k, l)
                            directional_rhs = inv_ds*(flux_left - flux_right)
                            ! Advection source from per-fluid rows for HLL and HLLD and from the shared adv%beg row otherwise
                            if (i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) then
                                if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
                                    source_left = flux_src_n_vf%vf(i)%sf(j - 1, k, l)
                                    source_right = flux_src_n_vf%vf(i)%sf(j, k, l)
                                    directional_rhs = directional_rhs + inv_ds*q_source_vf(eqn_idx%cont%end + 1)%sf(j, k, &
                                        & l)*(source_left - source_right)
                                else
                                    $:fused_boundary_alpha()
                                    source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j - 1, k, l)
                                    source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                    directional_rhs = directional_rhs + inv_ds*${SOURCE_ALPHA}$*(source_right - source_left)
                                end if
                            end if
                            ! Non-conservative pressure work of the 6-equation model
                            if (model_eqns == model_eqns_6eq .and. i >= eqn_idx%int_en%beg .and. i <= eqn_idx%int_en%end) then
                                source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j - 1, k, l)
                                source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                directional_rhs = directional_rhs - inv_ds*q_source_vf(i - eqn_idx%int_en%beg &
                                    & + eqn_idx%adv%beg)%sf(j, k, l)*q_source_vf(eqn_idx%E)%sf(j, k, l)*(source_right - source_left)
                            end if
                            #:if MFC_CASE_OPTIMIZATION and hllc_viscous_stress
                                ! The kernel wrote the viscous divergence of the momentum and energy rows of these
                                ! cells; the other rows take the directional RHS directly
                                if (i >= eqn_idx%mom%beg .and. i <= eqn_idx%E) then
                                    rhs_vf(i)%sf(j, k, l) = real(real(rhs_vf(i)%sf(j, k, l), wp) + directional_rhs, stp)
                                else
                                    rhs_vf(i)%sf(j, k, l) = directional_rhs
                                end if
                            #:else
                                rhs_vf(i)%sf(j, k, l) = directional_rhs
                            #:endif
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        case (2)
            $:GPU_PARALLEL_LOOP(collapse=4, &
                                & private='[i, j, k, l, inv_ds, flux_left, flux_right, source_left, source_right' &
                                & + (', source_alpha, kapila_term' if KAPILA else '') + ']')
            do i = 1, sys_size
                do l = llo, lhi
                    do k = normal_beg, normal_end
                        do j = jlo, jhi
                            inv_ds = 1._wp/dy(k)
                            flux_left = flux_n(2)%vf(i)%sf(j, k - 1, l)
                            flux_right = flux_n(2)%vf(i)%sf(j, k, l)
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + inv_ds*(flux_left - flux_right)
                            ! Same three-way split as the x branch above; see the comment there.
                            if (i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) then
                                if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
                                    source_left = flux_src_n_vf%vf(i)%sf(j, k - 1, l)
                                    source_right = flux_src_n_vf%vf(i)%sf(j, k, l)
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                           & l) + inv_ds*q_source_vf(eqn_idx%cont%end + 2)%sf(j, k, l)*(source_left - source_right)
                                else
                                    $:fused_boundary_alpha()
                                    source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k - 1, l)
                                    source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                           & l) + inv_ds*${SOURCE_ALPHA}$*(source_right - source_left)
                                    #:if KAPILA and cyl_coord
                                        if (i == eqn_idx%adv%beg .or. i == eqn_idx%adv%end) then
                                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + merge(kapila_term, -kapila_term, &
                                                   & i == eqn_idx%adv%beg)/(2._wp*y_cc(k))*(source_right + source_left)
                                        end if
                                    #:endif
                                end if
                            end if
                            if (model_eqns == model_eqns_6eq .and. i >= eqn_idx%int_en%beg .and. i <= eqn_idx%int_en%end) then
                                source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k - 1, l)
                                source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                       & l) - inv_ds*q_source_vf(i - eqn_idx%int_en%beg + eqn_idx%adv%beg)%sf(j, k, &
                                       & l)*q_source_vf(eqn_idx%E)%sf(j, k, l)*(source_right - source_left)
                                if (cyl_coord) then
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                           & l) - 5.e-1_wp/y_cc(k)*q_source_vf(i - eqn_idx%int_en%beg + eqn_idx%adv%beg)%sf(j, k, &
                                           & l)*q_source_vf(eqn_idx%E)%sf(j, k, l)*(source_right + source_left)
                                end if
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        case default
            $:GPU_PARALLEL_LOOP(collapse=4, &
                                & private='[i, j, k, l, inv_ds, flux_left, flux_right, source_left, source_right' &
                                & + (', source_alpha, kapila_term' if KAPILA else '') + ']')
            do i = 1, sys_size
                do l = normal_beg, normal_end
                    do k = klo, khi
                        do j = jlo, jhi
                            inv_ds = 1._wp/dz(l)
                            flux_left = flux_n(3)%vf(i)%sf(j, k, l - 1)
                            flux_right = flux_n(3)%vf(i)%sf(j, k, l)
                            rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + inv_ds*(flux_left - flux_right)
                            ! Same three-way split as the x branch above; see the comment there.
                            if (i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) then
                                if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
                                    source_left = flux_src_n_vf%vf(i)%sf(j, k, l - 1)
                                    source_right = flux_src_n_vf%vf(i)%sf(j, k, l)
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                           & l) + inv_ds*q_source_vf(eqn_idx%cont%end + 3)%sf(j, k, l)*(source_left - source_right)
                                else
                                    $:fused_boundary_alpha()
                                    source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l - 1)
                                    source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                           & l) + inv_ds*${SOURCE_ALPHA}$*(source_right - source_left)
                                end if
                            end if
                            if (model_eqns == model_eqns_6eq .and. i >= eqn_idx%int_en%beg .and. i <= eqn_idx%int_en%end) then
                                source_left = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l - 1)
                                source_right = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(j, k, l)
                                rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, &
                                       & l) - inv_ds*q_source_vf(i - eqn_idx%int_en%beg + eqn_idx%adv%beg)%sf(j, k, &
                                       & l)*q_source_vf(eqn_idx%E)%sf(j, k, l)*(source_right - source_left)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end select

    end subroutine s_update_fused_boundary_band
#endif

    !> Accumulate advection source contributions from a given coordinate direction into the RHS
    subroutine s_compute_advection_source_term(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, is_hat_L)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(vector_field), intent(inout) :: q_cons_vf
        type(vector_field), intent(inout) :: q_prim_vf
        type(vector_field), intent(inout) :: flux_src_n_vf
        logical, intent(in) :: is_hat_L
        integer :: j, k, l, q              !< Loop iterators from original, meaning varies
        integer :: k_loop, l_loop, q_loop  !< Standardized spatial loop iterators 0:m, 0:n, 0:p
        integer :: i_fluid_loop
        real(wp) :: inv_ds, inv_y, flux_face1, flux_face2
        real(wp) :: advected_qty_val, pressure_val
        real(wp) :: G1_eff, G2_eff, pres_K, alpha_K, alpha_rho_K, blkmod_K

        G1_eff = 0._wp
        G2_eff = 0._wp
        if (hypoelasticity) then
            G1_eff = fluid_pp(1)%G
            G2_eff = fluid_pp(2)%G
        end if

        if (alt_soundspeed) then
            $:GPU_PARALLEL_LOOP(private='[k_loop, l_loop, q_loop, pres_K, alpha_K, alpha_rho_K, blkmod_K]', collapse=3)
            do q_loop = 0, p
                do l_loop = 0, n
                    do k_loop = 0, m
                        ! Scalars in, scalar out: an element of a device-resident array passed by reference to a
                        ! device routine is read or written at the wrong address on Cray OpenACC.
                        pres_K = q_prim_vf%vf(eqn_idx%E)%sf(k_loop, l_loop, q_loop)
                        alpha_K = q_prim_vf%vf(eqn_idx%adv%beg)%sf(k_loop, l_loop, q_loop)
                        alpha_rho_K = q_prim_vf%vf(eqn_idx%cont%beg)%sf(k_loop, l_loop, q_loop)
                        call s_phase_bulk_modulus(pres_K, alpha_K, alpha_rho_K, 1, blkmod_K)
                        blkmod1(k_loop, l_loop, q_loop) = blkmod_K + (4._wp/3._wp)*G1_eff
                        alpha_K = q_prim_vf%vf(eqn_idx%adv%end)%sf(k_loop, l_loop, q_loop)
                        alpha_rho_K = q_prim_vf%vf(eqn_idx%cont%end)%sf(k_loop, l_loop, q_loop)
                        call s_phase_bulk_modulus(pres_K, alpha_K, alpha_rho_K, 2, blkmod_K)
                        blkmod2(k_loop, l_loop, q_loop) = blkmod_K + (4._wp/3._wp)*G2_eff
                        alpha1(k_loop, l_loop, q_loop) = q_cons_vf%vf(eqn_idx%adv%beg)%sf(k_loop, l_loop, q_loop)

                        if (bubbles_euler) then
                            alpha2(k_loop, l_loop, q_loop) = q_cons_vf%vf(eqn_idx%alf - 1)%sf(k_loop, l_loop, q_loop)
                        else
                            alpha2(k_loop, l_loop, q_loop) = q_cons_vf%vf(eqn_idx%adv%end)%sf(k_loop, l_loop, q_loop)
                        end if

                        Kterm(k_loop, l_loop, q_loop) = alpha1(k_loop, l_loop, q_loop)*alpha2(k_loop, l_loop, &
                              & q_loop)*(blkmod2(k_loop, l_loop, q_loop) - blkmod1(k_loop, l_loop, q_loop))/(alpha1(k_loop, &
                              & l_loop, q_loop)*blkmod2(k_loop, l_loop, q_loop) + alpha2(k_loop, l_loop, q_loop)*blkmod1(k_loop, &
                              & l_loop, q_loop))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        select case (idir)
        case (1)  ! x-direction
            if (bc_x%beg <= BC_CHAR_SLIP_WALL .and. bc_x%beg >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
            end if
            if (bc_x%end <= BC_CHAR_SLIP_WALL .and. bc_x%end >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
            end if

            if (hypo_nc_mode /= hypo_nc_mode_dual_pass) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face1, flux_face2]')
                do q_loop = 0, p
                    do l_loop = 0, n
                        do k_loop = 0, m
                            inv_ds = 1._wp/dx(k_loop)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_n(1)%vf(j)%sf(k_loop - 1, l_loop, q_loop)
                                flux_face2 = flux_n(1)%vf(j)%sf(k_loop, l_loop, q_loop)
                                rhs_vf(j)%sf(k_loop, l_loop, q_loop) = inv_ds*(flux_face1 - flux_face2)
#ifdef MFC_DIRECT_RIEMANN_FLUX
                                if (model_eqns == model_eqns_5eq .and. .not. cyl_coord .and. .not. alt_soundspeed &
                                    & .and. j >= eqn_idx%adv%beg .and. j <= eqn_idx%adv%end) then
                                    ! HLLC and HLL Method 2 share one face-normal velocity in flux_src(adv%beg)
                                    if (riemann_solver == riemann_solver_hllc .or. (riemann_solver == riemann_solver_hll &
                                        & .and. hll_u_interface)) then
                                        flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(k_loop, l_loop, q_loop)
                                        flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(k_loop - 1, l_loop, q_loop)
                                        rhs_vf(j)%sf(k_loop, l_loop, q_loop) = rhs_vf(j)%sf(k_loop, l_loop, &
                                               & q_loop) + inv_ds*q_cons_vf%vf(j)%sf(k_loop, l_loop, &
                                               & q_loop)*(flux_face1 - flux_face2)
                                    else if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
                                        flux_face1 = flux_src_n_vf%vf(j)%sf(k_loop - 1, l_loop, q_loop)
                                        flux_face2 = flux_src_n_vf%vf(j)%sf(k_loop, l_loop, q_loop)
                                        rhs_vf(j)%sf(k_loop, l_loop, q_loop) = rhs_vf(j)%sf(k_loop, l_loop, &
                                               & q_loop) + inv_ds*q_prim_vf%vf(eqn_idx%cont%end + 1)%sf(k_loop, l_loop, &
                                               & q_loop)*(flux_face1 - flux_face2)
                                    end if
                                end if
#endif
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (is_hat_L) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face2]')
                do q_loop = 0, p
                    do l_loop = 0, n
                        do k_loop = 0, m
                            inv_ds = 1._wp/dx(k_loop)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face2 = flux_n(1)%vf(j)%sf(k_loop, l_loop, q_loop)
                                rhs_vf(j)%sf(k_loop, l_loop, q_loop) = -inv_ds*flux_face2
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! is_hat_R
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face1]')
                do q_loop = 0, p
                    do l_loop = 0, n
                        do k_loop = 0, m
                            inv_ds = 1._wp/dx(k_loop)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_n(1)%vf(j)%sf(k_loop - 1, l_loop, q_loop)
                                rhs_vf(j)%sf(k_loop, l_loop, q_loop) = inv_ds*flux_face1
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (model_eqns == model_eqns_6eq) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[i_fluid_loop, k_loop, l_loop, q_loop, inv_ds, advected_qty_val, &
                                    & pressure_val, flux_face1, flux_face2]')
                do q_loop = 0, p
                    do l_loop = 0, n
                        do k_loop = 0, m
                            do i_fluid_loop = 1, num_fluids
                                inv_ds = 1._wp/dx(k_loop)
                                advected_qty_val = q_cons_vf%vf(i_fluid_loop + eqn_idx%adv%beg - 1)%sf(k_loop, l_loop, q_loop)
                                pressure_val = q_prim_vf%vf(eqn_idx%E)%sf(k_loop, l_loop, q_loop)
                                flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(k_loop, l_loop, q_loop)
                                flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(k_loop - 1, l_loop, q_loop)
                                rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(k_loop, l_loop, &
                                       & q_loop) = rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(k_loop, l_loop, &
                                       & q_loop) - inv_ds*advected_qty_val*pressure_val*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

#ifdef MFC_DIRECT_RIEMANN_FLUX
            ! Skipped when the flux-difference loop above already folded the advection source
            if (model_eqns /= model_eqns_5eq .or. cyl_coord .or. alt_soundspeed .or. (riemann_solver /= riemann_solver_hllc &
                & .and. riemann_solver /= riemann_solver_hll .and. riemann_solver /= riemann_solver_hlld)) then
                call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
            end if
#else
            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
#endif
        case (2)  ! y-direction
            if (bc_y%beg <= BC_CHAR_SLIP_WALL .and. bc_y%beg >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
            end if
            if (bc_y%end <= BC_CHAR_SLIP_WALL .and. bc_y%end >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
            end if

            if (hypo_nc_mode /= hypo_nc_mode_dual_pass) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face1, flux_face2]')
                do l = 0, p
                    do k = 0, n
                        do q = 0, m
                            inv_ds = 1._wp/dy(k)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_n(2)%vf(j)%sf(q, k - 1, l)
                                flux_face2 = flux_n(2)%vf(j)%sf(q, k, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) + inv_ds*(flux_face1 - flux_face2)
#ifdef MFC_DIRECT_RIEMANN_FLUX
                                if (model_eqns == model_eqns_5eq .and. .not. cyl_coord .and. .not. alt_soundspeed &
                                    & .and. j >= eqn_idx%adv%beg .and. j <= eqn_idx%adv%end) then
                                    ! HLLC and HLL Method 2 share one face-normal velocity in flux_src(adv%beg)
                                    if (riemann_solver == riemann_solver_hllc .or. (riemann_solver == riemann_solver_hll &
                                        & .and. hll_u_interface)) then
                                        flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(q, k, l)
                                        flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(q, k - 1, l)
                                        rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) + inv_ds*q_cons_vf%vf(j)%sf(q, k, &
                                               & l)*(flux_face1 - flux_face2)
                                    else if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) then
                                        flux_face1 = flux_src_n_vf%vf(j)%sf(q, k - 1, l)
                                        flux_face2 = flux_src_n_vf%vf(j)%sf(q, k, l)
                                        rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, &
                                               & l) + inv_ds*q_prim_vf%vf(eqn_idx%cont%end + 2)%sf(q, k, &
                                               & l)*(flux_face1 - flux_face2)
                                    end if
                                end if
#endif
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (is_hat_L) then
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face2]')
                do l = 0, p
                    do k = 0, n
                        do q = 0, m
                            inv_ds = 1._wp/dy(k)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face2 = flux_n(2)%vf(j)%sf(q, k, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - inv_ds*flux_face2
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! is_hat_R
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face1]')
                do l = 0, p
                    do k = 0, n
                        do q = 0, m
                            inv_ds = 1._wp/dy(k)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_n(2)%vf(j)%sf(q, k - 1, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) + inv_ds*flux_face1
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (model_eqns == model_eqns_6eq) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[i_fluid_loop, k, l, q, inv_ds, advected_qty_val, pressure_val, &
                                    & flux_face1, flux_face2]')
                do l = 0, p
                    do k = 0, n
                        do q = 0, m
                            do i_fluid_loop = 1, num_fluids
                                inv_ds = 1._wp/dy(k)
                                advected_qty_val = q_cons_vf%vf(i_fluid_loop + eqn_idx%adv%beg - 1)%sf(q, k, l)
                                pressure_val = q_prim_vf%vf(eqn_idx%E)%sf(q, k, l)
                                flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(q, k, l)
                                flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(q, k - 1, l)
                                rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(q, k, &
                                       & l) = rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(q, k, &
                                       & l) - inv_ds*advected_qty_val*pressure_val*(flux_face1 - flux_face2)
                                if (cyl_coord) then
                                    rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(q, k, &
                                           & l) = rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(q, k, &
                                           & l) - 5.e-1_wp/y_cc(k)*advected_qty_val*pressure_val*(flux_face1 + flux_face2)
                                end if
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (cyl_coord) then
                if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
                    ! Cell-owned anchored geometry: the hat_L pass owns the cell's outer radial face and the hat_R pass its
                    ! inner face (the same anchor ownership as the one-sided flux differences above), so the summed RHS
                    ! carries -(H_outer + H_inner)/(2*y_cc) with both faces anchored to the cell itself. With midpoint cell
                    ! centers this equals the annular finite-volume divergence (r_+ F_+ - r_- F_-)/(y_cc*dy) for the rows
                    ! whose geometric flux is the full radial flux.
                    if (is_hat_L) then
                        $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, flux_face2]')
                        do j = 1, sys_size
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        flux_face2 = flux_gsrc_n(2)%vf(j)%sf(q, k, l)
                                        rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - 5.e-1_wp/y_cc(k)*flux_face2
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    else  ! is_hat_R
                        $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, flux_face1]')
                        do j = 1, sys_size
                            do l = 0, p
                                do k = 0, n
                                    do q = 0, m
                                        flux_face1 = flux_gsrc_n(2)%vf(j)%sf(q, k - 1, l)
                                        rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - 5.e-1_wp/y_cc(k)*flux_face1
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                else
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, flux_face1, flux_face2]')
                    do j = 1, sys_size
                        do l = 0, p
                            do k = 0, n
                                do q = 0, m
                                    flux_face1 = flux_gsrc_n(2)%vf(j)%sf(q, k - 1, l)
                                    flux_face2 = flux_gsrc_n(2)%vf(j)%sf(q, k, l)
                                    rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - 5.e-1_wp/y_cc(k)*(flux_face1 + flux_face2)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

#ifdef MFC_DIRECT_RIEMANN_FLUX
            ! Skipped when the flux-difference loop above already folded the advection source
            if (model_eqns /= model_eqns_5eq .or. cyl_coord .or. alt_soundspeed .or. (riemann_solver /= riemann_solver_hllc &
                & .and. riemann_solver /= riemann_solver_hll .and. riemann_solver /= riemann_solver_hlld)) then
                call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
            end if
#else
            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
#endif
        case (3)  ! z-direction
            if (bc_z%beg <= BC_CHAR_SLIP_WALL .and. bc_z%beg >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
            end if
            if (bc_z%end <= BC_CHAR_SLIP_WALL .and. bc_z%end >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
            end if

            if (grid_geometry == 3) then  ! Cylindrical Coordinates
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face1, flux_face2]')
                do k = 0, p
                    do q = 0, n
                        do l = 0, m
                            inv_ds = 1._wp/(dz(k)*y_cc(q))
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_n(3)%vf(j)%sf(l, q, k - 1)
                                flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_y, flux_face1, flux_face2]')
                do k = 0, p
                    do q = 0, n
                        do l = 0, m
                            inv_y = 5.e-1_wp/y_cc(q)
                            $:GPU_LOOP(parallelism='[seq]')
                            do j = 1, sys_size
                                flux_face1 = flux_gsrc_n(3)%vf(j)%sf(l, q, k - 1)
                                flux_face2 = flux_gsrc_n(3)%vf(j)%sf(l, q, k)
                                rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) - inv_y*(flux_face1 + flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! Cartesian Coordinates
                if (hypo_nc_mode /= hypo_nc_mode_dual_pass) then
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face1, flux_face2]')
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                inv_ds = 1._wp/dz(k)
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, sys_size
                                    flux_face1 = flux_n(3)%vf(j)%sf(l, q, k - 1)
                                    flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                    rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*(flux_face1 - flux_face2)
#ifdef MFC_DIRECT_RIEMANN_FLUX
                                    if (model_eqns == model_eqns_5eq .and. .not. cyl_coord .and. .not. alt_soundspeed &
                                        & .and. j >= eqn_idx%adv%beg .and. j <= eqn_idx%adv%end) then
                                        ! HLLC and HLL Method 2 share one face-normal velocity in flux_src(adv%beg)
                                        if (riemann_solver == riemann_solver_hllc .or. (riemann_solver == riemann_solver_hll &
                                            & .and. hll_u_interface)) then
                                            flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(l, q, k)
                                            flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(l, q, k - 1)
                                            rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*q_cons_vf%vf(j)%sf(l, q, &
                                                   & k)*(flux_face1 - flux_face2)
                                        else if (riemann_solver == riemann_solver_hll .or. riemann_solver == riemann_solver_hlld) &
                                                 & then
                                            flux_face1 = flux_src_n_vf%vf(j)%sf(l, q, k - 1)
                                            flux_face2 = flux_src_n_vf%vf(j)%sf(l, q, k)
                                            rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, &
                                                   & k) + inv_ds*q_prim_vf%vf(eqn_idx%cont%end + 3)%sf(l, q, &
                                                   & k)*(flux_face1 - flux_face2)
                                        end if
                                    end if
#endif
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (is_hat_L) then
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face2]')
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                inv_ds = 1._wp/dz(k)
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, sys_size
                                    flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                    rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) - inv_ds*flux_face2
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else  ! is_hat_R
                    $:GPU_PARALLEL_LOOP(collapse=3,private='[j, k, l, q, inv_ds, flux_face1]')
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                inv_ds = 1._wp/dz(k)
                                $:GPU_LOOP(parallelism='[seq]')
                                do j = 1, sys_size
                                    flux_face1 = flux_n(3)%vf(j)%sf(l, q, k - 1)
                                    rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*flux_face1
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

            if (model_eqns == model_eqns_6eq) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[i_fluid_loop, k, l, q, inv_ds, advected_qty_val, pressure_val, &
                                    & flux_face1, flux_face2]')
                do k = 0, p
                    do q = 0, n
                        do l = 0, m
                            do i_fluid_loop = 1, num_fluids
                                inv_ds = 1._wp/dz(k)
                                advected_qty_val = q_cons_vf%vf(i_fluid_loop + eqn_idx%adv%beg - 1)%sf(l, q, k)
                                pressure_val = q_prim_vf%vf(eqn_idx%E)%sf(l, q, k)
                                flux_face1 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(l, q, k)
                                flux_face2 = flux_src_n_vf%vf(eqn_idx%adv%beg)%sf(l, q, k - 1)
                                rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(l, q, &
                                       & k) = rhs_vf(i_fluid_loop + eqn_idx%int_en%beg - 1)%sf(l, q, &
                                       & k) - inv_ds*advected_qty_val*pressure_val*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

#ifdef MFC_DIRECT_RIEMANN_FLUX
            ! Skipped when the flux-difference loop above already folded the advection source
            if (model_eqns /= model_eqns_5eq .or. cyl_coord .or. alt_soundspeed .or. (riemann_solver /= riemann_solver_hllc &
                & .and. riemann_solver /= riemann_solver_hll .and. riemann_solver /= riemann_solver_hlld)) then
                call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
            end if
#else
            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
#endif
        end select

    contains

        !> Add the advection source flux-difference terms for a single coordinate direction to the RHS
        subroutine s_add_directional_advection_source_terms(current_idir, rhs_vf_arg, q_cons_vf_arg, q_prim_vf_arg, &
            & flux_src_n_vf_arg, Kterm_arg)
            integer, intent(in)                                    :: current_idir
            type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf_arg
            type(vector_field), intent(in)                         :: q_cons_vf_arg
            type(vector_field), intent(in)                         :: q_prim_vf_arg
            type(vector_field), intent(in)                         :: flux_src_n_vf_arg
            real(wp), allocatable, dimension(:,:,:), intent(in)    :: Kterm_arg
            integer                                                :: j_adv, k_idx, l_idx, q_idx
            real(wp)                                               :: local_inv_ds, local_term_coeff, local_flux1, local_flux2
            real(wp)                                               :: local_k_term_val

            ! Three mutually exclusive modes for the NC volume fraction advection source term:

            !   adv_src_mode_alpha_iface - HLL Method 1: flux_src carries per-fluid interface alpha_k

            !   adv_src_mode_vel_iface   - HLLC, HLL Method 2, Exact, LF: flux_src carries shared face-normal velocity

            !   adv_src_mode_none        - HLLD: MHD has no volume fractions; hypo uses dual-pass Riemann flux

            select case (current_idir)
            case (1)  ! x-direction
                if (adv_src_mode == adv_src_mode_alpha_iface) then
                    ! Alpha-interface: flux_src(j_adv) supplies interface alpha_k. RHS applies velocity * d(alpha_k)/dx.
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                        & local_flux1, local_flux2]')
                    do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                        do q_idx = 0, p  ! z_extent
                            do l_idx = 0, n  ! y_extent
                                do k_idx = 0, m  ! x_extent
                                    local_inv_ds = 1._wp/dx(k_idx)
                                    local_term_coeff = q_prim_vf_arg%vf(eqn_idx%cont%end + current_idir)%sf(k_idx, l_idx, q_idx)
                                    local_flux1 = flux_src_n_vf_arg%vf(j_adv)%sf(k_idx - 1, l_idx, q_idx)
                                    local_flux2 = flux_src_n_vf_arg%vf(j_adv)%sf(k_idx, l_idx, q_idx)
                                    rhs_vf_arg(j_adv)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(j_adv)%sf(k_idx, l_idx, &
                                               & q_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    if (alt_soundspeed) then  ! K*div(u) correction
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                            & local_flux1, local_flux2]')
                        do q_idx = 0, p; do l_idx = 0, n; do k_idx = 0, m
                            local_inv_ds = 1._wp/dx(k_idx)
                            local_k_term_val = Kterm_arg(k_idx, l_idx, q_idx)
                            local_flux1 = nc_iface_vel_n(1)%vf(1)%sf(k_idx, l_idx, q_idx)
                            local_flux2 = nc_iface_vel_n(1)%vf(1)%sf(k_idx - 1, l_idx, q_idx)
                            rhs_vf_arg(eqn_idx%adv%beg)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(k_idx, l_idx, &
                                       & q_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            rhs_vf_arg(eqn_idx%adv%end)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(k_idx, l_idx, &
                                       & q_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                        end do; end do; end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                else if (adv_src_mode == adv_src_mode_vel_iface) then
                    ! u-interface: flux_src(adv%beg) supplies one shared face-normal velocity. RHS applies alpha_k * du/dx.
                    $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                        & local_flux1, local_flux2]')
                    do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                        do q_idx = 0, p; do l_idx = 0, n; do k_idx = 0, m
                            local_inv_ds = 1._wp/dx(k_idx)
                            local_term_coeff = q_cons_vf_arg%vf(j_adv)%sf(k_idx, l_idx, q_idx)
                            local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(k_idx, l_idx, q_idx)
                            local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(k_idx - 1, l_idx, q_idx)
                            rhs_vf_arg(j_adv)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(j_adv)%sf(k_idx, l_idx, &
                                       & q_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                        end do; end do; end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    if (alt_soundspeed) then  ! K*div(u) correction
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                            & local_flux1, local_flux2]')
                        do q_idx = 0, p; do l_idx = 0, n; do k_idx = 0, m
                            local_inv_ds = 1._wp/dx(k_idx)
                            local_k_term_val = Kterm_arg(k_idx, l_idx, q_idx)
                            local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(k_idx, l_idx, q_idx)
                            local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(k_idx - 1, l_idx, q_idx)
                            rhs_vf_arg(eqn_idx%adv%beg)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(k_idx, l_idx, &
                                       & q_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            rhs_vf_arg(eqn_idx%adv%end)%sf(k_idx, l_idx, q_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(k_idx, l_idx, &
                                       & q_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                        end do; end do; end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                end if
            case (2)  ! y-direction
                if (adv_src_mode == adv_src_mode_alpha_iface) then
                    ! Alpha-interface: flux_src(j_adv) supplies interface alpha_k. RHS applies velocity * d(alpha_k)/dy.
                    $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                        & local_flux1, local_flux2]')
                    do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                        do l_idx = 0, p  ! z_extent
                            do k_idx = 0, n  ! y_extent
                                do q_idx = 0, m  ! x_extent
                                    local_inv_ds = 1._wp/dy(k_idx)
                                    local_term_coeff = q_prim_vf_arg%vf(eqn_idx%cont%end + current_idir)%sf(q_idx, k_idx, l_idx)
                                    local_flux1 = flux_src_n_vf_arg%vf(j_adv)%sf(q_idx, k_idx - 1, l_idx)
                                    local_flux2 = flux_src_n_vf_arg%vf(j_adv)%sf(q_idx, k_idx, l_idx)
                                    rhs_vf_arg(j_adv)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(j_adv)%sf(q_idx, k_idx, &
                                               & l_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    if (alt_soundspeed) then  ! K*div(u) correction
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                            & local_flux1, local_flux2]')
                        do l_idx = 0, p; do k_idx = 0, n; do q_idx = 0, m
                            local_inv_ds = 1._wp/dy(k_idx)
                            local_k_term_val = Kterm_arg(q_idx, k_idx, l_idx)
                            local_flux1 = nc_iface_vel_n(2)%vf(2)%sf(q_idx, k_idx, l_idx)
                            local_flux2 = nc_iface_vel_n(2)%vf(2)%sf(q_idx, k_idx - 1, l_idx)
                            rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, &
                                       & l_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, &
                                       & l_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                        end do; end do; end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (cyl_coord) then
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_k_term_val]')
                            do l_idx = 0, p; do k_idx = 0, n; do q_idx = 0, m
                                local_k_term_val = Kterm_arg(q_idx, k_idx, l_idx)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, &
                                           & k_idx, l_idx) + local_k_term_val*q_prim_vf_arg%vf(eqn_idx%mom%beg + 1)%sf(q_idx, &
                                           & k_idx, l_idx)/y_cc(k_idx)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, &
                                           & k_idx, l_idx) - local_k_term_val*q_prim_vf_arg%vf(eqn_idx%mom%beg + 1)%sf(q_idx, &
                                           & k_idx, l_idx)/y_cc(k_idx)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                else if (adv_src_mode == adv_src_mode_vel_iface) then
                    ! u-interface: flux_src(adv%beg) supplies one shared face-normal velocity. RHS applies alpha_k * du/dy.
                    $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                        & local_flux1, local_flux2]')
                    do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                        do l_idx = 0, p; do k_idx = 0, n; do q_idx = 0, m
                            local_inv_ds = 1._wp/dy(k_idx)
                            local_term_coeff = q_cons_vf_arg%vf(j_adv)%sf(q_idx, k_idx, l_idx)
                            local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx)
                            local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx - 1, l_idx)
                            rhs_vf_arg(j_adv)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(j_adv)%sf(q_idx, k_idx, &
                                       & l_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                        end do; end do; end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    if (alt_soundspeed) then  ! K*div(u) correction
                        $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                            & local_flux1, local_flux2]')
                        do l_idx = 0, p; do k_idx = 0, n; do q_idx = 0, m
                            local_inv_ds = 1._wp/dy(k_idx)
                            local_k_term_val = Kterm_arg(q_idx, k_idx, l_idx)
                            local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx)
                            local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx - 1, l_idx)
                            rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, &
                                       & l_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, &
                                       & l_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                        end do; end do; end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (cyl_coord) then
                            $:GPU_PARALLEL_LOOP(collapse=3, &
                                                & private='[k_idx, l_idx, q_idx, local_k_term_val, local_flux1, local_flux2]')
                            do l_idx = 0, p; do k_idx = 0, n; do q_idx = 0, m
                                local_k_term_val = Kterm_arg(q_idx, k_idx, l_idx)
                                local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx)
                                local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(q_idx, k_idx - 1, l_idx)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(q_idx, &
                                           & k_idx, l_idx) + (local_k_term_val/(2._wp*y_cc(k_idx)))*(local_flux1 + local_flux2)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, k_idx, l_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(q_idx, &
                                           & k_idx, l_idx) - (local_k_term_val/(2._wp*y_cc(k_idx)))*(local_flux1 + local_flux2)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            case (3)  ! z-direction
                if (grid_geometry == 3) then
                    if (adv_src_mode == adv_src_mode_alpha_iface) then
                        $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                            & local_flux1, local_flux2]')
                        do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                            do k_idx = 0, p  ! z_extent
                                do q_idx = 0, n  ! y_extent
                                    do l_idx = 0, m  ! x_extent
                                        local_inv_ds = 1._wp/(dz(k_idx)*y_cc(q_idx))
                                        local_term_coeff = q_prim_vf_arg%vf(eqn_idx%cont%end + current_idir)%sf(l_idx, q_idx, k_idx)
                                        local_flux1 = flux_src_n_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx - 1)
                                        local_flux2 = flux_src_n_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx)
                                        rhs_vf_arg(j_adv)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(j_adv)%sf(l_idx, q_idx, &
                                                   & k_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (alt_soundspeed) then  ! K*div(u) correction
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                                & local_flux1, local_flux2]')
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/(dz(k_idx)*y_cc(q_idx))
                                local_k_term_val = Kterm_arg(l_idx, q_idx, k_idx)
                                local_flux1 = nc_iface_vel_n(3)%vf(3)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = nc_iface_vel_n(3)%vf(3)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, &
                                           & q_idx, k_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, &
                                           & q_idx, k_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    else if (adv_src_mode == adv_src_mode_vel_iface) then
                        ! u-interface: flux_src(adv%beg) supplies one shared face-normal velocity. RHS applies alpha_k * du/dz.
                        $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                            & local_flux1, local_flux2]')
                        do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/(dz(k_idx)*y_cc(q_idx))
                                local_term_coeff = q_cons_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx)
                                local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(j_adv)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(j_adv)%sf(l_idx, q_idx, &
                                           & k_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                            end do; end do; end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (alt_soundspeed) then  ! K*div(u) correction
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                                & local_flux1, local_flux2]')
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/(dz(k_idx)*y_cc(q_idx))
                                local_k_term_val = Kterm_arg(l_idx, q_idx, k_idx)
                                local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, &
                                           & q_idx, k_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, &
                                           & q_idx, k_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                else  ! non-cylindrical z-direction
                    if (adv_src_mode == adv_src_mode_alpha_iface) then
                        ! Alpha-interface: flux_src(j_adv) supplies interface alpha_k. RHS applies velocity * d(alpha_k)/dz.
                        $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                            & local_flux1, local_flux2]')
                        do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                            do k_idx = 0, p  ! z_extent
                                do q_idx = 0, n  ! y_extent
                                    do l_idx = 0, m  ! x_extent
                                        local_inv_ds = 1._wp/dz(k_idx)
                                        local_term_coeff = q_prim_vf_arg%vf(eqn_idx%cont%end + current_idir)%sf(l_idx, q_idx, k_idx)
                                        local_flux1 = flux_src_n_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx - 1)
                                        local_flux2 = flux_src_n_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx)
                                        rhs_vf_arg(j_adv)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(j_adv)%sf(l_idx, q_idx, &
                                                   & k_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (alt_soundspeed) then  ! K*div(u) correction
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                                & local_flux1, local_flux2]')
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/dz(k_idx)
                                local_k_term_val = Kterm_arg(l_idx, q_idx, k_idx)
                                local_flux1 = nc_iface_vel_n(3)%vf(3)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = nc_iface_vel_n(3)%vf(3)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, &
                                           & q_idx, k_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, &
                                           & q_idx, k_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    else if (adv_src_mode == adv_src_mode_vel_iface) then
                        ! u-interface: flux_src(adv%beg) supplies one shared face-normal velocity. RHS applies alpha_k * du/dz.
                        $:GPU_PARALLEL_LOOP(collapse=4, private='[j_adv, k_idx, l_idx, q_idx, local_inv_ds, local_term_coeff, &
                                            & local_flux1, local_flux2]')
                        do j_adv = eqn_idx%adv%beg, eqn_idx%adv%end
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/dz(k_idx)
                                local_term_coeff = q_cons_vf_arg%vf(j_adv)%sf(l_idx, q_idx, k_idx)
                                local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(j_adv)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(j_adv)%sf(l_idx, q_idx, &
                                           & k_idx) + local_inv_ds*local_term_coeff*(local_flux1 - local_flux2)
                            end do; end do; end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                        if (alt_soundspeed) then  ! K*div(u) correction
                            $:GPU_PARALLEL_LOOP(collapse=3, private='[k_idx, l_idx, q_idx, local_inv_ds, local_k_term_val, &
                                                & local_flux1, local_flux2]')
                            do k_idx = 0, p; do q_idx = 0, n; do l_idx = 0, m
                                local_inv_ds = 1._wp/dz(k_idx)
                                local_k_term_val = Kterm_arg(l_idx, q_idx, k_idx)
                                local_flux1 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx)
                                local_flux2 = flux_src_n_vf_arg%vf(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx - 1)
                                rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%beg)%sf(l_idx, &
                                           & q_idx, k_idx) + local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                                rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, q_idx, k_idx) = rhs_vf_arg(eqn_idx%adv%end)%sf(l_idx, &
                                           & q_idx, k_idx) - local_k_term_val*local_inv_ds*(local_flux1 - local_flux2)
                            end do; end do; end do
                            $:END_GPU_PARALLEL_LOOP()
                        end if
                    end if
                end if
            end select

        end subroutine s_add_directional_advection_source_terms

    end subroutine s_compute_advection_source_term

    !> Add viscous, surface-tension, species-diffusion, and heat-conduction source flux contributions to the RHS for a direction
    subroutine s_compute_additional_physics_rhs(idir, q_prim_vf, q_T_sf, rhs_vf, flux_src_n_in, dq_prim_dx_vf, dq_prim_dy_vf, &
        & dq_prim_dz_vf, stress_in_kernel)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), intent(in)                         :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: flux_src_n_in
        type(scalar_field), dimension(sys_size), intent(in)    :: dq_prim_dx_vf, dq_prim_dy_vf, dq_prim_dz_vf
        !> True when the fused kernel assembled and differenced the viscous stress, leaving flux_src_n_in(mom..E) unset
        logical, intent(in) :: stress_in_kernel
        integer             :: i, j, k, l

        if (idir == 1) then  ! x-direction
            if (surface_tension) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(eqn_idx%c)%sf(j, k, l) = rhs_vf(eqn_idx%c)%sf(j, k, &
                                   & l) + 1._wp/dx(j)*q_prim_vf(eqn_idx%c)%sf(j, k, l)*(flux_src_n_in(eqn_idx%adv%beg)%sf(j, k, &
                                   & l) - flux_src_n_in(eqn_idx%adv%beg)%sf(j - 1, k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if ((.not. stress_in_kernel .and. (surface_tension .or. viscous)) .or. chem_params%diffusion .or. heat_conduction) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (surface_tension .or. (viscous .and. .not. stress_in_kernel) .or. heat_conduction) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dx(j)*(flux_src_n_in(i)%sf(j - 1, k, &
                                           & l) - flux_src_n_in(i)%sf(j, k, l))
                                end do
                            end if

                            if (chem_params%diffusion) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dx(j)*(flux_src_n_in(i)%sf(j - 1, k, &
                                           & l) - flux_src_n_in(i)%sf(j, k, l))
                                end do

                                if (.not. viscous .or. stress_in_kernel) then
                                    rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, &
                                           & l) + 1._wp/dx(j)*(flux_src_n_in(eqn_idx%E)%sf(j - 1, k, &
                                           & l) - flux_src_n_in(eqn_idx%E)%sf(j, k, l))
                                end if
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (idir == 2) then  ! y-direction
            if (surface_tension) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(eqn_idx%c)%sf(j, k, l) = rhs_vf(eqn_idx%c)%sf(j, k, &
                                   & l) + 1._wp/dy(k)*q_prim_vf(eqn_idx%c)%sf(j, k, l)*(flux_src_n_in(eqn_idx%adv%beg)%sf(j, k, &
                                   & l) - flux_src_n_in(eqn_idx%adv%beg)%sf(j, k - 1, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (cyl_coord .and. ((bc_y%beg == -2) .or. (bc_y%beg == -14))) then
                if (viscous) then
                    call s_compute_viscous_stress_cylindrical_boundary(q_prim_vf, dq_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                        & dq_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), dq_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                        & tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                end if

                if (heat_conduction) then
                    call s_compute_conduction_axis_source(q_prim_vf, q_T_sf, tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                end if

                if (viscous .or. heat_conduction) then
                    $:GPU_PARALLEL_LOOP(private='[i, j, l]', collapse=2)
                    do l = 0, p
                        do j = 0, m
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = eqn_idx%mom%beg, eqn_idx%E
                                rhs_vf(i)%sf(j, 0, l) = rhs_vf(i)%sf(j, 0, l) + 1._wp/(y_cc(1) - y_cc(-1))*(tau_Re_vf(i)%sf(j, &
                                       & -1, l) - tau_Re_vf(i)%sf(j, 1, l))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if

                if (.not. stress_in_kernel) then
                    $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                    do l = 0, p
                        do k = 1, n
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dy(k)*(flux_src_n_in(i)%sf(j, k - 1, &
                                           & l) - flux_src_n_in(i)%sf(j, k, l))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            else
                if ((.not. stress_in_kernel .and. (surface_tension .or. viscous)) .or. chem_params%diffusion .or. heat_conduction) &
                    & then
                    $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                if (surface_tension .or. (viscous .and. .not. stress_in_kernel) .or. heat_conduction) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%mom%beg, eqn_idx%E
                                        rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dy(k)*(flux_src_n_in(i)%sf(j, &
                                               & k - 1, l) - flux_src_n_in(i)%sf(j, k, l))
                                    end do
                                end if

                                if (chem_params%diffusion) then
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%species%beg, eqn_idx%species%end
                                        rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dy(k)*(flux_src_n_in(i)%sf(j, &
                                               & k - 1, l) - flux_src_n_in(i)%sf(j, k, l))
                                    end do
                                    if (.not. viscous .or. stress_in_kernel) then
                                        rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, &
                                               & l) + 1._wp/dy(k)*(flux_src_n_in(eqn_idx%E)%sf(j, k - 1, &
                                               & l) - flux_src_n_in(eqn_idx%E)%sf(j, k, l))
                                    end if
                                end if
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if

            ! Applying the geometrical viscous Riemann source fluxes calculated as average of values at cell boundaries
            if (cyl_coord) then
                if ((bc_y%beg == -2) .or. (bc_y%beg == -14)) then
                    if (.not. stress_in_kernel) then
                        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                        do l = 0, p
                            do k = 1, n
                                do j = 0, m
                                    $:GPU_LOOP(parallelism='[seq]')
                                    do i = eqn_idx%mom%beg, eqn_idx%E
                                        rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - 5.e-1_wp/y_cc(k)*(flux_src_n_in(i)%sf(j, &
                                               & k - 1, l) + flux_src_n_in(i)%sf(j, k, l))
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if

                    if (viscous .or. heat_conduction) then
                        $:GPU_PARALLEL_LOOP(private='[i, j, l]', collapse=2)
                        do l = 0, p
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, 0, l) = rhs_vf(i)%sf(j, 0, l) - 1._wp/y_cc(0)*tau_Re_vf(i)%sf(j, 0, l)
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
                    end if
                else if (.not. stress_in_kernel) then
                    $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) - 5.e-1_wp/y_cc(k)*(flux_src_n_in(i)%sf(j, &
                                           & k - 1, l) + flux_src_n_in(i)%sf(j, k, l))
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end if
            end if
        else if (idir == 3) then  ! z-direction
            if (surface_tension) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(eqn_idx%c)%sf(j, k, l) = rhs_vf(eqn_idx%c)%sf(j, k, &
                                   & l) + 1._wp/dz(l)*q_prim_vf(eqn_idx%c)%sf(j, k, l)*(flux_src_n_in(eqn_idx%adv%beg)%sf(j, k, &
                                   & l) - flux_src_n_in(eqn_idx%adv%beg)%sf(j, k, l - 1))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if ((.not. stress_in_kernel .and. (surface_tension .or. viscous)) .or. chem_params%diffusion .or. heat_conduction) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (surface_tension .or. (viscous .and. .not. stress_in_kernel) .or. heat_conduction) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%mom%beg, eqn_idx%E
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dz(l)*(flux_src_n_in(i)%sf(j, k, &
                                           & l - 1) - flux_src_n_in(i)%sf(j, k, l))
                                end do
                            end if

                            if (chem_params%diffusion) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = eqn_idx%species%beg, eqn_idx%species%end
                                    rhs_vf(i)%sf(j, k, l) = rhs_vf(i)%sf(j, k, l) + 1._wp/dz(l)*(flux_src_n_in(i)%sf(j, k, &
                                           & l - 1) - flux_src_n_in(i)%sf(j, k, l))
                                end do
                                if (.not. viscous .or. stress_in_kernel) then
                                    rhs_vf(eqn_idx%E)%sf(j, k, l) = rhs_vf(eqn_idx%E)%sf(j, k, &
                                           & l) + 1._wp/dz(l)*(flux_src_n_in(eqn_idx%E)%sf(j, k, &
                                           & l - 1) - flux_src_n_in(eqn_idx%E)%sf(j, k, l))
                                end if
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            ! With stress_in_kernel the fused column carries the cross-momentum geometric terms
            if (grid_geometry == 3 .and. .not. stress_in_kernel) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, l) = rhs_vf(eqn_idx%mom%beg + 1)%sf(j, k, &
                                   & l) + 5.e-1_wp*(flux_src_n_in(eqn_idx%mom%end)%sf(j, k, &
                                   & l - 1) + flux_src_n_in(eqn_idx%mom%end)%sf(j, k, l))

                            rhs_vf(eqn_idx%mom%end)%sf(j, k, l) = rhs_vf(eqn_idx%mom%end)%sf(j, k, &
                                   & l) - 5.e-1_wp*(flux_src_n_in(eqn_idx%mom%beg + 1)%sf(j, k, &
                                   & l - 1) + flux_src_n_in(eqn_idx%mom%beg + 1)%sf(j, k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_compute_additional_physics_rhs

#ifdef MFC_FUSED_RHS
    !> WENO-reconstruct the cell-centred velocity gradients onto the cell boundaries, giving the Riemann solver one-sided face
    !! gradients. It runs for weno_Re_flux; the default viscous path builds the face gradients directly.
    subroutine s_reconstruct_viscous_derivatives(norm_dir)

        integer, intent(in) :: norm_dir

        iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end

        call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dx_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, dqR_rsx_vf, norm_dir, &
            & dqL_prim_dx_n(norm_dir)%vf(iv%beg:iv%end), dqR_prim_dx_n(norm_dir)%vf(iv%beg:iv%end), idwbuff(1), idwbuff(2), &
            & idwbuff(3))
        if (n > 0) then
            call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dy_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, dqR_rsx_vf, &
                & norm_dir, dqL_prim_dy_n(norm_dir)%vf(iv%beg:iv%end), dqR_prim_dy_n(norm_dir)%vf(iv%beg:iv%end), idwbuff(1), &
                & idwbuff(2), idwbuff(3))
            if (p > 0) then
                call s_reconstruct_cell_boundary_values_visc_deriv(dq_prim_dz_qp(1)%vf(iv%beg:iv%end), dqL_rsx_vf, dqR_rsx_vf, &
                    & norm_dir, dqL_prim_dz_n(norm_dir)%vf(iv%beg:iv%end), dqR_prim_dz_n(norm_dir)%vf(iv%beg:iv%end), idwbuff(1), &
                    & idwbuff(2), idwbuff(3))
            end if
        end if

    end subroutine s_reconstruct_viscous_derivatives

    !> Average the two face-gradient families into dqR_prim_d?_n, the family the in-kernel stress reads
    subroutine s_average_fused_face_gradients(norm_dir)

        integer, intent(in) :: norm_dir
        integer             :: i, j, k, l, xo, yo, zo, xb, yb, zb
        ! Local bounds, independent of irx/iry/irz

        xo = 0; yo = 0; zo = 0
        xb = 0; yb = 0; zb = 0
        if (norm_dir == 1) then
            xo = 1; xb = -1
        else if (norm_dir == 2) then
            yo = 1; yb = -1
        else
            zo = 1; zb = -1
        end if

        #:for DIR, ARRL, ARRR in [(1, 'dqL_prim_dx_n', 'dqR_prim_dx_n'), (2, 'dqL_prim_dy_n', 'dqR_prim_dy_n'), &
    & (3, 'dqL_prim_dz_n', 'dqR_prim_dz_n')]
            if (num_dims >= ${DIR}$) then
                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = eqn_idx%mom%beg, eqn_idx%mom%end
                    do l = zb, p
                        do k = yb, n
                            do j = xb, m
                                ${ARRR}$(norm_dir)%vf(i)%sf(j, k, l) = 5.e-1_wp*(${ARRR}$(norm_dir)%vf(i)%sf(j, k, &
                                         & l) + ${ARRL}$(norm_dir)%vf(i)%sf(j + xo, k + yo, l + zo))
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        #:endfor

    end subroutine s_average_fused_face_gradients
#endif

    !> Reconstruct left and right cell-boundary values from cell-averaged variables
    subroutine s_reconstruct_cell_boundary_values(v_vf, vL_x, vR_x, norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_x
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vR_x
        integer, intent(in) :: norm_dir
        integer :: recon_dir  !< Coordinate direction of the reconstruction
        integer :: i, j, k, l

        #:for SCHEME, TYPE in [('weno','recon_type_weno'), ('muscl','recon_type_muscl')]
            if (recon_type == ${TYPE}$) then
                ! Reconstruction in s1-direction
                if (norm_dir == 1) then
                    is1 = idwbuff(1); is2 = idwbuff(2); is3 = idwbuff(3)
                    recon_dir = 1; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                else if (norm_dir == 2) then
                    is1 = idwbuff(2); is2 = idwbuff(1); is3 = idwbuff(3)
                    recon_dir = 2; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                else
                    is1 = idwbuff(3); is2 = idwbuff(2); is3 = idwbuff(1)
                    recon_dir = 3; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                end if

                call s_${SCHEME}$ (v_vf(iv%beg:iv%end), vL_x(:,:,:,iv%beg:iv%end), vR_x(:,:,:,iv%beg:iv%end), recon_dir, is1, &
                                   & is2, is3)
            end if
        #:endfor

    end subroutine s_reconstruct_cell_boundary_values

    !> Perform first-order (piecewise constant) reconstruction of left and right cell-boundary values
    subroutine s_reconstruct_cell_boundary_values_first_order(v_vf, vL_x, vR_x, norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in) :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vL_x
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: vR_x
        integer, intent(in) :: norm_dir
        integer :: recon_dir  !< Coordinate direction of the reconstruction
        integer :: i, j, k, l
        integer :: is1_b, is1_e, is2_b, is2_e, is3_b, is3_e, iv_b, iv_e
        ! Reconstruction in s1-direction

        #:for SCHEME, TYPE in [('weno','recon_type_weno'), ('muscl', 'recon_type_muscl')]
            if (recon_type == ${TYPE}$) then
                if (norm_dir == 1) then
                    is1 = idwbuff(1); is2 = idwbuff(2); is3 = idwbuff(3)
                    recon_dir = 1; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                else if (norm_dir == 2) then
                    is1 = idwbuff(2); is2 = idwbuff(1); is3 = idwbuff(3)
                    recon_dir = 2; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                else
                    is1 = idwbuff(3); is2 = idwbuff(2); is3 = idwbuff(1)
                    recon_dir = 3; is1%beg = is1%beg + ${SCHEME}$_polyn
                    is1%end = is1%end - ${SCHEME}$_polyn
                end if
            end if
        #:endfor
        ! Local copies so the kernels take their bounds as arguments
        is1_b = is1%beg; is1_e = is1%end
        is2_b = is2%beg; is2_e = is2%end
        is3_b = is3%beg; is3_e = is3%end
        iv_b = iv%beg; iv_e = iv%end

        if (recon_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv_b, iv_e
                do l = is3_b, is3_e
                    do k = is2_b, is2_e
                        do j = is1_b, is1_e
                            vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (recon_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv_b, iv_e
                do l = is3_b, is3_e
                    do j = is1_b, is1_e
                        do k = is2_b, is2_e
                            vL_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                            vR_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (recon_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv_b, iv_e
                do j = is1_b, is1_e
                    do k = is2_b, is2_e
                        do l = is3_b, is3_e
                            vL_x(l, k, j, i) = v_vf(i)%sf(l, k, j)
                            vR_x(l, k, j, i) = v_vf(i)%sf(l, k, j)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_reconstruct_cell_boundary_values_first_order

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_rhs_module

        integer :: i, j, l

#ifdef MFC_CUDA_GRAPHS
        call s_finalize_cuda_graph_module
#endif

        if (.not. igr) then
            do j = eqn_idx%cont%beg, eqn_idx%cont%end
                if (relativity) then
                    ! Cons and Prim densities are different for relativity
                    @:DEALLOCATE(q_cons_qp%vf(j)%sf)
                    @:DEALLOCATE(q_prim_qp%vf(j)%sf)
                else
                    nullify (q_prim_qp%vf(j)%sf)
                end if
            end do

            do j = eqn_idx%adv%beg, eqn_idx%adv%end
                nullify (q_prim_qp%vf(j)%sf)
            end do

            do j = eqn_idx%mom%beg, eqn_idx%E
                @:DEALLOCATE(q_cons_qp%vf(j)%sf)
                @:DEALLOCATE(q_prim_qp%vf(j)%sf)
            end do
        end if

#if defined(MFC_FUSED_RHS)
        if (.not. igr) then
            do j = 1, sys_size
                @:DEALLOCATE(q_prim_next%vf(j)%sf)
            end do
        end if
        @:DEALLOCATE(q_prim_next%vf)
#endif
        @:DEALLOCATE(q_cons_qp%vf, q_prim_qp%vf)

        if (.not. igr) then
            @:DEALLOCATE(qL_rsx_vf, qR_rsx_vf)

            if (viscous) then
                do l = eqn_idx%mom%beg, eqn_idx%mom%end
                    @:DEALLOCATE(dq_prim_dx_qp(1)%vf(l)%sf)
                end do

                if (n > 0) then
                    do l = eqn_idx%mom%beg, eqn_idx%mom%end
                        @:DEALLOCATE(dq_prim_dy_qp(1)%vf(l)%sf)
                    end do

                    if (p > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:DEALLOCATE(dq_prim_dz_qp(1)%vf(l)%sf)
                        end do
                    end if
                end if

                @:DEALLOCATE(dq_prim_dx_qp(1)%vf)
                @:DEALLOCATE(dq_prim_dy_qp(1)%vf)
                @:DEALLOCATE(dq_prim_dz_qp(1)%vf)

                do i = num_dims, 1, -1
                    do l = eqn_idx%mom%beg, eqn_idx%mom%end
                        @:DEALLOCATE(dqL_prim_dx_n(i)%vf(l)%sf)
                        @:DEALLOCATE(dqR_prim_dx_n(i)%vf(l)%sf)
                    end do

                    if (n > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:DEALLOCATE(dqL_prim_dy_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dy_n(i)%vf(l)%sf)
                        end do
                    end if

                    if (p > 0) then
                        do l = eqn_idx%mom%beg, eqn_idx%mom%end
                            @:DEALLOCATE(dqL_prim_dz_n(i)%vf(l)%sf)
                            @:DEALLOCATE(dqR_prim_dz_n(i)%vf(l)%sf)
                        end do
                    end if

                    @:DEALLOCATE(dqL_prim_dx_n(i)%vf)
                    @:DEALLOCATE(dqL_prim_dy_n(i)%vf)
                    @:DEALLOCATE(dqL_prim_dz_n(i)%vf)
                    @:DEALLOCATE(dqR_prim_dx_n(i)%vf)
                    @:DEALLOCATE(dqR_prim_dy_n(i)%vf)
                    @:DEALLOCATE(dqR_prim_dz_n(i)%vf)
                end do

                if (weno_Re_flux) then
                    @:DEALLOCATE(dqL_rsx_vf, dqR_rsx_vf)
                end if
            end if

            if (viscous .or. heat_conduction) then
                do i = 1, num_dims
                    @:DEALLOCATE(tau_Re_vf(eqn_idx%cont%end + i)%sf)
                end do
                @:DEALLOCATE(tau_Re_vf(eqn_idx%E)%sf)
                @:DEALLOCATE(tau_Re_vf)
            end if

            @:DEALLOCATE(dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n)
            @:DEALLOCATE(dqR_prim_dx_n, dqR_prim_dy_n, dqR_prim_dz_n)
        end if

        if (mpp_lim .and. bubbles_euler) then
            @:DEALLOCATE(alf_sum%sf)
        end if

        if (.not. igr) then
            do i = num_dims, 1, -1
                if (i /= 1) then
                    do l = 1, sys_size
                        nullify (flux_n(i)%vf(l)%sf)
                        nullify (flux_src_n(i)%vf(l)%sf)
                        @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                    end do
                else
                    do l = 1, sys_size
#ifdef MFC_DIRECT_RIEMANN_FLUX
                        nullify (flux_n(i)%vf(l)%sf)
#else
                        @:DEALLOCATE(flux_n(i)%vf(l)%sf)
#endif
                        @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                    end do

                    if (viscous .or. surface_tension .or. heat_conduction) then
                        do l = eqn_idx%mom%beg, eqn_idx%E
                            @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        end do
                    end if

                    if (chem_params%diffusion .and. .not. (viscous .or. surface_tension .or. heat_conduction)) then
                        @:DEALLOCATE(flux_src_n(i)%vf(eqn_idx%E)%sf)
                    end if

#ifdef MFC_DIRECT_RIEMANN_FLUX
                    do l = eqn_idx%adv%beg, eqn_idx%adv%end
                        nullify (flux_src_n(i)%vf(l)%sf)
                    end do
#else
                    if (adv_src_mode == adv_src_mode_alpha_iface .or. adv_src_mode == adv_src_mode_none) then
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        end do
                    else
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            nullify (flux_src_n(i)%vf(l)%sf)
                        end do
                    end if

                    @:DEALLOCATE(flux_src_n(i)%vf(eqn_idx%adv%beg)%sf)
#endif
                end if

                @:DEALLOCATE(flux_n(i)%vf, flux_src_n(i)%vf, flux_gsrc_n(i)%vf)
            end do

            @:DEALLOCATE(flux_n, flux_src_n, flux_gsrc_n)
            do i = 1, num_dims
                do l = eqn_idx%mom%beg, eqn_idx%mom%end
                    @:DEALLOCATE(qL_prim(i)%vf(l)%sf)
                    @:DEALLOCATE(qR_prim(i)%vf(l)%sf)
                end do
                @:DEALLOCATE(qL_prim(i)%vf, qR_prim(i)%vf)
            end do
            @:DEALLOCATE(qL_prim, qR_prim)
        end if

        if (alt_soundspeed) then
            @:DEALLOCATE(blkmod1, blkmod2, alpha1, alpha2, Kterm)
        end if

        if (qbmm) then
            do i = 0, 2
                do j = 0, 2
                    do l = 1, nb
                        @:DEALLOCATE(mom_3d(i, j, l)%sf)
                    end do
                end do
            end do
            do i = 1, nmomsp
                @:DEALLOCATE(mom_sp(i)%sf)
            end do
            @:DEALLOCATE(mom_sp, mom_3d)
        end if

        if (use_nc_iface_vel) then
            do i = 1, num_dims
                do l = 1, num_dims
                    @:DEALLOCATE(nc_iface_vel_n(i)%vf(l)%sf)
                end do
                @:DEALLOCATE(nc_iface_vel_n(i)%vf)
            end do
            @:DEALLOCATE(nc_iface_vel_n)
            if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
                do i = 1, num_dims
                    do l = 1, num_dims
                        @:DEALLOCATE(nc_iface_vel_hatR_n(i)%vf(l)%sf)
                    end do
                    @:DEALLOCATE(nc_iface_vel_hatR_n(i)%vf)
                end do
                @:DEALLOCATE(nc_iface_vel_hatR_n)
            end if
        end if

        if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
            do i = 1, sys_size
                @:DEALLOCATE(rhs_hatL_vf(i)%sf)
                @:DEALLOCATE(rhs_hatR_vf(i)%sf)
            end do
            @:DEALLOCATE(rhs_hatL_vf)
            @:DEALLOCATE(rhs_hatR_vf)
        end if

        if (bf_spatial_support) then
            @:DEALLOCATE(spbf_source_x, spbf_source_y)
        end if

    end subroutine s_finalize_rhs_module

end module m_rhs
