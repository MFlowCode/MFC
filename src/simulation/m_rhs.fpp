!>
!! @file
!! @brief Contains module m_rhs

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Assembles the right-hand side of the governing equations using finite-volume flux differencing, Riemann solvers, and
!! physical source terms
module m_rhs

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_weno
    use m_constants, only: riemann_solver_hll, riemann_solver_hlld, model_eqns_6eq, int_comp_mthinc, recon_type_weno, &
        & recon_type_muscl
    use m_muscl
    use m_riemann_solvers
    use m_cbc
    use m_bubbles_EE
    use m_bubbles_EL
    use m_qbmm
    use m_hypoelastic
    use m_hyperelastic
    use m_acoustic_src
    use m_viscous
    use m_ibm
    use m_nvtx
    use m_boundary_common
    use m_helper
    use m_surface_tension
    use m_body_forces
    use m_chemistry
    use m_igr
    use m_thinc
    use m_pressure_relaxation

    implicit none

    private; public :: s_initialize_rhs_module, s_compute_rhs, s_finalize_rhs_module

    type(vector_field) :: q_cons_qp  !< WENO-reconstructed cell-average conservative variables at quadrature points
    $:GPU_DECLARE(create='[q_cons_qp]')

    type(vector_field) :: q_prim_qp  !< Primitive variables at cell-interior quadrature points
    $:GPU_DECLARE(create='[q_prim_qp]')

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
    $:GPU_DECLARE(create='[nc_iface_vel_n]')

    !> hat_R-pass interface velocities for fused dual-pass HLLD: the axisymmetric geometric source consumes both passes' face
    !! velocities after the direction loop, so the hat_R values need their own field set (hat_L uses nc_iface_vel_n).
    type(vector_field), allocatable, dimension(:) :: nc_iface_vel_hatR_n
    $:GPU_DECLARE(create='[nc_iface_vel_hatR_n]')

    type(vector_field), allocatable, dimension(:) :: qL_prim, qR_prim
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[qL_prim, qR_prim]')
#endif

    type(int_bounds_info) :: iv  !< Vector field indical bounds
    $:GPU_DECLARE(create='[iv]')

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(int_bounds_info) :: irx, iry, irz
    $:GPU_DECLARE(create='[irx, iry, irz]')

    type(int_bounds_info) :: is1, is2, is3
    !> @}
    $:GPU_DECLARE(create='[is1, is2, is3]')

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

    integer :: iglob
    $:GPU_DECLARE(create='[iglob]')

    !> @name Partial right-hand sides of the two anchored passes of dual-pass HLLD
    !> @{
    type(scalar_field), allocatable, dimension(:) :: rhs_hatL_vf, rhs_hatR_vf
    !> @}
    $:GPU_DECLARE(create='[rhs_hatL_vf, rhs_hatR_vf]')

contains

    !> Initialize the RHS module
    impure subroutine s_initialize_rhs_module

        integer :: i, j, k, l, id  !< Generic loop iterators

        $:GPU_ENTER_DATA(copyin='[idwbuff]')
        $:GPU_UPDATE(device='[idwbuff]')

        @:ALLOCATE(q_cons_qp%vf(1:sys_size))
        @:ALLOCATE(q_prim_qp%vf(1:sys_size))

        if (.not. igr) then
            do l = 1, sys_size
                @:ALLOCATE(q_cons_qp%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            do l = eqn_idx%mom%beg, eqn_idx%E
                @:ALLOCATE(q_prim_qp%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        end if

        if (surface_tension) then
            do l = eqn_idx%adv%end + 1, eqn_idx%c - 1
                @:ALLOCATE(q_prim_qp%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        else
            do l = eqn_idx%adv%end + 1, sys_size
                @:ALLOCATE(q_prim_qp%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
        end if

        if (.not. igr) then
            @:ACC_SETUP_VFs(q_cons_qp, q_prim_qp)

            do l = 1, eqn_idx%cont%end
                if (relativity) then
                    ! Cons and Prim densities are different for relativity
                    @:ALLOCATE(q_prim_qp%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
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
                    do l = 1, sys_size
                        @:ALLOCATE(flux_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                        @:ALLOCATE(flux_gsrc_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                    end do

                    if (viscous .or. surface_tension) then
                        do l = eqn_idx%mom%beg, eqn_idx%E
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if

                    ! flux_src at adv%beg:adv%end is overloaded by design:
                    !
                    ! adv_src_alpha_iface: flux_src(adv%beg:adv%end) = per-fluid interface alpha values. Each fluid gets a separate
                    ! 3D array.
                    !
                    ! adv_src_vel_iface: flux_src(adv%beg) = one shared face-normal interface velocity. adv%beg+1:adv%end are
                    ! pointer-aliased to adv%beg so loops over adv%beg:adv%end can keep fluid indexing while still reading one
                    ! value. Saves (num_fluids - 1) fields.
                    !
                    ! adv_src_none: flux_src(adv%beg:adv%end) = zeros. No separate NC advection source term. Allocated for
                    ! structural consistency with s_finalize_riemann_solver.
                    @:ALLOCATE(flux_src_n(i)%vf(eqn_idx%adv%beg)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))

                    if (adv_src_alpha_iface .or. adv_src_none) then
                        ! Alpha-interface needs separate per-fluid arrays. HLLD (adv_src_none) allocates for structural consistency
                        ! with s_finalize_riemann_solver.
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                    end if

                    if (chemistry) then
                        do l = eqn_idx%species%beg, eqn_idx%species%end
                            @:ALLOCATE(flux_src_n(i)%vf(l)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                       & idwbuff(3)%beg:idwbuff(3)%end))
                        end do
                        if (chem_params%diffusion .and. .not. viscous) then
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
                    if (adv_src_vel_iface) then
                        ! u-interface: flux_src(adv%beg) holds one shared face-normal velocity. Pointer-alias adv%beg+1:adv%end to
                        ! the same memory so loops over adv%beg:adv%end can keep fluid indexing while still reading one value. This
                        ! saves (num_fluids - 1) 3D field allocations.
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            flux_src_n(i)%vf(l)%sf => flux_src_n(i)%vf(eqn_idx%adv%beg)%sf
                            $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
                        end do
                    end if
                else
                    do l = 1, sys_size
                        flux_n(i)%vf(l)%sf => flux_n(1)%vf(l)%sf
                        $:GPU_ENTER_DATA(attach='[flux_n(i)%vf(l)%sf]')
                        flux_src_n(i)%vf(l)%sf => flux_src_n(1)%vf(l)%sf
                        $:GPU_ENTER_DATA(attach='[flux_src_n(i)%vf(l)%sf]')
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

            if (viscous) then
                @:ALLOCATE(tau_Re_vf(1:sys_size))
                do i = 1, num_dims
                    @:ALLOCATE(tau_Re_vf(eqn_idx%cont%end + i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                               & idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(tau_Re_vf(eqn_idx%cont%end + i))
                end do
                @:ALLOCATE(tau_Re_vf(eqn_idx%E)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(tau_Re_vf(eqn_idx%E))

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
                        @:ALLOCATE(mom_3d(i, j, k)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                                   & idwbuff(3)%beg:idwbuff(3)%end))
                        @:ACC_SETUP_SFs(mom_3d(i, j, k))
                    end do
                end do
            end do

            do i = 1, nmomsp
                @:ALLOCATE(mom_sp(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(mom_sp(i))
            end do
        end if

        if (mpp_lim .and. bubbles_euler) then
            @:ALLOCATE(alf_sum%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
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
            if (hypo_nc_dual_pass) then
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

        if (hypo_nc_dual_pass) then
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

        call s_initialize_pressure_relaxation_module

    end subroutine s_initialize_rhs_module

    !> Compute the right-hand side of the semi-discrete governing equations for a single time stage
    impure subroutine s_compute_rhs(q_cons_vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &

        & time_avg, stage)

        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_cons_vf
        type(scalar_field), intent(inout)                                                          :: q_T_sf
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                 :: bc_type
        type(scalar_field), dimension(sys_size), intent(inout)                                     :: rhs_vf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in

        ! TODO :: I think these other two variables need to be stp as well, but it doesn't compile like that right now
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: rhs_pb
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: mv_in
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: rhs_mv
        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg
        integer, intent(in) :: stage
        real(wp) :: t_start, t_finish
        integer :: id
        integer(kind=8) :: i, j, k, l, q  !< Generic loop iterators

        ! RHS: halo exchange -> reconstruct -> Riemann solve -> flux difference -> source terms

        call nvtxStartRange("COMPUTE-RHS")

        call cpu_time(t_start)

        if (.not. igr) then
            ! Association/Population of Working Variables
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

        if (igr) then
            call nvtxStartRange("RHS-COMMUNICATION")
            call s_populate_variables_buffers(bc_type, q_cons_vf, pb_in, mv_in, q_T_sf)
            call nvtxEndRange
        end if
        if (.not. igr) then
            call nvtxStartRange("RHS-CONVERT")
            call s_convert_conservative_to_primitive_variables(q_cons_qp%vf, q_T_sf, q_prim_qp%vf, idwint)
            call nvtxEndRange

            call nvtxStartRange("RHS-COMMUNICATION")
            call s_populate_variables_buffers(bc_type, q_prim_qp%vf, pb_in, mv_in, q_T_sf)
            call nvtxEndRange
        end if

        call nvtxStartRange("RHS-ELASTIC")
        if (hyperelasticity) call s_hyperelastic_rmt_stress_update(q_cons_qp%vf, q_prim_qp%vf)
        call nvtxEndRange

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        if (qbmm) call s_mom_inv(q_cons_qp%vf, q_prim_qp%vf, mom_sp, mom_3d, pb_in, rhs_pb, mv_in, rhs_mv, idwbuff(1), &
            & idwbuff(2), idwbuff(3))

        if ((viscous .and. .not. igr)) then
            call nvtxStartRange("RHS-VISCOUS")
            call s_get_viscous(qL_rsx_vf, dqL_prim_dx_n, dqL_prim_dy_n, dqL_prim_dz_n, qL_prim, qR_rsx_vf, dqR_prim_dx_n, &
                               & dqR_prim_dy_n, dqR_prim_dz_n, qR_prim, q_prim_qp, dq_prim_dx_qp, dq_prim_dy_qp, dq_prim_dz_qp, &
                               & idwbuff(1), idwbuff(2), idwbuff(3))
            call nvtxEndRange
        end if

        if (surface_tension) then
            call nvtxStartRange("RHS-SURFACE-TENSION")
            call s_get_capillary(q_prim_qp%vf, bc_type)
            call nvtxEndRange
        end if

        if (int_comp == int_comp_mthinc .and. n > 0) then
            call nvtxStartRange("RHS-COMPRESSION-NORMALS")
            call s_compute_mthinc_normals(q_prim_qp%vf)
            call nvtxEndRange
        end if

        if (.not. hypo_nc_dual_pass) then
            do id = 1, num_dims
                if (igr) then
                    if (id == 1) then
                        $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                        do l = -1, p + 1
                            do k = -1, n + 1
                                do j = -1, m + 1
                                    do i = 1, sys_size
                                        rhs_vf(i)%sf(j, k, l) = 0._stp
                                    end do
                                end do
                            end do
                        end do
                        $:END_GPU_PARALLEL_LOOP()
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
                    call s_reconstruct_riemann_states(id)

                    call s_compute_directional_rhs(id, rhs_vf, .false.)

                    ! RHS additions for hypoelasticity
                    if (hypo_nc_finite_diff) then
                        call nvtxStartRange("RHS-HYPOELASTICITY-FD-PER-SWEEP")
                        call s_compute_hypoelastic_rhs_finite_diff_per_sweep(id, q_prim_qp%vf, rhs_vf)
                        call nvtxEndRange
                    end if

                    ! RHS for diffusion
                    if (chemistry .and. chem_params%diffusion) then
                        call nvtxStartRange("RHS-CHEM-DIFFUSION")
                        call s_compute_chemistry_diffusion_flux(id, q_prim_qp%vf, flux_src_n(id)%vf, irx, iry, irz, q_T_sf)
                        call nvtxEndRange
                    end if

                    ! Viscous stress contribution to RHS
                    if (viscous .or. surface_tension .or. chem_params%diffusion) then
                        call nvtxStartRange("RHS-ADD-PHYSICS")
                        call s_compute_additional_physics_rhs(id, q_prim_qp%vf, rhs_vf, flux_src_n(id)%vf, dq_prim_dx_qp(1)%vf, &
                                                              & dq_prim_dy_qp(1)%vf, dq_prim_dz_qp(1)%vf)
                        call nvtxEndRange
                    end if

                    ! Bubble dynamics source terms
                    if (bubbles_euler) then
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
            ! flux field arrays and assembled. The full RHS is the sum of the two anchored partial RHS's; the axisymmetric geometric
            ! source is applied per pass at half weight so the sum carries the symmetric average of the hat_L/hat_R face velocities.
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
            if (grid_geometry == 2) then
                call nvtxStartRange("RHS-HYPOELASTICITY-AXISYM-HLLD")
                call s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_qp%vf, rhs_hatL_vf, nc_iface_vel_n(1)%vf, &
                    & nc_iface_vel_n(2)%vf, 0.5_wp)
                call s_compute_hypoelastic_rhs_axisym_geom_iface(q_prim_qp%vf, rhs_hatR_vf, nc_iface_vel_hatR_n(1)%vf, &
                    & nc_iface_vel_hatR_n(2)%vf, 0.5_wp)
                call nvtxEndRange
            end if

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
        end if
        ! END: Dimensional Splitting Loop

        ! RHS additions for hypoelasticity (interface-consistent path, after all sweeps)
        if (hypo_nc_interface) then
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

        ! Additional Physics and Source Terms Additions for acoustic_source
        if (acoustic_source) then
            call nvtxStartRange("RHS-ACOUSTIC-SRC")
            call s_acoustic_src_calculations(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf)
            call nvtxEndRange
        end if

        ! Add bubbles source term
        if (bubbles_euler .and. (.not. adap_dt) .and. (.not. qbmm)) then
            call nvtxStartRange("RHS-BUBBLES-SRC")
            call s_compute_bubble_EE_source(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf, divu)
            call nvtxEndRange
        end if

        if (bubbles_lagrange) then
            ! RHS additions for sub-grid bubbles_lagrange
            call nvtxStartRange("RHS-EL-BUBBLES-SRC")
            call s_compute_bubbles_EL_source(q_cons_qp%vf(1:sys_size), q_prim_qp%vf(1:sys_size), rhs_vf)
            call nvtxEndRange
            ! Compute bubble dynamics
            if (.not. adap_dt) then
                call nvtxStartRange("RHS-EL-BUBBLES-DYN")
                call s_compute_bubble_EL_dynamics(q_prim_qp%vf(1:sys_size), stage)
                call nvtxEndRange
            end if
        end if

        if (chemistry .and. chem_params%reactions) then
            call nvtxStartRange("RHS-CHEM-REACTIONS")
            call s_compute_chemistry_reaction_flux(rhs_vf, q_cons_qp%vf, q_T_sf, q_prim_qp%vf, idwint)
            call nvtxEndRange
        end if

        if (cont_damage) call s_compute_damage_state(q_cons_qp%vf, rhs_vf)

        ! END: Additional physics and source terms

        if (run_time_info .or. probe_wrt .or. ib .or. bubbles_lagrange) then
            if (.not. igr) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do i = 1, sys_size
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                q_prim_vf(i)%sf(j, k, l) = q_prim_qp%vf(i)%sf(j, k, l)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        call cpu_time(t_finish)

        if (t_step >= 2) then
            time_avg = (abs(t_finish - t_start) + (t_step - 2)*time_avg)/(t_step - 1)
        else
            time_avg = 0._wp
        end if

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

                iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end; iglob = id
                $:GPU_UPDATE(device='[iv, iglob]')

                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = iv%beg, iv%end
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                qL_rsx_vf(j, k, l, i) = qL_prim(iglob)%vf(i)%sf(j, k, l)
                                qR_rsx_vf(j, k, l, i) = qR_prim(iglob)%vf(i)%sf(j, k, l)
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

                iv%beg = eqn_idx%mom%beg; iv%end = eqn_idx%mom%end; iglob = id
                $:GPU_UPDATE(device='[iv, iglob]')

                $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
                do i = iv%beg, iv%end
                    do l = idwbuff(3)%beg, idwbuff(3)%end
                        do k = idwbuff(2)%beg, idwbuff(2)%end
                            do j = idwbuff(1)%beg, idwbuff(1)%end
                                qL_rsx_vf(j, k, l, i) = qL_prim(iglob)%vf(i)%sf(j, k, l)
                                qR_rsx_vf(j, k, l, i) = qR_prim(iglob)%vf(i)%sf(j, k, l)
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
        real(wp) :: inv_ds, flux_face1, flux_face2
        real(wp) :: advected_qty_val, pressure_val, velocity_val
        real(wp) :: geom_fac
        real(wp) :: G1_eff, G2_eff

        G1_eff = 0._wp
        G2_eff = 0._wp
        if (hypoelasticity) then
            G1_eff = fluid_pp(1)%G
            G2_eff = fluid_pp(2)%G
        end if

        if (alt_soundspeed) then
            $:GPU_PARALLEL_LOOP(private='[k_loop, l_loop, q_loop]', collapse=3)
            do q_loop = 0, p
                do l_loop = 0, n
                    do k_loop = 0, m
                        blkmod1(k_loop, l_loop, q_loop) = ((gammas(1) + 1._wp)*q_prim_vf%vf(eqn_idx%E)%sf(k_loop, l_loop, &
                                & q_loop) + pi_infs(1))/gammas(1) + (4._wp/3._wp)*G1_eff
                        blkmod2(k_loop, l_loop, q_loop) = ((gammas(2) + 1._wp)*q_prim_vf%vf(eqn_idx%E)%sf(k_loop, l_loop, &
                                & q_loop) + pi_infs(2))/gammas(2) + (4._wp/3._wp)*G2_eff
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

            if (.not. hypo_nc_dual_pass) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face1, flux_face2]')
                do j = 1, sys_size
                    do q_loop = 0, p
                        do l_loop = 0, n
                            do k_loop = 0, m
                                inv_ds = 1._wp/dx(k_loop)
                                flux_face1 = flux_n(1)%vf(j)%sf(k_loop - 1, l_loop, q_loop)
                                flux_face2 = flux_n(1)%vf(j)%sf(k_loop, l_loop, q_loop)
                                rhs_vf(j)%sf(k_loop, l_loop, q_loop) = inv_ds*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (is_hat_L) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face2]')
                do j = 1, sys_size
                    do q_loop = 0, p
                        do l_loop = 0, n
                            do k_loop = 0, m
                                inv_ds = 1._wp/dx(k_loop)
                                flux_face2 = flux_n(1)%vf(j)%sf(k_loop, l_loop, q_loop)
                                rhs_vf(j)%sf(k_loop, l_loop, q_loop) = -inv_ds*flux_face2
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! is_hat_R
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k_loop, l_loop, q_loop, inv_ds, flux_face1]')
                do j = 1, sys_size
                    do q_loop = 0, p
                        do l_loop = 0, n
                            do k_loop = 0, m
                                inv_ds = 1._wp/dx(k_loop)
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

            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
        case (2)  ! y-direction
            if (bc_y%beg <= BC_CHAR_SLIP_WALL .and. bc_y%beg >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
            end if
            if (bc_y%end <= BC_CHAR_SLIP_WALL .and. bc_y%end >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
            end if

            if (.not. hypo_nc_dual_pass) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face1, flux_face2]')
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                inv_ds = 1._wp/dy(k)
                                flux_face1 = flux_n(2)%vf(j)%sf(q, k - 1, l)
                                flux_face2 = flux_n(2)%vf(j)%sf(q, k, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) + inv_ds*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else if (is_hat_L) then
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face2]')
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                inv_ds = 1._wp/dy(k)
                                flux_face2 = flux_n(2)%vf(j)%sf(q, k, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - inv_ds*flux_face2
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! is_hat_R
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face1]')
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                inv_ds = 1._wp/dy(k)
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
                geom_fac = merge(0.25_wp, 5.e-1_wp, hypo_nc_dual_pass)
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, flux_face1, flux_face2]')
                do j = 1, sys_size
                    do l = 0, p
                        do k = 0, n
                            do q = 0, m
                                flux_face1 = flux_gsrc_n(2)%vf(j)%sf(q, k - 1, l)
                                flux_face2 = flux_gsrc_n(2)%vf(j)%sf(q, k, l)
                                rhs_vf(j)%sf(q, k, l) = rhs_vf(j)%sf(q, k, l) - geom_fac/y_cc(k)*(flux_face1 + flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
        case (3)  ! z-direction
            if (bc_z%beg <= BC_CHAR_SLIP_WALL .and. bc_z%beg >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, -1, irx, iry, irz)
            end if
            if (bc_z%end <= BC_CHAR_SLIP_WALL .and. bc_z%end >= BC_CHAR_SUP_OUTFLOW) then
                call s_cbc(q_prim_vf%vf, flux_n(idir)%vf, flux_src_n_vf%vf, idir, 1, irx, iry, irz)
            end if

            if (grid_geometry == 3) then  ! Cylindrical Coordinates
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, velocity_val, flux_face1, flux_face2]')
                do j = 1, sys_size
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                inv_ds = 1._wp/(dz(k)*y_cc(q))
                                velocity_val = q_prim_vf%vf(eqn_idx%cont%end + idir)%sf(l, q, k)
                                flux_face1 = flux_n(3)%vf(j)%sf(l, q, k - 1)
                                flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*velocity_val*(flux_face1 - flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, flux_face1, flux_face2]')
                do j = 1, sys_size
                    do k = 0, p
                        do q = 0, n
                            do l = 0, m
                                flux_face1 = flux_gsrc_n(3)%vf(j)%sf(l, q, k - 1)
                                flux_face2 = flux_gsrc_n(3)%vf(j)%sf(l, q, k)
                                rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) - 5.e-1_wp/y_cc(q)*(flux_face1 + flux_face2)
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            else  ! Cartesian Coordinates
                if (.not. hypo_nc_dual_pass) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face1, flux_face2]')
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    inv_ds = 1._wp/dz(k)
                                    flux_face1 = flux_n(3)%vf(j)%sf(l, q, k - 1)
                                    flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                    rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) + inv_ds*(flux_face1 - flux_face2)
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else if (is_hat_L) then
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face2]')
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    inv_ds = 1._wp/dz(k)
                                    flux_face2 = flux_n(3)%vf(j)%sf(l, q, k)
                                    rhs_vf(j)%sf(l, q, k) = rhs_vf(j)%sf(l, q, k) - inv_ds*flux_face2
                                end do
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                else  ! is_hat_R
                    $:GPU_PARALLEL_LOOP(collapse=4,private='[j, k, l, q, inv_ds, flux_face1]')
                    do j = 1, sys_size
                        do k = 0, p
                            do q = 0, n
                                do l = 0, m
                                    inv_ds = 1._wp/dz(k)
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

            call s_add_directional_advection_source_terms(idir, rhs_vf, q_cons_vf, q_prim_vf, flux_src_n_vf, Kterm)
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

            !   adv_src_alpha_iface - HLL Method 1: flux_src carries per-fluid interface alpha_k

            !   adv_src_vel_iface   - HLLC, HLL Method 2, Exact, LF: flux_src carries shared face-normal velocity

            !   adv_src_none        - HLLD: MHD has no volume fractions; hypo uses dual-pass Riemann flux

            select case (current_idir)
            case (1)  ! x-direction
                if (adv_src_alpha_iface) then
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
                else if (adv_src_vel_iface) then
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
                if (adv_src_alpha_iface) then
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
                else if (adv_src_vel_iface) then
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
                    if (adv_src_alpha_iface) then
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
                    else if (adv_src_vel_iface) then
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
                else  ! non-cylindrical z-direction
                    if (adv_src_alpha_iface) then
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
                    else if (adv_src_vel_iface) then
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

    !> Add viscous, surface-tension, and species-diffusion source flux contributions to the RHS for a given direction
    subroutine s_compute_additional_physics_rhs(idir, q_prim_vf, rhs_vf, flux_src_n_in, dq_prim_dx_vf, dq_prim_dy_vf, dq_prim_dz_vf)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), dimension(sys_size), intent(in)    :: flux_src_n_in
        type(scalar_field), dimension(sys_size), intent(in)    :: dq_prim_dx_vf, dq_prim_dy_vf, dq_prim_dz_vf
        integer                                                :: i, j, k, l

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

            if ((surface_tension .or. viscous) .or. chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (surface_tension .or. viscous) then
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

                                if (.not. viscous) then
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
                    if (p > 0) then
                        call s_compute_viscous_stress_cylindrical_boundary(q_prim_vf, &
                            & dq_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), dq_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                            & dq_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                    else
                        call s_compute_viscous_stress_cylindrical_boundary(q_prim_vf, &
                            & dq_prim_dx_vf(eqn_idx%mom%beg:eqn_idx%mom%end), dq_prim_dy_vf(eqn_idx%mom%beg:eqn_idx%mom%end), &
                            & dq_prim_dz_vf(eqn_idx%mom%beg:eqn_idx%mom%end), tau_Re_vf, idwbuff(1), idwbuff(2), idwbuff(3))
                    end if

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
            else
                if ((surface_tension .or. viscous) .or. chem_params%diffusion) then
                    $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                if (surface_tension .or. viscous) then
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
                                    if (.not. viscous) then
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

                    if (viscous) then
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
                else
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

            if ((surface_tension .or. viscous) .or. chem_params%diffusion) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (surface_tension .or. viscous) then
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
                                if (.not. viscous) then
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

            if (grid_geometry == 3) then
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

                $:GPU_UPDATE(device='[is1, is2, is3, iv]')
            end if
        #:endfor

        if (recon_dir == 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv%beg, iv%end
                do l = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = is1%beg, is1%end
                            vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                            vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (recon_dir == 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv%beg, iv%end
                do l = is3%beg, is3%end
                    do j = is1%beg, is1%end
                        do k = is2%beg, is2%end
                            vL_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                            vR_x(k, j, l, i) = v_vf(i)%sf(k, j, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (recon_dir == 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = iv%beg, iv%end
                do j = is1%beg, is1%end
                    do k = is2%beg, is2%end
                        do l = is3%beg, is3%end
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

        call s_finalize_pressure_relaxation_module

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
            $:GPU_EXIT_DATA(delete='[alf_sum%sf]')
            deallocate (alf_sum%sf)
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
                        @:DEALLOCATE(flux_n(i)%vf(l)%sf)
                        @:DEALLOCATE(flux_gsrc_n(i)%vf(l)%sf)
                    end do

                    if (viscous) then
                        do l = eqn_idx%mom%beg, eqn_idx%E
                            @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        end do
                    end if

                    if (chem_params%diffusion .and. .not. viscous) then
                        @:DEALLOCATE(flux_src_n(i)%vf(eqn_idx%E)%sf)
                    end if

                    if (adv_src_alpha_iface .or. adv_src_none) then
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            @:DEALLOCATE(flux_src_n(i)%vf(l)%sf)
                        end do
                    else
                        do l = eqn_idx%adv%beg + 1, eqn_idx%adv%end
                            nullify (flux_src_n(i)%vf(l)%sf)
                        end do
                    end if

                    @:DEALLOCATE(flux_src_n(i)%vf(eqn_idx%adv%beg)%sf)
                end if

                @:DEALLOCATE(flux_n(i)%vf, flux_src_n(i)%vf, flux_gsrc_n(i)%vf)
            end do

            @:DEALLOCATE(flux_n, flux_src_n, flux_gsrc_n)
        end if

        if (use_nc_iface_vel) then
            do i = 1, num_dims
                do l = 1, num_dims
                    @:DEALLOCATE(nc_iface_vel_n(i)%vf(l)%sf)
                end do
                @:DEALLOCATE(nc_iface_vel_n(i)%vf)
            end do
            @:DEALLOCATE(nc_iface_vel_n)
            if (hypo_nc_dual_pass) then
                do i = 1, num_dims
                    do l = 1, num_dims
                        @:DEALLOCATE(nc_iface_vel_hatR_n(i)%vf(l)%sf)
                    end do
                    @:DEALLOCATE(nc_iface_vel_hatR_n(i)%vf)
                end do
                @:DEALLOCATE(nc_iface_vel_hatR_n)
            end if
        end if

        if (hypo_nc_dual_pass) then
            do i = 1, sys_size
                @:DEALLOCATE(rhs_hatL_vf(i)%sf)
                @:DEALLOCATE(rhs_hatR_vf(i)%sf)
            end do
            @:DEALLOCATE(rhs_hatL_vf)
            @:DEALLOCATE(rhs_hatR_vf)
        end if

    end subroutine s_finalize_rhs_module

end module m_rhs
