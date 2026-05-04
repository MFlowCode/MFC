!>
!! @file
!! @brief Contains module m_initial_condition

!> @brief Assembles initial conditions by layering prioritized patches via constructive solid geometry
module m_initial_condition

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_variables_conversion
    use m_icpp_patches
    use m_assign_variables
    use m_perturbation
    use m_chemistry
    use m_boundary_conditions

    implicit none

    ! NOTE: Abstract interface enables dynamic dispatch without repeated model_eqns checks
    type(scalar_field), allocatable, dimension(:)    :: q_prim_vf  !< primitive variables
    type(scalar_field), allocatable, dimension(:)    :: q_cons_vf  !< conservative variables
    type(scalar_field)                               :: q_T_sf     !< Temperature field
    type(integer_field), dimension(:,:), allocatable :: bc_type    !< bc_type fields
    !> @cond
#ifdef MFC_MIXED_PRECISION
    integer(kind=1), allocatable, dimension(:,:,:) :: patch_id_fp
#else
    !> @endcond
    integer, allocatable, dimension(:,:,:) :: patch_id_fp
    !> @cond
#endif
    !> @endcond

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_initial_condition_module

        integer :: i, j, k, l

        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
            allocate (q_cons_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do

        if (chemistry) then
            allocate (q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end if

        allocate (patch_id_fp(0:m,0:n,0:p))

        if (qbmm .and. .not. polytropic) then
            allocate (pb%sf(0:m,0:n,0:p,1:nnode,1:nb))
            allocate (mv%sf(0:m,0:n,0:p,1:nnode,1:nb))
        end if

        do i = 1, sys_size
            q_cons_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp) ! TODO :: remove this magic number
            q_prim_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp)
        end do

        allocate (bc_type(1:num_dims,1:2))

        allocate (bc_type(1, 1)%sf(0:0,0:n,0:p))
        allocate (bc_type(1, 2)%sf(0:0,0:n,0:p))

        do l = 0, p
            do k = 0, n
                bc_type(1, 1)%sf(0, k, l) = int(min(bc_x%beg, 0), kind=1)
                bc_type(1, 2)%sf(0, k, l) = int(min(bc_x%end, 0), kind=1)
            end do
        end do

        if (n > 0) then
            allocate (bc_type(2, 1)%sf(-buff_size:m + buff_size,0:0,0:p))
            allocate (bc_type(2, 2)%sf(-buff_size:m + buff_size,0:0,0:p))

            do l = 0, p
                do j = -buff_size, m + buff_size
                    bc_type(2, 1)%sf(j, 0, l) = int(min(bc_y%beg, 0), kind=1)
                    bc_type(2, 2)%sf(j, 0, l) = int(min(bc_y%end, 0), kind=1)
                end do
            end do

            if (p > 0) then
                allocate (bc_type(3, 1)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (bc_type(3, 2)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))

                do k = -buff_size, n + buff_size
                    do j = -buff_size, m + buff_size
                        bc_type(3, 1)%sf(j, k, 0) = int(min(bc_z%beg, 0), kind=1)
                        bc_type(3, 2)%sf(j, k, 0) = int(min(bc_z%end, 0), kind=1)
                    end do
                end do
            end if
        end if

        ! Initial damage state is always zero
        if (cont_damage) then
            q_cons_vf(eqn_idx%damage)%sf = 0._wp
            q_prim_vf(eqn_idx%damage)%sf = 0._wp
        end if

        ! Initial hyper_cleaning state is always zero TODO more general
        if (hyper_cleaning) then
            q_cons_vf(eqn_idx%psi)%sf = 0._wp
            q_prim_vf(eqn_idx%psi)%sf = 0._wp
        end if

        ! Setting default values for patch identities bookkeeping variable. This is necessary to avoid any confusion in the
        ! assessment of the extent of application that the overwrite permissions give a patch when it is being applied in the
        ! domain.
        patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module

    !> Iterate over patches and, depending on the geometry type, call the related subroutine to setup the said geometry on the grid
    !! using the primitive variables included with the patch parameters. The subroutine is complete once the primitive variables are
    !! converted to conservative ones.
    impure subroutine s_generate_initial_condition

        integer :: i

        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwbuff)
        end if

        call s_apply_icpp_patches(patch_id_fp, q_prim_vf)

        if (num_bc_patches > 0) call s_apply_boundary_patches(q_prim_vf, bc_type)

        if (perturb_flow) call s_perturb_surrounding_flow(q_prim_vf)
        if (perturb_sph) call s_perturb_sphere(q_prim_vf)
        if (mixlayer_perturb) call s_perturb_mixlayer(q_prim_vf)
        if (simplex_perturb) call s_perturb_simplex(q_prim_vf)
        if (chemistry) call s_compute_T_from_primitives(q_T_sf, q_prim_vf, idwint)

        if (elliptic_smoothing .and. chemistry) then
            call s_elliptic_smoothing(q_prim_vf, bc_type, q_T_sf)
            call s_compute_T_from_primitives(q_T_sf, q_prim_vf, idwint)
        else if (elliptic_smoothing) then
            call s_elliptic_smoothing(q_prim_vf, bc_type)
        end if

        call s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        if (qbmm .and. .not. polytropic) then
            call s_initialize_mv(q_cons_vf, mv%sf)
            call s_initialize_pb(q_cons_vf, mv%sf, pb%sf)
        end if

    end subroutine s_generate_initial_condition

    !> Deallocation procedures for the module
    impure subroutine s_finalize_initial_condition_module

        integer :: i

        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        if (chemistry) then
            deallocate (q_T_sf%sf)
        end if

        deallocate (patch_id_fp)

        deallocate (bc_type(1, 1)%sf)
        deallocate (bc_type(1, 2)%sf)

        if (n > 0) then
            deallocate (bc_type(2, 1)%sf)
            deallocate (bc_type(2, 2)%sf)
        end if

        if (p > 0) then
            deallocate (bc_type(3, 1)%sf)
            deallocate (bc_type(3, 2)%sf)
        end if

        deallocate (bc_type)

    end subroutine s_finalize_initial_condition_module

end module m_initial_condition
