!>
!! @file m_initial_condition.f90
!! @brief Contains module m_initial_condition

!> @brief This module provides a platform that is analogous to constructive
!!              solid geometry techniques and in this way allows for the creation
!!              of a wide variety of initial conditions. Several 1D, 2D and 3D
!!              fundamental geometries are included that may further be combined
!!              into more complex shapes. This is achieved by carefully setting
!!              up the order in which the patches are laid out in the domain and
!!              specifying the priority that each patch has over the preceding
!!              ones. The resulting shapes may be identified both by the values
!!              of their primitive variables and the associated patch identities.
!!              Note that the user may choose to read in and modify a preexisting
!!              initial condition. The module m_start_up.f90 is responsible for
!!             reading in the relevant data files.
module m_initial_condition

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_helper

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    use m_patches

    use m_assign_variables

    use m_perturbation          ! Subroutines to perturb initial flow fields

    use m_chemistry

    use m_boundary_conditions

    implicit none

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables

    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables

    type(scalar_field) :: q_T_sf !< Temperature field

    type(integer_field), dimension(:, :), allocatable :: bc_type !< bc_type fields

    integer, allocatable, dimension(:, :, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

    type(integer_field) :: ib_markers !<
    !! Bookkepping variable used to track whether a given cell is within an
    !! immersed boundary. The default is 0, otherwise the value is assigned
    !! to the patch ID of the immersed boundary.

    type(levelset_field) :: levelset
    type(levelset_norm_field) :: levelset_norm

contains

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    impure subroutine s_initialize_initial_condition_module

        integer :: i, j, k, l !< generic loop iterators

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end))
            allocate (q_cons_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end))
        end do

        if (chemistry) then
            allocate (q_T_sf%sf(0:m, 0:n, 0:p))
        end if

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        if (ib) then
            allocate (ib_markers%sf(0:m, 0:n, 0:p))
            allocate (levelset%sf(0:m, 0:n, 0:p, 1:num_ibs))
            allocate (levelset_norm%sf(0:m, 0:n, 0:p, 1:num_ibs, 1:3))
            ib_markers%sf = 0
        end if

        if (qbmm .and. .not. polytropic) then
            !Allocate bubble pressure pb and vapor mass mv for non-polytropic qbmm at all quad nodes and R0 bins
            allocate (pb%sf(0:m, &
                            0:n, &
                            0:p, 1:nnode, 1:nb))
            allocate (mv%sf(0:m, &
                            0:n, &
                            0:p, 1:nnode, 1:nb))
        end if

        ! Setting default values for conservative and primitive variables so
        ! that in the case that the initial condition is wrongly laid out on
        ! the grid the simulation component will catch the problem on start-
        ! up. The conservative variables do not need to be similarly treated
        ! since they are computed directly from the primitive variables.
        do i = 1, sys_size
            q_cons_vf(i)%sf = dflt_real
            q_prim_vf(i)%sf = dflt_real
        end do

        ! Allocating arrays to store the bc types
        allocate (bc_type(1:num_dims, -1:1))

        allocate (bc_type(1, -1)%sf(0:0, 0:n, 0:p))
        allocate (bc_type(1, 1)%sf(0:0, 0:n, 0:p))

        do l = 0, p
            do k = 0, n
                bc_type(1, -1)%sf(0, k, l) = bc_x%beg
                bc_type(1, 1)%sf(0, k, l) = bc_x%end
            end do
        end do

        if (n > 0) then
            allocate (bc_type(2, -1)%sf(-buff_size:m + buff_size, 0:0, 0:p))
            allocate (bc_type(2, 1)%sf(-buff_size:m + buff_size, 0:0, 0:p))

            do l = 0, p
                do j = -buff_size, m + buff_size
                    bc_type(2, -1)%sf(j, 0, l) = bc_y%beg
                    bc_type(2, 1)%sf(j, 0, l) = bc_y%end
                end do
            end do

            if (p > 0) then
                allocate (bc_type(3, -1)%sf(-buff_size:m + buff_size, -buff_size:n + buff_size, 0:0))
                allocate (bc_type(3, 1)%sf(-buff_size:m + buff_size, -buff_size:n + buff_size, 0:0))

                do k = -buff_size, n + buff_size
                    do j = -buff_size, m + buff_size
                        bc_type(3, -1)%sf(j, k, 0) = bc_z%beg
                        bc_type(3, 1)%sf(j, k, 0) = bc_z%end
                    end do
                end do
            end if
        end if

        ! Initial damage state is always zero
        if (cont_damage) then
            q_cons_vf(damage_idx)%sf = 0._wp
            q_prim_vf(damage_idx)%sf = 0._wp
        end if

        ! Setting default values for patch identities bookkeeping variable.
        ! This is necessary to avoid any confusion in the assessment of the
        ! extent of application that the overwrite permissions give a patch
        ! when it is being applied in the domain.
        patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module

    !>  This subroutine peruses the patches and depending on the
        !!              type of geometry associated with a particular patch, it
        !!              calls the related subroutine to setup the said geometry
        !!              on the grid using the primitive variables included with
        !!              the patch parameters. The subroutine is complete once the
        !!              primitive variables are converted to conservative ones.
    impure subroutine s_generate_initial_condition

        ! Converting the conservative variables to the primitive ones given
        ! preexisting initial condition data files were read in on start-up
        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                               q_T_sf, &
                                                               q_prim_vf, &
                                                               idwbuff)
        end if

        if (ib) then
            call s_apply_domain_patches(patch_id_fp, q_prim_vf, ib_markers%sf, levelset, levelset_norm)
        else
            call s_apply_domain_patches(patch_id_fp, q_prim_vf)
        end if

        if (num_bc_patches > 0) call s_apply_boundary_patches(q_prim_vf, bc_type)

        if (perturb_flow) call s_perturb_surrounding_flow(q_prim_vf)
        if (perturb_sph) call s_perturb_sphere(q_prim_vf)
        if (mixlayer_perturb) call s_perturb_mixlayer(q_prim_vf)
        if (elliptic_smoothing) call s_elliptic_smoothing(q_prim_vf, bc_type)

        ! Converting the primitive variables to the conservative ones
        call s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        if (chemistry) call s_compute_T_from_primitives(q_T_sf, q_prim_vf, idwint)

        if (qbmm .and. .not. polytropic) then
            !Initialize pb and mv
            call s_initialize_mv(q_cons_vf, mv%sf)
            call s_initialize_pb(q_cons_vf, mv%sf, pb%sf)
        end if

    end subroutine s_generate_initial_condition

    !>  Deallocation procedures for the module
    impure subroutine s_finalize_initial_condition_module

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        if (chemistry) then
            deallocate (q_T_sf%sf)
        end if

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)

        if (ib) then
            deallocate (ib_markers%sf, levelset%sf, levelset_norm%sf)
        end if

    end subroutine s_finalize_initial_condition_module

end module m_initial_condition
