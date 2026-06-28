!>
!!@file
!!@brief Contains module m_active_box

#:include 'macros.fpp'

!> @brief Causal-envelope active-box restriction of the RHS compute window.
module m_active_box

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy

    implicit none

    private
    public :: s_initialize_active_box_module, s_finalize_active_box_module, s_initialize_active_box, s_grow_active_box, &
        & s_check_active_box_envelope, ab_x, ab_y, ab_z, ab_active, ab_ambient

    type(int_bounds_info) :: ab_x, ab_y, ab_z    !< Active-box interior cell ranges
    logical               :: ab_active           !< Whether the optimization is engaged
    real(wp), allocatable :: ab_ambient(:)       !< Uniform ambient conserved state
    real(wp), parameter   :: tol_ab = 1.e-10_wp  !< Ambient-deviation threshold

    $:GPU_DECLARE(create='[ab_x, ab_y, ab_z, ab_active]')

contains

    impure subroutine s_initialize_active_box_module

        @:ALLOCATE(ab_ambient(1:sys_size))
        ab_x%beg = 0; ab_x%end = m
        ab_y%beg = 0; ab_y%end = n
        ab_z%beg = 0; ab_z%end = p
        ab_active = .false.
        $:GPU_UPDATE(device='[ab_x, ab_y, ab_z, ab_active]')

    end subroutine s_initialize_active_box_module

    impure subroutine s_finalize_active_box_module

        @:DEALLOCATE(ab_ambient)

    end subroutine s_finalize_active_box_module

    !> Stub in Task 1: leaves ab_active=.false. (full-domain). Filled in Task 2.
    impure subroutine s_initialize_active_box(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

        ab_active = .false.
        $:GPU_UPDATE(device='[ab_active]')

    end subroutine s_initialize_active_box

    !> Stub in Task 1: no-op. Filled in Task 3.
    impure subroutine s_grow_active_box(dt_in)

        real(wp), intent(in) :: dt_in

    end subroutine s_grow_active_box

    !> Stub in Task 1: no-op. Filled in Task 6.
    impure subroutine s_check_active_box_envelope(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

    end subroutine s_check_active_box_envelope

end module m_active_box
