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

    type(int_bounds_info) :: ab_x, ab_y, ab_z     !< Active-box interior cell ranges
    logical               :: ab_active = .false.  !< Whether the optimization is engaged
    real(wp), allocatable :: ab_ambient(:)        !< Uniform ambient conserved state
    real(wp), parameter   :: tol_ab = 1.e-10_wp   !< Ambient-deviation threshold

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

    !> Detect the ambient state and set the initial active-box bounds from the IC support.
    impure subroutine s_initialize_active_box(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer                                             :: i, j, k, l
        integer                                             :: ib, ie, jb, je, kb, ke
        logical                                             :: deviates

        if (.not. active_box .or. num_procs /= 1) then
            ab_active = .false.
            $:GPU_UPDATE(device='[ab_active]')
            return
        end if

        ! Ambient = the (0,0,0) interior corner cell, assumed in the undisturbed region.
        do i = 1, sys_size
            ab_ambient(i) = q_cons_vf(i)%sf(0, 0, 0)
        end do

        ! Bounding box of cells deviating from ambient.
        ib = m + 1; ie = -1; jb = n + 1; je = -1; kb = p + 1; ke = -1
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    deviates = .false.
                    do i = 1, sys_size
                        if (abs(q_cons_vf(i)%sf(j, k, l) - ab_ambient(i)) > tol_ab) then
                            deviates = .true.; exit
                        end if
                    end do
                    if (deviates) then
                        ib = min(ib, j); ie = max(ie, j)
                        jb = min(jb, k); je = max(je, k)
                        kb = min(kb, l); ke = max(ke, l)
                    end if
                end do
            end do
        end do

        ! Empty deviation set -> nothing to do; disable.
        if (ie < ib) then
            ab_active = .false.
            $:GPU_UPDATE(device='[ab_active]')
            return
        end if

        ! Dilate by the reconstruction stencil and clamp to the interior.
        ab_x%beg = max(0, ib - buff_size); ab_x%end = min(m, ie + buff_size)
        ab_y%beg = max(0, jb - buff_size); ab_y%end = min(n, je + buff_size)
        ab_z%beg = max(0, kb - buff_size); ab_z%end = min(p, ke + buff_size)

        ! If the box already covers the whole domain there is no benefit; disable.
        ab_active = .not. (ab_x%beg == 0 .and. ab_x%end == m .and. ab_y%beg == 0 .and. ab_y%end == n .and. ab_z%beg == 0 &
                           & .and. ab_z%end == p)

        $:GPU_UPDATE(device='[ab_x, ab_y, ab_z, ab_active]')

#ifdef MFC_DEBUG
        if (ab_active) then
            print *, '[active_box] init box x[', ab_x%beg, ':', ab_x%end, '] y[', ab_y%beg, ':', ab_y%end, '] z[', ab_z%beg, ':', &
                & ab_z%end, ']'
        end if
#endif

    end subroutine s_initialize_active_box

    !> Grow the active box by one light-cone step and refresh the device copy.
    impure subroutine s_grow_active_box()

        integer :: g

        if (.not. ab_active) return

        ! Growth by buff_size cells/step outruns the physical front: a stable run has CFL <= ~1.4
        ! cells/step, so the box edge progressively leads the physical disturbance and sits in the
        ! exponentially-decaying, sub-tolerance numerical precursor. This gives agreement to
        ! round-off (~1e-14), not bit-identical, which is what the spec requires.
        ! Caveat: this is a finite-horizon round-off guarantee, not a strict numerical-light-cone
        ! containment. A bit-identical variant would grow by nstage*(weno_polyn+1) = ~9 cells/step
        ! (the full numerical domain of dependence) at the cost of a looser box and lower speedup -
        ! left as a future option. Under-growth (CFL > buff_size) implies an already-diverging run.
        g = buff_size

        ab_x%beg = max(0, ab_x%beg - g); ab_x%end = min(m, ab_x%end + g)
        ab_y%beg = max(0, ab_y%beg - g); ab_y%end = min(n, ab_y%end + g)
        ab_z%beg = max(0, ab_z%beg - g); ab_z%end = min(p, ab_z%end + g)

        ! Once the box fills the domain, disable the optimization (full-domain is correct and avoids the bookkeeping).
        if (ab_x%beg == 0 .and. ab_x%end == m .and. ab_y%beg == 0 .and. ab_y%end == n .and. ab_z%beg == 0 .and. ab_z%end == p) then
            ab_active = .false.
        end if

        $:GPU_UPDATE(device='[ab_x, ab_y, ab_z, ab_active]')

    end subroutine s_grow_active_box

    !> Abort in debug builds if the disturbance has reached the active-box boundary (under-growth).
    impure subroutine s_check_active_box_envelope(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

#ifdef MFC_DEBUG
        integer :: i, j, k, l
        logical :: in_margin
        if (.not. ab_active) return
#ifdef MFC_GPU
        ! Refresh host copy before reading conserved fields on CPU.
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_vf(i)%sf]')
        end do
#endif
        ! Check the buff_size-thick layer just INSIDE each box face. Those cells are updated by
        ! the RK loop, so if the disturbance front reaches them the reconstruction margin is
        ! compromised (box under-grown). The exterior layer is frozen-ambient by construction and
        ! cannot detect under-growth. Degenerate dimensions (n=0/p=0) contribute no faces.
        do l = ab_z%beg, ab_z%end
            do k = ab_y%beg, ab_y%end
                do j = ab_x%beg, ab_x%end
                    in_margin = (j <= ab_x%beg + buff_size - 1 .or. j >= ab_x%end - buff_size + 1)
                    if (n > 0) in_margin = in_margin .or. (k <= ab_y%beg + buff_size - 1 .or. k >= ab_y%end - buff_size + 1)
                    if (p > 0) in_margin = in_margin .or. (l <= ab_z%beg + buff_size - 1 .or. l >= ab_z%end - buff_size + 1)
                    if (.not. in_margin) cycle
                    do i = 1, sys_size
                        @:ASSERT(abs(q_cons_vf(i)%sf(j, k, l) - ab_ambient(i)) <= tol_ab, &
                                 & "active_box: disturbance reached the box boundary (under-grown)")
                    end do
                end do
            end do
        end do
#endif

    end subroutine s_check_active_box_envelope

end module m_active_box
