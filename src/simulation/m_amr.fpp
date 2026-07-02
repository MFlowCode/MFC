!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief AMR hierarchy (SP1: one static, inert refined level-1 patch alongside the base solve).
module m_amr

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters

    implicit none

    private
    public :: t_level, amr_fine, s_initialize_amr_module, s_populate_amr_fine, s_finalize_amr_module

    !> One refined level: its own grid + conservative fields. Host-only (no GPU in SP1).
    type t_level
        integer                         :: ref_ratio
        type(t_box)                     :: region   !< patch extent in parent (level-0) cell indices
        integer                         :: m, n, p  !< this level's interior extents
        integer                         :: buff_size
        type(int_bounds_info)           :: idwbuff(3)
        real(wp), allocatable           :: x_cb(:), x_cc(:), dx(:)
        real(wp), allocatable           :: y_cb(:), y_cc(:), dy(:)
        real(wp), allocatable           :: z_cb(:), z_cc(:), dz(:)
        type(scalar_field), allocatable :: q_cons(:)
    end type t_level

    type(t_level) :: amr_fine

contains

    !> Build the static refined level-1 patch. No-op unless amr. Called after level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set).
    impure subroutine s_initialize_amr_module()

        integer :: i

        if (.not. amr) return

        amr_fine%ref_ratio = 2
        amr_fine%region%lo = amr_patch_beg
        amr_fine%region%hi = amr_patch_end
        amr_fine%buff_size = buff_size

        ! interior extents: ref_ratio*(coarse cells in patch) - 1 per active dim; collapsed dims stay 0
        amr_fine%m = amr_fine%ref_ratio*(amr_patch_end(1) - amr_patch_beg(1) + 1) - 1
        amr_fine%n = 0; amr_fine%p = 0
        if (n_glb > 0) amr_fine%n = amr_fine%ref_ratio*(amr_patch_end(2) - amr_patch_beg(2) + 1) - 1
        if (p_glb > 0) amr_fine%p = amr_fine%ref_ratio*(amr_patch_end(3) - amr_patch_beg(3) + 1) - 1

        amr_fine%idwbuff(1)%beg = -buff_size; amr_fine%idwbuff(1)%end = amr_fine%m + buff_size
        amr_fine%idwbuff(2)%beg = -buff_size; amr_fine%idwbuff(2)%end = amr_fine%n + buff_size
        amr_fine%idwbuff(3)%beg = -buff_size; amr_fine%idwbuff(3)%end = amr_fine%p + buff_size

        ! level-1 coordinates by bisecting each covered coarse cell (handles stretched grids)
        call s_build_level_coords(x_cb, lbound(x_cb, 1), amr_patch_beg(1), amr_fine%m, amr_fine%x_cb, amr_fine%x_cc, amr_fine%dx)
        if (n_glb > 0) call s_build_level_coords(y_cb, lbound(y_cb, 1), amr_patch_beg(2), amr_fine%n, amr_fine%y_cb, &
            & amr_fine%y_cc, amr_fine%dy)
        if (p_glb > 0) call s_build_level_coords(z_cb, lbound(z_cb, 1), amr_patch_beg(3), amr_fine%p, amr_fine%z_cb, &
            & amr_fine%z_cc, amr_fine%dz)

        allocate (amr_fine%q_cons(1:sys_size))
        do i = 1, sys_size
            allocate (amr_fine%q_cons(i)%sf(amr_fine%idwbuff(1)%beg:amr_fine%idwbuff(1)%end, &
                      & amr_fine%idwbuff(2)%beg:amr_fine%idwbuff(2)%end,amr_fine%idwbuff(3)%beg:amr_fine%idwbuff(3)%end))
        end do

    end subroutine s_initialize_amr_module

    !> Fill level-1 fcb/fcc/fdx by bisecting parent cells; pcb_lb is lbound(parent_cb, 1). Passing pcb as assumed-shape resets
    !! lbound to 1; pcb_lb + idx_offset recovers original indexing.
    subroutine s_build_level_coords(pcb, pcb_lb, lo, nfine, fcb, fcc, fdx)

        real(wp), intent(in)               :: pcb(:)
        integer, intent(in)                :: pcb_lb, lo, nfine
        real(wp), allocatable, intent(out) :: fcb(:), fcc(:), fdx(:)
        integer                            :: fi, c, idx_offset
        real(wp)                           :: xl, xr, xm
        ! pcb(k) = parent_cb(k + pcb_lb - 1); to access parent_cb(j): k = j - pcb_lb + 1

        idx_offset = 1 - pcb_lb
        allocate (fcb(-1:nfine), fcc(0:nfine), fdx(0:nfine))
        ! fine cell fi (0..nfine) bisects coarse cell c = lo + fi/2
        do fi = 0, nfine
            c = lo + fi/2
            xl = pcb(c - 1 + idx_offset)  ! left boundary of coarse cell c
            xr = pcb(c + idx_offset)  ! right boundary of coarse cell c
            xm = 0.5_wp*(xl + xr)
            if (mod(fi, 2) == 0) then
                fcb(fi - 1) = xl
                fcb(fi) = xm
            else
                fcb(fi) = xr
            end if
        end do
        do fi = 0, nfine
            fdx(fi) = fcb(fi) - fcb(fi - 1)
            fcc(fi) = 0.5_wp*(fcb(fi - 1) + fcb(fi))
        end do

    end subroutine s_build_level_coords

    !> Piecewise-constant injection: each fine interior cell copies its covering coarse cell (level-0 read-only). Level-0 is never
    !! modified; the level-1 patch is inert and never fed back into the base solve.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        integer                                             :: i, fi, fj, fk, ci, cj, ck

        if (.not. amr) return
        do i = 1, sys_size
            do fk = 0, amr_fine%p
                ck = amr_patch_beg(3) + fk/amr_fine%ref_ratio
                if (p_glb == 0) ck = 0
                do fj = 0, amr_fine%n
                    cj = amr_patch_beg(2) + fj/amr_fine%ref_ratio
                    if (n_glb == 0) cj = 0
                    do fi = 0, amr_fine%m
                        ci = amr_patch_beg(1) + fi/amr_fine%ref_ratio
                        amr_fine%q_cons(i)%sf(fi, fj, fk) = q_cons_base(i)%sf(ci, cj, ck)
                    end do
                end do
            end do
        end do

    end subroutine s_populate_amr_fine

    impure subroutine s_finalize_amr_module()

        integer :: i

        if (.not. amr) return
        do i = 1, sys_size
            if (associated(amr_fine%q_cons(i)%sf)) deallocate (amr_fine%q_cons(i)%sf)
        end do
        deallocate (amr_fine%q_cons)
        if (allocated(amr_fine%x_cb)) deallocate (amr_fine%x_cb, amr_fine%x_cc, amr_fine%dx)
        if (allocated(amr_fine%y_cb)) deallocate (amr_fine%y_cb, amr_fine%y_cc, amr_fine%dy)
        if (allocated(amr_fine%z_cb)) deallocate (amr_fine%z_cb, amr_fine%z_cc, amr_fine%dz)

    end subroutine s_finalize_amr_module

end module m_amr
