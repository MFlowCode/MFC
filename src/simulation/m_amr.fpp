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
    public :: t_level, amr_fine, s_initialize_amr_module, s_populate_amr_fine, s_interpolate_coarse_to_fine, &
        & s_restrict_fine_to_coarse, s_amr_conservation_check, s_finalize_amr_module

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

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the ref_ratio^d children average to the coarse value.
    impure subroutine s_interpolate_coarse_to_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        integer                                             :: i, fi, fj, fk, ci, cj, ck
        real(wp)                                            :: u0, sx, sy, sz, xix, xiy, xiz

        do i = 1, sys_size
            do fk = 0, amr_fine%p
                ck = amr_patch_beg(3) + fk/amr_fine%ref_ratio; if (p_glb == 0) ck = 0
                xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fj = 0, amr_fine%n
                    cj = amr_patch_beg(2) + fj/amr_fine%ref_ratio; if (n_glb == 0) cj = 0
                    xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    do fi = 0, amr_fine%m
                        ci = amr_patch_beg(1) + fi/amr_fine%ref_ratio
                        xix = (real(mod(fi, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                        u0 = real(q_cons_base(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(q_cons_base(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_cons_base(i)%sf(ci - 1, cj, &
                                    & ck), wp))
                        sy = 0._wp
                        if (n_glb > 0) sy = minmod(real(q_cons_base(i)%sf(ci, cj + 1, ck), wp) - u0, &
                            & u0 - real(q_cons_base(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0) sz = minmod(real(q_cons_base(i)%sf(ci, cj, ck + 1), wp) - u0, &
                            & u0 - real(q_cons_base(i)%sf(ci, cj, ck - 1), wp))
                        amr_fine%q_cons(i)%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                    end do
                end do
            end do
        end do

    end subroutine s_interpolate_coarse_to_fine

    !> Dispatch prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base

        if (.not. amr) return
        call s_interpolate_coarse_to_fine(q_cons_base)

    end subroutine s_populate_amr_fine

    !> Volume-weighted restriction: each covered coarse cell = average of its ref_ratio^d fine children. Writes ONLY the caller's
    !! coarse target (SP2: a scratch buffer) - never level-0.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer                                                :: i, ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild
        real(wp)                                               :: acc

        nchild = amr_fine%ref_ratio
        if (n_glb > 0) nchild = nchild*amr_fine%ref_ratio
        if (p_glb > 0) nchild = nchild*amr_fine%ref_ratio
        do i = 1, sys_size
            do ck = amr_patch_beg(3), merge(amr_patch_end(3), amr_patch_beg(3), p_glb > 0)
                fk0 = (ck - amr_patch_beg(3))*amr_fine%ref_ratio
                do cj = amr_patch_beg(2), merge(amr_patch_end(2), amr_patch_beg(2), n_glb > 0)
                    fj0 = (cj - amr_patch_beg(2))*amr_fine%ref_ratio
                    do ci = amr_patch_beg(1), amr_patch_end(1)
                        fi0 = (ci - amr_patch_beg(1))*amr_fine%ref_ratio
                        acc = 0._wp
                        do ddk = 0, merge(amr_fine%ref_ratio - 1, 0, p_glb > 0)
                            do ddj = 0, merge(amr_fine%ref_ratio - 1, 0, n_glb > 0)
                                acc = acc + real(amr_fine%q_cons(i)%sf(fi0, fj0 + ddj, fk0 + ddk), &
                                                 & wp) + real(amr_fine%q_cons(i)%sf(fi0 + 1, fj0 + ddj, fk0 + ddk), wp)
                            end do
                        end do
                        coarse_tgt(i)%sf(ci, cj, ck) = acc/real(nchild, wp)
                    end do
                end do
            end do
        end do

    end subroutine s_restrict_fine_to_coarse

    !> SP2 gate: restrict(prolong(coarse)) must reproduce coarse over the patch interior (conservation). Init-only diagnostic;
    !! allocates a scratch coarse target, never touches level-0 or the solve.
    impure subroutine s_amr_conservation_check(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        type(scalar_field), dimension(:), allocatable       :: scratch
        integer                                             :: i, ci, cj, ck
        real(wp)                                            :: err, e

        if (.not. amr) return
        allocate (scratch(1:sys_size))
        do i = 1, sys_size
            allocate (scratch(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do
        call s_restrict_fine_to_coarse(scratch)
        err = 0._wp
        do i = 1, sys_size
            do ck = amr_patch_beg(3), merge(amr_patch_end(3), amr_patch_beg(3), p_glb > 0)
                do cj = amr_patch_beg(2), merge(amr_patch_end(2), amr_patch_beg(2), n_glb > 0)
                    do ci = amr_patch_beg(1), amr_patch_end(1)
                        e = abs(real(scratch(i)%sf(ci, cj, ck), wp) - real(q_cons_base(i)%sf(ci, cj, ck), wp))
                        if (e > err) err = e
                    end do
                end do
            end do
        end do
        if (proc_rank == 0) print '(A,ES12.4)', ' [amr] restrict-prolong conservation err = ', err
        do i = 1, sys_size
            deallocate (scratch(i)%sf)
        end do
        deallocate (scratch)

    end subroutine s_amr_conservation_check

    !> minmod slope limiter: 0 if a,b differ in sign, else the smaller-magnitude argument.
    pure elemental function minmod(a, b) result(m)

        real(wp), intent(in) :: a, b
        real(wp)             :: m

        if (a*b <= 0._wp) then
            m = 0._wp
        else if (abs(a) < abs(b)) then
            m = a
        else
            m = b
        end if

    end function minmod

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
