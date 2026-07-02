!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief AMR hierarchy (SP1: one static, inert refined level-1 patch alongside the base solve).
module m_amr

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort

    implicit none

    private
    public :: t_level, amr_fine, s_initialize_amr_module, s_populate_amr_fine, s_interpolate_coarse_to_fine, &
        & s_restrict_fine_to_coarse, s_amr_conservation_check, s_finalize_amr_module, s_amr_swap_to_fine, s_amr_restore_coarse, &
        & s_amr_fill_fine_ghosts, s_amr_operator_checks

    !> One refined level: its own grid + conservative fields. Host-only (no GPU in SP1).
    type t_level
        integer                         :: ref_ratio
        type(t_box)                     :: region          !< patch extent in parent (level-0) cell indices
        integer                         :: m, n, p         !< this level's interior extents
        integer                         :: buff_size
        type(int_bounds_info)           :: idwbuff(3)
        real(wp), allocatable           :: x_cb(:), x_cc(:), dx(:)
        real(wp), allocatable           :: y_cb(:), y_cc(:), dy(:)
        real(wp), allocatable           :: z_cb(:), z_cc(:), dz(:)
        type(scalar_field), allocatable :: q_cons(:)
        type(scalar_field), allocatable :: q_cons_stor(:)  !< stage storage (fine advance, same bounds as q_cons)
        type(scalar_field), allocatable :: q_prim(:)       !< primitive stage (fine advance, same bounds as q_cons)
        type(scalar_field), allocatable :: rhs(:)          !< RHS (fine interior only: 0:m, 0:n, 0:p)
    end type t_level

    type(t_level) :: amr_fine

    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)

contains

    !> Build the static refined level-1 patch. No-op unless amr. Called after level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set).
    impure subroutine s_initialize_amr_module()

        integer :: i

        if (.not. amr) return

        ! Runtime margin abort: patch must be at least buff_size coarse cells from every boundary.
        ! buff_size is not available at checker time, so this check must be here.
        if (amr_patch_beg(1) < buff_size .or. amr_patch_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_patch_beg(2) &
            & < buff_size .or. amr_patch_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_patch_beg(3) < buff_size &
            & .or. amr_patch_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr patch must be at least buff_size cells from every domain boundary')
        end if

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

        ! Allocate fine conservative fields (idwbuff shape for q_cons/q_cons_stor/q_prim; interior for rhs)
        allocate (amr_fine%q_cons(1:sys_size))
        allocate (amr_fine%q_cons_stor(1:sys_size))
        allocate (amr_fine%q_prim(1:sys_size))
        allocate (amr_fine%rhs(1:sys_size))
        do i = 1, sys_size
            allocate (amr_fine%q_cons(i)%sf(amr_fine%idwbuff(1)%beg:amr_fine%idwbuff(1)%end, &
                      & amr_fine%idwbuff(2)%beg:amr_fine%idwbuff(2)%end,amr_fine%idwbuff(3)%beg:amr_fine%idwbuff(3)%end))
            allocate (amr_fine%q_cons_stor(i)%sf(amr_fine%idwbuff(1)%beg:amr_fine%idwbuff(1)%end, &
                      & amr_fine%idwbuff(2)%beg:amr_fine%idwbuff(2)%end,amr_fine%idwbuff(3)%beg:amr_fine%idwbuff(3)%end))
            allocate (amr_fine%q_prim(i)%sf(amr_fine%idwbuff(1)%beg:amr_fine%idwbuff(1)%end, &
                      & amr_fine%idwbuff(2)%beg:amr_fine%idwbuff(2)%end,amr_fine%idwbuff(3)%beg:amr_fine%idwbuff(3)%end))
            allocate (amr_fine%rhs(i)%sf(0:amr_fine%m,0:amr_fine%n,0:amr_fine%p))
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

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, qf)

        type(scalar_field), intent(in)    :: qc
        type(scalar_field), intent(inout) :: qf
        integer                           :: fi, fj, fk, ci, cj, ck
        real(wp)                          :: u0, sx, sy, sz, xix, xiy, xiz

        do fk = 0, amr_fine%p
            ck = amr_patch_beg(3) + fk/amr_fine%ref_ratio; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_fine%n
                cj = amr_patch_beg(2) + fj/amr_fine%ref_ratio; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_fine%m
                    ci = amr_patch_beg(1) + fi/amr_fine%ref_ratio
                    xix = (real(mod(fi, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (n_glb > 0) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (p_glb > 0) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    qf%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                end do
            end do
        end do

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the ref_ratio^d children average to the coarse value.
    impure subroutine s_interpolate_coarse_to_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        integer                                             :: i

        do i = 1, sys_size
            call s_prolong_one_var(q_cons_base(i), amr_fine%q_cons(i))
        end do

    end subroutine s_interpolate_coarse_to_fine

    !> Dispatch prolongation. Guard: no-op unless amr.
    impure subroutine s_populate_amr_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base

        if (.not. amr) return
        call s_interpolate_coarse_to_fine(q_cons_base)

    end subroutine s_populate_amr_fine

    !> Volume-weighted restriction for a single variable pair. Reads from qf (fine, must include interior 0:amr_fine%m etc.); writes
    !! to qc (coarse, over the patch).
    impure subroutine s_restrict_one_var(qf, qc)

        type(scalar_field), intent(in)    :: qf
        type(scalar_field), intent(inout) :: qc
        integer                           :: ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild
        real(wp)                          :: acc

        nchild = amr_fine%ref_ratio
        if (n_glb > 0) nchild = nchild*amr_fine%ref_ratio
        if (p_glb > 0) nchild = nchild*amr_fine%ref_ratio
        do ck = amr_patch_beg(3), merge(amr_patch_end(3), amr_patch_beg(3), p_glb > 0)
            fk0 = (ck - amr_patch_beg(3))*amr_fine%ref_ratio
            do cj = amr_patch_beg(2), merge(amr_patch_end(2), amr_patch_beg(2), n_glb > 0)
                fj0 = (cj - amr_patch_beg(2))*amr_fine%ref_ratio
                do ci = amr_patch_beg(1), amr_patch_end(1)
                    fi0 = (ci - amr_patch_beg(1))*amr_fine%ref_ratio
                    acc = 0._wp
                    do ddk = 0, merge(amr_fine%ref_ratio - 1, 0, p_glb > 0)
                        do ddj = 0, merge(amr_fine%ref_ratio - 1, 0, n_glb > 0)
                            acc = acc + real(qf%sf(fi0, fj0 + ddj, fk0 + ddk), wp) + real(qf%sf(fi0 + 1, fj0 + ddj, fk0 + ddk), wp)
                        end do
                    end do
                    qc%sf(ci, cj, ck) = acc/real(nchild, wp)
                end do
            end do
        end do

    end subroutine s_restrict_one_var

    !> Volume-weighted restriction: each covered coarse cell = average of its ref_ratio^d fine children. Writes ONLY the caller's
    !! coarse target (SP2: a scratch buffer) - never level-0.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer                                                :: i

        do i = 1, sys_size
            call s_restrict_one_var(amr_fine%q_cons(i), coarse_tgt(i))
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

    !> Swap the global grid state to the fine patch. MUST be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        sw_m = m; sw_n = n; sw_p = p
        sw_idwint = idwint; sw_idwbuff = idwbuff
        m = amr_fine%m; n = amr_fine%n; p = amr_fine%p
        idwint(1)%beg = 0; idwint(1)%end = m
        idwint(2)%beg = 0; idwint(2)%end = n
        idwint(3)%beg = 0; idwint(3)%end = p
        idwbuff = amr_fine%idwbuff
        call s_swap_coords()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse()

        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        call s_swap_coords()

    end subroutine s_amr_restore_coarse

    !> Exchange global x/y/z coordinate arrays with amr_fine's via move_alloc (symmetric swap). Strategy: move_alloc. Pointer scan
    !! shows only local-scope pointers (s_cb in m_weno, s_cc in m_ibm); these are re-associated at each call entry and do not
    !! persist across subroutine calls. t_level coord arrays carry the target attribute to match the global declarations for
    !! portability.
    impure subroutine s_swap_coords()

        real(wp), allocatable :: tmp(:)

        call move_alloc(x_cb, tmp); call move_alloc(amr_fine%x_cb, x_cb); call move_alloc(tmp, amr_fine%x_cb)
        call move_alloc(x_cc, tmp); call move_alloc(amr_fine%x_cc, x_cc); call move_alloc(tmp, amr_fine%x_cc)
        call move_alloc(dx, tmp); call move_alloc(amr_fine%dx, dx); call move_alloc(tmp, amr_fine%dx)
        if (n_glb > 0) then
            call move_alloc(y_cb, tmp); call move_alloc(amr_fine%y_cb, y_cb); call move_alloc(tmp, amr_fine%y_cb)
            call move_alloc(y_cc, tmp); call move_alloc(amr_fine%y_cc, y_cc); call move_alloc(tmp, amr_fine%y_cc)
            call move_alloc(dy, tmp); call move_alloc(amr_fine%dy, dy); call move_alloc(tmp, amr_fine%dy)
        end if
        if (p_glb > 0) then
            call move_alloc(z_cb, tmp); call move_alloc(amr_fine%z_cb, z_cb); call move_alloc(tmp, amr_fine%z_cb)
            call move_alloc(z_cc, tmp); call move_alloc(amr_fine%z_cc, z_cc); call move_alloc(tmp, amr_fine%z_cc)
            call move_alloc(dz, tmp); call move_alloc(amr_fine%dz, dz); call move_alloc(tmp, amr_fine%dz)
        end if

    end subroutine s_swap_coords

    !> Fill the fine ghost shell of q_fine by conservative-linear prolongation from q_coarse. floor/modulo mapping is valid for
    !! negative fine indices (ghosts). Interior untouched.
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, q_fine)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_coarse
        type(scalar_field), dimension(sys_size), intent(inout) :: q_fine
        integer                                                :: i, fi, fj, fk, ci, cj, ck
        real(wp)                                               :: u0, sx, sy, sz, xix, xiy, xiz

        do i = 1, sys_size
            do fk = amr_fine%idwbuff(3)%beg, amr_fine%idwbuff(3)%end
                ck = 0; xiz = 0._wp
                if (p_glb > 0) then
                    ck = amr_patch_beg(3) + floor(real(fk, wp)/real(amr_fine%ref_ratio, wp))
                    xiz = (real(modulo(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                end if
                do fj = amr_fine%idwbuff(2)%beg, amr_fine%idwbuff(2)%end
                    cj = 0; xiy = 0._wp
                    if (n_glb > 0) then
                        cj = amr_patch_beg(2) + floor(real(fj, wp)/real(amr_fine%ref_ratio, wp))
                        xiy = (real(modulo(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    end if
                    do fi = amr_fine%idwbuff(1)%beg, amr_fine%idwbuff(1)%end
                        ! skip the interior: only the ghost shell is filled
                        if (fi >= 0 .and. fi <= amr_fine%m .and. fj >= 0 .and. fj <= amr_fine%n .and. fk >= 0 &
                            & .and. fk <= amr_fine%p) cycle
                        ci = amr_patch_beg(1) + floor(real(fi, wp)/real(amr_fine%ref_ratio, wp))
                        xix = (real(modulo(fi, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                        u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (n_glb > 0) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                            & cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                            & cj, ck - 1), wp))
                        q_fine(i)%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                    end do
                end do
            end do
        end do

    end subroutine s_amr_fill_fine_ghosts

    !> Init-time operator verification: (b) linear reproduction, (c) restriction of an independent field. Uses amr_fine%q_cons(1) as
    !! scratch; called before s_populate_amr_fine overwrites it.
    impure subroutine s_amr_operator_checks()

        type(scalar_field), allocatable :: cscr(:)
        integer                         :: fi, fj, fk, ci, cj, ck
        real(wp)                        :: e, errb, errc, si_f, si_c, dvf, dvc, want

        if (.not. amr) return

        ! (b) fill a coarse scratch with an exactly-linear field, prolong, compare pointwise
        allocate (cscr(1:1))
        allocate (cscr(1)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        do ck = idwbuff(3)%beg, idwbuff(3)%end
            do cj = idwbuff(2)%beg, idwbuff(2)%end
                do ci = idwbuff(1)%beg, idwbuff(1)%end
                    cscr(1)%sf(ci, cj, ck) = 1._wp + 2._wp*x_cc(ci)
                    if (n_glb > 0) cscr(1)%sf(ci, cj, ck) = cscr(1)%sf(ci, cj, ck) + 3._wp*y_cc(cj)
                    if (p_glb > 0) cscr(1)%sf(ci, cj, ck) = cscr(1)%sf(ci, cj, ck) + 4._wp*z_cc(ck)
                end do
            end do
        end do
        call s_prolong_one_var(cscr(1), amr_fine%q_cons(1))
        errb = 0._wp
        do fk = 0, amr_fine%p
            do fj = 0, amr_fine%n
                do fi = 0, amr_fine%m
                    want = 1._wp + 2._wp*amr_fine%x_cc(fi)
                    if (n_glb > 0) want = want + 3._wp*amr_fine%y_cc(fj)
                    if (p_glb > 0) want = want + 4._wp*amr_fine%z_cc(fk)
                    e = abs(real(amr_fine%q_cons(1)%sf(fi, fj, fk), wp) - want)
                    if (e > errb) errb = e
                end do
            end do
        end do

        ! (c) fill the fine patch with a quadratic (NOT from prolongation), restrict, compare integrals
        do fk = 0, amr_fine%p
            do fj = 0, amr_fine%n
                do fi = 0, amr_fine%m
                    amr_fine%q_cons(1)%sf(fi, fj, fk) = amr_fine%x_cc(fi)**2
                    if (n_glb > 0) amr_fine%q_cons(1)%sf(fi, fj, fk) = amr_fine%q_cons(1)%sf(fi, fj, fk) + amr_fine%y_cc(fj)**2
                    if (p_glb > 0) amr_fine%q_cons(1)%sf(fi, fj, fk) = amr_fine%q_cons(1)%sf(fi, fj, fk) + amr_fine%z_cc(fk)**2
                end do
            end do
        end do
        call s_restrict_one_var(amr_fine%q_cons(1), cscr(1))
        si_f = 0._wp; si_c = 0._wp
        do fk = 0, amr_fine%p
            do fj = 0, amr_fine%n
                do fi = 0, amr_fine%m
                    dvf = amr_fine%dx(fi)
                    if (n_glb > 0) dvf = dvf*amr_fine%dy(fj)
                    if (p_glb > 0) dvf = dvf*amr_fine%dz(fk)
                    si_f = si_f + dvf*real(amr_fine%q_cons(1)%sf(fi, fj, fk), wp)
                end do
            end do
        end do
        do ck = amr_patch_beg(3), merge(amr_patch_end(3), amr_patch_beg(3), p_glb > 0)
            do cj = amr_patch_beg(2), merge(amr_patch_end(2), amr_patch_beg(2), n_glb > 0)
                do ci = amr_patch_beg(1), amr_patch_end(1)
                    dvc = dx(ci)
                    if (n_glb > 0) dvc = dvc*dy(cj)
                    if (p_glb > 0) dvc = dvc*dz(ck)
                    si_c = si_c + dvc*real(cscr(1)%sf(ci, cj, ck), wp)
                end do
            end do
        end do
        errc = abs(si_f - si_c)/max(abs(si_f), 1.e-30_wp)
        if (proc_rank == 0) then
            print '(A,ES12.4)', ' [amr] prolong linear-reproduction err = ', errb
            print '(A,ES12.4)', ' [amr] restrict independent-integral err = ', errc
        end if
        deallocate (cscr(1)%sf); deallocate (cscr)

    end subroutine s_amr_operator_checks

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
            if (associated(amr_fine%q_cons_stor(i)%sf)) deallocate (amr_fine%q_cons_stor(i)%sf)
            if (associated(amr_fine%q_prim(i)%sf)) deallocate (amr_fine%q_prim(i)%sf)
            if (associated(amr_fine%rhs(i)%sf)) deallocate (amr_fine%rhs(i)%sf)
        end do
        deallocate (amr_fine%q_cons)
        deallocate (amr_fine%q_cons_stor)
        deallocate (amr_fine%q_prim)
        deallocate (amr_fine%rhs)
        if (allocated(amr_fine%x_cb)) deallocate (amr_fine%x_cb, amr_fine%x_cc, amr_fine%dx)
        if (allocated(amr_fine%y_cb)) deallocate (amr_fine%y_cb, amr_fine%y_cc, amr_fine%dy)
        if (allocated(amr_fine%z_cb)) deallocate (amr_fine%z_cb, amr_fine%z_cc, amr_fine%dz)

    end subroutine s_finalize_amr_module

end module m_amr
