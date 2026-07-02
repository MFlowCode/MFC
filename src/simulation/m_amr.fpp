!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief AMR hierarchy (SP1: one static, inert refined level-1 patch alongside the base solve).
module m_amr

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum
    use m_rhs, only: s_compute_rhs
    use m_amr_registers, only: s_amr_zero_fine_registers

    implicit none

    private
    public :: t_level, amr_fine, amr_maxc, amr_dt_fine, s_initialize_amr_module, s_populate_amr_fine, &
        & s_interpolate_coarse_to_fine, s_restrict_fine_to_coarse, s_amr_conservation_check, s_finalize_amr_module, &
        & s_amr_swap_to_fine, s_amr_restore_coarse, s_amr_fill_fine_ghosts, s_amr_operator_checks, s_advance_amr_fine_stage, &
        & s_advance_amr_fine_substeps, s_amr_conservation_defect, s_set_amr_fine_geometry, s_amr_regrid, amr_owner_win_lo, &
        & amr_owner_win_hi

    !> Fine-level time step for subcycling (= 0.5*dt after init; 0 when amr is off).
    real(wp) :: amr_dt_fine = 0._wp

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
        type(scalar_field), allocatable :: q_ghost_a(:)    !< subcycle ghost-lerp source at coarse t^n (ghost shell only)
        type(scalar_field), allocatable :: q_ghost_b(:)    !< subcycle ghost-lerp source at coarse t^{n+1} (ghost shell only)
    end type t_level

    type(t_level) :: amr_fine
    integer       :: amr_maxc(3)  !< max coarse patch cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

    !> Owner containment window in GLOBAL coarse indices (fixed for the run, known on all ranks; 0 in collapsed dims): a patch is
    !! owner-contained iff amr_owner_win_lo(d) <= beg(d) and end(d) <= amr_owner_win_hi(d) per active dim.
    integer :: amr_owner_win_lo(3) = 0, amr_owner_win_hi(3) = 0
    logical :: amr_win_clamp_warned = .false.  !< one-shot regrid warning (SP7a mobility limit)

    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)

    !> Conservation-defect baselines (level-0 interior integrals at init)
    real(wp) :: amr_mass0 = 0._wp, amr_energy0 = 0._wp

contains

    !> Build the static refined level-1 patch. No-op unless amr. Called after level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Preallocates all fine arrays at max size so regrid only needs to call s_set_amr_fine_geometry.
    impure subroutine s_initialize_amr_module()

        integer :: i, d, max_f1, max_f2, max_f3
        integer :: mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi
        integer :: sidx(3), ext(3), owner_count, bad_loc, bad_glb

        if (.not. amr) return

        amr_dt_fine = 0.5_wp*dt

        ! Containment (owner model): the patch + a buff_size shell must lie within THIS rank's subdomain.
        ! (np=1: start_idx=0, m=m_glb - degenerates to the old domain-margin rule.)
        ! buff_size is not available at checker time, so this check must be here.
        sidx = 0
        sidx(1) = start_idx(1)
        if (n_glb > 0) sidx(2) = start_idx(2)
        if (p_glb > 0) sidx(3) = start_idx(3)
        amr_rank_owns_patch = amr_patch_beg(1) >= sidx(1) + buff_size .and. amr_patch_end(1) <= sidx(1) + m - buff_size
        if (n_glb > 0) amr_rank_owns_patch = amr_rank_owns_patch .and. amr_patch_beg(2) >= sidx(2) + buff_size &
            & .and. amr_patch_end(2) <= sidx(2) + n - buff_size
        if (p_glb > 0) amr_rank_owns_patch = amr_rank_owns_patch .and. amr_patch_beg(3) >= sidx(3) + buff_size &
            & .and. amr_patch_end(3) <= sidx(3) + p - buff_size
        call s_mpi_allreduce_integer_max(merge(1, 0, amr_rank_owns_patch), owner_count)
        if (owner_count == 0) then
            call s_mpi_abort('amr patch must lie within a single rank subdomain (>= buff_size from its edges); ' &
                             & // 'move the patch, use fewer ranks, or await SP7b')
        end if

        ! The fine advance reuses the solver scratch (m_rhs/WENO/Riemann work arrays), which is sized to the
        ! owner's LOCAL grid, so every fine extent must fit the owner's local extent. (np=1: the checker's
        ! 2*patch-1 <= m_glb bound makes this a no-op.)
        bad_loc = 0
        if (amr_rank_owns_patch) then
            if (2*(amr_patch_end(1) - amr_patch_beg(1) + 1) - 1 > m) bad_loc = 1
            if (n_glb > 0 .and. 2*(amr_patch_end(2) - amr_patch_beg(2) + 1) - 1 > n) bad_loc = 1
            if (p_glb > 0 .and. 2*(amr_patch_end(3) - amr_patch_beg(3) + 1) - 1 > p) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        if (bad_glb == 1) then
            call s_mpi_abort('amr fine extent exceeds the owner rank local grid (solver scratch is local-sized); ' &
                             & // 'shrink the patch or use fewer ranks')
        end if

        ! Owner containment window in global indices, distributed to all ranks (fixed for the run; 0 in collapsed dims).
        ext(1) = m; ext(2) = n; ext(3) = p
        amr_owner_win_lo = 0; amr_owner_win_hi = 0
        do d = 1, num_dims
            call s_mpi_allreduce_integer_max(merge(sidx(d) + buff_size, -huge(1), amr_rank_owns_patch), amr_owner_win_lo(d))
            call s_mpi_allreduce_integer_min(merge(sidx(d) + ext(d) - buff_size, huge(1), amr_rank_owns_patch), amr_owner_win_hi(d))
        end do

        amr_fine%ref_ratio = 2
        amr_fine%buff_size = buff_size

        ! max coarse patch cells per dim (upper bound for any future regrid box); 1 for collapsed dims
        amr_maxc(1) = (m_glb + 1)/2
        amr_maxc(2) = 1; amr_maxc(3) = 1
        if (n_glb > 0) amr_maxc(2) = (n_glb + 1)/2
        if (p_glb > 0) amr_maxc(3) = (p_glb + 1)/2

        ! max fine extents and buffered bounds for preallocation
        max_f1 = 2*amr_maxc(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = 2*amr_maxc(2) - 1
        if (p_glb > 0) max_f3 = 2*amr_maxc(3) - 1
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if

        ! owner-only storage: non-owner ranks run pure coarse and never touch the fine level
        if (amr_rank_owns_patch) then
            ! preallocate coordinates at max fine extents
            allocate (amr_fine%x_cb(-1:max_f1), amr_fine%x_cc(0:max_f1), amr_fine%dx(0:max_f1))
            if (n_glb > 0) allocate (amr_fine%y_cb(-1:max_f2), amr_fine%y_cc(0:max_f2), amr_fine%dy(0:max_f2))
            if (p_glb > 0) allocate (amr_fine%z_cb(-1:max_f3), amr_fine%z_cc(0:max_f3), amr_fine%dz(0:max_f3))

            ! preallocate fine fields at max buffered extents; loops always use current amr_fine%m/n/p
            allocate (amr_fine%q_cons(1:sys_size))
            allocate (amr_fine%q_cons_stor(1:sys_size))
            allocate (amr_fine%q_prim(1:sys_size))
            allocate (amr_fine%rhs(1:sys_size))
            allocate (amr_fine%q_ghost_a(1:sys_size))
            allocate (amr_fine%q_ghost_b(1:sys_size))
            do i = 1, sys_size
                allocate (amr_fine%q_cons(i)%sf(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi))
                allocate (amr_fine%q_cons_stor(i)%sf(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi))
                allocate (amr_fine%q_prim(i)%sf(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi))
                allocate (amr_fine%rhs(i)%sf(0:max_f1,0:max_f2,0:max_f3))
                allocate (amr_fine%q_ghost_a(i)%sf(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi))
                allocate (amr_fine%q_ghost_b(i)%sf(mbuf1_lo:mbuf1_hi,mbuf2_lo:mbuf2_hi,mbuf3_lo:mbuf3_hi))
            end do
        end if

        ! set geometry (region, m/n/p, idwbuff, coordinates) for the initial patch
        call s_set_amr_fine_geometry(amr_patch_beg, amr_patch_end)

    end subroutine s_initialize_amr_module

    !> Fill level-1 fcb/fcc/fdx by bisecting parent cells; pcb_lb is lbound(parent_cb, 1). Passing pcb as assumed-shape resets
    !! lbound to 1; pcb_lb + idx_offset recovers original indexing. Arrays must be preallocated at max size; only 0..nfine filled.
    subroutine s_build_level_coords(pcb, pcb_lb, lo, nfine, fcb, fcc, fdx)

        real(wp), intent(in)                 :: pcb(:)
        integer, intent(in)                  :: pcb_lb, lo, nfine
        real(wp), allocatable, intent(inout) :: fcb(:), fcc(:), fdx(:)
        integer                              :: fi, c, idx_offset
        real(wp)                             :: xl, xr, xm
        ! pcb(k) = parent_cb(k + pcb_lb - 1); to access parent_cb(j): k = j - pcb_lb + 1

        idx_offset = 1 - pcb_lb
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

    !> Set the fine level's geometry (region, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated at max size;
    !! this only updates metadata and refills coords.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)

        amr_fine%region%lo = lo; amr_fine%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        amr_fine%m = 2*(hi(1) - lo(1) + 1) - 1
        amr_fine%n = 0; amr_fine%p = 0
        if (n_glb > 0) amr_fine%n = 2*(hi(2) - lo(2) + 1) - 1
        if (p_glb > 0) amr_fine%p = 2*(hi(3) - lo(3) + 1) - 1
        amr_fine%idwbuff(1)%beg = -buff_size; amr_fine%idwbuff(1)%end = amr_fine%m + buff_size
        amr_fine%idwbuff(2)%beg = 0; amr_fine%idwbuff(2)%end = 0
        amr_fine%idwbuff(3)%beg = 0; amr_fine%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_fine%idwbuff(2)%beg = -buff_size; amr_fine%idwbuff(2)%end = amr_fine%n + buff_size
        end if
        if (p_glb > 0) then
            amr_fine%idwbuff(3)%beg = -buff_size; amr_fine%idwbuff(3)%end = amr_fine%p + buff_size
        end if
        ! coord building is owner-only (non-owner coord arrays are unallocated); the parent origin is
        ! converted to LOCAL indexing so the bisection reads the owner's local x_cb slice
        if (amr_rank_owns_patch) then
            call s_build_level_coords(x_cb, lbound(x_cb, 1), lo(1) - start_idx(1), amr_fine%m, amr_fine%x_cb, amr_fine%x_cc, &
                                      & amr_fine%dx)
            if (n_glb > 0) call s_build_level_coords(y_cb, lbound(y_cb, 1), lo(2) - start_idx(2), amr_fine%n, amr_fine%y_cb, &
                & amr_fine%y_cc, amr_fine%dy)
            if (p_glb > 0) call s_build_level_coords(z_cb, lbound(z_cb, 1), lo(3) - start_idx(3), amr_fine%p, amr_fine%z_cb, &
                & amr_fine%z_cc, amr_fine%dz)
        end if

    end subroutine s_set_amr_fine_geometry

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, qf)

        type(scalar_field), intent(in)    :: qc
        type(scalar_field), intent(inout) :: qf
        integer                           :: fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                          :: u0, sx, sy, sz, xix, xiy, xiz

        ! region indices are GLOBAL; the coarse source qc is rank-LOCAL (identical at np=1: start_idx=0)

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do fk = 0, amr_fine%p
            ck = amr_fine%region%lo(3) + fk/amr_fine%ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_fine%n
                cj = amr_fine%region%lo(2) + fj/amr_fine%ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_fine%m
                    ci = amr_fine%region%lo(1) + fi/amr_fine%ref_ratio - ox
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
        if (.not. amr_rank_owns_patch) return
        call s_interpolate_coarse_to_fine(q_cons_base)

    end subroutine s_populate_amr_fine

    !> Volume-weighted restriction for a single variable pair. Reads from qf (fine, must include interior 0:amr_fine%m etc.); writes
    !! to qc (coarse, over the patch).
    impure subroutine s_restrict_one_var(qf, qc)

        type(scalar_field), intent(in)    :: qf
        type(scalar_field), intent(inout) :: qc
        integer                           :: ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild, ox, oy, oz
        real(wp)                          :: acc

        ! region loop indices are GLOBAL; the coarse target qc is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        nchild = amr_fine%ref_ratio
        if (n_glb > 0) nchild = nchild*amr_fine%ref_ratio
        if (p_glb > 0) nchild = nchild*amr_fine%ref_ratio
        do ck = amr_fine%region%lo(3), merge(amr_fine%region%hi(3), amr_fine%region%lo(3), p_glb > 0)
            fk0 = (ck - amr_fine%region%lo(3))*amr_fine%ref_ratio
            do cj = amr_fine%region%lo(2), merge(amr_fine%region%hi(2), amr_fine%region%lo(2), n_glb > 0)
                fj0 = (cj - amr_fine%region%lo(2))*amr_fine%ref_ratio
                do ci = amr_fine%region%lo(1), amr_fine%region%hi(1)
                    fi0 = (ci - amr_fine%region%lo(1))*amr_fine%ref_ratio
                    acc = 0._wp
                    do ddk = 0, merge(amr_fine%ref_ratio - 1, 0, p_glb > 0)
                        do ddj = 0, merge(amr_fine%ref_ratio - 1, 0, n_glb > 0)
                            acc = acc + real(qf%sf(fi0, fj0 + ddj, fk0 + ddk), wp) + real(qf%sf(fi0 + 1, fj0 + ddj, fk0 + ddk), wp)
                        end do
                    end do
                    qc%sf(ci - ox, cj - oy, ck - oz) = acc/real(nchild, wp)
                end do
            end do
        end do

    end subroutine s_restrict_one_var

    !> Volume-weighted restriction: each covered coarse cell = average of its ref_ratio^d fine children. Writes ONLY the caller's
    !! coarse target (SP2: a scratch buffer) - never level-0.
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt
        integer                                                :: i

        if (.not. amr_rank_owns_patch) return
        do i = 1, sys_size
            call s_restrict_one_var(amr_fine%q_cons(i), coarse_tgt(i))
        end do

    end subroutine s_restrict_fine_to_coarse

    !> SP2 gate: restrict(prolong(coarse)) must reproduce coarse over the patch interior (conservation). Init-only diagnostic;
    !! allocates a scratch coarse target, never touches level-0 or the solve.
    impure subroutine s_amr_conservation_check(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        type(scalar_field), dimension(:), allocatable       :: scratch
        integer                                             :: i, ci, cj, ck, ox, oy, oz
        real(wp)                                            :: err, e

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        allocate (scratch(1:sys_size))
        do i = 1, sys_size
            allocate (scratch(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do
        call s_restrict_fine_to_coarse(scratch)
        err = 0._wp
        do i = 1, sys_size
            do ck = amr_fine%region%lo(3), merge(amr_fine%region%hi(3), amr_fine%region%lo(3), p_glb > 0)
                do cj = amr_fine%region%lo(2), merge(amr_fine%region%hi(2), amr_fine%region%lo(2), n_glb > 0)
                    do ci = amr_fine%region%lo(1), amr_fine%region%hi(1)
                        e = abs(real(scratch(i)%sf(ci - ox, cj - oy, ck - oz), wp) - real(q_cons_base(i)%sf(ci - ox, cj - oy, &
                                & ck - oz), wp))
                        if (e > err) err = e
                    end do
                end do
            end do
        end do
        print '(A,ES12.4)', ' [amr] restrict-prolong conservation err = ', err  ! owner rank only
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
        integer                                                :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                                               :: u0, sx, sy, sz, xix, xiy, xiz

        ! region indices are GLOBAL; the coarse source q_coarse is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do i = 1, sys_size
            do fk = amr_fine%idwbuff(3)%beg, amr_fine%idwbuff(3)%end
                ck = 0; xiz = 0._wp
                if (p_glb > 0) then
                    ck = amr_fine%region%lo(3) + floor(real(fk, wp)/real(amr_fine%ref_ratio, wp)) - oz
                    xiz = (real(modulo(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                end if
                do fj = amr_fine%idwbuff(2)%beg, amr_fine%idwbuff(2)%end
                    cj = 0; xiy = 0._wp
                    if (n_glb > 0) then
                        cj = amr_fine%region%lo(2) + floor(real(fj, wp)/real(amr_fine%ref_ratio, wp)) - oy
                        xiy = (real(modulo(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    end if
                    do fi = amr_fine%idwbuff(1)%beg, amr_fine%idwbuff(1)%end
                        ! skip the interior: only the ghost shell is filled
                        if (fi >= 0 .and. fi <= amr_fine%m .and. fj >= 0 .and. fj <= amr_fine%n .and. fk >= 0 &
                            & .and. fk <= amr_fine%p) cycle
                        ci = amr_fine%region%lo(1) + floor(real(fi, wp)/real(amr_fine%ref_ratio, wp)) - ox
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

    !> Advance the fine level through RK stage s (same dt as level-0, no subcycling). Called BETWEEN the coarse RHS and the coarse
    !! RK update, so q_cons_coarse is the coarse STAGE-ENTRY state. Fine prim ghosts are obtained by widening idwint to the fine
    !! buffer so the cons->prim conversion inside s_compute_rhs covers the (prolonged) cons ghost shell; the BC populate is skipped
    !! via amr_in_fine_advance. The update mirrors the coarse non-IGR rk_coef form in s_tvd_rk.
    impure subroutine s_advance_amr_fine_stage(s, coefs, q_cons_coarse, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
        & time_avg)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(scalar_field), dimension(sys_size), intent(in)        :: q_cons_coarse
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        real(wp), intent(inout)                                    :: time_avg
        integer                                                    :: i, fi, fj, fk

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return

        ! ghost prolongation from the coarse stage-entry conservative state
        call s_amr_fill_fine_ghosts(q_cons_coarse, amr_fine%q_cons)

        ! step-entry backup for the SSP-RK combination
        if (s == 1) then
            do i = 1, sys_size
                amr_fine%q_cons_stor(i)%sf(:,:,:) = amr_fine%q_cons(i)%sf(:,:,:)
            end do
        end if

        amr_in_fine_advance = .true.
        call s_amr_swap_to_fine()
        idwint = amr_fine%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        call s_compute_rhs(amr_fine%q_cons, q_T_sf, amr_fine%q_prim, bc_type, amr_fine%rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                           & time_avg, s)
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update (mirror of the coarse non-IGR form; compute in wp, store stp)
        do i = 1, sys_size
            do fk = 0, amr_fine%p
                do fj = 0, amr_fine%n
                    do fi = 0, amr_fine%m
                        amr_fine%q_cons(i)%sf(fi, fj, fk) = (coefs(1)*real(amr_fine%q_cons(i)%sf(fi, fj, fk), &
                                        & wp) + coefs(2)*real(amr_fine%q_cons_stor(i)%sf(fi, fj, fk), &
                                        & wp) + coefs(3)*dt*real(amr_fine%rhs(i)%sf(fi, fj, fk), wp))/coefs(4)
                    end do
                end do
            end do
        end do

    end subroutine s_advance_amr_fine_stage

    !> Subcycled fine advance (amr_subcycle): two dt/2 SSP-RK3 substeps AFTER the coarse step. q_old/q_new are the coarse t^n and
    !! t^{n+1} states; each stage's ghosts are the linear time interpolation at the stage time theta = (substep-1 + c_s)/2 with
    !! SSP-RK3 abscissae c = [0, 1, 1/2]. Fine flux registers are zeroed here and accumulate over all six stages (0.5*rk3_w each) so
    !! the end-of-step state reflux sees the time-averaged effective fine flux.
    impure subroutine s_advance_amr_fine_substeps(q_old, q_new, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
        & time_avg)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_old, q_new
        real(wp), dimension(:,:), intent(in)                       :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp), intent(inout)                                    :: time_avg
        real(wp), parameter                                        :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer                                                    :: sub, s, i, fi, fj, fk
        real(wp)                                                   :: th

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return

        ! fill both lerp sources once: ghost shells prolonged from coarse t^n and t^{n+1}
        call s_amr_fill_fine_ghosts(q_old, amr_fine%q_ghost_a)
        call s_amr_fill_fine_ghosts(q_new, amr_fine%q_ghost_b)
        call s_amr_zero_fine_registers()

        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp

                ! lerp the ghost shell into q_cons at the stage time (interior untouched)
                do i = 1, sys_size
                    do fk = amr_fine%idwbuff(3)%beg, amr_fine%idwbuff(3)%end
                        do fj = amr_fine%idwbuff(2)%beg, amr_fine%idwbuff(2)%end
                            do fi = amr_fine%idwbuff(1)%beg, amr_fine%idwbuff(1)%end
                                if (fi >= 0 .and. fi <= amr_fine%m .and. fj >= 0 .and. fj <= amr_fine%n .and. fk >= 0 &
                                    & .and. fk <= amr_fine%p) cycle
                                amr_fine%q_cons(i)%sf(fi, fj, fk) = (1._wp - th)*real(amr_fine%q_ghost_a(i)%sf(fi, fj, fk), &
                                                & wp) + th*real(amr_fine%q_ghost_b(i)%sf(fi, fj, fk), wp)
                            end do
                        end do
                    end do
                end do

                ! substep-entry backup for the SSP-RK combination
                if (s == 1) then
                    do i = 1, sys_size
                        amr_fine%q_cons_stor(i)%sf(0:amr_fine%m,0:amr_fine%n,0:amr_fine%p) = amr_fine%q_cons(i)%sf(0:amr_fine%m, &
                                             & 0:amr_fine%n,0:amr_fine%p)
                    end do
                end if

                amr_in_fine_advance = .true.
                call s_amr_swap_to_fine()
                idwint = amr_fine%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
                call s_compute_rhs(amr_fine%q_cons, q_T_sf, amr_fine%q_prim, bc_type, amr_fine%rhs, pb_in, rhs_pb, mv_in, rhs_mv, &
                                   & t_step, time_avg, s)
                call s_amr_restore_coarse()
                amr_in_fine_advance = .false.

                ! RK stage update at the FINE time step (compute in wp, store stp)
                do i = 1, sys_size
                    do fk = 0, amr_fine%p
                        do fj = 0, amr_fine%n
                            do fi = 0, amr_fine%m
                                amr_fine%q_cons(i)%sf(fi, fj, fk) = (coefs(s, 1)*real(amr_fine%q_cons(i)%sf(fi, fj, fk), &
                                                & wp) + coefs(s, 2)*real(amr_fine%q_cons_stor(i)%sf(fi, fj, fk), wp) + coefs(s, &
                                                & 3)*amr_dt_fine*real(amr_fine%rhs(i)%sf(fi, fj, fk), wp))/coefs(s, 4)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_advance_amr_fine_substeps

    !> Regrid: tag by relative density gradient, form the padded/clamped bounding box, rebuild the fine level (copy old-fine on
    !! overlap via q_cons_stor bounce; prolong new cells). Called between steps only. No-op if nothing is tagged or the box is
    !! unchanged.
    impure subroutine s_amr_regrid(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        integer                                             :: lo(3), hi(3), old_lo(3), old_m1, old_n1, old_p1
        integer                                             :: ci, cj, ck, fi, fj, fk, ofi, ofj, ofk, sh(3), i, d
        integer                                             :: sidx(3), lo_r, hi_r, maxc_eff
        real(wp)                                            :: r0, g
        logical                                             :: tagged, win_cut

        ! 1) tag + bounding box over the LOCAL interior (sweep away from local edges; ghosts not needed), emitting
        !    GLOBAL indices. No rank-dependent branching before the reductions: every rank reaches both allreduces.

        sidx = 0
        sidx(1) = start_idx(1)
        if (n_glb > 0) sidx(2) = start_idx(2)
        if (p_glb > 0) sidx(3) = start_idx(3)
        lo = huge(1); hi = -huge(1); tagged = .false.
        do ck = merge(1, 0, p_glb > 0), merge(p - 1, 0, p_glb > 0)
            do cj = merge(1, 0, n_glb > 0), merge(n - 1, 0, n_glb > 0)
                do ci = 1, m - 1
                    r0 = max(abs(real(q_cons_base(1)%sf(ci, cj, ck), wp)), 1.e-30_wp)
                    g = abs(real(q_cons_base(1)%sf(ci + 1, cj, ck), wp) - real(q_cons_base(1)%sf(ci - 1, cj, ck), wp))
                    if (n_glb > 0) g = max(g, abs(real(q_cons_base(1)%sf(ci, cj + 1, ck), wp) - real(q_cons_base(1)%sf(ci, &
                        & cj - 1, ck), wp)))
                    if (p_glb > 0) g = max(g, abs(real(q_cons_base(1)%sf(ci, cj, ck + 1), wp) - real(q_cons_base(1)%sf(ci, cj, &
                        & ck - 1), wp)))
                    if (g/(2._wp*r0) > amr_tag_eps) then
                        tagged = .true.
                        lo(1) = min(lo(1), ci + sidx(1)); hi(1) = max(hi(1), ci + sidx(1))
                        lo(2) = min(lo(2), cj + sidx(2)); hi(2) = max(hi(2), cj + sidx(2))
                        lo(3) = min(lo(3), ck + sidx(3)); hi(3) = max(hi(3), ck + sidx(3))
                    end if
                end do
            end do
        end do
        ! reduce the box per active dim; the result is identical on all ranks
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min(merge(lo(d), huge(1), tagged), lo_r)
            call s_mpi_allreduce_integer_max(merge(hi(d), -huge(1), tagged), hi_r)
            lo(d) = lo_r; hi(d) = hi_r
        end do
        if (hi(1) < lo(1)) return  ! nothing tagged on any rank (all ranks agree)

        ! 2) pad + clamp: domain margin, owner containment window (the SP7a mobility limit: the patch cannot leave its
        !    owner rank; at np=1 the window equals the domain margin), and max extent (also bounded so the fine grid
        !    fits the owner's LOCAL solver scratch, whose extent = window + its two buff_size margins); collapsed dims
        !    pinned to 0. All inputs are globally identical, so every rank computes the same box.
        win_cut = .false.
        lo(1) = max(lo(1) - amr_buf, buff_size); hi(1) = min(hi(1) + amr_buf, m_glb - buff_size)
        if (lo(1) < amr_owner_win_lo(1) .or. hi(1) > amr_owner_win_hi(1)) win_cut = .true.
        lo(1) = max(lo(1), amr_owner_win_lo(1)); hi(1) = min(hi(1), amr_owner_win_hi(1))
        maxc_eff = min(amr_maxc(1), (amr_owner_win_hi(1) - amr_owner_win_lo(1) + 2*buff_size + 1)/2)
        if (hi(1) - lo(1) + 1 > maxc_eff) then
            if (amr_rank_owns_patch) print '(A)', ' [amr] WARNING: tagged region exceeds max patch; clamping'
            hi(1) = lo(1) + maxc_eff - 1
        end if
        if (n_glb > 0) then
            lo(2) = max(lo(2) - amr_buf, buff_size); hi(2) = min(hi(2) + amr_buf, n_glb - buff_size)
            if (lo(2) < amr_owner_win_lo(2) .or. hi(2) > amr_owner_win_hi(2)) win_cut = .true.
            lo(2) = max(lo(2), amr_owner_win_lo(2)); hi(2) = min(hi(2), amr_owner_win_hi(2))
            maxc_eff = min(amr_maxc(2), (amr_owner_win_hi(2) - amr_owner_win_lo(2) + 2*buff_size + 1)/2)
            if (hi(2) - lo(2) + 1 > maxc_eff) hi(2) = lo(2) + maxc_eff - 1
        else
            lo(2) = 0; hi(2) = 0
        end if
        if (p_glb > 0) then
            lo(3) = max(lo(3) - amr_buf, buff_size); hi(3) = min(hi(3) + amr_buf, p_glb - buff_size)
            if (lo(3) < amr_owner_win_lo(3) .or. hi(3) > amr_owner_win_hi(3)) win_cut = .true.
            lo(3) = max(lo(3), amr_owner_win_lo(3)); hi(3) = min(hi(3), amr_owner_win_hi(3))
            maxc_eff = min(amr_maxc(3), (amr_owner_win_hi(3) - amr_owner_win_lo(3) + 2*buff_size + 1)/2)
            if (hi(3) - lo(3) + 1 > maxc_eff) hi(3) = lo(3) + maxc_eff - 1
        else
            lo(3) = 0; hi(3) = 0
        end if
        if (win_cut .and. amr_rank_owns_patch .and. .not. amr_win_clamp_warned) then
            print '(A)', ' [amr] WARNING: regrid box clamped to the owner rank window (SP7a: the patch cannot migrate ranks)'
            amr_win_clamp_warned = .true.
        end if
        if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return  ! tag box fully outside the window; keep the old region
        if (all(lo == amr_fine%region%lo) .and. all(hi == amr_fine%region%hi)) return

        ! 3) owner-only: stash old fine interior in the (dead-between-steps) RK bounce buffer; non-owner field arrays
        !    are unallocated (metadata below is still updated on ALL ranks)
        old_lo = amr_fine%region%lo
        old_m1 = amr_fine%m; old_n1 = amr_fine%n; old_p1 = amr_fine%p
        if (amr_rank_owns_patch) then
            do i = 1, sys_size
                amr_fine%q_cons_stor(i)%sf(0:old_m1,0:old_n1,0:old_p1) = amr_fine%q_cons(i)%sf(0:old_m1,0:old_n1,0:old_p1)
            end do
        end if

        ! 4) rebuild geometry (region/mirrors/extents on all ranks; coord fill owner-guarded inside), then owner-only:
        !    prolong the whole new patch and overwrite the overlap with old fine data
        call s_set_amr_fine_geometry(lo, hi)
        if (.not. amr_rank_owns_patch) return
        call s_interpolate_coarse_to_fine(q_cons_base)
        sh = 2*(lo - old_lo)  ! old fine index = new fine index + sh (per dim; collapsed dims sh=0)
        do i = 1, sys_size
            do fk = 0, amr_fine%p
                ofk = fk + sh(3)
                if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_p1)) cycle
                do fj = 0, amr_fine%n
                    ofj = fj + sh(2)
                    if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_n1)) cycle
                    do fi = 0, amr_fine%m
                        ofi = fi + sh(1)
                        if (ofi < 0 .or. ofi > old_m1) cycle
                        amr_fine%q_cons(i)%sf(fi, fj, fk) = amr_fine%q_cons_stor(i)%sf(ofi, ofj, ofk)
                    end do
                end do
            end do
        end do
        ! owner-only code path (the unique owner prints)
        print '(A,I0,A,I0,A,I0,A)', ' [amr] regrid: box x ', lo(1), ':', hi(1), ' (', (hi(1) - lo(1) + 1), ' coarse cells)'

    end subroutine s_amr_regrid

    !> Global Sum(dV*U) for continuity (var 1) and energy (eqn_idx%E) over the level-0 interior. First call (finalize_report=F)
    !! stores the baseline; the finalize call prints the relative drift (expected small nonzero until SP4 refluxing).
    impure subroutine s_amr_conservation_defect(q_cons_base, finalize_report)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        logical, intent(in)                                 :: finalize_report
        real(wp)                                            :: sm, se, dv, sm_glb, se_glb
        integer                                             :: ci, cj, ck

        if (.not. amr) return
        sm = 0._wp; se = 0._wp
        do ck = 0, p
            do cj = 0, n
                do ci = 0, m
                    dv = dx(ci)
                    if (n_glb > 0) dv = dv*dy(cj)
                    if (p_glb > 0) dv = dv*dz(ck)
                    sm = sm + dv*real(q_cons_base(1)%sf(ci, cj, ck), wp)
                    se = se + dv*real(q_cons_base(eqn_idx%E)%sf(ci, cj, ck), wp)
                end do
            end do
        end do
        if (num_procs > 1) then
            call s_mpi_allreduce_sum(sm, sm_glb); sm = sm_glb
            call s_mpi_allreduce_sum(se, se_glb); se = se_glb
        end if
        if (.not. finalize_report) then
            amr_mass0 = sm; amr_energy0 = se
        else if (proc_rank == 0) then
            print '(A,ES12.4,A,ES12.4)', ' [amr] conservation defect: mass drift = ', abs(sm - amr_mass0)/max(abs(amr_mass0), &
                & 1.e-30_wp), '  energy drift = ', abs(se - amr_energy0)/max(abs(amr_energy0), 1.e-30_wp)
        end if

    end subroutine s_amr_conservation_defect

    !> Init-time operator verification: (b) linear reproduction, (c) restriction of an independent field. Uses amr_fine%q_cons(1) as
    !! scratch; called before s_populate_amr_fine overwrites it.
    impure subroutine s_amr_operator_checks()

        type(scalar_field), allocatable :: cscr(:)
        integer                         :: fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                        :: e, errb, errc, si_f, si_c, dvf, dvc, want

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)

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
        do ck = amr_fine%region%lo(3), merge(amr_fine%region%hi(3), amr_fine%region%lo(3), p_glb > 0)
            do cj = amr_fine%region%lo(2), merge(amr_fine%region%hi(2), amr_fine%region%lo(2), n_glb > 0)
                do ci = amr_fine%region%lo(1), amr_fine%region%hi(1)
                    dvc = dx(ci - ox)
                    if (n_glb > 0) dvc = dvc*dy(cj - oy)
                    if (p_glb > 0) dvc = dvc*dz(ck - oz)
                    si_c = si_c + dvc*real(cscr(1)%sf(ci - ox, cj - oy, ck - oz), wp)
                end do
            end do
        end do
        errc = abs(si_f - si_c)/max(abs(si_f), 1.e-30_wp)
        ! owner rank only
        print '(A,ES12.4)', ' [amr] prolong linear-reproduction err = ', errb
        print '(A,ES12.4)', ' [amr] restrict independent-integral err = ', errc
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
        if (.not. amr_rank_owns_patch) return  ! non-owner ranks never allocated the fine level
        do i = 1, sys_size
            if (associated(amr_fine%q_cons(i)%sf)) deallocate (amr_fine%q_cons(i)%sf)
            if (associated(amr_fine%q_cons_stor(i)%sf)) deallocate (amr_fine%q_cons_stor(i)%sf)
            if (associated(amr_fine%q_prim(i)%sf)) deallocate (amr_fine%q_prim(i)%sf)
            if (associated(amr_fine%rhs(i)%sf)) deallocate (amr_fine%rhs(i)%sf)
            if (associated(amr_fine%q_ghost_a(i)%sf)) deallocate (amr_fine%q_ghost_a(i)%sf)
            if (associated(amr_fine%q_ghost_b(i)%sf)) deallocate (amr_fine%q_ghost_b(i)%sf)
        end do
        deallocate (amr_fine%q_cons)
        deallocate (amr_fine%q_cons_stor)
        deallocate (amr_fine%q_prim)
        deallocate (amr_fine%rhs)
        deallocate (amr_fine%q_ghost_a)
        deallocate (amr_fine%q_ghost_b)
        if (allocated(amr_fine%x_cb)) deallocate (amr_fine%x_cb, amr_fine%x_cc, amr_fine%dx)
        if (allocated(amr_fine%y_cb)) deallocate (amr_fine%y_cb, amr_fine%y_cc, amr_fine%dy)
        if (allocated(amr_fine%z_cb)) deallocate (amr_fine%z_cb, amr_fine%z_cc, amr_fine%dz)

    end subroutine s_finalize_amr_module

end module m_amr
