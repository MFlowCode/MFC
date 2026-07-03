!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief AMR hierarchy (SP1: one static, inert refined level-1 patch alongside the base solve).
module m_amr

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort, s_initialize_amr_mpi_buffers, s_mpi_sendrecv_amr_fine_halo
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum, &
        & s_mpi_sendrecv_variables_buffers
    use m_rhs, only: s_compute_rhs
    use m_amr_registers, only: s_amr_zero_fine_registers

    implicit none

    private
    public :: t_level, amr_fine, amr_maxc, amr_dt_fine, s_initialize_amr_module, s_populate_amr_fine, &
        & s_interpolate_coarse_to_fine, s_restrict_fine_to_coarse, s_amr_conservation_check, s_finalize_amr_module, &
        & s_amr_swap_to_fine, s_amr_restore_coarse, s_amr_fill_fine_ghosts, s_amr_operator_checks, s_advance_amr_fine_stage, &
        & s_advance_amr_fine_substeps, s_amr_conservation_defect, s_set_amr_fine_geometry, s_amr_regrid

    !> Fine-level time step for subcycling (= 0.5*dt after init; 0 when amr is off).
    real(wp) :: amr_dt_fine = 0._wp

    !> One refined level: its own grid + conservative fields. Field arrays are device-resident (@:ALLOCATE); coords/metadata
    !! host-only.
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

    !> Regrid box size cap per dim (fixed for the run, identical on all ranks; 1 in collapsed dims): a box of at most min over ranks
    !! of (local extent + 1)/2 cells intersects EVERY rank in at most (its extent + 1)/2 cells, so the per-rank scratch constraint
    !! 2*(isect cells) - 1 <= local extent holds by construction. Equals amr_maxc at np=1.
    integer :: amr_maxc_fit(3) = 1

    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)
    real(wp), allocatable :: sw_x_cb(:), sw_x_cc(:), sw_dx(:)
    real(wp), allocatable :: sw_y_cb(:), sw_y_cc(:), sw_dy(:)
    real(wp), allocatable :: sw_z_cb(:), sw_z_cc(:), sw_dz(:)

    !> Conservation-defect baselines (level-0 interior integrals at init)
    real(wp) :: amr_mass0 = 0._wp, amr_energy0 = 0._wp

    !> True (identically on all ranks) iff some rank's fine ghost-fill stencil reads its coarse GHOST cells - the solver populates
    !! only PRIM ghosts, so the CONS ghosts the fill prolongs from must be halo-exchanged first. Never true at np=1 (patch faces sit
    !! >= buff_size inside the domain).
    logical :: amr_xchg_coarse_ghosts = .false.

contains

    !> Build the static refined level-1 patch. No-op unless amr. Called after level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Preallocates all fine arrays at max size so regrid only needs to call s_set_amr_fine_geometry.
    impure subroutine s_initialize_amr_module()

        integer :: i, d, max_f1, max_f2, max_f3
        integer :: mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi
        integer :: sidx(3), ext(3), maxc_loc(3), bad_loc, bad_glb, fit_d

        if (.not. amr) return

        amr_dt_fine = 0.5_wp*dt

        ! Mirror decomposition: each rank holds the fine cells covering patch /\ its own subdomain
        ! (np=1: the intersection is the whole patch). buff_size is not available at checker time,
        ! so the geometric aborts below must live here.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        call s_amr_compute_isect(amr_patch_beg, amr_patch_end)

        ! the fine ghost shell and the reflux outside cells must stay inside the global domain (identical
        ! inputs on all ranks; every rank takes the same branch)
        if (amr_patch_beg(1) < buff_size .or. amr_patch_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_patch_beg(2) &
            & < buff_size .or. amr_patch_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_patch_beg(3) < buff_size &
            & .or. amr_patch_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr patch must lie at least buff_size cells inside the domain boundaries')
        end if

        ! Scratch constraint: the fine advance reuses the solver scratch (m_rhs/WENO/Riemann work arrays),
        ! which is sized to THIS rank's local grid, so each rank's fine extent must fit its local extent.
        ! (np=1: the checker's 2*patch-1 <= m_glb bound makes this a no-op.)
        bad_loc = 0
        if (amr_rank_owns_patch) then
            if (2*(amr_isect_hi(1) - amr_isect_lo(1) + 1) - 1 > m) bad_loc = 1
            if (n_glb > 0 .and. 2*(amr_isect_hi(2) - amr_isect_lo(2) + 1) - 1 > n) bad_loc = 1
            if (p_glb > 0 .and. 2*(amr_isect_hi(3) - amr_isect_lo(3) + 1) - 1 > p) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        if (bad_glb == 1) then
            call s_mpi_abort('amr fine extent exceeds a rank local grid (solver scratch is local-sized): the patch ' &
                             & // 'may cover at most about half of any rank subdomain per dimension; shrink the patch ' &
                             & // 'or use fewer ranks')
        end if

        amr_fine%ref_ratio = 2
        amr_fine%buff_size = buff_size

        ! max coarse patch cells per dim (upper bound for any future regrid box); 1 for collapsed dims
        amr_maxc(1) = (m_glb + 1)/2
        amr_maxc(2) = 1; amr_maxc(3) = 1
        if (n_glb > 0) amr_maxc(2) = (n_glb + 1)/2
        if (p_glb > 0) amr_maxc(3) = (p_glb + 1)/2

        ! regrid size cap: min over ranks of the local half-extent (see the declaration; = amr_maxc at np=1),
        ! so any clamped box satisfies every rank's scratch constraint and can move freely across ranks
        amr_maxc_fit = amr_maxc
        do d = 1, num_dims
            call s_mpi_allreduce_integer_min((ext(d) + 1)/2, fit_d)
            amr_maxc_fit(d) = min(amr_maxc(d), fit_d)
        end do

        ! preallocation cap for MY fine arrays: the largest intersection any patch box can have with this
        ! rank's subdomain (the scratch constraint caps it at about half the local extent; = amr_maxc at np=1)
        maxc_loc(1) = min(amr_maxc(1), (m + 1)/2)
        maxc_loc(2) = 1; maxc_loc(3) = 1
        if (n_glb > 0) maxc_loc(2) = min(amr_maxc(2), (n + 1)/2)
        if (p_glb > 0) maxc_loc(3) = min(amr_maxc(3), (p + 1)/2)

        ! max fine extents and buffered bounds for preallocation
        max_f1 = 2*maxc_loc(1) - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = 2*maxc_loc(2) - 1
        if (p_glb > 0) max_f3 = 2*maxc_loc(3) - 1

        ! MPI exchange buffers for the fine halo (all ranks; no-op without MFC_MPI)
        call s_initialize_amr_mpi_buffers(max_f1, max_f2, max_f3)
        mbuf1_lo = -buff_size; mbuf1_hi = max_f1 + buff_size
        mbuf2_lo = 0; mbuf2_hi = 0; mbuf3_lo = 0; mbuf3_hi = 0
        if (n_glb > 0) then; mbuf2_lo = -buff_size; mbuf2_hi = max_f2 + buff_size; end if
        if (p_glb > 0) then; mbuf3_lo = -buff_size; mbuf3_hi = max_f3 + buff_size; end if

        ! fine-level storage on ALL ranks (pool-sized to MY max intersection): regrid can move the patch onto
        ! any rank, so allocation state never changes with amr_rank_owns_patch flips; ranks without a current
        ! intersection have amr_fine%m/n/p = -1 and all fine loops no-op
        ! preallocate coordinates at max fine extents
        allocate (amr_fine%x_cb(-1:max_f1), amr_fine%x_cc(0:max_f1), amr_fine%dx(0:max_f1))
        if (n_glb > 0) allocate (amr_fine%y_cb(-1:max_f2), amr_fine%y_cc(0:max_f2), amr_fine%dy(0:max_f2))
        if (p_glb > 0) allocate (amr_fine%z_cb(-1:max_f3), amr_fine%z_cc(0:max_f3), amr_fine%dz(0:max_f3))

        ! bounce buffers for copy-based coord swap (GPU-safe; same bounds as the base-level global arrays)
        allocate (sw_x_cb(-1 - buff_size:m + buff_size))
        allocate (sw_x_cc(-buff_size:m + buff_size))
        allocate (sw_dx(-buff_size:m + buff_size))
        if (n_glb > 0) then
            allocate (sw_y_cb(-1 - buff_size:n + buff_size))
            allocate (sw_y_cc(-buff_size:n + buff_size))
            allocate (sw_dy(-buff_size:n + buff_size))
        end if
        if (p_glb > 0) then
            allocate (sw_z_cb(-1 - buff_size:p + buff_size))
            allocate (sw_z_cc(-buff_size:p + buff_size))
            allocate (sw_dz(-buff_size:p + buff_size))
        end if

        ! preallocate fine fields at max buffered extents; loops always use current amr_fine%m/n/p.
        ! Device-resident (@:ALLOCATE) so s_compute_rhs kernels can run on the fine patch;
        ! q_cons/q_prim/rhs get the Cray pointer setup mirroring m_time_steppers' field vectors.
        @:ALLOCATE(amr_fine%q_cons(1:sys_size))
        @:ALLOCATE(amr_fine%q_cons_stor(1:sys_size))
        @:ALLOCATE(amr_fine%q_prim(1:sys_size))
        @:ALLOCATE(amr_fine%rhs(1:sys_size))
        @:ALLOCATE(amr_fine%q_ghost_a(1:sys_size))
        @:ALLOCATE(amr_fine%q_ghost_b(1:sys_size))
        do i = 1, sys_size
            @:ALLOCATE(amr_fine%q_cons(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_fine%q_cons_stor(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_fine%q_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_fine%rhs(i)%sf(0:max_f1, 0:max_f2, 0:max_f3))
            @:ALLOCATE(amr_fine%q_ghost_a(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ALLOCATE(amr_fine%q_ghost_b(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
            @:ACC_SETUP_SFs(amr_fine%q_cons(i))
            @:ACC_SETUP_SFs(amr_fine%q_prim(i))
            @:ACC_SETUP_SFs(amr_fine%rhs(i))
        end do

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

    !> Compute this rank's per-dim intersection of the box lo:hi with its subdomain (GLOBAL indices, mirrored to amr_isect_lo/hi)
    !! and whether it holds fine cells (amr_rank_owns_patch: nonempty in all active dims). Must be called with the COARSE grid state
    !! in m/n/p (never from inside the fine advance).
    impure subroutine s_amr_compute_isect(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), d

        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        do d = 1, 3
            amr_isect_lo(d) = max(lo(d), sidx(d))
            amr_isect_hi(d) = min(hi(d), sidx(d) + ext(d))
        end do
        amr_rank_owns_patch = amr_isect_lo(1) <= amr_isect_hi(1)
        if (n_glb > 0) amr_rank_owns_patch = amr_rank_owns_patch .and. amr_isect_lo(2) <= amr_isect_hi(2)
        if (p_glb > 0) amr_rank_owns_patch = amr_rank_owns_patch .and. amr_isect_lo(3) <= amr_isect_hi(3)

    end subroutine s_amr_compute_isect

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated
    !! at max size; this only updates metadata and refills coords. Collective: ALL ranks must call together (init and regrid do) -
    !! it also refreshes the allreduced amr_xchg_coarse_ghosts flag for the new box.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), nmar, bad_loc, bad_glb

        call s_amr_compute_isect(lo, hi)
        amr_fine%region%lo = lo; amr_fine%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        ! LOCAL fine extents cover this rank's intersection (the whole patch at np=1); -1 in a dim whose
        ! intersection is empty, so 0..extent loops are no-ops on ranks without fine cells
        amr_fine%m = 2*max(amr_isect_hi(1) - amr_isect_lo(1) + 1, 0) - 1
        amr_fine%n = 0; amr_fine%p = 0
        if (n_glb > 0) amr_fine%n = 2*max(amr_isect_hi(2) - amr_isect_lo(2) + 1, 0) - 1
        if (p_glb > 0) amr_fine%p = 2*max(amr_isect_hi(3) - amr_isect_lo(3) + 1, 0) - 1
        amr_fine%idwbuff(1)%beg = -buff_size; amr_fine%idwbuff(1)%end = amr_fine%m + buff_size
        amr_fine%idwbuff(2)%beg = 0; amr_fine%idwbuff(2)%end = 0
        amr_fine%idwbuff(3)%beg = 0; amr_fine%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_fine%idwbuff(2)%beg = -buff_size; amr_fine%idwbuff(2)%end = amr_fine%n + buff_size
        end if
        if (p_glb > 0) then
            amr_fine%idwbuff(3)%beg = -buff_size; amr_fine%idwbuff(3)%end = amr_fine%p + buff_size
        end if
        ! coord building only on ranks with fine cells (others never read their coord arrays); the parent origin
        ! is this rank's INTERSECTION start, converted to LOCAL indexing so the bisection reads its x_cb slice
        if (amr_rank_owns_patch) then
            call s_build_level_coords(x_cb, lbound(x_cb, 1), amr_isect_lo(1) - start_idx(1), amr_fine%m, amr_fine%x_cb, &
                                      & amr_fine%x_cc, amr_fine%dx)
            if (n_glb > 0) call s_build_level_coords(y_cb, lbound(y_cb, 1), amr_isect_lo(2) - start_idx(2), amr_fine%n, &
                & amr_fine%y_cb, amr_fine%y_cc, amr_fine%dy)
            if (p_glb > 0) call s_build_level_coords(z_cb, lbound(z_cb, 1), amr_isect_lo(3) - start_idx(3), amr_fine%p, &
                & amr_fine%z_cb, amr_fine%z_cc, amr_fine%dz)
        end if

        ! Fine ghost prolongation reads up to nmar coarse cells past each face of the intersection; if that
        ! stencil leaves ANY rank's interior (patch near/at/across a rank boundary), the coarse CONS ghosts it
        ! reads must be halo-exchanged before every fill (the solver populates only PRIM ghosts). All ranks
        ! agree on the flag, so the pairwise exchanges are called consistently.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        nmar = (buff_size + 1)/2 + 1
        bad_loc = 0
        if (amr_rank_owns_patch) then
            if (amr_isect_lo(1) - sidx(1) < nmar .or. sidx(1) + ext(1) - amr_isect_hi(1) < nmar) bad_loc = 1
            if (n_glb > 0 .and. (amr_isect_lo(2) - sidx(2) < nmar .or. sidx(2) + ext(2) - amr_isect_hi(2) < nmar)) bad_loc = 1
            if (p_glb > 0 .and. (amr_isect_lo(3) - sidx(3) < nmar .or. sidx(3) + ext(3) - amr_isect_hi(3) < nmar)) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1

    end subroutine s_set_amr_fine_geometry

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, qf)

        type(scalar_field), intent(in)    :: qc
        type(scalar_field), intent(inout) :: qf
        integer                           :: fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                          :: u0, sx, sy, sz, xix, xiy, xiz

        ! fine indices are LOCAL to this rank's intersection (the whole patch at np=1); amr_isect_lo is
        ! GLOBAL; the coarse source qc is rank-LOCAL (identical at np=1: isect = patch, start_idx = 0)

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do fk = 0, amr_fine%p
            ck = amr_isect_lo(3) + fk/amr_fine%ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_fine%n
                cj = amr_isect_lo(2) + fj/amr_fine%ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_fine%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_fine%m
                    ci = amr_isect_lo(1) + fi/amr_fine%ref_ratio - ox
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

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: i

        if (.not. amr) return
        ! interior prolongation slopes read one coarse cell past the intersection - valid CONS ghosts first
        ! (ALL ranks call: pairwise halo). This runs BEFORE s_initialize_gpu_vars, and the halo packs the
        ! DEVICE state, so push the host ICs first (redundant on CPU; repeated later by the GPU-vars init).
        if (amr_xchg_coarse_ghosts) then
            do i = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_base(i)%sf]')
            end do
            call s_amr_exchange_coarse_cons_halo(q_cons_base)
            ! the exchange unpacks on the device; the init prolongation below runs on the HOST
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
            end do
        end if
        if (.not. amr_rank_owns_patch) return
        call s_interpolate_coarse_to_fine(q_cons_base)
        ! prolonged fine state to the device (host prolongation; device consumers: ghost-fill/RK/RHS kernels)
        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_fine%q_cons(i)%sf]')
        end do

    end subroutine s_populate_amr_fine

    !> Volume-weighted restriction for a single variable pair. Reads from qf (fine, must include interior 0:amr_fine%m etc.); writes
    !! to qc (coarse, over the patch).
    impure subroutine s_restrict_one_var(qf, qc)

        type(scalar_field), intent(in)    :: qf
        type(scalar_field), intent(inout) :: qc
        integer                           :: ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild, ox, oy, oz
        real(wp)                          :: acc

        ! isect loop indices are GLOBAL (this rank's covered coarse cells); the coarse target qc is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        nchild = amr_fine%ref_ratio
        if (n_glb > 0) nchild = nchild*amr_fine%ref_ratio
        if (p_glb > 0) nchild = nchild*amr_fine%ref_ratio
        do ck = amr_isect_lo(3), merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
            fk0 = (ck - amr_isect_lo(3))*amr_fine%ref_ratio
            do cj = amr_isect_lo(2), merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
                fj0 = (cj - amr_isect_lo(2))*amr_fine%ref_ratio
                do ci = amr_isect_lo(1), amr_isect_hi(1)
                    fi0 = (ci - amr_isect_lo(1))*amr_fine%ref_ratio
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
    !! coarse target (SP2: a scratch buffer) - never level-0. Device kernel (per-cell arithmetic identical to s_restrict_one_var,
    !! which remains the host path for the init-time diagnostics).
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt

        if (.not. amr_rank_owns_patch) return
        call s_restrict_all_vars(amr_fine%q_cons, coarse_tgt)

    end subroutine s_restrict_fine_to_coarse

    !> Device restriction kernel over all sys_size variable pairs (fine source qf -> coarse target qc).
    impure subroutine s_restrict_all_vars(qf, qc)

        type(scalar_field), dimension(sys_size), intent(in)    :: qf
        type(scalar_field), dimension(sys_size), intent(inout) :: qc
        integer                                                :: i, ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild, ox, oy, oz, rr
        integer                                                :: c1lo, c1hi, c2lo, c2hi, c3lo, c3hi, dj_hi, dk_hi
        real(wp)                                               :: acc

        ! isect loop indices are GLOBAL (this rank's covered coarse cells); the coarse target is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        rr = amr_fine%ref_ratio
        nchild = rr
        if (n_glb > 0) nchild = nchild*rr
        if (p_glb > 0) nchild = nchild*rr
        c1lo = amr_isect_lo(1); c1hi = amr_isect_hi(1)
        c2lo = amr_isect_lo(2); c2hi = merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
        c3lo = amr_isect_lo(3); c3hi = merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        $:GPU_PARALLEL_LOOP(collapse=4, private='[fi0, fj0, fk0, acc, ddj, ddk]')
        do i = 1, sys_size
            do ck = c3lo, c3hi
                do cj = c2lo, c2hi
                    do ci = c1lo, c1hi
                        fi0 = (ci - c1lo)*rr; fj0 = (cj - c2lo)*rr; fk0 = (ck - c3lo)*rr
                        acc = 0._wp
                        do ddk = 0, dk_hi
                            do ddj = 0, dj_hi
                                acc = acc + real(qf(i)%sf(fi0, fj0 + ddj, fk0 + ddk), wp) + real(qf(i)%sf(fi0 + 1, fj0 + ddj, &
                                                 & fk0 + ddk), wp)
                            end do
                        end do
                        qc(i)%sf(ci - ox, cj - oy, ck - oz) = acc/real(nchild, wp)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_restrict_all_vars

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
        ! host restriction path: the scratch target is host-only and the fine state is host-current at init
        do i = 1, sys_size
            call s_restrict_one_var(amr_fine%q_cons(i), scratch(i))
        end do
        err = 0._wp
        do i = 1, sys_size
            do ck = amr_isect_lo(3), merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
                do cj = amr_isect_lo(2), merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
                    do ci = amr_isect_lo(1), amr_isect_hi(1)
                        e = abs(real(scratch(i)%sf(ci - ox, cj - oy, ck - oz), wp) - real(q_cons_base(i)%sf(ci - ox, cj - oy, &
                                & ck - oz), wp))
                        if (e > err) err = e
                    end do
                end do
            end do
        end do
        print '(A,ES12.4)', ' [amr] restrict-prolong conservation err = ', err  ! every rank with fine cells prints
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
        ! save coarse coords to bounce buffers, then copy fine coords into global arrays
        sw_x_cb = x_cb; sw_x_cc = x_cc; sw_dx = dx
        if (n_glb > 0) then; sw_y_cb = y_cb; sw_y_cc = y_cc; sw_dy = dy; end if
        if (p_glb > 0) then; sw_z_cb = z_cb; sw_z_cc = z_cc; sw_dz = dz; end if
        x_cb(-1:amr_fine%m) = amr_fine%x_cb(-1:amr_fine%m)
        x_cc(0:amr_fine%m) = amr_fine%x_cc(0:amr_fine%m)
        dx(0:amr_fine%m) = amr_fine%dx(0:amr_fine%m)
        if (n_glb > 0) then
            y_cb(-1:amr_fine%n) = amr_fine%y_cb(-1:amr_fine%n)
            y_cc(0:amr_fine%n) = amr_fine%y_cc(0:amr_fine%n)
            dy(0:amr_fine%n) = amr_fine%dy(0:amr_fine%n)
        end if
        if (p_glb > 0) then
            z_cb(-1:amr_fine%p) = amr_fine%z_cb(-1:amr_fine%p)
            z_cc(0:amr_fine%p) = amr_fine%z_cc(0:amr_fine%p)
            dz(0:amr_fine%p) = amr_fine%dz(0:amr_fine%p)
        end if
        ! sync the swapped extents/bounds/coordinates to the device: RHS kernels read the
        ! device copies of these GPU_DECLARE'd globals (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse()

        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        ! restore full coarse coords from bounce buffers
        x_cb = sw_x_cb; x_cc = sw_x_cc; dx = sw_dx
        if (n_glb > 0) then; y_cb = sw_y_cb; y_cc = sw_y_cc; dy = sw_dy; end if
        if (p_glb > 0) then; z_cb = sw_z_cb; z_cc = sw_z_cc; dz = sw_dz; end if
        ! sync the restored coarse extents/bounds/coordinates back to the device
        call s_amr_sync_grid_state_to_device()

    end subroutine s_amr_restore_coarse

    !> Push the (host-side) global grid state to its device copies after a swap/restore. m/n/p, idwint/idwbuff, and the coordinate
    !! arrays are GPU_DECLARE'd; kernels read the device copies. No-op on CPU builds.
    impure subroutine s_amr_sync_grid_state_to_device()

        $:GPU_UPDATE(device='[m, n, p, idwint, idwbuff]')
        $:GPU_UPDATE(device='[x_cb, x_cc, dx]')
        if (n_glb > 0) then
            $:GPU_UPDATE(device='[y_cb, y_cc, dy]')
        end if
        if (p_glb > 0) then
            $:GPU_UPDATE(device='[z_cb, z_cc, dz]')
        end if

    end subroutine s_amr_sync_grid_state_to_device

    !> Fill the fine ghost shell of q_fine by conservative-linear prolongation from q_coarse (device kernel: reads the coarse source
    !! and writes the fine target in device memory). floor/modulo mapping is valid for negative fine indices (ghosts). Interior
    !! untouched.
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, q_fine)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_coarse
        type(scalar_field), dimension(sys_size), intent(inout) :: q_fine
        integer                                                :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        integer                                                :: rr, lo1, lo2, lo3, fm, fn, fp, b1, e1, b2, e2, b3, e3
        logical                                                :: d2, d3
        real(wp)                                               :: u0, sx, sy, sz, xix, xiy, xiz

        ! fine indices are LOCAL to this rank's intersection; amr_isect_lo is GLOBAL; the coarse source q_coarse is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_fine%ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        fm = amr_fine%m; fn = amr_fine%n; fp = amr_fine%p
        b1 = amr_fine%idwbuff(1)%beg; e1 = amr_fine%idwbuff(1)%end
        b2 = amr_fine%idwbuff(2)%beg; e2 = amr_fine%idwbuff(2)%end
        b3 = amr_fine%idwbuff(3)%beg; e3 = amr_fine%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=4, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        ! skip the interior: only the ghost shell is filled
                        if (.not. (fi >= 0 .and. fi <= fm .and. fj >= 0 .and. fj <= fn .and. fk >= 0 .and. fk <= fp)) then
                            ck = 0; xiz = 0._wp
                            if (d3) then
                                ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                                xiz = (real(modulo(fk, rr), wp) - 0.5_wp)*0.5_wp
                            end if
                            cj = 0; xiy = 0._wp
                            if (d2) then
                                cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                                xiy = (real(modulo(fj, rr), wp) - 0.5_wp)*0.5_wp
                            end if
                            ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                            xix = (real(modulo(fi, rr), wp) - 0.5_wp)*0.5_wp
                            u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                            sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), &
                                        & wp))
                            sy = 0._wp
                            if (d2) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                & cj - 1, ck), wp))
                            sz = 0._wp
                            if (d3) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, cj, &
                                & ck - 1), wp))
                            q_fine(i)%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_fine_ghosts

    !> Lerp the fine ghost shell of q_tgt between q_a (coarse t^n) and q_b (coarse t^{n+1}) at time fraction th (device kernel).
    !! Interior untouched.
    impure subroutine s_amr_lerp_fine_ghosts(q_a, q_b, q_tgt, th)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_a, q_b
        type(scalar_field), dimension(sys_size), intent(inout) :: q_tgt
        real(wp), intent(in)                                   :: th
        integer                                                :: i, fi, fj, fk, fm, fn, fp, b1, e1, b2, e2, b3, e3

        fm = amr_fine%m; fn = amr_fine%n; fp = amr_fine%p
        b1 = amr_fine%idwbuff(1)%beg; e1 = amr_fine%idwbuff(1)%end
        b2 = amr_fine%idwbuff(2)%beg; e2 = amr_fine%idwbuff(2)%end
        b3 = amr_fine%idwbuff(3)%beg; e3 = amr_fine%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        if (.not. (fi >= 0 .and. fi <= fm .and. fj >= 0 .and. fj <= fn .and. fk >= 0 .and. fk <= fp)) then
                            q_tgt(i)%sf(fi, fj, fk) = (1._wp - th)*real(q_a(i)%sf(fi, fj, fk), wp) + th*real(q_b(i)%sf(fi, fj, &
                                  & fk), wp)
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_lerp_fine_ghosts

    !> Device copy q_src -> q_dst over [b1:e1, b2:e2, b3:e3] for all sys_size fields (RK step-entry backup).
    impure subroutine s_amr_copy_fine_fields(q_src, q_dst, b1, e1, b2, e2, b3, e3)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_src
        type(scalar_field), dimension(sys_size), intent(inout) :: q_dst
        integer, intent(in)                                    :: b1, e1, b2, e2, b3, e3
        integer                                                :: i, fi, fj, fk

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        q_dst(i)%sf(fi, fj, fk) = q_src(i)%sf(fi, fj, fk)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_copy_fine_fields

    !> Device RK stage update over the fine interior: q = (c1*q + c2*q_stor + c3*dt_in*rhs)/c4 (compute in wp, store stp). Mirrors
    !! the coarse non-IGR rk_coef form in s_tvd_rk.
    impure subroutine s_amr_fine_rk_update(q_upd, q_stor, q_rhs, c1, c2, c3, c4, dt_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_upd
        type(scalar_field), dimension(sys_size), intent(in)    :: q_stor, q_rhs
        real(wp), intent(in)                                   :: c1, c2, c3, c4, dt_in
        integer                                                :: i, fi, fj, fk, fm, fn, fp

        fm = amr_fine%m; fn = amr_fine%n; fp = amr_fine%p
        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = 1, sys_size
            do fk = 0, fp
                do fj = 0, fn
                    do fi = 0, fm
                        q_upd(i)%sf(fi, fj, fk) = (c1*real(q_upd(i)%sf(fi, fj, fk), wp) + c2*real(q_stor(i)%sf(fi, fj, fk), &
                              & wp) + c3*dt_in*real(q_rhs(i)%sf(fi, fj, fk), wp))/c4
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update

    !> Exchange the coarse conservative ghost layers at internal rank boundaries (physical-boundary ghosts untouched; per direction
    !! beg then end, mirroring s_populate_variables_buffers' dispatch). The solver never fills CONS ghosts (only prim), so ranks
    !! whose fine ghost-fill or prolongation stencil leaves their interior need this first. ALL ranks must call together (pairwise
    !! exchange per internal neighbor).
    impure subroutine s_amr_exchange_coarse_cons_halo(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons

        if (bc_x%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, -1, sys_size)
        if (bc_x%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 1, 1, sys_size)
        if (n_glb > 0) then
            if (bc_y%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, -1, sys_size)
            if (bc_y%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 2, 1, sys_size)
        end if
        if (p_glb > 0) then
            if (bc_z%beg >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, -1, sys_size)
            if (bc_z%end >= 0) call s_mpi_sendrecv_variables_buffers(q_cons, 3, 1, sys_size)
        end if

    end subroutine s_amr_exchange_coarse_cons_halo

    !> Advance the fine level through RK stage s (same dt as level-0, no subcycling). Called BETWEEN the coarse RHS and the coarse
    !! RK update, so q_cons_coarse is the coarse STAGE-ENTRY state. Fine prim ghosts are obtained by widening idwint to the fine
    !! buffer so the cons->prim conversion inside s_compute_rhs covers the (prolonged) cons ghost shell; the BC populate is skipped
    !! via amr_in_fine_advance. The update mirrors the coarse non-IGR rk_coef form in s_tvd_rk.
    impure subroutine s_advance_amr_fine_stage(s, coefs, q_cons_coarse, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
        & time_avg)

        integer, intent(in)                                        :: s, t_step
        real(wp), intent(in)                                       :: coefs(4)
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons_coarse
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        real(wp), intent(inout)                                    :: time_avg

        if (.not. amr) return
        ! valid coarse CONS ghosts for the ghost prolongation below (ALL ranks call: pairwise halo)
        if (amr_xchg_coarse_ghosts) call s_amr_exchange_coarse_cons_halo(q_cons_coarse)
        if (.not. amr_rank_owns_patch) return

        ! ghost prolongation from the coarse stage-entry conservative state (device kernel reads the
        ! device-current coarse stage-entry state directly)
        call s_amr_fill_fine_ghosts(q_cons_coarse, amr_fine%q_cons)

        ! continuation faces (the patch spans a rank boundary there): overwrite the prolonged ghosts with the
        ! neighbor's true fine data (no exchange fires at np=1 / for fully-contained patches)
        call s_mpi_sendrecv_amr_fine_halo(amr_fine%q_cons, amr_fine%m, amr_fine%n, amr_fine%p)

        ! step-entry backup for the SSP-RK combination (device copy over the current buffered extents)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_fine%q_cons, amr_fine%q_cons_stor, amr_fine%idwbuff(1)%beg, amr_fine%idwbuff(1)%end, &
                                        & amr_fine%idwbuff(2)%beg, amr_fine%idwbuff(2)%end, amr_fine%idwbuff(3)%beg, &
                                        & amr_fine%idwbuff(3)%end)
        end if

        amr_in_fine_advance = .true.
        call s_amr_swap_to_fine()
        idwint = amr_fine%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        $:GPU_UPDATE(device='[idwint]')
        call s_compute_rhs(amr_fine%q_cons, q_T_sf, amr_fine%q_prim, bc_type, amr_fine%rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
                           & time_avg, s)
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update (device kernel; mirror of the coarse non-IGR form)
        call s_amr_fine_rk_update(amr_fine%q_cons, amr_fine%q_cons_stor, amr_fine%rhs, coefs(1), coefs(2), coefs(3), coefs(4), dt)

    end subroutine s_advance_amr_fine_stage

    !> Subcycled fine advance (amr_subcycle): two dt/2 SSP-RK3 substeps AFTER the coarse step. q_old/q_new are the coarse t^n and
    !! t^{n+1} states; each stage's ghosts are the linear time interpolation at the stage time theta = (substep-1 + c_s)/2 with
    !! SSP-RK3 abscissae c = [0, 1, 1/2]. Fine flux registers are zeroed here and accumulate over all six stages (0.5*rk3_w each) so
    !! the end-of-step state reflux sees the time-averaged effective fine flux.
    impure subroutine s_advance_amr_fine_substeps(q_old, q_new, coefs, bc_type, q_T_sf, pb_in, rhs_pb, mv_in, rhs_mv, t_step, &
        & time_avg)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_old, q_new
        real(wp), dimension(:,:), intent(in)                       :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), intent(inout)                          :: q_T_sf
        real(stp), dimension(:,:,:,:,:), intent(inout)             :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)              :: rhs_pb, rhs_mv
        integer, intent(in)                                        :: t_step
        real(wp), intent(inout)                                    :: time_avg
        real(wp), parameter                                        :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer                                                    :: sub, s
        real(wp)                                                   :: th

        if (.not. amr) return
        ! valid coarse CONS ghosts on both lerp sources (ALL ranks call: pairwise halo); the exchanged t^n /
        ! t^{n+1} ghost layers make the prolonged patch-boundary ghosts correct even at rank boundaries
        if (amr_xchg_coarse_ghosts) then
            call s_amr_exchange_coarse_cons_halo(q_old)
            call s_amr_exchange_coarse_cons_halo(q_new)
        end if
        if (.not. amr_rank_owns_patch) return

        ! fill both lerp sources once: ghost shells prolonged from coarse t^n and t^{n+1} (device kernels
        ! read the device-current coarse states directly)
        call s_amr_fill_fine_ghosts(q_old, amr_fine%q_ghost_a)
        call s_amr_fill_fine_ghosts(q_new, amr_fine%q_ghost_b)
        call s_amr_zero_fine_registers()

        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp

                ! lerp the ghost shell into q_cons at the stage time (device kernel; interior untouched)
                call s_amr_lerp_fine_ghosts(amr_fine%q_ghost_a, amr_fine%q_ghost_b, amr_fine%q_cons, th)

                ! continuation-face ghosts AFTER the lerp: the substeps run in lockstep across ranks, so the
                ! neighbor's current fine q_cons is the same-time data (the lerp stays patch-boundary-only)
                call s_mpi_sendrecv_amr_fine_halo(amr_fine%q_cons, amr_fine%m, amr_fine%n, amr_fine%p)

                ! substep-entry backup for the SSP-RK combination (device copy, interior only)
                if (s == 1) then
                    call s_amr_copy_fine_fields(amr_fine%q_cons, amr_fine%q_cons_stor, 0, amr_fine%m, 0, amr_fine%n, 0, amr_fine%p)
                end if

                amr_in_fine_advance = .true.
                call s_amr_swap_to_fine()
                idwint = amr_fine%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
                $:GPU_UPDATE(device='[idwint]')
                call s_compute_rhs(amr_fine%q_cons, q_T_sf, amr_fine%q_prim, bc_type, amr_fine%rhs, pb_in, rhs_pb, mv_in, rhs_mv, &
                                   & t_step, time_avg, s)
                call s_amr_restore_coarse()
                amr_in_fine_advance = .false.

                ! RK stage update at the FINE time step (device kernel)
                call s_amr_fine_rk_update(amr_fine%q_cons, amr_fine%q_cons_stor, amr_fine%rhs, coefs(s, 1), coefs(s, 2), coefs(s, &
                                          & 3), coefs(s, 4), amr_dt_fine)
            end do
        end do

    end subroutine s_advance_amr_fine_substeps

    !> Regrid: tag by relative density gradient, form the padded/clamped bounding box, then every rank rebuilds its LOCAL piece
    !! (copy old-fine on overlap via q_cons_stor bounce; prolong new cells). Fine data never migrates: fine cells for coarse cell c
    !! always live on rank(c), so the overlap copy is rank-local by construction. Called between steps only. No-op if nothing is
    !! tagged or the box is unchanged.
    impure subroutine s_amr_regrid(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
        integer                                                :: lo(3), hi(3), old_ilo(3), old_m1, old_n1, old_p1
        integer                                                :: ci, cj, ck, fi, fj, fk, ofi, ofj, ofk, sh(3), i, d
        integer                                                :: sidx(3), lo_r, hi_r, tg_lo(3), tg_hi(3)
        real(wp)                                               :: r0, g
        logical                                                :: tagged

        ! host consumer: regrid (tag sweep + rebuild prolongation run on host every regrid_int steps;
        ! ALL ranks tag their local interior, so this is unguarded)

        ! valid coarse CONS ghosts at internal rank boundaries: the tag sweep reads +/-1 across seams and the
        ! rebuild prolongation reads past the new intersection (ALL ranks call: pairwise per-direction
        ! exchange; complete no-op at np=1). The exchange unpacks on the device; the host pull below carries
        ! the fresh ghosts to the host consumers.

        call s_amr_exchange_coarse_cons_halo(q_cons_base)
        do i = 1, sys_size
            $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
        end do

        ! 1) tag + bounding box, emitting GLOBAL indices: sweep every cell except GLOBAL cells 0 and m_glb etc.
        !    (the np=1 bounds - the swept union is decomposition-invariant, so the reduced box matches np=1
        !    exactly; seam-adjacent cells read the exchanged ghosts). No rank-dependent branching before the
        !    reductions: every rank reaches both allreduces.

        sidx = 0
        sidx(1) = start_idx(1)
        if (n_glb > 0) sidx(2) = start_idx(2)
        if (p_glb > 0) sidx(3) = start_idx(3)
        tg_lo = 0; tg_hi = 0
        tg_lo(1) = merge(1, 0, sidx(1) == 0); tg_hi(1) = merge(m - 1, m, sidx(1) + m == m_glb)
        if (n_glb > 0) then; tg_lo(2) = merge(1, 0, sidx(2) == 0); tg_hi(2) = merge(n - 1, n, sidx(2) + n == n_glb); end if
        if (p_glb > 0) then; tg_lo(3) = merge(1, 0, sidx(3) == 0); tg_hi(3) = merge(p - 1, p, sidx(3) + p == p_glb); end if
        lo = huge(1); hi = -huge(1); tagged = .false.
        do ck = tg_lo(3), tg_hi(3)
            do cj = tg_lo(2), tg_hi(2)
                do ci = tg_lo(1), tg_hi(1)
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

        ! 2) pad + clamp: domain margin and the rank-fit size cap (amr_maxc_fit: any box that small intersects
        !    EVERY rank within its scratch constraint, so the patch moves freely across rank boundaries);
        !    collapsed dims pinned to 0. All inputs are globally identical, so every rank computes the same box.
        lo(1) = max(lo(1) - amr_buf, buff_size); hi(1) = min(hi(1) + amr_buf, m_glb - buff_size)
        if (hi(1) - lo(1) + 1 > amr_maxc_fit(1)) then
            if (proc_rank == 0) print '(A)', ' [amr] WARNING: tagged region exceeds max patch; clamping'
            hi(1) = lo(1) + amr_maxc_fit(1) - 1
        end if
        if (n_glb > 0) then
            lo(2) = max(lo(2) - amr_buf, buff_size); hi(2) = min(hi(2) + amr_buf, n_glb - buff_size)
            if (hi(2) - lo(2) + 1 > amr_maxc_fit(2)) hi(2) = lo(2) + amr_maxc_fit(2) - 1
        else
            lo(2) = 0; hi(2) = 0
        end if
        if (p_glb > 0) then
            lo(3) = max(lo(3) - amr_buf, buff_size); hi(3) = min(hi(3) + amr_buf, p_glb - buff_size)
            if (hi(3) - lo(3) + 1 > amr_maxc_fit(3)) hi(3) = lo(3) + amr_maxc_fit(3) - 1
        else
            lo(3) = 0; hi(3) = 0
        end if
        if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return  ! tag box confined to the domain margin; keep the old region
        if (all(lo == amr_fine%region%lo) .and. all(hi == amr_fine%region%hi)) return

        ! 3) ranks with OLD fine cells stash their old fine interior in the (dead-between-steps) RK bounce buffer
        !    (amr_rank_owns_patch still holds the pre-rebuild value here; ranks without get old_m1 = -1 => no-op)
        old_ilo = amr_isect_lo
        old_m1 = amr_fine%m; old_n1 = amr_fine%n; old_p1 = amr_fine%p
        if (amr_rank_owns_patch) then
            ! host consumer: regrid stash + rebuild run on host; the rebuilt state is pushed back below
            do i = 1, sys_size
                $:GPU_UPDATE(host='[amr_fine%q_cons(i)%sf]')
            end do
            do i = 1, sys_size
                amr_fine%q_cons_stor(i)%sf(0:old_m1,0:old_n1,0:old_p1) = amr_fine%q_cons(i)%sf(0:old_m1,0:old_n1,0:old_p1)
            end do
        end if

        ! 4) rebuild geometry (region/mirrors/extents/xchg flag on all ranks; coord fill owner-guarded inside), then
        !    every rank with a NEW intersection prolongs its local piece and overwrites the overlap with its own old
        !    fine data (rank-local by construction)
        call s_set_amr_fine_geometry(lo, hi)
        if (proc_rank == 0) then
            print '(A,I0,A,I0,A,I0,A)', ' [amr] regrid: box x ', lo(1), ':', hi(1), ' (', (hi(1) - lo(1) + 1), ' coarse cells)'
        end if
        if (.not. amr_rank_owns_patch) return
        call s_interpolate_coarse_to_fine(q_cons_base)
        sh = 2*(amr_isect_lo - old_ilo)  ! old LOCAL fine index = new LOCAL fine index + sh (per dim; collapsed dims sh=0)
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
        ! rebuilt fine state back to the device for the stepping kernels
        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_fine%q_cons(i)%sf]')
        end do
        ! continuation-face ghosts from the neighbors' rebuilt interiors (fresh fine ghosts for the next step)
        call s_mpi_sendrecv_amr_fine_halo(amr_fine%q_cons, amr_fine%m, amr_fine%n, amr_fine%p)

    end subroutine s_amr_regrid

    !> Global Sum(dV*U) for continuity (var 1) and energy (eqn_idx%E) over the level-0 interior. First call (finalize_report=F)
    !! stores the baseline; the finalize call prints the relative drift (expected small nonzero until SP4 refluxing).
    impure subroutine s_amr_conservation_defect(q_cons_base, finalize_report)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        logical, intent(in)                                 :: finalize_report
        real(wp)                                            :: sm, se, dv, sm_glb, se_glb
        integer                                             :: ci, cj, ck

        if (.not. amr) return
        ! host consumer: diagnostics (host sum over exactly the two summed fields). The init baseline call
        ! runs BEFORE s_initialize_gpu_vars pushes the ICs to the device, so it must NOT pull the
        ! (uninitialized) device copies.
        if (finalize_report) then
            $:GPU_UPDATE(host='[q_cons_base(1)%sf, q_cons_base(eqn_idx%E)%sf]')
        end if
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
        do ck = amr_isect_lo(3), merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
            do cj = amr_isect_lo(2), merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
                do ci = amr_isect_lo(1), amr_isect_hi(1)
                    dvc = dx(ci - ox)
                    if (n_glb > 0) dvc = dvc*dy(cj - oy)
                    if (p_glb > 0) dvc = dvc*dz(ck - oz)
                    si_c = si_c + dvc*real(cscr(1)%sf(ci - ox, cj - oy, ck - oz), wp)
                end do
            end do
        end do
        errc = abs(si_f - si_c)/max(abs(si_f), 1.e-30_wp)
        ! every rank with fine cells prints
        print '(A,ES12.4)', ' [amr] prolong linear-reproduction err = ', errb
        print '(A,ES12.4)', ' [amr] restrict independent-integral err = ', errc
        deallocate (cscr(1)%sf); deallocate (cscr)

    end subroutine s_amr_operator_checks

    !> minmod slope limiter: 0 if a,b differ in sign, else the smaller-magnitude argument.
    pure elemental function minmod(a, b) result(m)

        $:GPU_ROUTINE(parallelism='[seq]')
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
            @:DEALLOCATE(amr_fine%q_cons(i)%sf)
            @:DEALLOCATE(amr_fine%q_cons_stor(i)%sf)
            @:DEALLOCATE(amr_fine%q_prim(i)%sf)
            @:DEALLOCATE(amr_fine%rhs(i)%sf)
            @:DEALLOCATE(amr_fine%q_ghost_a(i)%sf)
            @:DEALLOCATE(amr_fine%q_ghost_b(i)%sf)
        end do
        @:DEALLOCATE(amr_fine%q_cons)
        @:DEALLOCATE(amr_fine%q_cons_stor)
        @:DEALLOCATE(amr_fine%q_prim)
        @:DEALLOCATE(amr_fine%rhs)
        @:DEALLOCATE(amr_fine%q_ghost_a)
        @:DEALLOCATE(amr_fine%q_ghost_b)
        if (allocated(amr_fine%x_cb)) deallocate (amr_fine%x_cb, amr_fine%x_cc, amr_fine%dx)
        if (allocated(amr_fine%y_cb)) deallocate (amr_fine%y_cb, amr_fine%y_cc, amr_fine%dy)
        if (allocated(amr_fine%z_cb)) deallocate (amr_fine%z_cb, amr_fine%z_cc, amr_fine%dz)
        if (allocated(sw_x_cb)) deallocate (sw_x_cb, sw_x_cc, sw_dx)
        if (allocated(sw_y_cb)) deallocate (sw_y_cb, sw_y_cc, sw_dy)
        if (allocated(sw_z_cb)) deallocate (sw_z_cb, sw_z_cc, sw_dz)

    end subroutine s_finalize_amr_module

end module m_amr
