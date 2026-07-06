!>
!!@file
!!@brief Contains module m_amr

#:include 'macros.fpp'

!> @brief Block-structured AMR: up to amr_max_blocks 2:1 refined level-1 blocks, advanced with the shared solver via grid-state
!! swap, conservatively coupled to level 0 (ghost prolongation, Berger-Colella flux reflux, restriction), with optional dt/2
!! subcycling and dynamic regrid.
module m_amr

#ifdef MFC_MPI
    use mpi  !< MPI-IO for the parallel_io AMR restart file
#endif

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_constants, only: num_fluids_max, model_eqns_6eq
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_mpi_proxy, only: s_mpi_abort, s_initialize_amr_mpi_buffers, s_mpi_sendrecv_amr_fine_halo
    use m_mpi_common, only: s_mpi_allreduce_integer_min, s_mpi_allreduce_integer_max, s_mpi_allreduce_sum, &
        & s_mpi_allreduce_integer_sum, s_mpi_sendrecv_variables_buffers
    use m_rhs, only: s_compute_rhs
    use m_phase_change, only: s_infinite_relaxation_k
    use m_amr_registers, only: s_amr_zero_fine_registers
    use m_rank_timing, only: s_rank_time_tic, s_rank_time_toc
    use m_ibm, only: s_ibm_alloc_fine, s_ibm_setup_fine, s_ibm_swap_to_fine, s_ibm_restore_from_fine, s_ibm_correct_state, &
        & s_update_mib, moving_immersed_boundary_flag, num_gps
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs
    use m_weno, only: s_compute_weno_coefficients
    use m_acoustic_src, only: acoustic_supp_lo, acoustic_supp_hi

    implicit none

    private
    public :: t_level, amr_maxc, amr_dt_fine, s_initialize_amr_module, s_populate_amr_fine, s_interpolate_coarse_to_fine, &
        & s_restrict_fine_to_coarse, s_amr_conservation_check, s_finalize_amr_module, s_amr_swap_to_fine, s_amr_restore_coarse, &
        & s_amr_fill_fine_ghosts, s_amr_operator_checks, s_advance_amr_fine_stage, s_advance_amr_fine_substeps, &
        & s_amr_conservation_defect, s_set_amr_fine_geometry, s_amr_regrid, s_write_amr_restart, s_read_amr_restart, &
        & s_amr_relax_fine, s_amr_setup_ib

    !> Fine-level time step for subcycling (= 0.5*dt after init; 0 when amr is off).
    real(wp) :: amr_dt_fine = 0._wp

    !> Realizability floor for prolonged Euler-Euler bubble POSITIVE moments (radius nR, and non-polytropic partial pressure npb /
    !! vapor mass nmv): a positive fraction of the coarse parent so the derived R = nR/n, pb, mv stay >= 0. Minmod already keeps a
    !! positive field positive, so this fires only under floating-point edge cases (conservation defect ~0 otherwise).
    real(wp), parameter :: bub_pos_frac = 1.0e-10_wp

    !> One refined level: its own grid + conservative fields. Field arrays are device-resident (@:ALLOCATE); coords/metadata
    !! host-only.
    type t_level
        integer                         :: ref_ratio
        type(t_box)                     :: region          !< block extent in parent (level-0) cell indices
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
        !> non-polytropic QBMM quadrature side-state on the block (nnode x nb per cell). pb/mv evolve cell-locally (their rhs reads
        !! only the local cell + the block's own moment fluxes), so the fine treatment is prolong -> advance -> restrict with no
        !! reflux; ghosts feed the widened- idwint conversions and are prolonged piecewise-constant (CHyQMOM realizability, like the
        !! moments).
        real(stp), allocatable :: pb_f(:,:,:,:,:), mv_f(:,:,:,:,:)        !< fine pb/mv (ghost-inclusive)
        real(stp), allocatable :: pb_stor(:,:,:,:,:), mv_stor(:,:,:,:,:)  !< SSP-RK step-entry backup (also the regrid bounce)
        !> subcycle ghost-lerp sources at coarse t^n / t^{n+1} (ghost shell only): ghost pb feeds the mixture pressure in the
        !! widened conversion, so it needs the same time fidelity as q_cons
        real(stp), allocatable :: pb_ghost_a(:,:,:,:,:), mv_ghost_a(:,:,:,:,:)
        real(stp), allocatable :: pb_ghost_b(:,:,:,:,:), mv_ghost_b(:,:,:,:,:)
        real(wp), allocatable  :: rhs_pb_f(:,:,:,:,:), rhs_mv_f(:,:,:,:,:)  !< fine rhs (ghost-inclusive like rhs)
    end type t_level

    !> Fixed pool of refined-block slots (at init one slot is active; dynamic regrid activates up to amr_max_blocks). The working
    !! slot amr_cur (m_global_parameters) selects which slot every per-block routine operates on.
    type(t_level), allocatable :: amr_slots(:)
    integer                    :: amr_maxc(3)  !< max coarse block cells per dim: (m_glb+1)/2 etc.; 1 for collapsed dims

    !> Regrid box size cap per dim (fixed for the run, identical on all ranks; 1 in collapsed dims): a box of at most min over ranks
    !! of (local extent + 1)/2 cells intersects EVERY rank in at most (its extent + 1)/2 cells, so the per-rank scratch constraint
    !! 2*(isect cells) - 1 <= local extent holds by construction. Equals amr_maxc at np=1.
    integer :: amr_maxc_fit(3) = 1

    !> Saved coarse-level global state for swap/restore
    integer               :: sw_m, sw_n, sw_p
    type(int_bounds_info) :: sw_idwint(3), sw_idwbuff(3)
    logical               :: sw_acoustic_source
    logical               :: amr_swapped = .false.  !< paired-swap guard: a nested swap would corrupt the bounce buffers silently
    !> True when the coarse grid is nonuniform in an allowed way (the 2D-axisymmetric half-width axis cell): the WENO reconstruction
    !! coefficients are then per-cell, so the fine advance must recompute them for the block's own (uniform) grid on every swap and
    !! restore the coarse ones after. False on fully uniform grids - the recompute is skipped and behavior is bit-identical.
    logical               :: amr_weno_coef_recompute = .false.
    real(wp), allocatable :: sw_x_cb(:), sw_x_cc(:), sw_dx(:)
    real(wp), allocatable :: sw_y_cb(:), sw_y_cc(:), sw_dy(:)
    real(wp), allocatable :: sw_z_cb(:), sw_z_cc(:), sw_dz(:)

    !> Conservation-defect baselines (level-0 interior integrals at init; per-fluid masses + energy)
    real(wp) :: amr_mass0(num_fluids_max) = 0._wp, amr_energy0 = 0._wp

    !> True (identically on all ranks) iff some rank's fine ghost-fill stencil reads its coarse GHOST cells - the solver populates
    !! only PRIM ghosts, so the CONS ghosts the fill prolongs from must be halo-exchanged first. Never true at np=1 (block faces sit
    !! >= buff_size inside the domain).
    logical :: amr_xchg_coarse_ghosts = .false.

contains

    !> Build the static refined level-1 block. No-op unless amr. Called after level-0 grid (x_cb/dx ready) and time-steppers
    !! (sys_size/buff_size set). Preallocates all fine arrays at max size so regrid only needs to call s_set_amr_fine_geometry.
    impure subroutine s_initialize_amr_module()

        integer :: i, d, max_f1, max_f2, max_f3, islot
        integer :: mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi
        integer :: sidx(3), ext(3), maxc_loc(3), bad_loc, bad_glb, fit_d
        integer :: blk_lo(3), blk_hi(3)

        if (.not. amr) return

        amr_dt_fine = 0.5_wp*dt

        ! fixed pool of amr_max_blocks slots; init activates exactly one (slot amr_cur = 1); dynamic regrid clusters into up to
        ! amr_max_blocks slots
        allocate (amr_slots(1:amr_max_blocks))
        allocate (amr_region_lo_all(3, amr_max_blocks), amr_region_hi_all(3, amr_max_blocks))
        allocate (amr_isect_lo_all(3, amr_max_blocks), amr_isect_hi_all(3, amr_max_blocks))
        allocate (amr_owns_all(amr_max_blocks))
        amr_region_lo_all = 0; amr_region_hi_all = 0; amr_isect_lo_all = 0; amr_isect_hi_all = 0; amr_owns_all = .false.
        amr_num_blocks = 1
        amr_cur = 1

        ! Mirror decomposition: each rank holds the fine cells covering block /\ its own subdomain
        ! (np=1: the intersection is the whole block). buff_size is not available at checker time,
        ! so the geometric aborts below must live here.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        call s_amr_compute_isect(amr_block_beg, amr_block_end)

        ! the fine ghost shell and the reflux outside cells must stay inside the global domain (identical
        ! inputs on all ranks; every rank takes the same branch)
        if (amr_block_beg(1) < buff_size .or. amr_block_end(1) > m_glb - buff_size .or. (n_glb > 0 .and. (amr_block_beg(2) &
            & < buff_size .or. amr_block_end(2) > n_glb - buff_size)) .or. (p_glb > 0 .and. (amr_block_beg(3) < buff_size &
            & .or. amr_block_end(3) > p_glb - buff_size))) then
            call s_mpi_abort('amr block must lie at least buff_size cells inside the domain boundaries')
        end if

        ! Scratch constraint: the fine advance reuses the solver scratch (m_rhs/WENO/Riemann work arrays),
        ! which is sized to THIS rank's local grid, so each rank's fine extent must fit its local extent.
        ! (np=1: the checker's 2*block-1 <= m_glb bound makes this a no-op.)
        bad_loc = 0
        if (amr_rank_owns_block) then
            if (2*(amr_isect_hi(1) - amr_isect_lo(1) + 1) - 1 > m) bad_loc = 1
            if (n_glb > 0 .and. 2*(amr_isect_hi(2) - amr_isect_lo(2) + 1) - 1 > n) bad_loc = 1
            if (p_glb > 0 .and. 2*(amr_isect_hi(3) - amr_isect_lo(3) + 1) - 1 > p) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        if (bad_glb == 1) then
            call s_mpi_abort('amr fine extent exceeds a rank local grid (solver scratch is local-sized): the block ' &
                             & // 'may cover at most about half of any rank subdomain per dimension; shrink the block ' &
                             & // 'or use fewer ranks')
        end if

        ! max coarse block cells per dim (upper bound for any future regrid box); 1 for collapsed dims
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

        ! preallocation cap for MY fine arrays: the largest intersection any block box can have with this
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

        ! Grid uniformity policy. The two spacing-uniformity consumers are both handled exactly:
        ! the fine-block ghost-shell coordinates extend by exact parent-cell bisection (reads
        ! sw_*_cb), and the spacing-dependent WENO reconstruction coefficients are recomputed for
        ! the active grid on every swap/restore when the grid is nonuniform anywhere (stretched
        ! grids, or 2D-axisymmetric's half-width axis cell dy(0) = dy/2). On fully uniform grids
        ! the flag stays false and behavior is bit-identical to the reuse path. The stretch_*
        ! flags are pre_process-only, so the grid itself is checked (this also catches
        ! externally generated grids).
        if (maxval(dx(0:m)) - minval(dx(0:m)) > 1.e-12_wp*maxval(dx(0:m))) then
            amr_weno_coef_recompute = .true.
        end if
        if (n_glb > 0) then
            if (maxval(dy(0:n)) - minval(dy(0:n)) > 1.e-12_wp*maxval(dy(0:n))) then
                amr_weno_coef_recompute = .true.
            end if
        end if
        if (p_glb > 0) then
            if (maxval(dz(0:p)) - minval(dz(0:p)) > 1.e-12_wp*maxval(dz(0:p))) then
                amr_weno_coef_recompute = .true.
            end if
        end if
        if (weno_order == 1) amr_weno_coef_recompute = .false.  ! order 1 has no grid-dependent coefficients

        ! preallocate every slot at max buffered extents (each sized exactly as the single legacy block): coords (host) +
        ! fields (device-resident @:ALLOCATE so s_compute_rhs kernels can run on the fine block; q_cons/q_prim/rhs get the
        ! Cray pointer setup mirroring m_time_steppers' field vectors). Loops always use the current slot's m/n/p.
        do islot = 1, amr_max_blocks
            amr_slots(islot)%ref_ratio = 2
            amr_slots(islot)%buff_size = buff_size
            allocate (amr_slots(islot)%x_cb(-1:max_f1), amr_slots(islot)%x_cc(0:max_f1), amr_slots(islot)%dx(0:max_f1))
            if (n_glb > 0) allocate (amr_slots(islot)%y_cb(-1:max_f2), amr_slots(islot)%y_cc(0:max_f2), &
                & amr_slots(islot)%dy(0:max_f2))
            if (p_glb > 0) allocate (amr_slots(islot)%z_cb(-1:max_f3), amr_slots(islot)%z_cc(0:max_f3), &
                & amr_slots(islot)%dz(0:max_f3))
            @:ALLOCATE(amr_slots(islot)%q_cons(1:sys_size))
            @:ALLOCATE(amr_slots(islot)%q_cons_stor(1:sys_size))
            @:ALLOCATE(amr_slots(islot)%q_prim(1:sys_size))
            @:ALLOCATE(amr_slots(islot)%rhs(1:sys_size))
            @:ALLOCATE(amr_slots(islot)%q_ghost_a(1:sys_size))
            @:ALLOCATE(amr_slots(islot)%q_ghost_b(1:sys_size))
            do i = 1, sys_size
                @:ALLOCATE(amr_slots(islot)%q_cons(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ALLOCATE(amr_slots(islot)%q_cons_stor(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ALLOCATE(amr_slots(islot)%q_prim(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                ! ghost-inclusive (mbuf) like q_cons/q_prim: the fine advance widens idwint to the buffer, so the cell-local
                ! chemistry reaction source writes rhs over the ghost shell too (the RK update and reflux read only 0:m interior)
                @:ALLOCATE(amr_slots(islot)%rhs(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ALLOCATE(amr_slots(islot)%q_ghost_a(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ALLOCATE(amr_slots(islot)%q_ghost_b(i)%sf(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi))
                @:ACC_SETUP_SFs(amr_slots(islot)%q_cons(i))
                @:ACC_SETUP_SFs(amr_slots(islot)%q_prim(i))
                @:ACC_SETUP_SFs(amr_slots(islot)%rhs(i))
                @:ACC_SETUP_SFs(amr_slots(islot)%q_cons_stor(i))
                @:ACC_SETUP_SFs(amr_slots(islot)%q_ghost_a(i))
                @:ACC_SETUP_SFs(amr_slots(islot)%q_ghost_b(i))
            end do
            if (qbmm .and. .not. polytropic) then
                @:ALLOCATE(amr_slots(islot)%pb_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ALLOCATE(amr_slots(islot)%mv_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ALLOCATE(amr_slots(islot)%pb_stor(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ALLOCATE(amr_slots(islot)%mv_stor(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ALLOCATE(amr_slots(islot)%rhs_pb_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                @:ALLOCATE(amr_slots(islot)%rhs_mv_f(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                if (amr_subcycle) then
                    @:ALLOCATE(amr_slots(islot)%pb_ghost_a(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                    @:ALLOCATE(amr_slots(islot)%mv_ghost_a(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                    @:ALLOCATE(amr_slots(islot)%pb_ghost_b(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                    @:ALLOCATE(amr_slots(islot)%mv_ghost_b(mbuf1_lo:mbuf1_hi, mbuf2_lo:mbuf2_hi, mbuf3_lo:mbuf3_hi, 1:nnode, 1:nb))
                end if
            end if
        end do

        ! per-slot fine-grid IB marker fields (static-body AMR); sized to the same max buffered fine extents
        ! as q_cons so the fine IB pipeline can resolve the body on the block
        if (ib) call s_ibm_alloc_fine(amr_max_blocks, mbuf1_lo, mbuf1_hi, mbuf2_lo, mbuf2_hi, mbuf3_lo, mbuf3_hi)

        ! set geometry (region, m/n/p, idwbuff, coordinates) for the initial block (slot amr_cur = 1).
        ! Under dynamic regrid with bodies the initial block gets the same body-containment
        ! expansion regrid boxes get (the moving-body containment guard requires it from step 1);
        ! for a static block (amr_regrid_int = 0) the user's placement is authoritative.
        blk_lo = amr_block_beg; blk_hi = amr_block_end
        if (ib .and. amr_regrid_int > 0) call s_amr_expand_box_over_bodies(blk_lo, blk_hi)
        call s_set_amr_fine_geometry(blk_lo, blk_hi)

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
    !! and whether it holds fine cells (amr_rank_owns_block: nonempty in all active dims). Must be called with the COARSE grid state
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
        amr_rank_owns_block = amr_isect_lo(1) <= amr_isect_hi(1)
        if (n_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(2) <= amr_isect_hi(2)
        if (p_glb > 0) amr_rank_owns_block = amr_rank_owns_block .and. amr_isect_lo(3) <= amr_isect_hi(3)

        ! non-polytropic QBMM: the fine rank-seam halo exchanges q_cons only, so a block spanning
        ! ranks would advance with unexchanged pb/mv seam ghosts - a silent wrong answer. Keep the
        ! block inside one rank subdomain (fewer ranks or reposition) until the halo carries pb/mv.
        if (qbmm .and. (.not. polytropic) .and. amr_rank_owns_block) then
            if (amr_isect_lo(1) /= lo(1) .or. amr_isect_hi(1) /= hi(1) .or. (n_glb > 0 .and. (amr_isect_lo(2) /= lo(2) &
                & .or. amr_isect_hi(2) /= hi(2))) .or. (p_glb > 0 .and. (amr_isect_lo(3) /= lo(3) .or. amr_isect_hi(3) /= hi(3)))) &
                & then
                call s_mpi_abort('amr with non-polytropic qbmm: the fine block spans a rank boundary, but the ' &
                                 & // 'fine seam halo does not carry the pb/mv side-state; keep the block within ' &
                                 & // 'a single rank subdomain (use fewer ranks or reposition the block)')
            end if
        end if

    end subroutine s_amr_compute_isect

    !> Set the fine level's geometry (region, intersection, extents, bounds, coordinates) for the box lo:hi. Arrays are preallocated
    !! at max size; this only updates metadata and refills coords. Collective: ALL ranks must call together (init and regrid do) -
    !! it also refreshes the allreduced amr_xchg_coarse_ghosts flag for the new box.
    impure subroutine s_set_amr_fine_geometry(lo, hi)

        integer, intent(in) :: lo(3), hi(3)
        integer             :: sidx(3), ext(3), nmar, bad_loc, bad_glb

        call s_amr_compute_isect(lo, hi)
        amr_slots(amr_cur)%region%lo = lo; amr_slots(amr_cur)%region%hi = hi
        amr_region_lo = lo; amr_region_hi = hi  ! global mirror for m_amr_registers (no use-cycle)
        ! stash this slot's mirrors so s_amr_select_slot can revisit each block during the per-block advance and the coarse capture
        amr_region_lo_all(:,amr_cur) = lo; amr_region_hi_all(:,amr_cur) = hi
        amr_isect_lo_all(:,amr_cur) = amr_isect_lo; amr_isect_hi_all(:,amr_cur) = amr_isect_hi
        amr_owns_all(amr_cur) = amr_rank_owns_block
        ! LOCAL fine extents cover this rank's intersection (the whole block at np=1); -1 in a dim whose
        ! intersection is empty, so 0..extent loops are no-ops on ranks without fine cells
        amr_slots(amr_cur)%m = 2*max(amr_isect_hi(1) - amr_isect_lo(1) + 1, 0) - 1
        amr_slots(amr_cur)%n = 0; amr_slots(amr_cur)%p = 0
        if (n_glb > 0) amr_slots(amr_cur)%n = 2*max(amr_isect_hi(2) - amr_isect_lo(2) + 1, 0) - 1
        if (p_glb > 0) amr_slots(amr_cur)%p = 2*max(amr_isect_hi(3) - amr_isect_lo(3) + 1, 0) - 1
        amr_slots(amr_cur)%idwbuff(1)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(1)%end = amr_slots(amr_cur)%m + buff_size
        amr_slots(amr_cur)%idwbuff(2)%beg = 0; amr_slots(amr_cur)%idwbuff(2)%end = 0
        amr_slots(amr_cur)%idwbuff(3)%beg = 0; amr_slots(amr_cur)%idwbuff(3)%end = 0
        if (n_glb > 0) then
            amr_slots(amr_cur)%idwbuff(2)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(2)%end = amr_slots(amr_cur)%n + buff_size
        end if
        if (p_glb > 0) then
            amr_slots(amr_cur)%idwbuff(3)%beg = -buff_size; amr_slots(amr_cur)%idwbuff(3)%end = amr_slots(amr_cur)%p + buff_size
        end if
        ! coord building only on ranks with fine cells (others never read their coord arrays); the parent origin
        ! is this rank's INTERSECTION start, converted to LOCAL indexing so the bisection reads its x_cb slice
        if (amr_rank_owns_block) then
            call s_build_level_coords(x_cb, lbound(x_cb, 1), amr_isect_lo(1) - start_idx(1), amr_slots(amr_cur)%m, &
                                      & amr_slots(amr_cur)%x_cb, amr_slots(amr_cur)%x_cc, amr_slots(amr_cur)%dx)
            if (n_glb > 0) call s_build_level_coords(y_cb, lbound(y_cb, 1), amr_isect_lo(2) - start_idx(2), amr_slots(amr_cur)%n, &
                & amr_slots(amr_cur)%y_cb, amr_slots(amr_cur)%y_cc, amr_slots(amr_cur)%dy)
            if (p_glb > 0) call s_build_level_coords(z_cb, lbound(z_cb, 1), amr_isect_lo(3) - start_idx(3), amr_slots(amr_cur)%p, &
                & amr_slots(amr_cur)%z_cb, amr_slots(amr_cur)%z_cc, amr_slots(amr_cur)%dz)
        end if

        ! Fine ghost prolongation reads up to nmar coarse cells past each face of the intersection; if that
        ! stencil leaves ANY rank's interior (block near/at/across a rank boundary), the coarse CONS ghosts it
        ! reads must be halo-exchanged before every fill (the solver populates only PRIM ghosts). All ranks
        ! agree on the flag, so the pairwise exchanges are called consistently.
        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        nmar = (buff_size + 1)/2 + 1
        bad_loc = 0
        if (amr_rank_owns_block) then
            if (amr_isect_lo(1) - sidx(1) < nmar .or. sidx(1) + ext(1) - amr_isect_hi(1) < nmar) bad_loc = 1
            if (n_glb > 0 .and. (amr_isect_lo(2) - sidx(2) < nmar .or. sidx(2) + ext(2) - amr_isect_hi(2) < nmar)) bad_loc = 1
            if (p_glb > 0 .and. (amr_isect_lo(3) - sidx(3) < nmar .or. sidx(3) + ext(3) - amr_isect_hi(3) < nmar)) bad_loc = 1
        end if
        call s_mpi_allreduce_integer_max(bad_loc, bad_glb)
        amr_xchg_coarse_ghosts = bad_glb == 1

    end subroutine s_set_amr_fine_geometry

    !> Conservative-linear prolongation for a single variable pair. Reads coarse interior/ghost from qc; writes fine interior to qf.
    !! Minmod-limited slopes.
    impure subroutine s_prolong_one_var(qc, qf, pos, inject)

        type(scalar_field), intent(in)    :: qc
        type(scalar_field), intent(inout) :: qf
        logical, optional, intent(in)     :: pos     !< floor the child at bub_pos_frac*u0 (bubble radius-moment realizability)
        logical, optional, intent(in)     :: inject  !< piecewise-constant (child = u0): QBMM moment realizability preservation
        integer                           :: fi, fj, fk, ci, cj, ck, ox, oy, oz
        real(wp)                          :: u0, sx, sy, sz, xix, xiy, xiz, child
        logical                           :: floor_pos, pw_const

        floor_pos = .false.; if (present(pos)) floor_pos = pos
        pw_const = .false.; if (present(inject)) pw_const = inject

        ! fine indices are LOCAL to this rank's intersection (the whole block at np=1); amr_isect_lo is
        ! GLOBAL; the coarse source qc is rank-LOCAL (identical at np=1: isect = block, start_idx = 0)

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    u0 = real(qc%sf(ci, cj, ck), wp)
                    sx = minmod(real(qc%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc%sf(ci - 1, cj, ck), wp))
                    sy = 0._wp
                    if (n_glb > 0) sy = minmod(real(qc%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc%sf(ci, cj - 1, ck), wp))
                    sz = 0._wp
                    if (p_glb > 0) sz = minmod(real(qc%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc%sf(ci, cj, ck - 1), wp))
                    if (pw_const) then
                        sx = 0._wp; sy = 0._wp; sz = 0._wp
                    end if
                    child = u0 + sx*xix + sy*xiy + sz*xiz
                    if (floor_pos) child = max(child, bub_pos_frac*u0)
                    qf%sf(fi, fj, fk) = child
                end do
            end do
        end do

    end subroutine s_prolong_one_var

    !> Conservative-linear prolongation: fill amr_fine interior from coarse (level-0), minmod-limited. Symmetric child offsets
    !! (+/-1/4 of a coarse cell) => the ref_ratio^d children average to the coarse value. Multi-fluid volume fractions take the
    !! sum-preserving closure path instead (single-fluid runs never branch, so their prolongation is untouched).
    impure subroutine s_interpolate_coarse_to_fine(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        integer                                             :: i, bstride

        bstride = 1
        if (bubbles_euler) bstride = (eqn_idx%bub%end - eqn_idx%bub%beg + 1)/nb
        do i = 1, sys_size
            if (num_fluids > 1 .and. i >= eqn_idx%adv%beg .and. i <= eqn_idx%adv%end) cycle
            if (chemistry .and. i >= eqn_idx%species%beg .and. i <= eqn_idx%species%end) cycle  ! sum/positivity closure below
            ! QBMM carries a bivariate 6-moment set per R0 bin whose CHyQMOM inversion requires realizability
            ! (variance c20 = m20/m00 - (m10/m00)^2 > 0); per-component minmod prolongation can break that joint
            ! constraint, so the whole bub block is injected piecewise-constant (each child inherits the coarse
            ! cell's realizable moment set exactly). Non-QBMM Euler-Euler bubbles instead floor their POSITIVE
            ! moments (radius nR, non-polytropic partial pressure npb / vapor mass nmv); the signed velocity moment
            ! nV (offset 1 in each bin's stride) prolongs freely.
            call s_prolong_one_var(q_cons_base(i), amr_slots(amr_cur)%q_cons(i), &
                                   & pos=bubbles_euler .and. .not. qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end &
                                   & .and. mod(i - eqn_idx%bub%beg, bstride) /= 1, &
                                   & inject=qbmm .and. i >= eqn_idx%bub%beg .and. i <= eqn_idx%bub%end)
        end do
        if (num_fluids > 1) call s_prolong_alphas_closure(q_cons_base, amr_slots(amr_cur)%q_cons)
        if (chemistry) call s_prolong_species_closure(q_cons_base, amr_slots(amr_cur)%q_cons)

    end subroutine s_interpolate_coarse_to_fine

    !> Sum-preserving volume-fraction prolongation (num_fluids > 1): fluids adv%beg..adv%end-1 are interpolated with minmod slopes
    !! under a SHARED per-cell limiter switch (a sign change for ANY fluid in a dim zeroes that dim's slope for ALL fluids, so the
    !! closure fluid's effective slope is limited consistently) and clamped to [0,1]; the last fluid closes alpha_n = 1 -
    !! sum(others), so sum(alpha) = 1 on the fine level by construction. For two fluids the closure is also in [0,1]; for >2 fluids
    !! any residual closure undershoot is handled by mpp_lim (required by the checker). Same fine/coarse index mapping as
    !! s_prolong_one_var.
    impure subroutine s_prolong_alphas_closure(qc, qf)

        type(scalar_field), dimension(sys_size), intent(in)    :: qc
        type(scalar_field), dimension(sys_size), intent(inout) :: qf
        integer                                                :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        real(wp)                                               :: xix, xiy, xiz, u0, sx, sy, sz, av, asum
        logical                                                :: shx, shy, shz

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    call s_alpha_shared_switch(qc, ci, cj, ck, shx, shy, shz)
                    asum = 0._wp
                    do i = eqn_idx%adv%beg, eqn_idx%adv%end - 1
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = 0._wp
                        if (shx) sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (n_glb > 0 .and. shy) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, &
                            & cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0 .and. shz) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, &
                            & ck - 1), wp))
                        av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                        qf(i)%sf(fi, fj, fk) = av
                        asum = asum + av
                    end do
                    qf(eqn_idx%adv%end)%sf(fi, fj, fk) = 1._wp - asum
                end do
            end do
        end do

    end subroutine s_prolong_alphas_closure

    !> Species mass-fraction prolongation closure (chemistry): each partial density rho*Y_k is minmod-prolonged and clamped
    !! non-negative, then all species are rescaled so sum_k(rho*Y_k) equals the (already prolonged) continuity density at the fine
    !! cell. This keeps the fine species realizable (Y_k >= 0, and sum(Y_k) = 1 exactly under the cons->prim recovery rho = sum
    !! rho*Y_k) and consistent with the continuity variable the reaction source reads. Same fine/coarse index mapping as
    !! s_prolong_one_var; cont is prolonged in the main loop before this runs.
    impure subroutine s_prolong_species_closure(qc, qf)

        type(scalar_field), dimension(sys_size), intent(in)    :: qc
        type(scalar_field), dimension(sys_size), intent(inout) :: qf
        integer                                                :: fi, fj, fk, ci, cj, ck, ox, oy, oz, i
        real(wp)                                               :: xix, xiy, xiz, u0, sx, sy, sz, av, rsum, rscale

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        do fk = 0, amr_slots(amr_cur)%p
            ck = amr_isect_lo(3) + fk/amr_slots(amr_cur)%ref_ratio - oz; if (p_glb == 0) ck = 0
            xiz = 0._wp; if (p_glb > 0) xiz = (real(mod(fk, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
            do fj = 0, amr_slots(amr_cur)%n
                cj = amr_isect_lo(2) + fj/amr_slots(amr_cur)%ref_ratio - oy; if (n_glb == 0) cj = 0
                xiy = 0._wp; if (n_glb > 0) xiy = (real(mod(fj, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                do fi = 0, amr_slots(amr_cur)%m
                    ci = amr_isect_lo(1) + fi/amr_slots(amr_cur)%ref_ratio - ox
                    xix = (real(mod(fi, amr_slots(amr_cur)%ref_ratio), wp) - 0.5_wp)*0.5_wp
                    rsum = 0._wp
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        u0 = real(qc(i)%sf(ci, cj, ck), wp)
                        sx = minmod(real(qc(i)%sf(ci + 1, cj, ck), wp) - u0, u0 - real(qc(i)%sf(ci - 1, cj, ck), wp))
                        sy = 0._wp
                        if (n_glb > 0) sy = minmod(real(qc(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(qc(i)%sf(ci, cj - 1, ck), wp))
                        sz = 0._wp
                        if (p_glb > 0) sz = minmod(real(qc(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(qc(i)%sf(ci, cj, ck - 1), wp))
                        av = max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp)
                        qf(i)%sf(fi, fj, fk) = av
                        rsum = rsum + av
                    end do
                    rscale = real(qf(eqn_idx%cont%end)%sf(fi, fj, fk), wp)/max(rsum, 1.e-30_wp)
                    do i = eqn_idx%species%beg, eqn_idx%species%end
                        qf(i)%sf(fi, fj, fk) = real(qf(i)%sf(fi, fj, fk), wp)*rscale
                    end do
                end do
            end do
        end do

    end subroutine s_prolong_species_closure

    !> Shared per-cell limiter switch for the volume-fraction closure prolongation: per dim, slopes stay on only if NO fluid's
    !! centered differences change sign there (symmetric in the fluids, including the closure fluid).
    pure subroutine s_alpha_shared_switch(qc, ci, cj, ck, shx, shy, shz)

        type(scalar_field), dimension(sys_size), intent(in) :: qc
        integer, intent(in)                                 :: ci, cj, ck
        logical, intent(out)                                :: shx, shy, shz
        integer                                             :: i
        real(wp)                                            :: u0

        shx = .true.; shy = n_glb > 0; shz = p_glb > 0
        do i = eqn_idx%adv%beg, eqn_idx%adv%end
            u0 = real(qc(i)%sf(ci, cj, ck), wp)
            if ((real(qc(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci - 1, cj, ck), wp)) <= 0._wp) shx = .false.
            if (n_glb > 0) then
                if ((real(qc(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj - 1, ck), wp)) <= 0._wp) shy = .false.
            end if
            if (p_glb > 0) then
                if ((real(qc(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(qc(i)%sf(ci, cj, ck - 1), wp)) <= 0._wp) shz = .false.
            end if
        end do

    end subroutine s_alpha_shared_switch

    !> Disblock prolongation. Guard: no-op unless amr.
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
        if (.not. amr_rank_owns_block) return
        call s_interpolate_coarse_to_fine(q_cons_base)
        ! prolonged fine state to the device (host prolongation; device consumers: ghost-fill/RK/RHS kernels)
        do i = 1, sys_size
            $:GPU_UPDATE(device='[amr_slots(amr_cur)%q_cons(i)%sf]')
        end do
        ! non-polytropic QBMM: seed the block's quadrature side-state from the coarse fields
        if (qbmm .and. .not. polytropic) call s_amr_prolong_pbmv_host()

    end subroutine s_populate_amr_fine

    !> Volume-weighted restriction for a single variable pair. Reads from qf (fine, must include interior 0:amr_slots(amr_cur)%m
    !! etc.); writes to qc (coarse, over the block).
    impure subroutine s_restrict_one_var(qf, qc)

        type(scalar_field), intent(in)    :: qf
        type(scalar_field), intent(inout) :: qc
        integer                           :: ci, cj, ck, fi0, fj0, fk0, ddj, ddk, nchild, ox, oy, oz
        real(wp)                          :: acc

        ! isect loop indices are GLOBAL (this rank's covered coarse cells); the coarse target qc is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        nchild = amr_slots(amr_cur)%ref_ratio
        if (n_glb > 0) nchild = nchild*amr_slots(amr_cur)%ref_ratio
        if (p_glb > 0) nchild = nchild*amr_slots(amr_cur)%ref_ratio
        do ck = amr_isect_lo(3), merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
            fk0 = (ck - amr_isect_lo(3))*amr_slots(amr_cur)%ref_ratio
            do cj = amr_isect_lo(2), merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
                fj0 = (cj - amr_isect_lo(2))*amr_slots(amr_cur)%ref_ratio
                do ci = amr_isect_lo(1), amr_isect_hi(1)
                    fi0 = (ci - amr_isect_lo(1))*amr_slots(amr_cur)%ref_ratio
                    acc = 0._wp
                    do ddk = 0, merge(amr_slots(amr_cur)%ref_ratio - 1, 0, p_glb > 0)
                        do ddj = 0, merge(amr_slots(amr_cur)%ref_ratio - 1, 0, n_glb > 0)
                            acc = acc + real(qf%sf(fi0, fj0 + ddj, fk0 + ddk), wp) + real(qf%sf(fi0 + 1, fj0 + ddj, fk0 + ddk), wp)
                        end do
                    end do
                    qc%sf(ci - ox, cj - oy, ck - oz) = acc/real(nchild, wp)
                end do
            end do
        end do

    end subroutine s_restrict_one_var

    !> Volume-weighted restriction: each covered coarse cell = average of its ref_ratio^d fine children. Writes the caller's coarse
    !! target - in production the level-0 state q_cons_ts(1)%vf (the deliberate fold-back of fine data each step, plus coarse pb/mv
    !! for non-polytropic QBMM); the init-time diagnostics pass a scratch buffer instead. Device kernel (per-cell arithmetic
    !! identical to s_restrict_one_var, which remains the host path for the diagnostics).
    impure subroutine s_restrict_fine_to_coarse(coarse_tgt)

        type(scalar_field), dimension(sys_size), intent(inout) :: coarse_tgt

        if (.not. amr_rank_owns_block) return
        if (rank_time_wrt) call s_rank_time_tic()
        call s_restrict_all_vars(amr_slots(amr_cur)%q_cons, coarse_tgt)
        ! non-polytropic QBMM: the block's side-state restricts alongside the moments so the coarse
        ! cells under the block carry the fine-resolved quadrature state into the next coarse step
        if (qbmm .and. .not. polytropic) call s_restrict_pbmv(pb_ts(1)%sf, mv_ts(1)%sf)
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_restrict_fine_to_coarse

    !> Apply phase-change relaxation (relax) to the current fine block's interior, BEFORE restriction. Relaxation is a cell-local,
    !! mass/energy-conserving equilibration (no stencil, no ghosts), so it needs no coarse/fine coupling: it just runs over the fine
    !! interior. Swaps m/n/p to the fine extents so s_infinite_relaxation_k's 0:m,0:n,0:p loop covers this block. Matches the coarse
    !! timing (once per full step; the coarse relax runs once after s_tvd_rk on q_cons_ts(1)) but on the fine solution so the fine
    !! dynamics equilibrate at fine resolution rather than only the restricted coarse average.
    impure subroutine s_amr_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_infinite_relaxation_k(amr_slots(amr_cur)%q_cons)
        call s_amr_restore_coarse()

    end subroutine s_amr_relax_fine

    !> 6-equation model: apply the per-stage pressure relaxation to the fine block's interior (cell-local equilibration, no
    !! stencil), mirroring the coarse per-stage call. Swaps the grid so the routine's 0:m,0:n,0:p loop covers this block.
    impure subroutine s_amr_pressure_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_pressure_relaxation_procedure(amr_slots(amr_cur)%q_cons)
        call s_amr_restore_coarse()

    end subroutine s_amr_pressure_relax_fine

    !> Compute the fine-grid IB state (markers/ghost points/levelset) for every active block from the body geometry (static-body
    !! AMR). Called once after the coarse IB setup at init (regrid+IB is gated). Per slot with fine cells: swap the grid to the fine
    !! block, swap the IB globals to the slot store, run the fine IB pipeline (writing into the slot store), restore. No-op unless
    !! amr .and. ib.
    impure subroutine s_amr_setup_ib()

        integer         :: islot, save_cur
        integer(kind=8) :: my_ib_gps, nrank_ib

        if (.not. amr .or. .not. ib) return
        save_cur = amr_cur
        my_ib_gps = 0_8
        do islot = 1, amr_num_blocks
            call s_amr_select_slot(islot)
            if (.not. amr_rank_owns_block) cycle
            call s_amr_swap_to_fine()
            call s_ibm_swap_to_fine(islot, gps_on_device=.false.)
            call s_ibm_setup_fine()
            my_ib_gps = my_ib_gps + int(num_gps, 8)
            call s_ibm_restore_from_fine(islot)
            call s_amr_restore_coarse()
        end do
        call s_amr_select_slot(save_cur)

        ! The fine-IB image-point stencil is not decomposition-exact across a rank seam.
        ! If the body's fine ghost points appear on more than one rank (i.e. the body
        ! straddles a coarse/fine rank boundary), abort rather than return a wrong
        ! body-surface state. A body wholly within one rank is decomposition-exact.
        call s_mpi_allreduce_integer_sum(merge(1_8, 0_8, my_ib_gps > 0_8), nrank_ib)
        if (nrank_ib > 1_8) then
            call s_mpi_abort('amr with ib: the immersed body straddles a rank boundary, where the ' &
                             & // 'fine-IB image-point stencil is not yet decomposition-exact; keep the ' &
                             & // 'body within a single rank subdomain (use fewer ranks or reposition it).')
        end if

    end subroutine s_amr_setup_ib

    !> Apply the IB state correction on the current fine block after its RK update (static-body AMR). Mirrors the coarse per-stage
    !! s_ibm_correct_state: swap the grid + IB globals to the fine block, correct q_cons/q_prim at the fine body/ghost cells,
    !! restore. amr_cur / amr_rank_owns_block are set by the caller (the per-block advance loop). No-op unless ib.
    impure subroutine s_amr_ib_correct_fine()

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        if (qbmm .and. .not. polytropic) then
            ! mirror the coarse correct-state: non-polytropic QBMM also corrects the block's own
            ! pb/mv side-state at the body ghost points (bounds match the swapped fine idwbuff)
            call s_ibm_correct_state(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_prim, amr_slots(amr_cur)%pb_f, &
                                     & amr_slots(amr_cur)%mv_f)
        else
            call s_ibm_correct_state(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_prim)
        end if
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_ib_correct_fine

    !> Rebuild the current fine block's IB state (markers/ghost points/image points) from the moving body's position (prescribed
    !! motion, moving_ibm==1). Reuses the coarse s_update_mib recompute on the swapped-in fine slot (grid + IB globals swapped to
    !! the fine block, recompute writes into the slot store, restore). For the subcycled advance pass th in [0,1], the fine
    !! substep's fraction of the coarse step: s_update_mib snapshots the body to the linear time interpolation between the coarse
    !! t^n and t^{n+1} positions, the same time-interpolation the subcycle applies to the fluid ghost shell. Pass th < 0 for the
    !! non-subcycled lockstep stage (uses the body's current position). No-op unless ib. Must precede s_amr_ib_correct_fine.
    impure subroutine s_amr_update_mib_fine(th)

        real(wp), intent(in) :: th
        integer              :: i, blo(3), bhi(3)
        logical              :: ovl, inside

        if (.not. ib) return
        if (.not. amr_rank_owns_block) return
        ! Between regrids a moving body must stay inside its block: the regrid expansion
        ! contained it with margin max(amr_buf,4), and here we require the body plus the
        ! image-point stencil reach (2 coarse cells) to remain contained. Consecutive
        ! contained positions keep every sub-time interpolate contained (axis-aligned
        ! boxes are convex in the linearly-moving corners). A body that overlaps the
        ! block edge would get silently clipped forcing - abort instead.
        if (amr_regrid_int > 0) then
            do i = 1, num_ibs
                if (patch_ib(i)%moving_ibm == 0) cycle
                call s_amr_body_bbox(i, 2, blo, bhi)
                ovl = blo(1) <= amr_slots(amr_cur)%region%hi(1) .and. bhi(1) >= amr_slots(amr_cur)%region%lo(1)
                if (n_glb > 0) ovl = ovl .and. blo(2) <= amr_slots(amr_cur)%region%hi(2) .and. bhi(2) &
                    & >= amr_slots(amr_cur)%region%lo(2)
                if (p_glb > 0) ovl = ovl .and. blo(3) <= amr_slots(amr_cur)%region%hi(3) .and. bhi(3) &
                    & >= amr_slots(amr_cur)%region%lo(3)
                if (.not. ovl) cycle
                inside = blo(1) >= amr_slots(amr_cur)%region%lo(1) .and. bhi(1) <= amr_slots(amr_cur)%region%hi(1)
                if (n_glb > 0) inside = inside .and. blo(2) >= amr_slots(amr_cur)%region%lo(2) .and. bhi(2) &
                    & <= amr_slots(amr_cur)%region%hi(2)
                if (p_glb > 0) inside = inside .and. blo(3) >= amr_slots(amr_cur)%region%lo(3) .and. bhi(3) &
                    & <= amr_slots(amr_cur)%region%hi(3)
                if (.not. inside) then
                    call s_mpi_abort('amr dynamic regrid with moving ib: the body reached the fine-block ' &
                                     & // 'boundary between regrids; reduce amr_regrid_int or increase amr_buf')
                end if
            end do
        end if
        call s_amr_swap_to_fine()
        call s_ibm_swap_to_fine(amr_cur, gps_on_device=.true.)
        call s_update_mib(num_ibs, th)
        call s_ibm_restore_from_fine(amr_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_update_mib_fine

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
        rr = amr_slots(amr_cur)%ref_ratio
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

    !> Non-polytropic QBMM: piecewise-constant HOST prolongation of the coarse pb/mv side-state onto the block's fine interior, then
    !! push to the device. Piecewise-constant preserves the per-cell realizable quadrature set (mirroring the moments' injection).
    !! Called at init populate and at restart (the fine side-state is re-prolonged from the restored coarse fields - a one-time
    !! smoothing; q_cons is restored exactly).
    impure subroutine s_amr_prolong_pbmv_host()

        integer :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, ox, oy, oz

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        rr = amr_slots(amr_cur)%ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, amr_slots(amr_cur)%p
                    ck = 0
                    if (p_glb > 0) ck = lo3 + fk/rr - oz
                    do fj = 0, amr_slots(amr_cur)%n
                        cj = 0
                        if (n_glb > 0) cj = lo2 + fj/rr - oy
                        do fi = 0, amr_slots(amr_cur)%m
                            ci = lo1 + fi/rr - ox
                            amr_slots(amr_cur)%pb_f(fi, fj, fk, q, ib_) = pb_ts(1)%sf(ci, cj, ck, q, ib_)
                            amr_slots(amr_cur)%mv_f(fi, fj, fk, q, ib_) = mv_ts(1)%sf(ci, cj, ck, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:GPU_UPDATE(device='[amr_slots(amr_cur)%pb_f, amr_slots(amr_cur)%mv_f]')

    end subroutine s_amr_prolong_pbmv_host

    !> Non-polytropic QBMM: piecewise-constant prolongation of the fine pb/mv GHOST shell from the coarse side-state (device kernel;
    !! interior untouched). The ghosts feed the widened-idwint conversions and the qbmm rhs over the shell, mirroring the q_cons
    !! ghost fill.
    impure subroutine s_amr_fill_fine_ghosts_pbmv(pb_c, mv_c, pb_tgt, mv_tgt)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_c, mv_c

        real(stp), dimension(amr_slots(amr_cur)%idwbuff(1)%beg:,amr_slots(amr_cur)%idwbuff(2)%beg:, &
             & amr_slots(amr_cur)%idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_tgt, mv_tgt
        integer :: fi, fj, fk, q, ib_, ci, cj, ck, rr, lo1, lo2, lo3, fm, fn, fp, b1, e1, b2, e2, b3, e3, ox, oy, oz
        logical :: d2, d3

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5, private='[ci, cj, ck]')
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = b3, e3
                    do fj = b2, e2
                        do fi = b1, e1
                            if (.not. (fi >= 0 .and. fi <= fm .and. fj >= 0 .and. fj <= fn .and. fk >= 0 .and. fk <= fp)) then
                                ck = 0
                                if (d3) ck = lo3 + floor(real(fk, wp)/real(rr, wp)) - oz
                                cj = 0
                                if (d2) cj = lo2 + floor(real(fj, wp)/real(rr, wp)) - oy
                                ci = lo1 + floor(real(fi, wp)/real(rr, wp)) - ox
                                pb_tgt(fi, fj, fk, q, ib_) = pb_c(ci, cj, ck, q, ib_)
                                mv_tgt(fi, fj, fk, q, ib_) = mv_c(ci, cj, ck, q, ib_)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fill_fine_ghosts_pbmv

    !> Non-polytropic QBMM: device copy of the block's pb/mv into the step-entry backup (SSP-RK).
    impure subroutine s_amr_backup_pbmv()

        integer :: fi, fj, fk, q, ib_, b1, e1, b2, e2, b3, e3

        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = b3, e3
                    do fj = b2, e2
                        do fi = b1, e1
                            amr_slots(amr_cur)%pb_stor(fi, fj, fk, q, ib_) = amr_slots(amr_cur)%pb_f(fi, fj, fk, q, ib_)
                            amr_slots(amr_cur)%mv_stor(fi, fj, fk, q, ib_) = amr_slots(amr_cur)%mv_f(fi, fj, fk, q, ib_)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_backup_pbmv

    !> Non-polytropic QBMM: SSP-RK stage update of the block's pb/mv (device kernel, interior only; mirror of the coarse pb_ts/mv_ts
    !! stage combination in m_time_steppers).
    impure subroutine s_amr_fine_rk_update_pbmv(c1, c2, c3, c4, dtl)

        real(wp), intent(in) :: c1, c2, c3, c4, dtl
        integer              :: fi, fj, fk, q, ib_, fm, fn, fp

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = 0, fp
                    do fj = 0, fn
                        do fi = 0, fm
                            amr_slots(amr_cur)%pb_f(fi, fj, fk, q, ib_) = (c1*amr_slots(amr_cur)%pb_f(fi, fj, fk, q, &
                                      & ib_) + c2*amr_slots(amr_cur)%pb_stor(fi, fj, fk, q, &
                                      & ib_) + c3*dtl*amr_slots(amr_cur)%rhs_pb_f(fi, fj, fk, q, ib_))/c4
                            amr_slots(amr_cur)%mv_f(fi, fj, fk, q, ib_) = (c1*amr_slots(amr_cur)%mv_f(fi, fj, fk, q, &
                                      & ib_) + c2*amr_slots(amr_cur)%mv_stor(fi, fj, fk, q, &
                                      & ib_) + c3*dtl*amr_slots(amr_cur)%rhs_mv_f(fi, fj, fk, q, ib_))/c4
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_fine_rk_update_pbmv

    !> Non-polytropic QBMM: volume-weighted restriction of the block's pb/mv onto the coarse side-state under the block (device
    !! kernel; mirror of s_restrict_all_vars' child average).
    impure subroutine s_restrict_pbmv(pb_c, mv_c)

        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_c, mv_c
        integer :: ci, cj, ck, q, ib_, fi0, fj0, fk0, ddj, ddk, nchild, ox, oy, oz, rr
        integer :: c1lo, c1hi, c2lo, c2hi, c3lo, c3hi, dj_hi, dk_hi
        real(wp) :: accp, accm

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        rr = amr_slots(amr_cur)%ref_ratio
        nchild = rr
        if (n_glb > 0) nchild = nchild*rr
        if (p_glb > 0) nchild = nchild*rr
        c1lo = amr_isect_lo(1); c1hi = amr_isect_hi(1)
        c2lo = amr_isect_lo(2); c2hi = merge(amr_isect_hi(2), amr_isect_lo(2), n_glb > 0)
        c3lo = amr_isect_lo(3); c3hi = merge(amr_isect_hi(3), amr_isect_lo(3), p_glb > 0)
        dj_hi = merge(rr - 1, 0, n_glb > 0); dk_hi = merge(rr - 1, 0, p_glb > 0)
        $:GPU_PARALLEL_LOOP(collapse=5, private='[fi0, fj0, fk0, accp, accm, ddj, ddk]')
        do ib_ = 1, nb
            do q = 1, nnode
                do ck = c3lo, c3hi
                    do cj = c2lo, c2hi
                        do ci = c1lo, c1hi
                            fi0 = (ci - c1lo)*rr; fj0 = (cj - c2lo)*rr; fk0 = (ck - c3lo)*rr
                            accp = 0._wp; accm = 0._wp
                            do ddk = 0, dk_hi
                                do ddj = 0, dj_hi
                                    accp = accp + real(amr_slots(amr_cur)%pb_f(fi0, fj0 + ddj, fk0 + ddk, q, ib_), &
                                                       & wp) + real(amr_slots(amr_cur)%pb_f(fi0 + 1, fj0 + ddj, fk0 + ddk, q, &
                                                       & ib_), wp)
                                    accm = accm + real(amr_slots(amr_cur)%mv_f(fi0, fj0 + ddj, fk0 + ddk, q, ib_), &
                                                       & wp) + real(amr_slots(amr_cur)%mv_f(fi0 + 1, fj0 + ddj, fk0 + ddk, q, &
                                                       & ib_), wp)
                                end do
                            end do
                            pb_c(ci - ox, cj - oy, ck - oz, q, ib_) = accp/real(nchild, wp)
                            mv_c(ci - ox, cj - oy, ck - oz, q, ib_) = accm/real(nchild, wp)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_restrict_pbmv

    !> SP2 gate: restrict(prolong(coarse)) must reproduce coarse over the block interior (conservation). Init-only diagnostic;
    !! allocates a scratch coarse target, never touches level-0 or the solve.
    impure subroutine s_amr_conservation_check(q_cons_base)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
        type(scalar_field), dimension(:), allocatable       :: scratch
        integer                                             :: i, ci, cj, ck, ox, oy, oz
        real(wp)                                            :: err, e

        if (.not. amr) return
        if (.not. amr_rank_owns_block) return
        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        allocate (scratch(1:sys_size))
        do i = 1, sys_size
            allocate (scratch(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do
        ! host restriction path: the scratch target is host-only and the fine state is host-current at init
        do i = 1, sys_size
            call s_restrict_one_var(amr_slots(amr_cur)%q_cons(i), scratch(i))
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

    !> Swap the global grid state to the fine block. MUST be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine()

        ! a nested swap would overwrite the sw_* bounce buffers with FINE state, and the eventual
        ! restore would install fine extents as the coarse grid - silent corruption of everything after
        @:ASSERT(.not. amr_swapped, "nested s_amr_swap_to_fine (swap/restore must pair)")
        amr_swapped = .true.
        sw_m = m; sw_n = n; sw_p = p
        sw_idwint = idwint; sw_idwbuff = idwbuff
        ! the acoustic source's precomputed spatials are coarse-grid cell indices: applying them on
        ! the fine block would inject at wrong cells (or out of bounds). The support is guaranteed
        ! not to overlap the block (checked at startup), so the fine RHS correctly skips the source.
        sw_acoustic_source = acoustic_source; acoustic_source = .false.
        m = amr_slots(amr_cur)%m; n = amr_slots(amr_cur)%n; p = amr_slots(amr_cur)%p
        idwint(1)%beg = 0; idwint(1)%end = m
        idwint(2)%beg = 0; idwint(2)%end = n
        idwint(3)%beg = 0; idwint(3)%end = p
        idwbuff = amr_slots(amr_cur)%idwbuff
        ! save coarse coords to bounce buffers, then copy fine coords into global arrays
        sw_x_cb = x_cb; sw_x_cc = x_cc; sw_dx = dx
        if (n_glb > 0) then; sw_y_cb = y_cb; sw_y_cc = y_cc; sw_dy = dy; end if
        if (p_glb > 0) then; sw_z_cb = z_cb; sw_z_cc = z_cc; sw_dz = dz; end if
        x_cb(-1:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cb(-1:amr_slots(amr_cur)%m)
        x_cc(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%x_cc(0:amr_slots(amr_cur)%m)
        dx(0:amr_slots(amr_cur)%m) = amr_slots(amr_cur)%dx(0:amr_slots(amr_cur)%m)
        if (n_glb > 0) then
            y_cb(-1:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cb(-1:amr_slots(amr_cur)%n)
            y_cc(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%y_cc(0:amr_slots(amr_cur)%n)
            dy(0:amr_slots(amr_cur)%n) = amr_slots(amr_cur)%dy(0:amr_slots(amr_cur)%n)
        end if
        if (p_glb > 0) then
            z_cb(-1:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cb(-1:amr_slots(amr_cur)%p)
            z_cc(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%z_cc(0:amr_slots(amr_cur)%p)
            dz(0:amr_slots(amr_cur)%p) = amr_slots(amr_cur)%dz(0:amr_slots(amr_cur)%p)
        end if
        ! Extend the fine grid into the ghost shell (s_build_level_coords only fills the interior 0:m).
        ! Ghost cells use the EXACT parent-cell bisection - the same formula as the interior, with
        ! floor division for negative indices. The parent boundaries beyond the block edge come from
        ! the coarse coordinates saved in the sw_* bounce buffers moments ago; blocks stay buff_size
        ! inside the domain, so every ghost parent is an interior (or rank-seam-exchanged) coarse
        ! cell with exact coordinates. This makes the ghost shell stretched-grid-exact (the old
        ! locally-uniform continuation was only exact on uniform bases).
        block
            integer :: jg, cl
            do jg = amr_slots(amr_cur)%m + 1, amr_slots(amr_cur)%m + buff_size
                cl = amr_isect_lo(1) + floor(real(jg, wp)/2._wp) - start_idx(1)
                if (mod(jg, 2) == 0) then
                    x_cb(jg) = 0.5_wp*(sw_x_cb(cl - 1) + sw_x_cb(cl))
                else
                    x_cb(jg) = sw_x_cb(cl)
                end if
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            ! unified boundary formula (matches the interior bisection): boundary k belongs to
            ! parent c = isect_lo + floor(k/2); even k -> parent midpoint, odd k -> parent right edge
            do jg = -1 - buff_size, -1
                cl = amr_isect_lo(1) + floor(real(jg, wp)/2._wp) - start_idx(1)
                if (mod(abs(jg), 2) == 0) then
                    x_cb(jg) = 0.5_wp*(sw_x_cb(cl - 1) + sw_x_cb(cl))
                else
                    x_cb(jg) = sw_x_cb(cl)
                end if
            end do
            do jg = -buff_size, -1
                dx(jg) = x_cb(jg) - x_cb(jg - 1); x_cc(jg) = 0.5_wp*(x_cb(jg - 1) + x_cb(jg))
            end do
            if (n_glb > 0) then
                do jg = amr_slots(amr_cur)%n + 1, amr_slots(amr_cur)%n + buff_size
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/2._wp) - start_idx(2)
                    if (mod(jg, 2) == 0) then
                        y_cb(jg) = 0.5_wp*(sw_y_cb(cl - 1) + sw_y_cb(cl))
                    else
                        y_cb(jg) = sw_y_cb(cl)
                    end if
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
                ! unified boundary formula (matches the interior bisection): boundary k belongs to
                ! parent c = isect_lo + floor(k/2); even k -> parent midpoint, odd k -> parent right edge
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(2) + floor(real(jg, wp)/2._wp) - start_idx(2)
                    if (mod(abs(jg), 2) == 0) then
                        y_cb(jg) = 0.5_wp*(sw_y_cb(cl - 1) + sw_y_cb(cl))
                    else
                        y_cb(jg) = sw_y_cb(cl)
                    end if
                end do
                do jg = -buff_size, -1
                    dy(jg) = y_cb(jg) - y_cb(jg - 1); y_cc(jg) = 0.5_wp*(y_cb(jg - 1) + y_cb(jg))
                end do
            end if
            if (p_glb > 0) then
                do jg = amr_slots(amr_cur)%p + 1, amr_slots(amr_cur)%p + buff_size
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/2._wp) - start_idx(3)
                    if (mod(jg, 2) == 0) then
                        z_cb(jg) = 0.5_wp*(sw_z_cb(cl - 1) + sw_z_cb(cl))
                    else
                        z_cb(jg) = sw_z_cb(cl)
                    end if
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
                ! unified boundary formula (matches the interior bisection): boundary k belongs to
                ! parent c = isect_lo + floor(k/2); even k -> parent midpoint, odd k -> parent right edge
                do jg = -1 - buff_size, -1
                    cl = amr_isect_lo(3) + floor(real(jg, wp)/2._wp) - start_idx(3)
                    if (mod(abs(jg), 2) == 0) then
                        z_cb(jg) = 0.5_wp*(sw_z_cb(cl - 1) + sw_z_cb(cl))
                    else
                        z_cb(jg) = sw_z_cb(cl)
                    end if
                end do
                do jg = -buff_size, -1
                    dz(jg) = z_cb(jg) - z_cb(jg - 1); z_cc(jg) = 0.5_wp*(z_cb(jg - 1) + z_cb(jg))
                end do
            end if
        end block
        ! sync the swapped extents/bounds/coordinates to the device: RHS kernels read the
        ! device copies of these GPU_DECLARE'd globals (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()
        ! hypoelastic stress sources use grid-spacing-dependent FD coefficients: recompute
        ! them from the (now fine) grid, else every fine velocity gradient is halved
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        ! nonuniform coarse grid (axisymmetric axis half-cell): the per-cell WENO coefficients
        ! must be rebuilt for the block's own uniform grid (no-op flag on uniform grids)
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()

    end subroutine s_amr_swap_to_fine

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse()

        @:ASSERT(amr_swapped, "s_amr_restore_coarse without a matching s_amr_swap_to_fine")
        amr_swapped = .false.
        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        acoustic_source = sw_acoustic_source
        ! restore full coarse coords from bounce buffers
        x_cb = sw_x_cb; x_cc = sw_x_cc; dx = sw_dx
        if (n_glb > 0) then; y_cb = sw_y_cb; y_cc = sw_y_cc; dy = sw_dy; end if
        if (p_glb > 0) then; z_cb = sw_z_cb; z_cc = sw_z_cc; dz = sw_dz; end if
        ! sync the restored coarse extents/bounds/coordinates back to the device
        call s_amr_sync_grid_state_to_device()
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        if (amr_weno_coef_recompute) call s_amr_recompute_weno_coefs()

    end subroutine s_amr_restore_coarse

    !> Recompute the WENO reconstruction coefficient arrays from the CURRENT grid globals (the fine block's after a swap, the coarse
    !! grid's after a restore). s_compute_weno_coefficients reads the live cell-boundary arrays, refreshes uniform_grid, and pushes
    !! its own device updates; the coefficient arrays were sized to the coarse local ranges at init, which the fine ranges never
    !! exceed (the scratch-fit abort guarantees it).
    impure subroutine s_amr_recompute_weno_coefs()

        type(int_bounds_info) :: is1, is2, is3

        is1%beg = -buff_size; is1%end = m + buff_size
        call s_compute_weno_coefficients(1, is1)
        if (n_glb > 0) then
            is2%beg = -buff_size; is2%end = n + buff_size
            call s_compute_weno_coefficients(2, is2)
        end if
        if (p_glb > 0) then
            is3%beg = -buff_size; is3%end = p + buff_size
            call s_compute_weno_coefficients(3, is3)
        end if

    end subroutine s_amr_recompute_weno_coefs

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
    !! untouched. Multi-fluid volume fractions get the same sum-preserving closure as the interior prolongation (second kernel).
    impure subroutine s_amr_fill_fine_ghosts(q_coarse, q_fine)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_coarse
        type(scalar_field), dimension(sys_size), intent(inout) :: q_fine
        integer                                                :: i, fi, fj, fk, ci, cj, ck, ox, oy, oz
        integer                                                :: rr, lo1, lo2, lo3, fm, fn, fp, b1, e1, b2, e2, b3, e3
        integer                                                :: advb, adve, bbeg, bend, bstride
        logical                                                :: d2, d3, multi, shx, shy, shz, bubEE
        real(wp)                                               :: u0, sx, sy, sz, xix, xiy, xiz, av, asum

        ! fine indices are LOCAL to this rank's intersection; amr_isect_lo is GLOBAL; the coarse source q_coarse is rank-LOCAL

        ox = start_idx(1); oy = 0; oz = 0
        if (n_glb > 0) oy = start_idx(2)
        if (p_glb > 0) oz = start_idx(3)
        d2 = n_glb > 0; d3 = p_glb > 0
        rr = amr_slots(amr_cur)%ref_ratio
        lo1 = amr_isect_lo(1); lo2 = amr_isect_lo(2); lo3 = amr_isect_lo(3)
        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        multi = num_fluids > 1
        advb = eqn_idx%adv%beg; adve = eqn_idx%adv%end
        bubEE = bubbles_euler; bbeg = eqn_idx%bub%beg; bend = eqn_idx%bub%end
        bstride = 1; if (bubEE) bstride = (bend - bbeg + 1)/nb
        $:GPU_PARALLEL_LOOP(collapse=4, private='[ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz]')
        do i = 1, sys_size
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
                        ! skip the interior (only the ghost shell is filled) and, multi-fluid, the volume
                        ! fractions (closure kernel below)
                        if (.not. (fi >= 0 .and. fi <= fm .and. fj >= 0 .and. fj <= fn .and. fk >= 0 .and. fk <= fp) &
                            & .and. .not. (multi .and. i >= advb .and. i <= adve)) then
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
                            ! QBMM: inject the bub block piecewise-constant (child = u0) so the ghost inherits the
                            ! coarse cell's realizable 6-moment set (CHyQMOM needs variance c20 > 0; per-component
                            ! minmod slopes would break that joint constraint). Non-QBMM Euler-Euler bubbles instead
                            ! floor their positive moments (nR / npb / nmv); the signed velocity moment nV (offset 1)
                            ! is skipped.
                            if (qbmm .and. i >= bbeg .and. i <= bend) then
                                sx = 0._wp; sy = 0._wp; sz = 0._wp
                            end if
                            q_fine(i)%sf(fi, fj, fk) = u0 + sx*xix + sy*xiy + sz*xiz
                            if (bubEE .and. .not. qbmm .and. i >= bbeg .and. i <= bend) then
                                if (mod(i - bbeg, bstride) /= 1) q_fine(i)%sf(fi, fj, fk) = max(real(q_fine(i)%sf(fi, fj, fk), &
                                    & wp), bub_pos_frac*u0)
                            end if
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! multi-fluid volume-fraction ghosts: per-cell closure mirroring s_prolong_alphas_closure (shared
        ! limiter switch over all fluids; interpolate + clamp fluids advb..adve-1; alpha_n = 1 - sum)
        if (multi) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, ci, cj, ck, xix, xiy, xiz, u0, sx, sy, sz, av, asum, shx, shy, shz]')
            do fk = b3, e3
                do fj = b2, e2
                    do fi = b1, e1
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
                            shx = .true.; shy = d2; shz = d3
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advb, adve
                                u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                                if ((real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), &
                                    & wp)) <= 0._wp) shx = .false.
                                if (d2) then
                                    if ((real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj - 1, &
                                        & ck), wp)) <= 0._wp) shy = .false.
                                end if
                                if (d3) then
                                    if ((real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0)*(u0 - real(q_coarse(i)%sf(ci, cj, &
                                        & ck - 1), wp)) <= 0._wp) shz = .false.
                                end if
                            end do
                            asum = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advb, adve - 1
                                u0 = real(q_coarse(i)%sf(ci, cj, ck), wp)
                                sx = 0._wp
                                if (shx) sx = minmod(real(q_coarse(i)%sf(ci + 1, cj, ck), wp) - u0, &
                                    & u0 - real(q_coarse(i)%sf(ci - 1, cj, ck), wp))
                                sy = 0._wp
                                if (shy) sy = minmod(real(q_coarse(i)%sf(ci, cj + 1, ck), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj - 1, ck), wp))
                                sz = 0._wp
                                if (shz) sz = minmod(real(q_coarse(i)%sf(ci, cj, ck + 1), wp) - u0, u0 - real(q_coarse(i)%sf(ci, &
                                    & cj, ck - 1), wp))
                                av = min(max(u0 + sx*xix + sy*xiy + sz*xiz, 0._wp), 1._wp)
                                q_fine(i)%sf(fi, fj, fk) = av
                                asum = asum + av
                            end do
                            q_fine(adve)%sf(fi, fj, fk) = 1._wp - asum
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_fill_fine_ghosts

    !> Lerp the fine ghost shell of q_tgt between q_a (coarse t^n) and q_b (coarse t^{n+1}) at time fraction th (device kernel).
    !! Interior untouched.
    impure subroutine s_amr_lerp_fine_ghosts(q_a, q_b, q_tgt, th)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_a, q_b
        type(scalar_field), dimension(sys_size), intent(inout) :: q_tgt
        real(wp), intent(in)                                   :: th
        integer                                                :: i, fi, fj, fk, fm, fn, fp, b1, e1, b2, e2, b3, e3

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
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

    !> Non-polytropic QBMM twin of s_amr_lerp_fine_ghosts: lerp the block's pb/mv ghost shell between the coarse t^n and t^{n+1}
    !! sources at the substage time (device kernel; interior untouched). Ghost pb feeds the mixture pressure in the widened
    !! conversion, so it gets the same time treatment as the conservative ghosts.
    impure subroutine s_amr_lerp_fine_ghosts_pbmv(th)

        real(wp), intent(in) :: th
        integer              :: fi, fj, fk, q, ib_, fm, fn, fp, b1, e1, b2, e2, b3, e3

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
        b1 = amr_slots(amr_cur)%idwbuff(1)%beg; e1 = amr_slots(amr_cur)%idwbuff(1)%end
        b2 = amr_slots(amr_cur)%idwbuff(2)%beg; e2 = amr_slots(amr_cur)%idwbuff(2)%end
        b3 = amr_slots(amr_cur)%idwbuff(3)%beg; e3 = amr_slots(amr_cur)%idwbuff(3)%end
        $:GPU_PARALLEL_LOOP(collapse=5)
        do ib_ = 1, nb
            do q = 1, nnode
                do fk = b3, e3
                    do fj = b2, e2
                        do fi = b1, e1
                            if (.not. (fi >= 0 .and. fi <= fm .and. fj >= 0 .and. fj <= fn .and. fk >= 0 .and. fk <= fp)) then
                                amr_slots(amr_cur)%pb_f(fi, fj, fk, q, ib_) = (1._wp - th)*real(amr_slots(amr_cur)%pb_ghost_a(fi, &
                                          & fj, fk, q, ib_), wp) + th*real(amr_slots(amr_cur)%pb_ghost_b(fi, fj, fk, q, ib_), wp)
                                amr_slots(amr_cur)%mv_f(fi, fj, fk, q, ib_) = (1._wp - th)*real(amr_slots(amr_cur)%mv_ghost_a(fi, &
                                          & fj, fk, q, ib_), wp) + th*real(amr_slots(amr_cur)%mv_ghost_b(fi, fj, fk, q, ib_), wp)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_amr_lerp_fine_ghosts_pbmv

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

        fm = amr_slots(amr_cur)%m; fn = amr_slots(amr_cur)%n; fp = amr_slots(amr_cur)%p
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
    !! beg then end, mirroring s_populate_variables_buffers' disblock). The solver never fills CONS ghosts (only prim), so ranks
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
        if (.not. amr_rank_owns_block) return

        ! ghost prolongation from the coarse stage-entry conservative state (device kernel reads the
        ! device-current coarse stage-entry state directly); rank_time brackets cover the fine-advance
        ! compute segments and pause across the MPI exchanges (the inner s_compute_rhs pair nests to a no-op)
        if (rank_time_wrt) call s_rank_time_tic()
        call s_amr_fill_fine_ghosts(q_cons_coarse, amr_slots(amr_cur)%q_cons)
        if (rank_time_wrt) call s_rank_time_toc()

        ! continuation faces (the block spans a rank boundary there): overwrite the prolonged ghosts with the
        ! neighbor's true fine data (no exchange fires at np=1 / for fully-contained blocks)
        call s_mpi_sendrecv_amr_fine_halo(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%m, amr_slots(amr_cur)%n, &
                                          & amr_slots(amr_cur)%p)
        if (rank_time_wrt) call s_rank_time_tic()

        ! non-polytropic QBMM: prolong the block's pb/mv ghost shell from the coarse stage-entry
        ! side-state (its rhs is cell-local, so ghosts only feed the widened-idwint conversions)
        if (qbmm .and. .not. polytropic) call s_amr_fill_fine_ghosts_pbmv(pb_in, mv_in, amr_slots(amr_cur)%pb_f, &
            & amr_slots(amr_cur)%mv_f)

        ! step-entry backup for the SSP-RK combination (device copy over the current buffered extents)
        if (s == 1) then
            call s_amr_copy_fine_fields(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, &
                                        & amr_slots(amr_cur)%idwbuff(1)%beg, amr_slots(amr_cur)%idwbuff(1)%end, &
                                        & amr_slots(amr_cur)%idwbuff(2)%beg, amr_slots(amr_cur)%idwbuff(2)%end, &
                                        & amr_slots(amr_cur)%idwbuff(3)%beg, amr_slots(amr_cur)%idwbuff(3)%end)
            if (qbmm .and. .not. polytropic) call s_amr_backup_pbmv()
        end if

        amr_in_fine_advance = .true.
        call s_amr_swap_to_fine()
        idwint = amr_slots(amr_cur)%idwbuff  ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
        $:GPU_UPDATE(device='[idwint]')
        if (qbmm .and. .not. polytropic) then
            ! the block's OWN side-state and rhs scratch: the coarse pb_in/rhs_pb must not be
            ! touched at fine indices (the coarse stage consumes them after this fine stage)
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & amr_slots(amr_cur)%pb_f, amr_slots(amr_cur)%rhs_pb_f, amr_slots(amr_cur)%mv_f, &
                               & amr_slots(amr_cur)%rhs_mv_f, t_step, time_avg, s)
        else
            call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, amr_slots(amr_cur)%rhs, &
                               & pb_in, rhs_pb, mv_in, rhs_mv, t_step, time_avg, s)
        end if
        call s_amr_restore_coarse()
        amr_in_fine_advance = .false.

        ! RK stage update (device kernel; mirror of the coarse non-IGR form)
        call s_amr_fine_rk_update(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, amr_slots(amr_cur)%rhs, coefs(1), &
                                  & coefs(2), coefs(3), coefs(4), dt)
        if (qbmm .and. .not. polytropic) call s_amr_fine_rk_update_pbmv(coefs(1), coefs(2), coefs(3), coefs(4), dt)
        ! 6-equation model: per-stage pressure relaxation on the block (before IB correct, coarse order)
        if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
        ! moving body: rebuild the fine-block IB state at the current (lockstep-stage) body position before the correct-state
        if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(-1._wp)
        ! IB state correction on the fine block (mirrors the coarse per-stage correct-state; no-op unless ib)
        call s_amr_ib_correct_fine()
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_advance_amr_fine_stage

    !> Subcycled fine advance (amr_subcycle): two dt/2 SSP-RK3 substeps AFTER the coarse step. q_old/q_new are the coarse t^n and
    !! t^{n+1} states; each stage's ghosts are the linear time interpolation at the stage time theta = (substep-1 + c_s)/2 with
    !! SSP-RK3 abscissae c = [0, 1, 1/2]. Fine flux registers are zeroed here and accumulate over all six stages (0.5*rk3_w each) so
    !! the end-of-step state reflux sees the time-averaged effective fine flux.
    impure subroutine s_advance_amr_fine_substeps(q_old, q_new, coefs, bc_type, q_T_sf, pb_old, mv_old, pb_in, rhs_pb, mv_in, &
        & rhs_mv, t_step, time_avg)

        type(scalar_field), dimension(sys_size), intent(inout)                                  :: q_old, q_new
        real(wp), dimension(:,:), intent(in)                                                    :: coefs  !< rk_coef(1:3, 1:4)
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                              :: bc_type
        type(scalar_field), intent(inout)                                                       :: q_T_sf
        real(stp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(in) :: pb_old, mv_old
        real(stp), dimension(:,:,:,:,:), intent(inout)                                          :: pb_in, mv_in
        real(wp), dimension(:,:,:,:,:), intent(inout)                                           :: rhs_pb, rhs_mv
        integer, intent(in)                                                                     :: t_step
        real(wp), intent(inout)                                                                 :: time_avg
        real(wp), parameter                                                                     :: c_abs(3) = [0._wp, 1._wp, 0.5_wp]
        integer                                                                                 :: sub, s
        real(wp)                                                                                :: th

        if (.not. amr) return
        ! valid coarse CONS ghosts on both lerp sources (ALL ranks call: pairwise halo); the exchanged t^n /
        ! t^{n+1} ghost layers make the prolonged block-boundary ghosts correct even at rank boundaries
        if (amr_xchg_coarse_ghosts) then
            call s_amr_exchange_coarse_cons_halo(q_old)
            call s_amr_exchange_coarse_cons_halo(q_new)
        end if
        if (.not. amr_rank_owns_block) return

        ! fill both lerp sources once: ghost shells prolonged from coarse t^n and t^{n+1} (device kernels
        ! read the device-current coarse states directly); rank_time brackets cover the fine-advance
        ! compute segments and pause across the MPI exchanges (the inner s_compute_rhs pair nests to a no-op)
        if (rank_time_wrt) call s_rank_time_tic()
        call s_amr_fill_fine_ghosts(q_old, amr_slots(amr_cur)%q_ghost_a)
        call s_amr_fill_fine_ghosts(q_new, amr_slots(amr_cur)%q_ghost_b)
        ! non-polytropic QBMM: the pb/mv ghost shell gets the same two-source time-lerp treatment
        if (qbmm .and. .not. polytropic) then
            call s_amr_fill_fine_ghosts_pbmv(pb_old, mv_old, amr_slots(amr_cur)%pb_ghost_a, amr_slots(amr_cur)%mv_ghost_a)
            call s_amr_fill_fine_ghosts_pbmv(pb_in, mv_in, amr_slots(amr_cur)%pb_ghost_b, amr_slots(amr_cur)%mv_ghost_b)
        end if
        call s_amr_zero_fine_registers()

        do sub = 1, 2
            do s = 1, 3
                th = (real(sub - 1, wp) + c_abs(s))*0.5_wp

                ! lerp the ghost shell into q_cons at the stage time (device kernel; interior untouched)
                call s_amr_lerp_fine_ghosts(amr_slots(amr_cur)%q_ghost_a, amr_slots(amr_cur)%q_ghost_b, &
                                            & amr_slots(amr_cur)%q_cons, th)
                if (qbmm .and. .not. polytropic) call s_amr_lerp_fine_ghosts_pbmv(th)
                if (rank_time_wrt) call s_rank_time_toc()

                ! continuation-face ghosts AFTER the lerp: the substeps run in lockstep across ranks, so the
                ! neighbor's current fine q_cons is the same-time data (the lerp stays block-boundary-only)
                call s_mpi_sendrecv_amr_fine_halo(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%m, amr_slots(amr_cur)%n, &
                                                  & amr_slots(amr_cur)%p)
                if (rank_time_wrt) call s_rank_time_tic()

                ! substep-entry backup for the SSP-RK combination (device copy, interior only)
                if (s == 1) then
                    call s_amr_copy_fine_fields(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, 0, &
                                                & amr_slots(amr_cur)%m, 0, amr_slots(amr_cur)%n, 0, amr_slots(amr_cur)%p)
                    if (qbmm .and. .not. polytropic) call s_amr_backup_pbmv()
                end if

                amr_in_fine_advance = .true.
                call s_amr_swap_to_fine()
                ! widen the conversion range to the ghost shell (restored by s_amr_restore_coarse)
                idwint = amr_slots(amr_cur)%idwbuff
                $:GPU_UPDATE(device='[idwint]')
                if (qbmm .and. .not. polytropic) then
                    ! the block's OWN side-state and rhs scratch (the coarse arrays stay untouched)
                    call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, &
                                       & amr_slots(amr_cur)%rhs, amr_slots(amr_cur)%pb_f, amr_slots(amr_cur)%rhs_pb_f, &
                                       & amr_slots(amr_cur)%mv_f, amr_slots(amr_cur)%rhs_mv_f, t_step, time_avg, s)
                else
                    call s_compute_rhs(amr_slots(amr_cur)%q_cons, q_T_sf, amr_slots(amr_cur)%q_prim, bc_type, &
                                       & amr_slots(amr_cur)%rhs, pb_in, rhs_pb, mv_in, rhs_mv, t_step, time_avg, s)
                end if
                call s_amr_restore_coarse()
                amr_in_fine_advance = .false.

                ! RK stage update at the FINE time step (device kernel)
                call s_amr_fine_rk_update(amr_slots(amr_cur)%q_cons, amr_slots(amr_cur)%q_cons_stor, amr_slots(amr_cur)%rhs, &
                                          & coefs(s, 1), coefs(s, 2), coefs(s, 3), coefs(s, 4), amr_dt_fine)
                if (qbmm .and. .not. polytropic) then
                    call s_amr_fine_rk_update_pbmv(coefs(s, 1), coefs(s, 2), coefs(s, 3), coefs(s, 4), amr_dt_fine)
                end if
                ! 6-equation model: per-substage pressure relaxation (instantaneous equilibration -
                ! per stage at fine dt is the same infinite-rate limit the coarse applies per stage)
                if (model_eqns == model_eqns_6eq .and. (.not. relax)) call s_amr_pressure_relax_fine()
                ! moving body: rebuild the fine-block IB state at the body's fine sub-time position (th matches the fluid-ghost
                ! lerp)
                if (moving_immersed_boundary_flag) call s_amr_update_mib_fine(th)
                ! IB state correction on the fine block after each substep RK update (no-op unless ib)
                call s_amr_ib_correct_fine()
            end do
        end do
        if (rank_time_wrt) call s_rank_time_toc()

    end subroutine s_advance_amr_fine_substeps

    !> Shrink box [blo:bhi] to the tight bounding box of the tagged (gtag==1) cells inside it. ok=.false. if no tagged cell.
    !! Collapsed dims (lo=hi=0) survive unchanged. Deterministic (integer scan of the identical global tag field).
    impure subroutine s_amr_trim_box(gtag, blo, bhi, ok)

        integer, intent(in)    :: gtag(0:,0:,0:)
        integer, intent(inout) :: blo(3), bhi(3)
        logical, intent(out)   :: ok
        integer                :: tlo(3), thi(3), i, j, k

        tlo = huge(1); thi = -huge(1)
        do k = blo(3), bhi(3)
            do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)
                    if (gtag(i, j, k) == 1) then
                        tlo(1) = min(tlo(1), i); thi(1) = max(thi(1), i)
                        tlo(2) = min(tlo(2), j); thi(2) = max(thi(2), j)
                        tlo(3) = min(tlo(3), k); thi(3) = max(thi(3), k)
                    end if
                end do
            end do
        end do
        ok = thi(1) >= tlo(1)
        if (ok) then; blo = tlo; bhi = thi; end if

    end subroutine s_amr_trim_box

    !> Berger-Rigoutsos bisection of one (already tagged-trimmed) candidate box on the global tag field: pick the longest splittable
    !! axis, prefer a zero-signature hole (widest interior run), else the strongest signature inflection (Laplacian sign change).
    !! ok=.false. if no axis admits a split leaving both children >= 2 cells. Integer-only => identical on all ranks.
    impure subroutine s_amr_find_split(gtag, blo, bhi, sax, spos, ok)

        integer, intent(in)  :: gtag(0:,0:,0:)
        integer, intent(in)  :: blo(3), bhi(3)
        integer, intent(out) :: sax, spos
        logical, intent(out) :: ok
        integer, parameter   :: min_child = 2
        integer              :: axord(3), ext(3), d, ax, t, u, v, s
        integer              :: run, run_start, best_run, best_start, lap, prevlap, bestmag, bestpos
        integer, allocatable :: sig(:)

        ok = .false.; sax = 0; spos = 0
        ext = bhi - blo + 1
        axord = [1, 2, 3]  ! sort axes by descending extent (deterministic bubble)
        do d = 1, 2
            do ax = 1, 3 - d
                if (ext(axord(ax)) < ext(axord(ax + 1))) then
                    s = axord(ax); axord(ax) = axord(ax + 1); axord(ax + 1) = s
                end if
            end do
        end do
        allocate (sig(0:maxval(bhi)))
        do d = 1, 3
            ax = axord(d)
            if (ax > num_dims) cycle
            if (ext(ax) < 2*min_child) cycle
            do t = blo(ax), bhi(ax)
                sig(t) = 0
                select case (ax)
                case (1); do v = blo(3), bhi(3); do u = blo(2), bhi(2); sig(t) = sig(t) + gtag(t, u, v); end do; end do
                case (2); do v = blo(3), bhi(3); do u = blo(1), bhi(1); sig(t) = sig(t) + gtag(u, t, v); end do; end do
                case (3); do v = blo(2), bhi(2); do u = blo(1), bhi(1); sig(t) = sig(t) + gtag(u, v, t); end do; end do
                end select
            end do
            ! (1) widest interior zero run (box is trimmed => sig(blo)>0 and sig(bhi)>0, so any run is interior)
            best_run = 0; best_start = -1; run = 0; run_start = -1
            do t = blo(ax), bhi(ax)
                if (sig(t) == 0) then
                    if (run == 0) run_start = t
                    run = run + 1
                else
                    if (run > best_run) then; best_run = run; best_start = run_start; end if
                    run = 0
                end if
            end do
            if (best_start > blo(ax)) then
                spos = best_start
                if (spos - blo(ax) >= min_child .and. bhi(ax) - spos + 1 >= min_child) then
                    sax = ax; ok = .true.; deallocate (sig); return
                end if
            end if
            ! (2) strongest inflection: Laplacian sign change with the largest jump
            bestmag = -1; bestpos = -1; prevlap = 0
            do t = blo(ax) + 1, bhi(ax) - 1
                lap = sig(t - 1) - 2*sig(t) + sig(t + 1)
                if (t > blo(ax) + 1) then
                    if (((lap < 0) .neqv. (prevlap < 0)) .and. abs(lap - prevlap) > bestmag .and. t - blo(ax) >= min_child &
                        & .and. bhi(ax) - t + 1 >= min_child) then
                        bestmag = abs(lap - prevlap); bestpos = t
                    end if
                end if
                prevlap = lap
            end do
            if (bestpos > 0) then
                sax = ax; spos = bestpos; ok = .true.; deallocate (sig); return
            end if
        end do
        deallocate (sig)

    end subroutine s_amr_find_split

    !> True iff global level-0 cell (gi, gj, gk) lies inside any acoustic source support bbox.
    pure logical function f_in_acoustic_support(gi, gj, gk) result(insup)

        integer, intent(in) :: gi, gj, gk
        integer             :: s

        insup = .false.
        do s = 1, num_source
            if (gi >= acoustic_supp_lo(1, s) .and. gi <= acoustic_supp_hi(1, s) .and. (n_glb == 0 .or. (gj >= acoustic_supp_lo(2, &
                & s) .and. gj <= acoustic_supp_hi(2, s))) .and. (p_glb == 0 .or. (gk >= acoustic_supp_lo(3, &
                & s) .and. gk <= acoustic_supp_hi(3, s)))) then
                insup = .true.; return
            end if
        end do

    end function f_in_acoustic_support

    !> Clip a candidate regrid box (global indices) so it does not overlap any acoustic source support bbox: per overlapping source,
    !! remove the overlap along the single axis/side that keeps the largest remaining extent (deterministic: lower axis, then begin
    !! side, wins ties). Clipping only shrinks the box and may empty it (hi < lo); the caller drops empties.
    impure subroutine s_amr_clip_box_from_sources(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: s, d, best_d, best_side, best_ext, ext_l, ext_r
        logical                :: ovl

        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return  ! emptied by an earlier clip
            ovl = lo(1) <= acoustic_supp_hi(1, s) .and. hi(1) >= acoustic_supp_lo(1, s)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= acoustic_supp_hi(2, s) .and. hi(2) >= acoustic_supp_lo(2, s)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= acoustic_supp_hi(3, s) .and. hi(3) >= acoustic_supp_lo(3, s)
            if (.not. ovl) cycle
            best_d = 1; best_side = 1; best_ext = -1
            do d = 1, num_dims
                ext_l = acoustic_supp_lo(d, s) - lo(d)  ! cells kept by [lo(d), supp_lo-1]
                ext_r = hi(d) - acoustic_supp_hi(d, s)  ! cells kept by [supp_hi+1, hi(d)]
                if (ext_l > best_ext) then; best_ext = ext_l; best_d = d; best_side = 1; end if
                if (ext_r > best_ext) then; best_ext = ext_r; best_d = d; best_side = 2; end if
            end do
            if (best_side == 1) then
                hi(best_d) = acoustic_supp_lo(best_d, s) - 1
            else
                lo(best_d) = acoustic_supp_hi(best_d, s) + 1
            end if
        end do
        ! safety net: clipping removed every overlap by construction - anything left is a bug
        do s = 1, num_source
            if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) return
            ovl = lo(1) <= acoustic_supp_hi(1, s) .and. hi(1) >= acoustic_supp_lo(1, s)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= acoustic_supp_hi(2, s) .and. hi(2) >= acoustic_supp_lo(2, s)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= acoustic_supp_hi(3, s) .and. hi(3) >= acoustic_supp_lo(3, s)
            if (ovl) call s_mpi_abort('amr regrid: acoustic source exclusion clip failed (internal error)')
        end do

    end subroutine s_amr_clip_box_from_sources

    !> Expand a candidate regrid box (global indices) to fully contain every immersed body it overlaps, with a buff_size margin (the
    !! IB image-point stencils need resolved surroundings). Expansion is re-clamped to the domain interior by the caller's own
    !! guards; a body too large for the per-rank block cap aborts with a named message. The bbox reads the live centroid, so a
    !! moving body's box tracks its current position; between regrids s_amr_update_mib_fine guards containment.
    !> Margin-padded global coarse-index bounding box of immersed body i (supported analytic geometries only; aborts on others).
    !! Reads the body's LIVE centroid, so a moving body's box tracks its current position.
    impure subroutine s_amr_body_bbox(i, mrg, blo, bhi)

        integer, intent(in)  :: i, mrg
        integer, intent(out) :: blo(3), bhi(3)
        real(wp)             :: c(3), half(3)

        c = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, patch_ib(i)%z_centroid]
        select case (patch_ib(i)%geometry)
        case (2, 8, 10)  ! circle, sphere, cylinder: radius-bounded (cylinder length adds below)
            half = patch_ib(i)%radius
            if (patch_ib(i)%geometry == 10) then
                half(1) = max(half(1), 0.5_wp*patch_ib(i)%length_x)
                half(2) = max(half(2), 0.5_wp*patch_ib(i)%length_y)
                half(3) = max(half(3), 0.5_wp*patch_ib(i)%length_z)
            end if
        case (3, 9)  ! rectangle, box
            half = 0.5_wp*[patch_ib(i)%length_x, patch_ib(i)%length_y, patch_ib(i)%length_z]
        case default
            call s_mpi_abort('amr dynamic regrid with ib: unsupported body geometry for the ' &
                             & // 'containment bounding box (supported: circle/rectangle/sphere/box/cylinder)')
        end select
        ! physical bbox -> global coarse indices (uniform spacing is enforced for amr; the
        ! axisymmetric half axis cell only shrinks dy(0), so the floor is still conservative)
        blo(1) = int((c(1) - half(1) - x_domain%beg)/dx(0)) - mrg
        bhi(1) = int((c(1) + half(1) - x_domain%beg)/dx(0)) + mrg
        blo(2) = 0; bhi(2) = 0; blo(3) = 0; bhi(3) = 0
        if (n_glb > 0) then
            blo(2) = int((c(2) - half(2) - y_domain%beg)/dy(min(1, n))) - mrg
            bhi(2) = int((c(2) + half(2) - y_domain%beg)/dy(min(1, n))) + mrg
        end if
        if (p_glb > 0) then
            blo(3) = int((c(3) - half(3) - z_domain%beg)/dz(0)) - mrg
            bhi(3) = int((c(3) + half(3) - z_domain%beg)/dz(0)) + mrg
        end if

    end subroutine s_amr_body_bbox

    impure subroutine s_amr_expand_box_over_bodies(lo, hi)

        integer, intent(inout) :: lo(3), hi(3)
        integer                :: i, d, blo(3), bhi(3), mrg
        logical                :: ovl

        ! containment margin: the IB image-point stencil reaches a few cells beyond the surface
        ! (the validated static-block goldens keep ~5); buff_size (floored to 10 by ib) would
        ! exceed the per-rank block cap for ordinary bodies

        mrg = max(amr_buf, 4)

        do i = 1, num_ibs
            call s_amr_body_bbox(i, mrg, blo, bhi)
            ! blocks must stay buff_size inside the domain: a body whose margin-padded bbox does
            ! not fit cannot be contained - fail with a named message instead of a clipped body
            if (blo(1) < buff_size .or. bhi(1) > m_glb - buff_size .or. (n_glb > 0 .and. (blo(2) < buff_size .or. bhi(2) > n_glb &
                & - buff_size)) .or. (p_glb > 0 .and. (blo(3) < buff_size .or. bhi(3) > p_glb - buff_size))) then
                call s_mpi_abort('amr dynamic regrid with ib: the immersed body plus its containment ' &
                                 & // 'margin does not fit inside the refinable domain interior (blocks stay buff_size off the edges)')
            end if
            ovl = lo(1) <= bhi(1) .and. hi(1) >= blo(1)
            if (n_glb > 0) ovl = ovl .and. lo(2) <= bhi(2) .and. hi(2) >= blo(2)
            if (p_glb > 0) ovl = ovl .and. lo(3) <= bhi(3) .and. hi(3) >= blo(3)
            if (.not. ovl) cycle
            do d = 1, num_dims
                lo(d) = min(lo(d), blo(d))
                hi(d) = max(hi(d), bhi(d))
            end do
            if (hi(1) - lo(1) + 1 > amr_maxc_fit(1) .or. (n_glb > 0 .and. hi(2) - lo(2) + 1 > amr_maxc_fit(2)) .or. (p_glb > 0 &
                & .and. hi(3) - lo(3) + 1 > amr_maxc_fit(3))) then
                call s_mpi_abort('amr dynamic regrid with ib: containing the immersed body plus margin ' &
                                 & // 'exceeds the per-rank block size cap; use fewer ranks or a larger amr_maxc_fit')
            end if
        end do

    end subroutine s_amr_expand_box_over_bodies

    !> Cluster the local per-cell tag field into a LIST of separated block boxes (global level-0 cell indices), identically on every
    !! rank. Gathers the tags into a global field (allreduce MAX), runs Berger-Rigoutsos recursive bisection until each box's tag
    !! efficiency reaches amr_cluster_eff (or it is atomic / the amr_max_blocks cap is reached), then merges any two boxes whose
    !! amr_buf-padded extents come within buff_size (guaranteeing no fine-fine adjacency: separated boxes stay >= buff_size apart,
    !! nearby ones collapse to a single box == the legacy bounding box). Boxes are the raw tagged extents; the caller pads, clamps
    !! and size-caps each one.
    impure subroutine s_amr_cluster(tag_grid, boxes, nboxes)

        logical, intent(in)                   :: tag_grid(0:,0:,0:)
        type(t_box), allocatable, intent(out) :: boxes(:)
        integer, intent(out)                  :: nboxes
        integer, allocatable                  :: gtag(:,:,:), slo(:,:), shi(:,:), alo(:,:), ahi(:,:)
        integer                               :: mg, ng, pg, sidx(3), ci, cj, ck, gi, gj, gk
        integer                               :: cap, nwork, nacc, i, j, d, sax, spos, thr, ntag, vol
        integer                               :: blo(3), bhi(3)
        logical                               :: ok, force, capped, changed, tooclose
        real(wp)                              :: eff

#ifdef MFC_MPI
        integer :: ierr
#endif

        nboxes = 0
        mg = m_glb; ng = 0; pg = 0
        if (n_glb > 0) ng = n_glb
        if (p_glb > 0) pg = p_glb
        allocate (gtag(0:mg,0:ng,0:pg)); gtag = 0
        sidx = 0; sidx(1) = start_idx(1)
        if (n_glb > 0) sidx(2) = start_idx(2)
        if (p_glb > 0) sidx(3) = start_idx(3)
        do ck = 0, p
            do cj = 0, n
                do ci = 0, m
                    if (tag_grid(ci, cj, ck)) then
                        gi = ci + sidx(1); gj = 0; gk = 0
                        if (n_glb > 0) gj = cj + sidx(2)
                        if (p_glb > 0) gk = ck + sidx(3)
                        gtag(gi, gj, gk) = 1
                    end if
                end do
            end do
        end do
#ifdef MFC_MPI
        ! every rank ORs in its local tags => an identical global tag field, so the bisection below is rank-invariant (SP7a)
        if (num_procs > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, gtag, (mg + 1)*(ng + 1)*(pg + 1), MPI_INTEGER, MPI_MAX, &
            & MPI_COMM_WORLD, ierr)
#endif
        if (sum(gtag) == 0) then; deallocate (gtag); return; end if

        cap = amr_max_blocks
        allocate (slo(3, 4*cap + 8), shi(3, 4*cap + 8), alo(3, cap), ahi(3, cap))
        nwork = 1; slo(:,1) = [0, 0, 0]; shi(:,1) = [mg, ng, pg]  ! first pop trims to the global tagged bbox
        nacc = 0; capped = .false.
        do while (nwork > 0)
            blo = slo(:,nwork); bhi = shi(:,nwork); nwork = nwork - 1
            call s_amr_trim_box(gtag, blo, bhi, ok)
            if (.not. ok) cycle
            ntag = 0
            do ck = blo(3), bhi(3); do cj = blo(2), bhi(2); do ci = blo(1), bhi(1)
                ntag = ntag + gtag(ci, cj, ck)
            end do; end do; end do
            vol = 1
            do d = 1, num_dims; vol = vol*(bhi(d) - blo(d) + 1); end do
                eff = real(ntag, wp)/real(max(vol, 1), wp)
                call s_amr_find_split(gtag, blo, bhi, sax, spos, ok)
                force = (nacc + nwork + 1 >= cap)  ! splitting now could overflow the amr_max_blocks cap
                if (eff >= amr_cluster_eff .or. .not. ok .or. force) then
                    if (nacc < cap) then; nacc = nacc + 1; alo(:,nacc) = blo; ahi(:,nacc) = bhi; end if
                    if (force .and. ok .and. eff < amr_cluster_eff) capped = .true.
                else
                    slo(:,nwork + 1) = blo; shi(:,nwork + 1) = bhi; shi(sax, nwork + 1) = spos - 1
                    slo(:,nwork + 2) = blo; shi(:,nwork + 2) = bhi; slo(sax, nwork + 2) = spos
                    nwork = nwork + 2
                end if
            end do
            deallocate (gtag)

            ! min-separation merge: two boxes are separated only if some active dim's gap reaches thr; else fuse to their bounding
            ! box
            thr = buff_size + 2*amr_buf
            changed = .true.
            do while (changed)
                changed = .false.
                outer: do i = 1, nacc - 1
                    do j = i + 1, nacc
                        tooclose = .true.
                        do d = 1, num_dims
                            if (max(alo(d, i), alo(d, j)) - min(ahi(d, i), ahi(d, j)) - 1 >= thr) tooclose = .false.
                        end do
                        if (tooclose) then
                            alo(:,i) = min(alo(:,i), alo(:,j)); ahi(:,i) = max(ahi(:,i), ahi(:,j))
                            alo(:,j) = alo(:,nacc); ahi(:,j) = ahi(:,nacc)
                            nacc = nacc - 1; changed = .true.
                            exit outer
                        end if
                    end do
                end do outer
            end do
            if (capped .and. proc_rank == 0) print '(A,I0)', ' [amr] WARNING: tag clustering capped at amr_max_blocks = ', cap

            nboxes = nacc
            allocate (boxes(nboxes))
            do i = 1, nboxes
                boxes(i)%lo = alo(:,i); boxes(i)%hi = ahi(:,i)
            end do
            deallocate (slo, shi, alo, ahi)

        end subroutine s_amr_cluster

        !> Regrid: tag by relative density gradient into a per-cell field, cluster (Berger-Rigoutsos + min-separation merge) into a
        !! list of separated boxes, pad/clamp/size-cap each, and rebuild every active slot. Each new box's slot prolongs from coarse
        !! then overwrites its overlap with whichever OLD slot(s) covered it (rank-local by construction; a split copies from one
        !! old slot, a merge from both). Called between steps only. No-op if nothing is tagged or the box set is unchanged.
        impure subroutine s_amr_regrid(q_cons_base)

            type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_base
            logical, allocatable                                   :: tag_grid(:,:,:)
            type(t_box), allocatable                               :: boxes(:)
            integer                                                :: lo(3), hi(3), sh(3), old_np, k, kk
            integer                                                :: old_ilo(3, amr_max_blocks), old_ext(3, amr_max_blocks)
            logical                                                :: old_owns(amr_max_blocks), any_xchg, same, merged
            integer                                                :: ci, cj, ck, fi, fj, fk, ofi, ofj, ofk, i
            integer                                                :: sidx(3), tg_lo(3), tg_hi(3), nboxes
            real(wp)                                               :: r0, g

            ! valid coarse CONS ghosts at internal rank boundaries: the tag sweep reads +/-1 across seams and the rebuild
            ! prolongation
            ! reads past the new intersection (ALL ranks call: pairwise per-direction exchange; complete no-op at np=1).

            call s_amr_exchange_coarse_cons_halo(q_cons_base)
            do i = 1, sys_size
                $:GPU_UPDATE(host='[q_cons_base(i)%sf]')
            end do

            ! 1) per-cell tag field (density-gradient criterion, unchanged), skipping the two global boundary cells per active dim
            sidx = 0
            sidx(1) = start_idx(1)
            if (n_glb > 0) sidx(2) = start_idx(2)
            if (p_glb > 0) sidx(3) = start_idx(3)
            tg_lo = 0; tg_hi = 0
            tg_lo(1) = merge(1, 0, sidx(1) == 0); tg_hi(1) = merge(m - 1, m, sidx(1) + m == m_glb)
            if (n_glb > 0) then; tg_lo(2) = merge(1, 0, sidx(2) == 0); tg_hi(2) = merge(n - 1, n, sidx(2) + n == n_glb); end if
            if (p_glb > 0) then; tg_lo(3) = merge(1, 0, sidx(3) == 0); tg_hi(3) = merge(p - 1, p, sidx(3) + p == p_glb); end if
            allocate (tag_grid(0:m,0:n,0:p)); tag_grid = .false.
            do ck = tg_lo(3), tg_hi(3)
                do cj = tg_lo(2), tg_hi(2)
                    do ci = tg_lo(1), tg_hi(1)
                        ! total density gradient (sum of the continuity variables): degenerates to the single-fluid tagger and is
                        ! immune to trace-fluid noise. Matched-density composition-only interfaces are invisible (documented limit).
                        r0 = max(abs(f_amr_rho_tot(q_cons_base, ci, cj, ck)), 1.e-30_wp)
                        g = abs(f_amr_rho_tot(q_cons_base, ci + 1, cj, ck) - f_amr_rho_tot(q_cons_base, ci - 1, cj, ck))
                        if (n_glb > 0) g = max(g, abs(f_amr_rho_tot(q_cons_base, ci, cj + 1, ck) - f_amr_rho_tot(q_cons_base, ci, &
                            & cj - 1, ck)))
                        if (p_glb > 0) g = max(g, abs(f_amr_rho_tot(q_cons_base, ci, cj, ck + 1) - f_amr_rho_tot(q_cons_base, ci, &
                            & cj, ck - 1)))
                        if (g/(2._wp*r0) > amr_tag_eps) tag_grid(ci, cj, ck) = .true.
                        ! the acoustic source support stays coarse (its spatials are coarse cell
                        ! indices): suppress tags there so the clusterer splits around the source
                        if (acoustic_source .and. tag_grid(ci, cj, ck)) then
                            if (f_in_acoustic_support(ci + sidx(1), cj + sidx(2), ck + sidx(3))) tag_grid(ci, cj, ck) = .false.
                        end if
                    end do
                end do
            end do

            ! 2) cluster into a list of separated boxes (deterministic on all ranks)
            call s_amr_cluster(tag_grid, boxes, nboxes)
            deallocate (tag_grid)
            if (nboxes == 0) return  ! nothing tagged on any rank; keep the current blocks

            ! 3) pad + clamp + size-cap each box (amr_maxc_fit lets each box move freely across rank boundaries); drop margin-only
            ! boxes
            k = 0
            do kk = 1, nboxes
                lo = boxes(kk)%lo; hi = boxes(kk)%hi
                lo(1) = max(lo(1) - amr_buf, buff_size); hi(1) = min(hi(1) + amr_buf, m_glb - buff_size)
                if (hi(1) - lo(1) + 1 > amr_maxc_fit(1)) hi(1) = lo(1) + amr_maxc_fit(1) - 1
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
                ! keep candidate boxes clear of every acoustic source support (the source acts on the
                ! coarse grid only); clipping only shrinks, so boxes stay disjoint - empties drop below
                if (acoustic_source) call s_amr_clip_box_from_sources(lo, hi)
                ! a fine block that PARTIALLY covers an immersed body is an untested regime (ghost
                ! prolongation through body-interior cells, refluxing across the body): any box that
                ! overlaps a body's bounding box is expanded to contain the whole body plus margin
                if (ib) call s_amr_expand_box_over_bodies(lo, hi)
                if (hi(1) < lo(1) .or. hi(2) < lo(2) .or. hi(3) < lo(3)) cycle  ! confined to the domain margin
                k = k + 1; boxes(k)%lo = lo; boxes(k)%hi = hi
            end do
            nboxes = k
            if (nboxes == 0) return

            if (ib) then
                ! body-containment expansion can make boxes overlap (bisection guaranteed disjoint
                ! boxes; two boxes near one body both grow over it): merge overlapping pairs to a
                ! bounding box until none remain - overlapping blocks would double-restrict/reflux
                merged = .true.
                do while (merged)
                    merged = .false.
                    outer: do k = 1, nboxes - 1
                        do kk = k + 1, nboxes
                            if (boxes(k)%lo(1) <= boxes(kk)%hi(1) .and. boxes(k)%hi(1) >= boxes(kk)%lo(1) .and. (n_glb == 0 &
                                & .or. (boxes(k)%lo(2) <= boxes(kk)%hi(2) .and. boxes(k)%hi(2) >= boxes(kk)%lo(2))) &
                                & .and. (p_glb == 0 .or. (boxes(k)%lo(3) <= boxes(kk)%hi(3) .and. boxes(k)%hi(3) &
                                & >= boxes(kk)%lo(3)))) then
                                boxes(k)%lo = min(boxes(k)%lo, boxes(kk)%lo)
                                boxes(k)%hi = max(boxes(k)%hi, boxes(kk)%hi)
                                boxes(kk) = boxes(nboxes); nboxes = nboxes - 1
                                if (boxes(k)%hi(1) - boxes(k)%lo(1) + 1 > amr_maxc_fit(1) .or. (n_glb > 0 .and. boxes(k)%hi(2) &
                                    & - boxes(k)%lo(2) + 1 > amr_maxc_fit(2)) .or. (p_glb > 0 .and. boxes(k)%hi(3) &
                                    & - boxes(k)%lo(3) + 1 > amr_maxc_fit(3))) then
                                    call s_mpi_abort('amr regrid: merging body-containing blocks exceeds ' &
                                                     & // 'the per-rank block size cap')
                                end if
                                merged = .true.
                                exit outer
                            end if
                        end do
                    end do outer
                end do
                ! the expansion may also have grown a box onto an acoustic source support: the two
                ! constraints (contain the body, exclude the source) cannot both hold - fail closed
                if (acoustic_source) then
                    do k = 1, nboxes
                        lo = boxes(k)%lo; hi = boxes(k)%hi
                        call s_amr_clip_box_from_sources(lo, hi)
                        if (any(lo /= boxes(k)%lo) .or. any(hi /= boxes(k)%hi)) then
                            call s_mpi_abort('amr regrid: a block must contain an immersed body AND stay ' &
                                             & // 'clear of an acoustic source support - the constraints conflict; ' &
                                             & // 'move the body or the source apart')
                        end if
                    end do
                end if
            end if

            ! 4) unchanged? (same count and boxes as the live slots -> keep them; a rebuild would reproduce them exactly anyway)
            if (nboxes == amr_num_blocks) then
                same = .true.
                do k = 1, nboxes
                    if (any(boxes(k)%lo /= amr_slots(k)%region%lo) .or. any(boxes(k)%hi /= amr_slots(k)%region%hi)) same = .false.
                end do
                if (same) return
            end if

            ! 5) stash every live slot's fine interior (dead-between-steps q_cons_stor bounce), keeping its old intersection origin
            old_np = amr_num_blocks
            do k = 1, old_np
                old_ilo(:,k) = amr_isect_lo_all(:,k)
                old_ext(:,k) = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                old_owns(k) = amr_owns_all(k)
                if (old_owns(k)) then
                    do i = 1, sys_size
                        $:GPU_UPDATE(host='[amr_slots(k)%q_cons(i)%sf]')
                    end do
                    do i = 1, sys_size
                        amr_slots(k)%q_cons_stor(i)%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, &
                                  & k)) = amr_slots(k)%q_cons(i)%sf(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k))
                    end do
                    ! non-polytropic QBMM: the side-state bounces through pb/mv_stor exactly like
                    ! q_cons (both stors are dead between steps)
                    if (qbmm .and. .not. polytropic) then
                        $:GPU_UPDATE(host='[amr_slots(k)%pb_f, amr_slots(k)%mv_f]')
                        amr_slots(k)%pb_stor(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                                  & :) = amr_slots(k)%pb_f(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                        amr_slots(k)%mv_stor(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:, &
                                  & :) = amr_slots(k)%mv_f(0:old_ext(1, k),0:old_ext(2, k),0:old_ext(3, k),:,:)
                    end if
                end if
            end do
            ! coarse pb/mv host-current for the per-block re-prolongation below
            if (qbmm .and. .not. polytropic) then
                $:GPU_UPDATE(host='[pb_ts(1)%sf, mv_ts(1)%sf]')
            end if

            ! 6) build each new slot: geometry (collective on all ranks), prolong, then overlap-copy from every covering old slot
            amr_num_blocks = nboxes
            any_xchg = .false.
            if (proc_rank == 0) print '(A,I0,A)', ' [amr] regrid: ', nboxes, ' block(s)'
            do k = 1, nboxes
                amr_cur = k
                call s_set_amr_fine_geometry(boxes(k)%lo, boxes(k)%hi)
                any_xchg = any_xchg .or. amr_xchg_coarse_ghosts
                if (proc_rank == 0) print '(A,I0,A,I0,A,I0,A,I0,A)', ' [amr]   block ', k, ': box x ', boxes(k)%lo(1), ':', &
                    & boxes(k)%hi(1), ' (', (boxes(k)%hi(1) - boxes(k)%lo(1) + 1), ' coarse cells)'
                if (.not. amr_rank_owns_block) cycle
                call s_interpolate_coarse_to_fine(q_cons_base)
                do kk = 1, old_np
                    if (.not. old_owns(kk)) cycle
                    sh = 2*(amr_isect_lo - old_ilo(:,kk))  ! old LOCAL fine index = new LOCAL fine index + sh (collapsed dims sh=0)
                    do i = 1, sys_size
                        do fk = 0, amr_slots(k)%p
                            ofk = fk + sh(3)
                            if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_ext(3, kk))) cycle
                            do fj = 0, amr_slots(k)%n
                                ofj = fj + sh(2)
                                if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_ext(2, kk))) cycle
                                do fi = 0, amr_slots(k)%m
                                    ofi = fi + sh(1)
                                    if (ofi < 0 .or. ofi > old_ext(1, kk)) cycle
                                    amr_slots(k)%q_cons(i)%sf(fi, fj, fk) = amr_slots(kk)%q_cons_stor(i)%sf(ofi, ofj, ofk)
                                end do
                            end do
                        end do
                    end do
                end do
                do i = 1, sys_size
                    $:GPU_UPDATE(device='[amr_slots(k)%q_cons(i)%sf]')
                end do
                ! non-polytropic QBMM: prolong the side-state from coarse (piecewise-constant),
                ! then overwrite the overlap with the old blocks' fine data (same index shift)
                if (qbmm .and. .not. polytropic) then
                    call s_amr_prolong_pbmv_host()
                    do kk = 1, old_np
                        if (.not. old_owns(kk)) cycle
                        sh = 2*(amr_isect_lo - old_ilo(:,kk))
                        do fk = 0, amr_slots(k)%p
                            ofk = fk + sh(3)
                            if (p_glb > 0 .and. (ofk < 0 .or. ofk > old_ext(3, kk))) cycle
                            do fj = 0, amr_slots(k)%n
                                ofj = fj + sh(2)
                                if (n_glb > 0 .and. (ofj < 0 .or. ofj > old_ext(2, kk))) cycle
                                do fi = 0, amr_slots(k)%m
                                    ofi = fi + sh(1)
                                    if (ofi < 0 .or. ofi > old_ext(1, kk)) cycle
                                    amr_slots(k)%pb_f(fi, fj, fk,:,:) = amr_slots(kk)%pb_stor(ofi, ofj, ofk,:,:)
                                    amr_slots(k)%mv_f(fi, fj, fk,:,:) = amr_slots(kk)%mv_stor(ofi, ofj, ofk,:,:)
                                end do
                            end do
                        end do
                    end do
                    $:GPU_UPDATE(device='[amr_slots(k)%pb_f, amr_slots(k)%mv_f]')
                end if
                call s_mpi_sendrecv_amr_fine_halo(amr_slots(k)%q_cons, amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p)
            end do
            amr_xchg_coarse_ghosts = any_xchg  ! coarse halo exchanged once per step if ANY block needs it
            ! rebuild every block's fine-grid IB state for the NEW geometry (markers/ghost points/
            ! image points recomputed from the body definitions; no state carries across regrids)
            if (ib) call s_amr_setup_ib()
            call s_amr_select_slot(1)

        end subroutine s_amr_regrid

        !> Write the fine-level restart file for save step t_step alongside the level-0 restart (whose format stays untouched): the
        !! writing rank count, the active-block count, and for EACH block its box + each rank's intersection-local fine conservative
        !! state. Serial mode: one unformatted file per rank inside its level-0 step directory. Parallel mode: one shared MPI-IO
        !! file (3-int global header [np, nboxes, sys_size], then per block a 6-int box header, a 3*np-int per-rank fine-extents
        !! record [m,n,p per rank, 0s for non-owners; validated on read], followed by the ranks' fine blocks concatenated in rank
        !! order). Same rank count + decomposition required to restart (enforced by the extents record).
        impure subroutine s_write_amr_restart(t_step)

            integer, intent(in)                  :: t_step
            character(LEN=path_len + 3*name_len) :: file_loc
            integer                              :: i, k

#ifdef MFC_MPI
            integer                             :: ifile, ierr, cnt, idx, fi, fj, fk, reg(6), ibytes, sbytes
            integer                             :: myext(3)
            integer, allocatable                :: wext(:)
            integer, dimension(MPI_STATUS_SIZE) :: status
            integer(kind=MPI_OFFSET_KIND)       :: my_cnt, my_off, tot_cnt, disp0, ddisp
            logical                             :: file_exist
            real(stp), allocatable              :: buf(:)
#endif

            if (.not. amr) return
            ! host consumer: the fine state is device-current during stepping (pull every owned slot)
            do k = 1, amr_num_blocks
                if (amr_owns_all(k)) then
                    do i = 1, sys_size
                        $:GPU_UPDATE(host='[amr_slots(k)%q_cons(i)%sf]')
                    end do
                end if
            end do

            if (.not. parallel_io) then
                ! per-rank file in the step directory freshly created by the level-0 serial write
                write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step, '/amr_fine.dat'
                open (2, FILE=trim(file_loc), form='unformatted', STATUS='new')
                write (2) num_procs, amr_num_blocks, sys_size
                do k = 1, amr_num_blocks
                    write (2) amr_slots(k)%region%lo, amr_slots(k)%region%hi, amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p
                    if (amr_owns_all(k)) then
                        do i = 1, sys_size
                            write (2) amr_slots(k)%q_cons(i)%sf(0:amr_slots(k)%m,0:amr_slots(k)%n,0:amr_slots(k)%p)
                        end do
                    end if
                end do
                close (2)
            else
#ifdef MFC_MPI
                ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
                write (file_loc, '(A,I0,A)') 'amr_', t_step, '.dat'
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
                inquire (FILE=trim(file_loc), EXIST=file_exist)
                if (file_exist .and. proc_rank == 0) then
                    call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
                end if
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
                ! MPI-IO file handles default to MPI_ERRORS_RETURN: failures are silent unless checked
                if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_OPEN failed for ' // trim(file_loc))
                if (proc_rank == 0) call MPI_FILE_WRITE_AT(ifile, int(0, MPI_OFFSET_KIND), [num_procs, amr_num_blocks, sys_size], &
                    & 3, MPI_INTEGER, status, ierr)
                disp0 = int(3*ibytes, MPI_OFFSET_KIND)  ! running byte offset past the 3-int global header
                do k = 1, amr_num_blocks
                    cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                    if (.not. amr_owns_all(k)) cnt = 0
                    my_cnt = int(cnt, MPI_OFFSET_KIND)
                    my_off = int(0, MPI_OFFSET_KIND)
                    call MPI_EXSCAN(my_cnt, my_off, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                    if (proc_rank == 0) my_off = int(0, MPI_OFFSET_KIND)
                    call MPI_ALLREDUCE(my_cnt, tot_cnt, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                    if (proc_rank == 0) then
                        reg(1:3) = amr_slots(k)%region%lo; reg(4:6) = amr_slots(k)%region%hi
                        call MPI_FILE_WRITE_AT(ifile, disp0, reg, 6, MPI_INTEGER, status, ierr)
                    end if
                    ! per-rank fine extents (0s for non-owning ranks): readers rebuild this vector
                    ! from their own decomposition and abort on mismatch - a different rank count,
                    ! ownership pattern, or load_balance split would otherwise silently misalign
                    ! the concatenated per-rank data slices below
                    if (.not. allocated(wext)) allocate (wext(3*num_procs))
                    myext = 0
                    if (amr_owns_all(k)) myext = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                    call MPI_ALLGATHER(myext, 3, MPI_INTEGER, wext, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                    if (proc_rank == 0) then
                        call MPI_FILE_WRITE_AT(ifile, disp0 + int(6*ibytes, MPI_OFFSET_KIND), wext, 3*num_procs, MPI_INTEGER, &
                                               & status, ierr)
                    end if
                    ddisp = disp0 + int((6 + 3*num_procs)*ibytes, MPI_OFFSET_KIND)
                    allocate (buf(max(cnt, 1)))
                    idx = 0
                    do i = 1, sys_size
                        do fk = 0, amr_slots(k)%p
                            do fj = 0, amr_slots(k)%n
                                do fi = 0, amr_slots(k)%m
                                    idx = idx + 1
                                    buf(idx) = amr_slots(k)%q_cons(i)%sf(fi, fj, fk)
                                end do
                            end do
                        end do
                    end do
                    call MPI_FILE_WRITE_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, &
                                               & mpi_io_p, status, ierr)
                    if (ierr /= MPI_SUCCESS) &
                        & call s_mpi_abort('amr restart write: data write failed (disk full/quota?); the file is unusable')
                    deallocate (buf)
                    disp0 = ddisp + tot_cnt*int(sbytes, MPI_OFFSET_KIND)
                end do
                ! the close is where buffered MPI-IO data flushes on many stacks - a failure here truncates the file
                call MPI_FILE_CLOSE(ifile, ierr)
                if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart write: MPI_FILE_CLOSE failed; the file may be truncated')
#endif
            end if

        end subroutine s_write_amr_restart

        !> Restore the fine level from the AMR restart file at t_step_start (n_start under cfl_dt), if one exists: for each saved
        !! block rebuild the box via s_set_amr_fine_geometry, then read each rank's intersection-local fine state (exact stp
        !! round-trip). Requires the rank count + decomposition that wrote the file. restored = false on a fresh start, or - with a
        !! one-line warning - on a legacy restart without the file; the caller then re-prolongs from coarse. Collective: ALL ranks
        !! call together.
        impure subroutine s_read_amr_restart(restored)

            logical, intent(out)                 :: restored
            character(LEN=path_len + 3*name_len) :: file_loc
            character(LEN=300)                   :: msg
            logical                              :: file_exist
            integer                              :: i, k, ts, have_loc, have_glb, ghdr(3), reg(6), rm, rn, rp

#ifdef MFC_MPI
            integer                             :: ifile, ierr, cnt, idx, fi, fj, fk, ibytes, sbytes
            integer                             :: myext(3)
            integer, allocatable                :: wext(:), rext(:)
            integer, dimension(MPI_STATUS_SIZE) :: status
            integer(kind=MPI_OFFSET_KIND)       :: my_cnt, my_off, tot_cnt, disp0, ddisp, fsz
            real(stp), allocatable              :: buf(:)
#endif

            restored = .false.
            if (.not. amr) return
            if (cfl_dt) then
                ts = n_start
            else
                ts = t_step_start
            end if
            if (ts == 0) return  ! fresh start: the fine level is prolonged from the pre_process ICs

            if (.not. parallel_io) then
                write (file_loc, '(A,I0,A,I0,A)') trim(case_dir) // '/p_all/p', proc_rank, '/', ts, '/amr_fine.dat'
            else
                write (file_loc, '(A,I0,A)') 'amr_', ts, '.dat'
                file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
            end if
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            have_loc = merge(1, 0, file_exist)
            call s_mpi_allreduce_integer_min(have_loc, have_glb)
            if (have_glb == 0) then
                if (proc_rank == 0) then
                    print '(A)', &
                        & ' [amr] WARNING: no AMR restart file at this step; the fine level is re-initialized by ' &
                        & // 'prolongation from coarse (fine-level accuracy is lost across this restart)'
                end if
                return
            end if

            if (.not. parallel_io) then
                open (2, FILE=trim(file_loc), form='unformatted', ACTION='read', STATUS='old')
                read (2) ghdr
                if (ghdr(1) /= num_procs) then
                    write (msg, '(A,I0,A,I0,A)') 'amr restart rank-count mismatch: the AMR restart file was written with ', &
                           & ghdr(1), ' ranks but this run has ', num_procs, '; restart with the same rank count'
                    call s_mpi_abort(trim(msg))
                end if
                if (ghdr(3) /= sys_size) then
                    write (msg, '(A,I0,A,I0,A)') 'amr restart sys_size mismatch: the AMR restart file has ', ghdr(3), &
                           & ' conserved variables but this run has ', sys_size, &
                           & '; the physics configuration ' &
                           & // '(num_fluids/model_eqns/bubbles/chemistry) must match the run that wrote the restart'
                    call s_mpi_abort(trim(msg))
                end if
                if (ghdr(2) < 1 .or. ghdr(2) > amr_max_blocks) then
                    call s_mpi_abort('amr restart: the file holds more fine blocks than amr_max_blocks ' &
                                     & // 'in this run; restart with amr_max_blocks at least the written block count')
                end if
                amr_num_blocks = ghdr(2)
                do k = 1, amr_num_blocks
                    amr_cur = k
                    read (2) reg, rm, rn, rp
                    ! corrupt/foreign-file guard: a box outside the global domain would drive the geometry
                    ! build and coordinate reads out of bounds silently in release builds
                    if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb &
                        & .or. reg(2) > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
                        call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
                    end if
                    call s_set_amr_fine_geometry(reg(1:3), reg(4:6))
                    if (rm /= amr_slots(k)%m .or. rn /= amr_slots(k)%n .or. rp /= amr_slots(k)%p) then
                        call s_mpi_abort('amr restart decomposition mismatch: this rank''s fine extents differ from the file' &
                                         & // ' (identical decomposition - rank count and load_balance settings - required)')
                    end if
                    if (amr_rank_owns_block) then
                        do i = 1, sys_size
                            read (2) amr_slots(k)%q_cons(i)%sf(0:amr_slots(k)%m,0:amr_slots(k)%n,0:amr_slots(k)%p)
                        end do
                    end if
                end do
                close (2)
            else
#ifdef MFC_MPI
                ibytes = storage_size(0)/8; sbytes = storage_size(0._stp)/8
                call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
                ! MPI-IO errors are silent by default (MPI_ERRORS_RETURN on file handles) and a read past EOF
                ! is not even an error - it returns short with an uninitialized tail. Grab the size up front;
                ! the exact expected byte count is compared after the layout records are consumed below.
                if (ierr /= MPI_SUCCESS) call s_mpi_abort('amr restart read: MPI_FILE_OPEN failed for ' // trim(file_loc))
                call MPI_FILE_GET_SIZE(ifile, fsz, ierr)
                call MPI_FILE_READ_AT_ALL(ifile, int(0, MPI_OFFSET_KIND), ghdr, 3, MPI_INTEGER, status, ierr)
                if (ghdr(1) /= num_procs) then
                    write (msg, '(A,I0,A,I0,A)') 'amr restart rank-count mismatch: the AMR restart file was written with ', &
                           & ghdr(1), ' ranks but this run has ', num_procs, '; restart with the same rank count'
                    call s_mpi_abort(trim(msg))
                end if
                if (ghdr(3) /= sys_size) then
                    write (msg, '(A,I0,A,I0,A)') 'amr restart sys_size mismatch: the AMR restart file has ', ghdr(3), &
                           & ' conserved variables but this run has ', sys_size, &
                           & '; the physics configuration ' &
                           & // '(num_fluids/model_eqns/bubbles/chemistry) must match the run that wrote the restart'
                    call s_mpi_abort(trim(msg))
                end if
                if (ghdr(2) < 1 .or. ghdr(2) > amr_max_blocks) then
                    call s_mpi_abort('amr restart: the file holds more fine blocks than amr_max_blocks ' &
                                     & // 'in this run; restart with amr_max_blocks at least the written block count')
                end if
                amr_num_blocks = ghdr(2)
                disp0 = int(3*ibytes, MPI_OFFSET_KIND)
                do k = 1, amr_num_blocks
                    amr_cur = k
                    call MPI_FILE_READ_AT_ALL(ifile, disp0, reg, 6, MPI_INTEGER, status, ierr)
                    ! corrupt/foreign-file guard: a box outside the global domain would drive the geometry
                    ! build and coordinate reads out of bounds silently in release builds
                    if (reg(1) < 0 .or. reg(4) > m_glb .or. reg(1) > reg(4) .or. (n_glb > 0 .and. (reg(2) < 0 .or. reg(5) > n_glb &
                        & .or. reg(2) > reg(5))) .or. (p_glb > 0 .and. (reg(3) < 0 .or. reg(6) > p_glb .or. reg(3) > reg(6)))) then
                        call s_mpi_abort('amr restart: corrupt block record (box outside the global domain)')
                    end if
                    if (.not. allocated(wext)) allocate (wext(3*num_procs), rext(3*num_procs))
                    call MPI_FILE_READ_AT_ALL(ifile, disp0 + int(6*ibytes, MPI_OFFSET_KIND), wext, 3*num_procs, MPI_INTEGER, &
                                              & status, ierr)
                    call s_set_amr_fine_geometry(reg(1:3), reg(4:6))
                    ! validate the writer's per-rank layout against this run's decomposition: a
                    ! mismatch (e.g. re-derived load_balance splits at restart) would silently
                    ! misalign every rank's slice of the concatenated data
                    myext = 0
                    if (amr_rank_owns_block) myext = [amr_slots(k)%m, amr_slots(k)%n, amr_slots(k)%p]
                    call MPI_ALLGATHER(myext, 3, MPI_INTEGER, rext, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                    if (any(rext /= wext)) then
                        call s_mpi_abort('amr restart: the per-rank fine-block layout in the file does not match ' &
                                         & // 'this run''s decomposition; the rank count, ownership, and (with ' &
                                         & // 'load_balance) the weighted splits must match the run that wrote the restart')
                    end if
                    cnt = sys_size*(amr_slots(k)%m + 1)*(amr_slots(k)%n + 1)*(amr_slots(k)%p + 1)
                    if (.not. amr_rank_owns_block) cnt = 0
                    my_cnt = int(cnt, MPI_OFFSET_KIND)
                    my_off = int(0, MPI_OFFSET_KIND)
                    call MPI_EXSCAN(my_cnt, my_off, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                    if (proc_rank == 0) my_off = int(0, MPI_OFFSET_KIND)
                    call MPI_ALLREDUCE(my_cnt, tot_cnt, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD, ierr)
                    ddisp = disp0 + int((6 + 3*num_procs)*ibytes, MPI_OFFSET_KIND)
                    allocate (buf(max(cnt, 1)))
                    call MPI_FILE_READ_AT_ALL(ifile, ddisp + my_off*int(sbytes, MPI_OFFSET_KIND), buf, cnt*mpi_io_type, mpi_io_p, &
                                              & status, ierr)
                    idx = 0
                    do i = 1, sys_size
                        do fk = 0, amr_slots(k)%p
                            do fj = 0, amr_slots(k)%n
                                do fi = 0, amr_slots(k)%m
                                    idx = idx + 1
                                    amr_slots(k)%q_cons(i)%sf(fi, fj, fk) = buf(idx)
                                end do
                            end do
                        end do
                    end do
                    deallocate (buf)
                    disp0 = ddisp + tot_cnt*int(sbytes, MPI_OFFSET_KIND)
                end do
                ! disp0 now equals the exact byte count a complete file must have: a truncated file (crashed
                ! writer, filesystem hiccup) passes every layout check above but returns short reads with
                ! garbage tails - fail closed instead of restoring uninitialized data as the fine level
                if (disp0 /= fsz) then
                    call s_mpi_abort('amr restart read: file size does not match the expected layout ' &
                                     & // '(truncated or corrupt amr restart file)')
                end if
                call MPI_FILE_CLOSE(ifile, ierr)
#endif
            end if

            ! restored fine state to the device (mirrors s_populate_amr_fine's push; host reads above)
            do k = 1, amr_num_blocks
                if (amr_owns_all(k)) then
                    do i = 1, sys_size
                        $:GPU_UPDATE(device='[amr_slots(k)%q_cons(i)%sf]')
                    end do
                end if
            end do
            ! non-polytropic QBMM: the restart file carries q_cons only; re-prolong each block's
            ! side-state from the restored coarse pb/mv (one-time piecewise-constant smoothing)
            if (qbmm .and. .not. polytropic) then
                do k = 1, amr_num_blocks
                    call s_amr_select_slot(k)
                    if (amr_owns_all(k)) call s_amr_prolong_pbmv_host()
                end do
            end if
            call s_amr_select_slot(1)
            restored = .true.
            if (proc_rank == 0) then
                print '(A,I0,A)', ' [amr] restart: restored fine level, ', amr_num_blocks, ' block(s)'
            end if

        end subroutine s_read_amr_restart

        !> Global Sum(dV*U) for the per-fluid masses (continuity variables) and energy (eqn_idx%E) over the level-0 interior. First
        !! call (finalize_report=F) stores the baselines; the finalize call prints the relative drifts (~roundoff with refluxing).
        impure subroutine s_amr_conservation_defect(q_cons_base, finalize_report)

            type(scalar_field), dimension(sys_size), intent(in) :: q_cons_base
            logical, intent(in)                                 :: finalize_report
            real(wp)                                            :: sm(num_fluids_max), se, dv, s_glb
            integer                                             :: ci, cj, ck, f

            if (.not. amr) return
            ! host consumer: diagnostics (host sum over exactly the summed fields). The init baseline call
            ! runs BEFORE s_initialize_gpu_vars pushes the ICs to the device, so it must NOT pull the
            ! (uninitialized) device copies.
            if (finalize_report) then
                do f = 1, num_fluids
                    $:GPU_UPDATE(host='[q_cons_base(f)%sf]')
                end do
                $:GPU_UPDATE(host='[q_cons_base(eqn_idx%E)%sf]')
            end if
            sm = 0._wp; se = 0._wp
            do ck = 0, p
                do cj = 0, n
                    do ci = 0, m
                        dv = dx(ci)
                        if (n_glb > 0) dv = dv*dy(cj)
                        if (p_glb > 0) dv = dv*dz(ck)
                        do f = 1, num_fluids
                            sm(f) = sm(f) + dv*real(q_cons_base(f)%sf(ci, cj, ck), wp)
                        end do
                        se = se + dv*real(q_cons_base(eqn_idx%E)%sf(ci, cj, ck), wp)
                    end do
                end do
            end do
            if (num_procs > 1) then
                do f = 1, num_fluids
                    call s_mpi_allreduce_sum(sm(f), s_glb); sm(f) = s_glb
                end do
                call s_mpi_allreduce_sum(se, s_glb); se = s_glb
            end if
            if (.not. finalize_report) then
                amr_mass0(1:num_fluids) = sm(1:num_fluids); amr_energy0 = se
            else if (proc_rank == 0) then
                do f = 1, num_fluids
                    print '(A,I0,A,ES12.4)', ' [amr] conservation defect: mass(', f, ') drift = ', &
                        & abs(sm(f) - amr_mass0(f))/max(abs(amr_mass0(f)), 1.e-30_wp)
                end do
                print '(A,ES12.4)', ' [amr] conservation defect: energy drift = ', abs(se - amr_energy0)/max(abs(amr_energy0), &
                    & 1.e-30_wp)
            end if

        end subroutine s_amr_conservation_defect

        !> Init-time operator verification: (b) linear reproduction, (c) restriction of an independent field. Uses
        !! amr_slots(amr_cur)%q_cons(1) as scratch; called before s_populate_amr_fine overwrites it.
        impure subroutine s_amr_operator_checks()

            type(scalar_field), allocatable :: cscr(:)
            integer                         :: fi, fj, fk, ci, cj, ck, ox, oy, oz
            real(wp)                        :: e, errb, errc, si_f, si_c, dvf, dvc, want

            if (.not. amr) return
            if (.not. amr_rank_owns_block) return
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
            call s_prolong_one_var(cscr(1), amr_slots(amr_cur)%q_cons(1))
            errb = 0._wp
            do fk = 0, amr_slots(amr_cur)%p
                do fj = 0, amr_slots(amr_cur)%n
                    do fi = 0, amr_slots(amr_cur)%m
                        want = 1._wp + 2._wp*amr_slots(amr_cur)%x_cc(fi)
                        if (n_glb > 0) want = want + 3._wp*amr_slots(amr_cur)%y_cc(fj)
                        if (p_glb > 0) want = want + 4._wp*amr_slots(amr_cur)%z_cc(fk)
                        e = abs(real(amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, fk), wp) - want)
                        if (e > errb) errb = e
                    end do
                end do
            end do

            ! (c) fill the fine block with a quadratic (NOT from prolongation), restrict, compare integrals
            do fk = 0, amr_slots(amr_cur)%p
                do fj = 0, amr_slots(amr_cur)%n
                    do fi = 0, amr_slots(amr_cur)%m
                        amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, fk) = amr_slots(amr_cur)%x_cc(fi)**2
                        if (n_glb > 0) amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, fk) = amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, &
                            & fk) + amr_slots(amr_cur)%y_cc(fj)**2
                        if (p_glb > 0) amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, fk) = amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, &
                            & fk) + amr_slots(amr_cur)%z_cc(fk)**2
                    end do
                end do
            end do
            call s_restrict_one_var(amr_slots(amr_cur)%q_cons(1), cscr(1))
            si_f = 0._wp; si_c = 0._wp
            do fk = 0, amr_slots(amr_cur)%p
                do fj = 0, amr_slots(amr_cur)%n
                    do fi = 0, amr_slots(amr_cur)%m
                        dvf = amr_slots(amr_cur)%dx(fi)
                        if (n_glb > 0) dvf = dvf*amr_slots(amr_cur)%dy(fj)
                        if (p_glb > 0) dvf = dvf*amr_slots(amr_cur)%dz(fk)
                        si_f = si_f + dvf*real(amr_slots(amr_cur)%q_cons(1)%sf(fi, fj, fk), wp)
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

        !> Total density (sum of the continuity variables) at one cell: the regrid tag field. Reduces to variable 1 for one fluid.
        pure function f_amr_rho_tot(q, ci, cj, ck) result(r)

            type(scalar_field), dimension(:), intent(in) :: q
            integer, intent(in)                          :: ci, cj, ck
            real(wp)                                     :: r
            integer                                      :: f

            r = 0._wp
            do f = eqn_idx%cont%beg, eqn_idx%cont%end
                r = r + real(q(f)%sf(ci, cj, ck), wp)
            end do

        end function f_amr_rho_tot

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

            integer :: i, islot

            if (.not. amr) return
            do islot = 1, amr_max_blocks
                do i = 1, sys_size
                    @:DEALLOCATE(amr_slots(islot)%q_cons(i)%sf)
                    @:DEALLOCATE(amr_slots(islot)%q_cons_stor(i)%sf)
                    @:DEALLOCATE(amr_slots(islot)%q_prim(i)%sf)
                    @:DEALLOCATE(amr_slots(islot)%rhs(i)%sf)
                    @:DEALLOCATE(amr_slots(islot)%q_ghost_a(i)%sf)
                    @:DEALLOCATE(amr_slots(islot)%q_ghost_b(i)%sf)
                end do
                @:DEALLOCATE(amr_slots(islot)%q_cons)
                @:DEALLOCATE(amr_slots(islot)%q_cons_stor)
                @:DEALLOCATE(amr_slots(islot)%q_prim)
                @:DEALLOCATE(amr_slots(islot)%rhs)
                @:DEALLOCATE(amr_slots(islot)%q_ghost_a)
                @:DEALLOCATE(amr_slots(islot)%q_ghost_b)
                if (qbmm .and. .not. polytropic) then
                    @:DEALLOCATE(amr_slots(islot)%pb_f)
                    @:DEALLOCATE(amr_slots(islot)%mv_f)
                    @:DEALLOCATE(amr_slots(islot)%pb_stor)
                    @:DEALLOCATE(amr_slots(islot)%mv_stor)
                    @:DEALLOCATE(amr_slots(islot)%rhs_pb_f)
                    @:DEALLOCATE(amr_slots(islot)%rhs_mv_f)
                    if (amr_subcycle) then
                        @:DEALLOCATE(amr_slots(islot)%pb_ghost_a)
                        @:DEALLOCATE(amr_slots(islot)%mv_ghost_a)
                        @:DEALLOCATE(amr_slots(islot)%pb_ghost_b)
                        @:DEALLOCATE(amr_slots(islot)%mv_ghost_b)
                    end if
                end if
                if (allocated(amr_slots(islot)%x_cb)) deallocate (amr_slots(islot)%x_cb, amr_slots(islot)%x_cc, amr_slots(islot)%dx)
                if (allocated(amr_slots(islot)%y_cb)) deallocate (amr_slots(islot)%y_cb, amr_slots(islot)%y_cc, amr_slots(islot)%dy)
                if (allocated(amr_slots(islot)%z_cb)) deallocate (amr_slots(islot)%z_cb, amr_slots(islot)%z_cc, amr_slots(islot)%dz)
            end do
            deallocate (amr_slots)
            deallocate (amr_region_lo_all, amr_region_hi_all, amr_isect_lo_all, amr_isect_hi_all, amr_owns_all)
            if (allocated(sw_x_cb)) deallocate (sw_x_cb, sw_x_cc, sw_dx)
            if (allocated(sw_y_cb)) deallocate (sw_y_cb, sw_y_cc, sw_dy)
            if (allocated(sw_z_cb)) deallocate (sw_z_cb, sw_z_cc, sw_dz)

        end subroutine s_finalize_amr_module

    end module m_amr
