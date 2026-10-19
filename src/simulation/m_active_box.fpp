!>
!!@file
!!@brief Contains module m_active_box

#:include 'macros.fpp'

#! Under CUDA graphs the scan runs asynchronously on the compute queue, and the deferred apply reads its row partials after an event
#! recorded behind it (mfc_scan_record/mfc_scan_wait).
#:set AB_ASYNC_SCAN = (defined('MFC_CUDA_GRAPHS') and MFC_CUDA_GRAPHS)
#:set AB_SCAN_QUEUE = (ACC_QUEUE if ACC_QUEUE is not None else 0)

!> @brief Per-rank active-box restriction of the fused-kernel compute window.
!!
!! Each rank restricts its fused-kernel compute to a window (local indices) around the cells of its subdomain that deviate from its
!! far-corner state, while the halo exchange stays full size. Skipped cells keep their values, so the solution agrees with an
!! unwindowed run to roundoff. The deviation test is relative, with a floor at a fraction of the largest ambient component so the
!! roundoff drift of computed quiescent cells does not register as activity.
!!
!! Every <steps> steps the window is re-derived from a scan of the window grown by one margin, or of the whole subdomain when the
!! window is empty or activity escapes the scanned region. The margin covers the stencil plus the largest front travel until the
!! next update (an explicit growth can only enlarge it). Bounds snap outward to multiples of ab_snap, grow at once, and tighten
!! only to a target repeated on two consecutive updates.
!!
!! When a window reaches a decomposition seam, the ranks across it add a band on their side (s_propagate_active_box_seams), so their
!! compute already runs when the front crosses. A periodic axis held by one rank is a seam with itself, and every physical face
!! whose boundary condition can move the cells next to it is a permanent band. A rank that has measured no activity keeps scanning
!! its full subdomain.
!!
!! The update runs in two phases, the scan on the cadence step and the apply (window derivation, seam exchange, trim and NVSHMEM
!! reconfigure) ab_defer steps later. Without CUDA graphs the scan is a synchronous reduction. With them each row writes its partial
!! box to page-locked memory and the apply waits on an event behind the scan, so the host never drains the queue. The apply uses two
!! allgathers and no other collective.
#ifdef MFC_ACTIVE_BOX
module m_active_box

    use iso_c_binding
    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_mpi_proxy
    use m_constants, only: BC_PERIODIC, BC_GHOST_EXTRAP, BC_RIEMANN_EXTRAP

    implicit none

    #:if AB_ASYNC_SCAN
        interface
            subroutine mfc_scan_record(queue) bind(c)

                import :: c_long_long
                integer(c_long_long), value :: queue

            end subroutine mfc_scan_record
            function mfc_scan_wait(queue) bind(c)

                import :: c_int, c_long_long
                integer(c_long_long), value :: queue
                integer(c_int)              :: mfc_scan_wait

            end function mfc_scan_wait
            function mfc_scan_part_init(lanes) bind(c)

                import :: c_int, c_long_long
                integer(c_long_long), value :: lanes
                integer(c_int)              :: mfc_scan_part_init

            end function mfc_scan_part_init
            function mfc_scan_part_ptr() bind(c)

                import :: c_ptr
                type(c_ptr) :: mfc_scan_part_ptr

            end function mfc_scan_part_ptr
        end interface

        !> Compute queue of the async scan.
        integer(c_long_long), parameter :: ab_scan_queue = ${AB_SCAN_QUEUE}$_c_long_long
        !> Per-row partial boxes, page-locked and device-mapped by the graph backend (mfc_scan_part_init). Slots 1/2 and 3/4 hold
        !! the x-extents of the row's first and last deviating run; a row with three or more runs stores its hull in both.
        integer(c_int), pointer :: ab_part(:,:,:) => null()
        !> The staged scan's cell range (x/y/z beg/end) whose rows the partials cover.
        integer :: ab_scan_rng(1:6) = 0
    #:endif

    private
    public :: s_initialize_active_box_module, s_finalize_active_box_module, s_update_active_box, s_flush_active_box, ab_x, ab_y, &
        & ab_z, ab_active, ab_nb, ab_box, ab_w, ab_cfl_recompute, ab_cfl_valid, ab_cfl_dt_out, ab_inner

    type(int_bounds_info) :: ab_x, ab_y, ab_z  !< Active compute window, local cell indices
    logical :: ab_active = .false.  !< Compute-window optimization engaged (this rank, or any rank via the cadence exchange)
    logical :: ab_local_engaged = .false.  !< This rank's own window has committed non-empty
    logical :: ab_measured = .false.  !< A window has been derived from the solution
    logical :: ab_settled = .false.  !< No rank's window can change again; updates stopped
    integer :: ab_steps = -1  !< Recompute cadence in steps (0 = fixed from ICs, -1 = disabled)
    integer :: ab_growth = 0  !< Window growth in cells, stencil margin plus front travel between updates
    integer :: ab_vstep = 2  !< Bound on a front's travel per step, in cells
    integer :: ab_travel = 0  !< Running max of the raw box's per-cadence edge travel, in cells
    logical :: ab_adaptive = .false.  !< The growth margin tracks the measured travel (no explicit growth)
    integer :: ab_raw(1:6) = 0  !< Raw (ungrown) box of the last measurement, x/y/z beg/end
    logical :: ab_raw_valid = .false.  !< ab_raw holds a measured box
    logical :: ab_seam_band = .false.  !< The last propagation planted a seam band
    logical :: ab_drive(1:3,1:2) = .false.  !< Physical face whose boundary condition can move frozen cells next to it
    integer :: ab_tprev(1:6) = 0  !< Previous cadence's re-derived window bounds (shrink hysteresis)
    logical :: ab_tprev_valid = .false.  !< ab_tprev holds a cadence's target
    integer :: ab_scan(1:6) = 0  !< Staged raw box from the scan phase, awaiting the apply
    logical :: ab_scan_first = .false.  !< The staged scan was the first measurement
    logical :: ab_scan_static = .false.  !< local_static at the staged scan (no scan ran)
    integer :: ab_scan_step = 0  !< The step the staged scan ran on
    logical :: ab_pending = .false.  !< An apply phase is owed at a later step's top
    integer :: ab_defer = 0  !< Steps between a scan and its deferred apply

    !> Multibox compute cover, up to ab_maxbox boxes inside the window that together contain every deviating cell. The fused-kernel
    !! column loops skip columns outside every box. The halo exchange, seam propagation, NVSHMEM trim and graph signature key on the
    !! bounding window. ab_nb = 1 with ab_box(:,1) equal to the window is the single-box case.
    integer, parameter :: AB_MAX_BOXES = 4
    integer            :: ab_maxbox = 1                   !< Box cap from the environment (1 = single box)
    integer            :: ab_nb = 1                       !< Committed box count (1..AB_MAX_BOXES)
    integer            :: ab_box(1:6,1:AB_MAX_BOXES) = 0  !< Committed boxes, x/y/z beg/end per box
    !> Device mirror of the window bounds (x/y/z beg/end) for device-side masks such as the CFL complement pass.
    integer :: ab_w(1:6) = 0
    !> Hull of the seam bands planted at the last apply (inverted when none), kept in the inner cover.
    integer :: ab_band(1:6) = (/0, -1, 0, -1, 0, -1/)

    !> MFC_ACTIVE_BOX_INNER=<k> rebuilds the single-box cover on the device every k steps from a scan of the window, with a margin
    !! of the stencil plus k steps of front travel, no snap, and the seam bands added. The window, launch bounds and graph signature
    !! follow the host cadence. Requires single-box mode.
    integer              :: ab_inner = 0     !< Rebuild cadence in steps (0 = off); device-resident (column-skip gate)
    integer              :: ab_imargin = 0   !< Cover margin, stencil plus one rebuild interval of travel
    integer, allocatable :: ab_ipart(:,:,:)  !< Per-row x-extent partials of the inner-cover scan
    integer, allocatable :: ab_iz(:,:)       !< Per-slice partial boxes of the inner-cover reduction

    !> Windowed CFL. The per-step CFL reduction covers the window, and the frozen cells outside it contribute a constant from a
    !! masked full pass that reruns after every apply.
    logical  :: ab_cfl_recompute = .true.   !< The out-of-window constant is stale (set again at every apply)
    logical  :: ab_cfl_valid = .false.      !< ab_cfl_dt_out holds a global minimum
    real(wp) :: ab_cfl_dt_out(1:4) = 0._wp  !< Global min dt candidates (ICFL/VCFL/CCFL/TCFL) outside the window

    !> Device resident so the column skip reads the current cover at every graph replay without a recapture.
    $:GPU_DECLARE(create='[ab_nb, ab_box, ab_w, ab_band, ab_inner, ab_ipart, ab_iz]')

    !> Relative deviation threshold, far below real waves and far above the roundoff drift of computed quiescent cells.
    real(wp), parameter :: tol_ab = 1.e-12_wp

    !> Floor fraction of the deviation threshold. A component with ambient value near zero borrows the largest ambient component's
    !! scale, giving a threshold of tol_ab*ab_floor_frac*|amb|_max. MFC_ACTIVE_BOX_FLOOR raises it for long runs whose transverse
    !! momentum roundoff grows past the default.
    real(wp) :: ab_floor_frac = 1.e-4_wp

    !> Bound quantization in cells, a multiple of every transverse lane stride of the fused column kernels. Coarse snapping lets
    !! most updates leave every bound in place.
    integer, parameter :: ab_snap = 16
    #:set AB_LANES = [MFC_FUSION_TRANSVERSE_LANES if defined('MFC_FUSION_TRANSVERSE_LANES') else 2, &
        & MFC_FUSION_TRANSVERSE_WIDE_LANES if defined('MFC_FUSION_TRANSVERSE_WIDE_LANES') else 2, &
        & MFC_FUSION_COMPACT_LANES if defined('MFC_FUSION_COMPACT_LANES') else 2] &
        & + list(MFC_FUSION_CHEM_LANES if defined('MFC_FUSION_CHEM_LANES') else (1, 2))
    #:assert all(16 % lanes == 0 for lanes in AB_LANES)

contains

    impure subroutine s_initialize_active_box_module

        character(len=64) :: spec
        integer           :: ierr, comma, comma2, i, j, bcode

        ab_active = .false.
        ab_local_engaged = .false.
        ab_measured = .false.
        ab_settled = .false.
        ab_steps = -1
        ab_growth = 0
        ab_vstep = 2
        ab_travel = 0
        ab_adaptive = .false.
        ab_raw_valid = .false.
        ab_tprev_valid = .false.
        ab_pending = .false.
        ab_seam_band = .false.
        ab_drive = .false.
        ab_maxbox = 1
        ab_nb = 1
        ab_box = 0
        ab_box(2, 1) = -1; ab_box(4, 1) = -1; ab_box(6, 1) = -1
        ab_w = 0
        ab_w(1) = m + 1; ab_w(2) = -1; ab_w(3) = n + 1; ab_w(4) = -1; ab_w(5) = p + 1; ab_w(6) = -1
        @:ALLOCATE(ab_ipart(1:2, 0:n, 0:p), ab_iz(1:6, 0:p))
        ab_band = (/m + 1, -1, n + 1, -1, p + 1, -1/)
        $:GPU_UPDATE(device='[ab_nb, ab_box, ab_w, ab_band, ab_inner]')
        ! Empty until the first measurement; the first window is the box of deviating cells grown and snapped outward.
        ab_x%beg = m + 1; ab_x%end = -1
        ab_y%beg = n + 1; ab_y%end = -1
        ab_z%beg = p + 1; ab_z%end = -1

        call get_environment_variable('MFC_ACTIVE_BOX', value=spec, status=ierr)
        ! Unset or empty disables the window.
        if (ierr /= 0 .or. len_trim(spec) == 0) return
#ifndef MFC_FUSED_RHS
        call s_mpi_abort('MFC_ACTIVE_BOX is set but this build has no fused RHS kernel (MFC_FUSED_RHS) to window.')
#endif

        ! Programmed burn and GRCBC boundaries change cells the deviation test sees as quiescent, and cylindrical geometric
        ! sources move a uniform flowing state, so they disable the window.
        if (cyl_coord .or. reactive_burn .or. bc_x%grcbc_in .or. bc_x%grcbc_out .or. bc_x%grcbc_vel_out .or. bc_y%grcbc_in &
            & .or. bc_y%grcbc_out .or. bc_y%grcbc_vel_out .or. bc_z%grcbc_in .or. bc_z%grcbc_out .or. bc_z%grcbc_vel_out) return

        ! Ghost extrapolation reproduces the frozen boundary cells, so it never moves them. Walls act only on flow that
        ! deviates before it arrives, so the scan covers them. Characteristic, Dirichlet, axis and patch boundaries drive
        ! the flow themselves, and their faces stay full-width bands of the window.
        do i = 1, num_dims
            do j = 1, 2
                if (i == 1) bcode = merge(bc_x%beg, bc_x%end, j == 1)
                if (i == 2) bcode = merge(bc_y%beg, bc_y%end, j == 1)
                if (i == 3) bcode = merge(bc_z%beg, bc_z%end, j == 1)
                ab_drive(i, &
                         & j) = bcode < 0 .and. bcode /= BC_PERIODIC .and. (num_bc_patches > 0 .or. (bcode /= BC_GHOST_EXTRAP &
                         & .and. bcode /= BC_RIEMANN_EXTRAP .and. bcode /= BC_REFLECTIVE .and. bcode /= BC_SLIP_WALL &
                         & .and. bcode /= BC_NO_SLIP_WALL))
            end do
        end do

        comma = index(spec, ',')
        if (comma > 0) then
            read (spec(1:comma - 1), *, iostat=ierr) ab_steps
            if (ierr == 0) then
                comma2 = index(spec(comma + 1:), ',')
                if (comma2 > 0) then
                    ! An empty growth field ("20,,4") keeps the default.
                    if (comma2 > 1) read (spec(comma + 1:comma + comma2 - 1), *, iostat=ierr) ab_growth
                    if (ierr == 0 .and. len_trim(spec(comma + comma2 + 1:)) > 0) read (spec(comma + comma2 + 1:), *, &
                        & iostat=ierr) ab_maxbox
                else if (len_trim(spec(comma + 1:)) > 0) then
                    read (spec(comma + 1:), *, iostat=ierr) ab_growth
                end if
            end if
        else
            read (spec, *, iostat=ierr) ab_steps
        end if
        if (ierr /= 0 .or. ab_steps < 0 .or. ab_growth < 0 .or. ab_maxbox < 1) then
            call s_mpi_abort('MFC_ACTIVE_BOX: expected <steps>[,<growth>[,<boxes>]] with nonnegative integers.')
        end if
        ab_maxbox = min(ab_maxbox, AB_MAX_BOXES)

        ! Optional deviation floor override (see ab_floor_frac).
        call get_environment_variable('MFC_ACTIVE_BOX_FLOOR', value=spec, status=ierr)
        if (ierr == 0 .and. len_trim(spec) > 0) then
            read (spec, *, iostat=ierr) ab_floor_frac
            if (ierr /= 0 .or. ab_floor_frac <= 0._wp .or. ab_floor_frac >= 1._wp) then
                call s_mpi_abort('MFC_ACTIVE_BOX_FLOOR: expected a real fraction in (0,1).')
            end if
        end if

        ! Inner cover rebuilt on the device every <k> steps (see ab_inner).
        call get_environment_variable('MFC_ACTIVE_BOX_INNER', value=spec, status=ierr)
        if (ierr == 0 .and. len_trim(spec) > 0) then
            read (spec, *, iostat=ierr) ab_inner
            if (ierr /= 0 .or. ab_inner < 1) call s_mpi_abort('MFC_ACTIVE_BOX_INNER: expected a positive step count.')
            if (ab_maxbox > 1) call s_mpi_abort('MFC_ACTIVE_BOX_INNER: single-box mode only (third MFC_ACTIVE_BOX field 1).')
            $:GPU_UPDATE(device='[ab_inner]')
        end if

        #:if AB_ASYNC_SCAN
            if (mfc_scan_part_init(int(4*(n + 1)*(p + 1), &
                & c_long_long)) /= 0) call s_mpi_abort('MFC: active-box scan partial buffer unavailable.')
            call c_f_pointer(mfc_scan_part_ptr(), ab_part, (/4, n + 1, p + 1/))
        #:endif

        ! Deferred-apply depth of half a cadence, leaving the host the other half to run ahead before the next scan.
        if (ab_steps > 0) ab_defer = min(ab_steps - 1, max(1, ab_steps/2))

        ! Growth is the stencil plus the front travel over a window's exposure (one cadence plus the deferral, or the whole run for
        ! a fixed window). Without an explicit growth the margin starts at one cell per step and then tracks the measured travel,
        ! so a slow front is not padded to the CFL-1 bound all run. An explicit growth replaces the default.
        if (cfl_dt .and. cfl_target > 0._wp) ab_vstep = int(cfl_target) + 1
        if (ab_steps == 0) then
            if (cfl_dt) call s_mpi_abort('MFC_ACTIVE_BOX: a fixed window (<steps> = 0) needs a fixed step count; use a cadence.')
            ab_growth = max(ab_growth, buff_size + ab_vstep*(t_step_stop - t_step_start))
        else if (ab_growth == 0) then
            ab_growth = buff_size + ab_steps + ab_defer
            if (cfl_dt .and. cfl_target > 0._wp) ab_growth = buff_size + int(cfl_target*real(ab_steps + ab_defer, wp)) + 1
            ab_adaptive = .true.
        end if
        ab_imargin = buff_size + ab_vstep*ab_inner

    end subroutine s_initialize_active_box_module

    impure subroutine s_finalize_active_box_module

        ab_active = .false.
        ab_local_engaged = .false.
        ab_measured = .false.
        ab_settled = .false.
        ab_steps = -1
        ab_travel = 0
        ab_adaptive = .false.
        ab_raw_valid = .false.
        ab_tprev_valid = .false.
        ab_pending = .false.
        ab_maxbox = 1
        ab_nb = 1
        ab_box = 0
        ab_inner = 0
        if (allocated(ab_ipart)) then
            @:DEALLOCATE(ab_ipart)
        end if
        if (allocated(ab_iz)) then
            @:DEALLOCATE(ab_iz)
        end if
        #:if AB_ASYNC_SCAN
            ! The partial buffer is owned and freed by the graph backend.
            if (associated(ab_part)) nullify (ab_part)
        #:endif

    end subroutine s_finalize_active_box_module

    !> Refresh the per-rank active window. Called once per time step before the stage loop and once at startup with the initial
    !! condition. <steps> = 0 keeps the window fixed at its initial value. With eligible false (the fused final-RK path is off) the
    !! window never engages.
    !!
    !! The scan runs on the cadence step and its apply (s_apply_active_box) ab_defer steps later, so the apply's collectives overlap
    !! the device work enqueued in between.
    impure subroutine s_update_active_box(q_cons_vf, t_step, eligible)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer, intent(in)                                 :: t_step
        logical, intent(in)                                 :: eligible
        integer                                             :: slx, shx, sly, shy, slz, shz
        logical                                             :: full, first, local_static

        if (ab_steps < 0 .or. .not. eligible) return
        if (ab_settled) return

        ! Apply the previous cadence's scan once the deferral has elapsed.
        if (ab_pending .and. t_step >= ab_scan_step + ab_defer) call s_flush_active_box(q_cons_vf)

        ! Rebuild the inner cover (MFC_ACTIVE_BOX_INNER) after any deferred apply, so a cadence's cover publish cannot overwrite it.
        if (ab_inner > 0 .and. ab_active .and. .not. ab_settled .and. ab_x%beg <= ab_x%end .and. mod(t_step - t_step_start, &
            & ab_inner) == 0) then
            call s_rebuild_inner_cover(q_cons_vf)
        end if

        if (ab_measured .and. (ab_steps == 0 .or. t_step == t_step_start)) return
        if (t_step > t_step_start .and. (ab_steps == 0 .or. mod(t_step - t_step_start, ab_steps) /= 0)) return

        ! A rank whose window spans its subdomain, whose window is fixed, or that runs alone with an empty window cannot change, so
        ! it skips the scan but still joins the exchange at the apply. With several ranks an empty-window rank keeps scanning, since
        ! activity reaches it through seam propagation.
        local_static = ab_measured .and. (ab_steps == 0 .or. (num_procs == 1 .and. ab_x%end < ab_x%beg) .or. (ab_x%beg == 0 &
                                          & .and. ab_x%end == m .and. ab_y%beg == 0 .and. ab_y%end == n .and. ab_z%beg == 0 &
                                          & .and. ab_z%end == p))

        first = .not. ab_measured
        if (.not. local_static) then
            ! Bounding box of the cells deviating from the far-corner state. The first measurement and an empty window scan the
            ! whole subdomain; otherwise the scan covers the window grown by one margin, and activity beyond the window triggers a
            ! full remeasure.
            full = first .or. ab_x%end < ab_x%beg
            if (full) then
                slx = 0; shx = m; sly = 0; shy = n; slz = 0; shz = p
            else
                slx = max(0, ab_x%beg - ab_growth); shx = min(m, ab_x%end + ab_growth)
                sly = max(0, ab_y%beg - ab_growth); shy = min(n, ab_y%end + ab_growth)
                slz = max(0, ab_z%beg - ab_growth); shz = min(p, ab_z%end + ab_growth)
            end if

            #:if AB_ASYNC_SCAN
                ! Each (k, l) row writes its partial box to page-locked memory, and an event recorded behind the region marks when
                ! the partials are visible. The escape check needs the box, so any full remeasure happens at the apply.
                ab_scan_rng = (/slx, shx, sly, shy, slz, shz/)
                call s_scan_active_box(q_cons_vf, slx, shx, sly, shy, slz, shz)
                call mfc_scan_record(ab_scan_queue)
            #:else
                do
                    call s_scan_active_box(q_cons_vf, slx, shx, sly, shy, slz, shz, ab_scan(1), ab_scan(2), ab_scan(3), &
                                           & ab_scan(4), ab_scan(5), ab_scan(6))
                    if (full) exit
                    if (ab_scan(1) >= ab_x%beg .and. ab_scan(2) <= ab_x%end .and. ab_scan(3) >= ab_y%beg .and. ab_scan(4) &
                        & <= ab_y%end .and. ab_scan(5) >= ab_z%beg .and. ab_scan(6) <= ab_z%end) exit
                    full = .true.
                    slx = 0; shx = m; sly = 0; shy = n; slz = 0; shz = p
                end do
            #:endif

            ! Activation latches in s_apply_active_box on the first committed non-empty window, since an empty window with the
            ! consumers engaged would freeze an at-rest initial condition.
            ab_measured = .true.
        end if

        ! Stage the measurement for the apply. A static rank stages no box; under the async scan the box arrives at the apply's
        ! event wait.
        ab_scan_first = first
        ab_scan_static = local_static
        ab_scan_step = t_step
        ab_pending = .true.

    end subroutine s_update_active_box

    #:if AB_ASYNC_SCAN
        !> Scan the given cell range for cells deviating from the far-corner state. Each (k, l) row stores the x-extents of its
        !! first and last deviating run in ab_part (inverted when none deviates); a row with three or more runs stores its hull in
        !! both slots. The rows run asynchronously on the compute queue, and s_combine_active_box_scan reads them after the event
        !! recorded behind the region fires.
        impure subroutine s_scan_active_box(q_cons_vf, slx, shx, sly, shy, slz, shz)

            type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
            integer, intent(in)                                 :: slx, shx, sly, shy, slz, shz
            integer                                             :: j, k, l, i
            integer                                             :: jb1, je1, jb2, je2, runs, curb
            logical                                             :: deviates
            real(wp)                                            :: ambient, val, amb_scale

            $:GPU_PARALLEL_LOOP(collapse=2, private='[j, i, jb1, je1, jb2, je2, runs, curb, deviates, ambient, val, amb_scale]')
            do l = slz, shz
                do k = sly, shy
                    jb1 = m + 1; je1 = -1; jb2 = m + 1; je2 = -1
                    runs = 0; curb = -1
                    do j = slx, shx
                        deviates = .false.
                        amb_scale = 0._wp
                        do i = 1, sys_size
                            amb_scale = max(amb_scale, abs(real(q_cons_vf(i)%sf(m, n, p), wp)))
                        end do
                        do i = 1, sys_size
                            ambient = real(q_cons_vf(i)%sf(m, n, p), wp)
                            val = real(q_cons_vf(i)%sf(j, k, l), wp)
                            if (abs(val - ambient) > tol_ab*max(abs(ambient), ab_floor_frac*amb_scale)) then
                                deviates = .true.
                                exit
                            end if
                        end do
                        if (deviates) then
                            if (curb < 0) then
                                curb = j
                                runs = runs + 1
                                if (runs == 1) jb1 = j
                            end if
                            jb2 = curb; je2 = j
                        else if (curb >= 0) then
                            if (runs == 1) je1 = j - 1
                            curb = -1
                        end if
                    end do
                    if (curb >= 0 .and. runs == 1) je1 = je2
                    if (runs > 2) then
                        je1 = je2
                        jb2 = jb1
                    end if
                    ab_part(1, k + 1, l + 1) = jb1
                    ab_part(2, k + 1, l + 1) = je1
                    ab_part(3, k + 1, l + 1) = jb2
                    ab_part(4, k + 1, l + 1) = je2
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_scan_active_box

        !> Combine the row partials over the given row range into the bounding box <jb>..<le> (inverted when nothing deviates), on
        !! the host after the scan's event has fired.
        subroutine s_combine_active_box_scan(sly, shy, slz, shz, jb, je, kb, ke, lb, le)

            integer, intent(in)  :: sly, shy, slz, shz
            integer, intent(out) :: jb, je, kb, ke, lb, le
            integer              :: k, l

            jb = m + 1; je = -1
            kb = n + 1; ke = -1
            lb = p + 1; le = -1
            do l = slz, shz
                do k = sly, shy
                    if (ab_part(2, k + 1, l + 1) >= ab_part(1, k + 1, l + 1)) then
                        ! The bounding box spans both runs of the row.
                        jb = min(jb, ab_part(1, k + 1, l + 1)); je = max(je, ab_part(2, k + 1, l + 1))
                        jb = min(jb, ab_part(3, k + 1, l + 1)); je = max(je, ab_part(4, k + 1, l + 1))
                        kb = min(kb, k); ke = max(ke, k)
                        lb = min(lb, l); le = max(le, l)
                    end if
                end do
            end do

        end subroutine s_combine_active_box_scan
    #:else
        !> Synchronous scan reducing the bounding box of the cells deviating from the far-corner state over the given cell range
        !! into <jb>..<le> (inverted when nothing deviates).
        impure subroutine s_scan_active_box(q_cons_vf, slx, shx, sly, shy, slz, shz, jb, je, kb, ke, lb, le)

            type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
            integer, intent(in)                                 :: slx, shx, sly, shy, slz, shz
            integer, intent(out)                                :: jb, je, kb, ke, lb, le
            integer                                             :: j, k, l, i
            logical                                             :: deviates
            real(wp)                                            :: ambient, val, amb_scale

            jb = m + 1; je = -1
            kb = n + 1; ke = -1
            lb = p + 1; le = -1

            $:GPU_PARALLEL_LOOP(collapse=3, private='[i, deviates, ambient, val, amb_scale]', reduction='[[jb, kb, lb], [je, ke, &
                                & le]]', reductionOp='[min, max]')
            do l = slz, shz
                do k = sly, shy
                    do j = slx, shx
                        deviates = .false.
                        amb_scale = 0._wp
                        do i = 1, sys_size
                            amb_scale = max(amb_scale, abs(real(q_cons_vf(i)%sf(m, n, p), wp)))
                        end do
                        do i = 1, sys_size
                            ambient = real(q_cons_vf(i)%sf(m, n, p), wp)
                            val = real(q_cons_vf(i)%sf(j, k, l), wp)
                            if (abs(val - ambient) > tol_ab*max(abs(ambient), ab_floor_frac*amb_scale)) then
                                deviates = .true.
                                exit
                            end if
                        end do
                        if (deviates) then
                            jb = min(jb, j); je = max(je, j)
                            kb = min(kb, k); ke = max(ke, k)
                            lb = min(lb, l); le = max(le, l)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()

        end subroutine s_scan_active_box
    #:endif

    !> Rebuild the single-box cover on the device. A row scan of the window fills ab_ipart, a per-slice pass reduces the window rows
    !! into ab_iz, and a single-gang pass grows the box by the inner margin, adds the seam bands, clamps it into the window and
    !! stores it in ab_box.
    impure subroutine s_rebuild_inner_cover(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer                                             :: j, k, l, i, dummy
        integer                                             :: jb, je, kb, ke, lb, le
        integer                                             :: jlo, jhi, klo, khi, llo, lhi, margin
        logical                                             :: deviates
        real(wp)                                            :: ambient, val, amb_scale

        jlo = max(0, ab_x%beg); jhi = min(m, ab_x%end)
        klo = max(0, ab_y%beg); khi = min(n, ab_y%end)
        llo = max(0, ab_z%beg); lhi = min(p, ab_z%end)
        margin = ab_imargin

        $:GPU_PARALLEL_LOOP(collapse=2, private='[j, i, jb, je, deviates, ambient, val, amb_scale]')
        do l = llo, lhi
            do k = klo, khi
                jb = m + 1; je = -1
                do j = jlo, jhi
                    deviates = .false.
                    amb_scale = 0._wp
                    do i = 1, sys_size
                        amb_scale = max(amb_scale, abs(real(q_cons_vf(i)%sf(m, n, p), wp)))
                    end do
                    do i = 1, sys_size
                        ambient = real(q_cons_vf(i)%sf(m, n, p), wp)
                        val = real(q_cons_vf(i)%sf(j, k, l), wp)
                        if (abs(val - ambient) > tol_ab*max(abs(ambient), ab_floor_frac*amb_scale)) then
                            deviates = .true.
                            exit
                        end if
                    end do
                    if (deviates) then
                        jb = min(jb, j); je = max(je, j)
                    end if
                end do
                ab_ipart(1, k + 1, l + 1) = jb
                ab_ipart(2, k + 1, l + 1) = je
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(private='[k, jb, je, kb, ke]', firstprivate='[klo, khi]')
        do l = llo, lhi
            jb = m + 1; je = -1; kb = n + 1; ke = -1
            do k = klo, khi
                if (ab_ipart(2, k + 1, l + 1) >= ab_ipart(1, k + 1, l + 1)) then
                    jb = min(jb, ab_ipart(1, k + 1, l + 1)); je = max(je, ab_ipart(2, k + 1, l + 1))
                    kb = min(kb, k); ke = max(ke, k)
                end if
            end do
            ab_iz(1, l) = jb; ab_iz(2, l) = je; ab_iz(3, l) = kb; ab_iz(4, l) = ke
            ab_iz(5, l) = merge(l, p + 1, je >= jb); ab_iz(6, l) = merge(l, -1, je >= jb)
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(parallelism='[gang]', private='[l, jb, je, kb, ke, lb, le]', firstprivate='[jlo, jhi, klo, khi, llo, &
                            & lhi, margin]')
        do dummy = 1, 1
            jb = m + 1; je = -1; kb = n + 1; ke = -1; lb = p + 1; le = -1
            $:GPU_LOOP(parallelism='[seq]')
            do l = llo, lhi
                if (ab_iz(2, l) >= ab_iz(1, l)) then
                    jb = min(jb, ab_iz(1, l)); je = max(je, ab_iz(2, l))
                    kb = min(kb, ab_iz(3, l)); ke = max(ke, ab_iz(4, l))
                    lb = min(lb, ab_iz(5, l)); le = max(le, ab_iz(6, l))
                end if
            end do
            ! An empty measurement keeps the current cover.
            if (je >= jb) then
                ab_box(1, 1) = max(jlo, jb - margin)
                ab_box(2, 1) = min(jhi, je + margin)
                ab_box(3, 1) = max(klo, kb - margin)
                ab_box(4, 1) = min(khi, ke + margin)
                ab_box(5, 1) = max(llo, lb - margin)
                ab_box(6, 1) = min(lhi, le + margin)
                ! The seam bands lie inside the window and carry activity no local cell shows yet.
                if (ab_band(2) >= ab_band(1)) then
                    ab_box(1, 1) = max(jlo, min(ab_box(1, 1), ab_band(1)))
                    ab_box(2, 1) = min(jhi, max(ab_box(2, 1), ab_band(2)))
                    ab_box(3, 1) = max(klo, min(ab_box(3, 1), ab_band(3)))
                    ab_box(4, 1) = min(khi, max(ab_box(4, 1), ab_band(4)))
                    ab_box(5, 1) = max(llo, min(ab_box(5, 1), ab_band(5)))
                    ab_box(6, 1) = min(lhi, max(ab_box(6, 1), ab_band(6)))
                end if
            end if
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_rebuild_inner_cover

    !> Apply a staged measurement. Re-derive the window, exchange it, propagate seam bands and re-trim the exchange. Under the async
    !! scan, first wait on the scan's event, combine the partials and remeasure the whole subdomain when the box escaped the scan
    !! region. The result depends only on the staged values and the module state.
    impure subroutine s_apply_active_box(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf
        integer                                             :: ierr
        integer                                             :: jb, je, kb, ke, lb, le
        integer                                             :: flags(2)
        integer                                             :: oxb, oxe, oyb, oye, ozb, oze, tprev(1:6)
        logical                                             :: tprev_valid
        logical                                             :: changed, clustered, any_engaged

        #:if AB_ASYNC_SCAN
            logical :: remeas
        #:endif
        type(int_bounds_info) :: wx, wy, wz
        integer               :: mine(1:18), share(1:18,1:num_procs)
        integer               :: share2(1:21,1:num_procs)

        if (.not. ab_pending) return
        ab_pending = .false.

        #:if AB_ASYNC_SCAN
            ! Wait on the event recorded behind the staged scan, not on the whole queue.
            remeas = .false.
            if (mfc_scan_wait(ab_scan_queue) /= 0) call s_mpi_abort('MFC: active-box scan wait failed.')
            call s_combine_active_box_scan(ab_scan_rng(3), ab_scan_rng(4), ab_scan_rng(5), ab_scan_rng(6), ab_scan(1), &
                                           & ab_scan(2), ab_scan(3), ab_scan(4), ab_scan(5), ab_scan(6))
            ! A non-first measurement whose box escaped the scan region remeasures the whole subdomain.
            if (.not. ab_scan_first .and. ab_x%beg <= ab_x%end .and. .not. (ab_scan(1) >= ab_x%beg .and. ab_scan(2) <= ab_x%end &
                & .and. ab_scan(3) >= ab_y%beg .and. ab_scan(4) <= ab_y%end .and. ab_scan(5) >= ab_z%beg .and. ab_scan(6) &
                & <= ab_z%end)) then
                call s_scan_active_box(q_cons_vf, 0, m, 0, n, 0, p)
                call mfc_scan_record(ab_scan_queue)
                if (mfc_scan_wait(ab_scan_queue) /= 0) call s_mpi_abort('MFC: active-box scan wait failed.')
                call s_combine_active_box_scan(0, n, 0, p, ab_scan(1), ab_scan(2), ab_scan(3), ab_scan(4), ab_scan(5), ab_scan(6))
                remeas = .true.
            end if
        #:endif

        ! Adapt the margin to the measured per-cadence edge travel, scaled to the window's exposure. The running max keeps it
        ! monotone, so an accelerating front widens the margin and a slower one does not narrow it below a travel already seen.
        if (.not. ab_scan_static .and. ab_raw_valid .and. ab_scan(2) >= ab_scan(1)) then
            ab_travel = max(ab_travel, abs(ab_scan(1) - ab_raw(1)), abs(ab_scan(2) - ab_raw(2)), abs(ab_scan(3) - ab_raw(3)), &
                            & abs(ab_scan(4) - ab_raw(4)), abs(ab_scan(5) - ab_raw(5)), abs(ab_scan(6) - ab_raw(6)))
            if (ab_adaptive) ab_growth = max(buff_size + 1, &
                & buff_size + (ab_travel*(ab_steps + ab_defer) + ab_steps - 1)/ab_steps + ab_snap/2)
        end if
        if (.not. ab_scan_static .and. ab_scan(2) >= ab_scan(1)) then
            ab_raw = ab_scan
            ab_raw_valid = .true.
        end if

        oxb = ab_x%beg; oxe = ab_x%end; oyb = ab_y%beg; oye = ab_y%end; ozb = ab_z%beg; oze = ab_z%end
        tprev = ab_tprev; tprev_valid = ab_tprev_valid
        jb = ab_scan(1); je = ab_scan(2); kb = ab_scan(3); ke = ab_scan(4); lb = ab_scan(5); le = ab_scan(6)

        if (.not. ab_scan_static) then
            if (je >= jb) then
                ! Re-derive the window as the measured box grown by the margin and snapped outward to ab_snap. Column kernels launch
                ! whole lane blocks, so an unaligned upper bound would let a partial block update cells outside the window.
                wx%beg = ab_snap*(max(0, jb - ab_growth)/ab_snap)
                wx%end = min(m, ab_snap*((min(m, je + ab_growth) + ab_snap)/ab_snap) - 1)
                wy%beg = ab_snap*(max(0, kb - ab_growth)/ab_snap)
                wy%end = min(n, ab_snap*((min(n, ke + ab_growth) + ab_snap)/ab_snap) - 1)
                wz%beg = ab_snap*(max(0, lb - ab_growth)/ab_snap)
                wz%end = min(p, ab_snap*((min(p, le + ab_growth) + ab_snap)/ab_snap) - 1)
                ab_x = wx; ab_y = wy; ab_z = wz
            else
                ! Nothing deviates in the scanned region, so the window is empty. Seam bands and the full-subdomain scan of the next
                ! cadence catch new activity.
                ab_x%beg = m + 1; ab_x%end = -1
                ab_y%beg = n + 1; ab_y%end = -1
                ab_z%beg = p + 1; ab_z%end = -1
            end if
        end if

        ! Exchange the windows every cadence so seam bands follow the neighbors' current activity.
        ab_seam_band = .false.
        ab_band = (/m + 1, -1, n + 1, -1, p + 1, -1/)
        mine = f_active_box_share()
#ifdef MFC_MPI
        if (num_procs > 1) then
            call MPI_ALLGATHER(mine, 18, MPI_INTEGER, share, 18, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        else
            share(:,1) = mine
        end if
#else
        share(:,1) = mine
#endif
        call s_propagate_active_box_seams(share)

        ! Multibox cover from the staged partials (async scan builds). A static or empty measurement keeps the single window, and
        ! any seam band falls back to the window, since a band carries activity no local row measured.
        clustered = .false.
        #:if AB_ASYNC_SCAN
            if (ab_maxbox > 1 .and. .not. ab_scan_static .and. ab_x%beg <= ab_x%end) then
                if (remeas) then
                    call s_derive_active_box_cover(0, n, 0, p, ab_seam_band)
                else
                    call s_derive_active_box_cover(ab_scan_rng(3), ab_scan_rng(4), ab_scan_rng(5), ab_scan_rng(6), ab_seam_band)
                end if
                clustered = .true.
            end if
        #:endif

        ! Store this cadence's target, then commit it against the entry window. Growth applies at once; a tighten applies to a
        ! target repeated from the previous cadence, so a box oscillating across a snap boundary does not move a bound every
        ! cadence.
        ab_tprev = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
        ab_tprev_valid = .true.
        ! An empty entry or target window commits as is.
        if (oxb <= oxe .and. ab_x%beg <= ab_x%end) then
            ab_x%beg = -f_ab_bound(-ab_x%beg, -oxb, -tprev(1), tprev_valid)
            ab_x%end = f_ab_bound(ab_x%end, oxe, tprev(2), tprev_valid)
            ab_y%beg = -f_ab_bound(-ab_y%beg, -oyb, -tprev(3), tprev_valid)
            ab_y%end = f_ab_bound(ab_y%end, oye, tprev(4), tprev_valid)
            ab_z%beg = -f_ab_bound(-ab_z%beg, -ozb, -tprev(5), tprev_valid)
            ab_z%end = f_ab_bound(ab_z%end, oze, tprev(6), tprev_valid)
        end if

        ! An unclustered cover is the committed window. A clustered cover derives from the pre-commit target, which the commit can
        ! only widen.
        if (.not. clustered) then
            ab_nb = 1
            ab_box(:,1) = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
        end if

        ! Activation latches on the first committed non-empty window of any rank (slot 21 of the exchange). Until then every rank
        ! computes its full subdomain; afterwards an empty-window rank computes nothing until a seam band reaches it.
        if (ab_x%beg <= ab_x%end .and. ab_y%beg <= ab_y%end .and. ab_z%beg <= ab_z%end) ab_local_engaged = .true.

        ! One exchange carries whether any window moved and whether any rank still scans. With nothing moved the trim and NVSHMEM
        ! reconfigure are skipped; with no rank scanning the updates stop. The cover is device data read live, so it does not enter
        ! this flag.
        changed = ab_scan_first .or. ab_x%beg /= oxb .or. ab_x%end /= oxe .or. ab_y%beg /= oyb .or. ab_y%end /= oye &
            & .or. ab_z%beg /= ozb .or. ab_z%end /= oze
#ifdef MFC_MPI
        if (num_procs > 1) then
            ! Re-exchange after propagation and hysteresis so every rank holds every final window for the trimmed exchange.
            mine = f_active_box_share()
            share2(1:18,proc_rank + 1) = mine
            share2(19, proc_rank + 1) = merge(1, 0, changed)
            share2(20, proc_rank + 1) = merge(1, 0, .not. ab_scan_static)
            share2(21, proc_rank + 1) = merge(1, 0, ab_local_engaged)
            call MPI_ALLGATHER(MPI_IN_PLACE, 21, MPI_INTEGER, share2, 21, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            flags(1) = merge(1, 0, any(share2(19,:) == 1))
            flags(2) = merge(1, 0, any(share2(20,:) == 1))
        else
            flags(1) = merge(1, 0, changed)
            flags(2) = merge(1, 0, .not. ab_scan_static)
        end if
#else
        flags(1) = merge(1, 0, changed)
        flags(2) = merge(1, 0, .not. ab_scan_static)
#endif
        any_engaged = ab_local_engaged
#ifdef MFC_MPI
        if (num_procs > 1) any_engaged = any_engaged .or. any(share2(21,:) == 1)
#endif
        if (any_engaged) ab_active = .true.

        if (flags(1) == 1) then
#ifdef MFC_MPI
            if (num_procs > 1) then
#if defined(MFC_NVSHMEM) && defined(__NVCOMPILER) && defined(MFC_HALO_SINGLE_ROUND)
                ! Distribute the windows so the NVSHMEM 26-neighbor exchange trims each partner's put to the sub-box its compute
                ! consumes (s_mpi_apply_26_trim).
                call s_mpi_apply_26_trim(share2)
#endif
            end if
#endif
        end if

        if (flags(1) == 0 .and. flags(2) == 0) ab_settled = .true.

        ! Publish the cover and window bounds to the device, ordered on the compute queue ahead of the next replay that reads them.
        ab_w = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
        if (ab_active) ab_cfl_recompute = .true.
        $:GPU_UPDATE(device='[ab_nb, ab_box, ab_w, ab_band]')

    end subroutine s_apply_active_box

    !> Apply any staged measurement immediately. The startup sequence calls it so the first window is in place before the time loop.
    impure subroutine s_flush_active_box(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf

        call s_apply_active_box(q_cons_vf)

    end subroutine s_flush_active_box

    !> Grow a raw cell range by the margin and snap it outward to ab_snap, the window's own bound arithmetic. The map is monotone,
    !! so a box derived from a sub-range stays inside the window.
    pure subroutine s_ab_grow_snap(lo, hi, lim, olo, ohi)

        integer, intent(in)  :: lo, hi, lim
        integer, intent(out) :: olo, ohi

        olo = ab_snap*(max(0, lo - ab_growth)/ab_snap)
        ohi = min(lim, ab_snap*((min(lim, hi + ab_growth) + ab_snap)/ab_snap) - 1)

    end subroutine s_ab_grow_snap

    #:if AB_ASYNC_SCAN
        !> Derive the multibox cover from the staged row partials over the given row range. A seam band planted this cadence
        !! (seam_widened) or a single-box request keeps the bounding window.
        impure subroutine s_derive_active_box_cover(sly, shy, slz, shz, seam_widened)

            integer, intent(in) :: sly, shy, slz, shz
            logical, intent(in) :: seam_widened

            if (seam_widened .or. ab_maxbox <= 1) then
                ab_nb = 1
                ab_box(:,1) = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
                return
            end if

            call s_derive_active_box_cover_br(sly, shy, slz, shz)

        end subroutine s_derive_active_box_cover

        !> Berger-Rigoutsos style clustering of the row partials. Starting from the activity hull, repeatedly split the leaf with
        !! the most covered inactive volume at the deepest signature dip or strongest inflection of an axis, whichever drops the
        !! most volume once both halves shrink to their activity hulls. The leaves always contain the measured activity, and the
        !! committed boxes grow and snap like the window, so each stays inside it.
        impure subroutine s_derive_active_box_cover_br(sly, shy, slz, shz)

            integer, intent(in)  :: sly, shy, slz, shz
            integer              :: cap, nrows, i, j, k, l, c, d, ia1, ia2, ib1, ib2, act_a, gc, gmin, gmax
            integer              :: nleaf, best, bestw, bgain, bia1, bia2, w, cv, ns, lo, hi, cov
            integer              :: cdip, cinf, cut, t1(1:6), t2(1:6), bt1(1:6), bt2(1:6)
            integer              :: leaf(1:6,AB_MAX_BOXES), lact(AB_MAX_BOXES), seg(1:2,AB_MAX_BOXES)
            logical              :: ldone(AB_MAX_BOXES)
            integer, allocatable :: rk(:), rl(:), ra1(:), ra2(:), rb1(:), rb2(:)
            integer              :: sx(0:m + 1), sy(0:n), sz(0:p)

            ! Active rows of the scan range; a row's two x-intervals merge when they touch or overlap.

            cap = (shy - sly + 1)*(shz - slz + 1)
            allocate (rk(0:cap - 1), rl(0:cap - 1), ra1(0:cap - 1), ra2(0:cap - 1), rb1(0:cap - 1), rb2(0:cap - 1))
            nrows = 0
            do l = slz, shz
                do k = sly, shy
                    ia1 = ab_part(1, k + 1, l + 1); ia2 = ab_part(2, k + 1, l + 1)
                    ib1 = ab_part(3, k + 1, l + 1); ib2 = ab_part(4, k + 1, l + 1)
                    if (ia2 < ia1) cycle
                    if (ib2 >= ib1 .and. ib1 <= ia2 + 1) then
                        ia2 = max(ia2, ib2)
                        ib1 = 1; ib2 = 0
                    end if
                    rk(nrows) = k; rl(nrows) = l
                    ra1(nrows) = ia1; ra2(nrows) = ia2
                    rb1(nrows) = ib1; rb2(nrows) = ib2
                    nrows = nrows + 1
                end do
            end do

            if (nrows == 0) then
                ! No deviating rows, so the window is the cover.
                ab_nb = 1
                ab_box(:,1) = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
                return
            end if

            leaf(:,1) = (/0, m, sly, shy, slz, shz/)
            call s_ab_br_hull(rk, rl, ra1, ra2, rb1, rb2, nrows, leaf(:,1), lact(1))
            nleaf = 1
            ldone = .false.

            do
                ! The not-yet-irreducible leaf with the most covered-but-inactive volume.
                best = 0
                bestw = -1
                do i = 1, nleaf
                    if (ldone(i)) cycle
                    w = f_ab_br_vol(leaf(:,i)) - lact(i)
                    if (w > bestw) then
                        bestw = w
                        best = i
                    end if
                end do
                if (best == 0 .or. nleaf >= ab_maxbox) exit

                ! Axis signatures of the leaf's rows, the active cells per x column (difference array, prefix-summed below) and per
                ! y and z slice.
                sx(leaf(1, best):leaf(2, best) + 1) = 0
                sy(leaf(3, best):leaf(4, best)) = 0
                sz(leaf(5, best):leaf(6, best)) = 0
                do i = 0, nrows - 1
                    if (rk(i) < leaf(3, best) .or. rk(i) > leaf(4, best) .or. rl(i) < leaf(5, best) .or. rl(i) > leaf(6, &
                        & best)) cycle
                    act_a = 0
                    ia1 = max(ra1(i), leaf(1, best)); ia2 = min(ra2(i), leaf(2, best))
                    if (ia2 >= ia1) then
                        sx(ia1) = sx(ia1) + 1
                        sx(ia2 + 1) = sx(ia2 + 1) - 1
                        act_a = ia2 - ia1 + 1
                    end if
                    ib1 = max(rb1(i), leaf(1, best)); ib2 = min(rb2(i), leaf(2, best))
                    if (ib2 >= ib1) then
                        sx(ib1) = sx(ib1) + 1
                        sx(ib2 + 1) = sx(ib2 + 1) - 1
                        act_a = act_a + ib2 - ib1 + 1
                    end if
                    if (act_a > 0) then
                        sy(rk(i)) = sy(rk(i)) + act_a
                        sz(rl(i)) = sz(rl(i)) + act_a
                    end if
                end do
                do j = leaf(1, best) + 1, leaf(2, best) + 1
                    sx(j) = sx(j) + sx(j - 1)
                end do

                ! Per axis, the deepest-dip and the strongest-inflection cut; the winner is
                ! the candidate whose tightened halves drop the most volume.
                bgain = 0
                do d = 1, 3
                    if (leaf(2*d, best) < leaf(2*d - 1, best) + 1) cycle
                    cdip = 0
                    cinf = 0
                    gmin = huge(1)
                    gmax = -huge(1)
                    do c = leaf(2*d - 1, best) + 1, leaf(2*d, best)
                        gc = f_ab_br_sig(d, sx, sy, sz, c - 1) + f_ab_br_sig(d, sx, sy, sz, c)
                        if (gc < gmin) then
                            gmin = gc
                            cdip = c
                        end if
                        if (c <= leaf(2*d, best) - 1) then
                            gc = f_ab_br_sig(d, sx, sy, sz, c - 1) + f_ab_br_sig(d, sx, sy, sz, c + 1) - 2*f_ab_br_sig(d, sx, sy, &
                                             & sz, c)
                            if (gc > gmax) then
                                gmax = gc
                                cinf = c
                            end if
                        end if
                    end do
                    do c = 1, 2
                        if (c == 1) then
                            cut = cdip
                        else
                            cut = cinf
                            if (cut == cdip) cycle
                        end if
                        if (cut == 0) cycle
                        t1 = leaf(:,best)
                        t2 = leaf(:,best)
                        t1(2*d) = cut - 1
                        t2(2*d - 1) = cut
                        call s_ab_br_hull(rk, rl, ra1, ra2, rb1, rb2, nrows, t1, ia1)
                        call s_ab_br_hull(rk, rl, ra1, ra2, rb1, rb2, nrows, t2, ia2)
                        gc = f_ab_br_vol(leaf(:,best)) - f_ab_br_vol(t1) - f_ab_br_vol(t2)
                        if (gc > bgain) then
                            bgain = gc
                            bt1 = t1
                            bt2 = t2
                            bia1 = ia1
                            bia2 = ia2
                        end if
                    end do
                end do

                if (bgain <= 0) then
                    ldone(best) = .true.
                    cycle
                end if

                ! An empty half is a pure tighten of the leaf, not a split.
                if (bia1 == 0) then
                    leaf(:,best) = bt2
                    lact(best) = bia2
                else if (bia2 == 0) then
                    leaf(:,best) = bt1
                    lact(best) = bia1
                else
                    leaf(:,best) = bt1
                    lact(best) = bia1
                    nleaf = nleaf + 1
                    leaf(:,nleaf) = bt2
                    lact(nleaf) = bia2
                    ldone(nleaf) = .false.
                end if
            end do

            ! Grow and snap like the window; the map is monotone, so every box lies inside the window.
            do i = 1, nleaf
                call s_ab_grow_snap(leaf(1, i), leaf(2, i), m, ab_box(1, i), ab_box(2, i))
                call s_ab_grow_snap(leaf(3, i), leaf(4, i), n, ab_box(3, i), ab_box(4, i))
                call s_ab_grow_snap(leaf(5, i), leaf(6, i), p, ab_box(5, i), ab_box(6, i))
            end do

            ! Abort unless the boxes contain every active row interval.
            do i = 0, nrows - 1
                do c = 1, 2
                    if (c == 1) then
                        lo = ra1(i); hi = ra2(i)
                    else
                        lo = rb1(i); hi = rb2(i)
                    end if
                    if (hi < lo) cycle
                    ! Segments of [lo, hi] covered by each box holding this row; walk the chain of furthest-reaching segments from
                    ! lo and require it to reach hi.
                    ns = 0
                    do j = 1, nleaf
                        if (rk(i) < ab_box(3, j) .or. rk(i) > ab_box(4, j) .or. rl(i) < ab_box(5, j) .or. rl(i) > ab_box(6, &
                            & j)) cycle
                        ia1 = max(ab_box(1, j), lo); ia2 = min(ab_box(2, j), hi)
                        if (ia2 >= ia1) then
                            ns = ns + 1
                            seg(1, ns) = ia1
                            seg(2, ns) = ia2
                        end if
                    end do
                    cov = lo
                    do while (cov <= hi)
                        cut = cov - 1
                        do j = 1, ns
                            if (seg(1, j) <= cov .and. seg(2, j) > cut) cut = seg(2, j)
                        end do
                        if (cut < cov) exit
                        cov = cut + 1
                    end do
                    if (cov <= hi) call s_mpi_abort('MFC: active-box cover violated activity containment.')
                end do
            end do

            ! A cover no tighter than the window is replaced by the single window box.
            cv = 0
            do i = 1, nleaf
                cv = cv + f_ab_br_vol(ab_box(:,i))
            end do
            if (real(cv, wp) >= 0.9_wp*real(max(0, ab_x%end - ab_x%beg + 1)*max(0, ab_y%end - ab_y%beg + 1)*max(0, &
                & ab_z%end - ab_z%beg + 1), wp)) then
                ab_nb = 1
                ab_box(:,1) = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end/)
                return
            end if
            ! Inactive slots are inverted, so they intersect no column.
            do i = nleaf + 1, AB_MAX_BOXES
                ab_box(:,i) = 0
                ab_box(2, i) = -1; ab_box(4, i) = -1; ab_box(6, i) = -1
            end do
            ab_nb = nleaf

        end subroutine s_derive_active_box_cover_br

        !> Tighten <bx> (x/y/z beg/end) to the hull of the rows' activity inside it; <act> is the active cell count, 0 with an
        !! inverted box when nothing inside is active.
        pure subroutine s_ab_br_hull(rk, rl, ra1, ra2, rb1, rb2, nrows, bx, act)

            integer, intent(in)    :: rk(0:), rl(0:), ra1(0:), ra2(0:), rb1(0:), rb2(0:), nrows
            integer, intent(inout) :: bx(1:6)
            integer, intent(out)   :: act
            integer                :: i, ia1, ia2, ib1, ib2, hx1, hx2, hy1, hy2, hz1, hz2

            hx1 = m + 1; hx2 = -1
            hy1 = n + 1; hy2 = -1
            hz1 = p + 1; hz2 = -1
            act = 0
            do i = 0, nrows - 1
                if (rk(i) < bx(3) .or. rk(i) > bx(4) .or. rl(i) < bx(5) .or. rl(i) > bx(6)) cycle
                ia1 = max(ra1(i), bx(1)); ia2 = min(ra2(i), bx(2))
                ib1 = max(rb1(i), bx(1)); ib2 = min(rb2(i), bx(2))
                if (ia2 < ia1 .and. ib2 < ib1) cycle
                if (ia2 >= ia1) then
                    hx1 = min(hx1, ia1); hx2 = max(hx2, ia2)
                    act = act + ia2 - ia1 + 1
                end if
                if (ib2 >= ib1) then
                    hx1 = min(hx1, ib1); hx2 = max(hx2, ib2)
                    act = act + ib2 - ib1 + 1
                end if
                hy1 = min(hy1, rk(i)); hy2 = max(hy2, rk(i))
                hz1 = min(hz1, rl(i)); hz2 = max(hz2, rl(i))
            end do
            if (hx2 >= hx1) then
                bx = (/hx1, hx2, hy1, hy2, hz1, hz2/)
            else
                bx = 0
                bx(2) = -1; bx(4) = -1; bx(6) = -1
                act = 0
            end if

        end subroutine s_ab_br_hull

        !> Cell count of a box, 0 when inverted.
        pure function f_ab_br_vol(bx) result(v)

            integer, intent(in) :: bx(1:6)
            integer             :: v

            v = max(0, bx(2) - bx(1) + 1)*max(0, bx(4) - bx(3) + 1)*max(0, bx(6) - bx(5) + 1)

        end function f_ab_br_vol

        !> Signature value along axis <d> (1: x column, 2: y slice, 3: z slice) at index <c>.
        pure function f_ab_br_sig(d, sx, sy, sz, c) result(v)

            integer, intent(in) :: d, c, sx(0:), sy(0:), sz(0:)
            integer             :: v

            if (d == 1) then
                v = sx(c)
            else if (d == 2) then
                v = sy(c)
            else
                v = sz(c)
            end if

        end function f_ab_br_sig
    #:endif

    !> Commit a re-derived upper window bound against the current one (lower bounds pass negated values). Growth applies at once; a
    !! tighten applies to a target repeated on the previous cadence.
    pure function f_ab_bound(new, old, prev, prev_valid) result(b)

        integer, intent(in) :: new, old, prev
        logical, intent(in) :: prev_valid
        integer             :: b

        if (new >= old .or. (prev_valid .and. new == prev)) then
            b = new
        else
            b = old
        end if

    end function f_ab_bound

    !> Window sharing payload, per axis the local window bounds, whether the window reaches within buff_size of the low and high
    !! subdomain edges, and the raw measured box (inverted when none was measured).
    pure function f_active_box_share() result(mine)

        integer :: mine(1:18)

        mine = (/ab_x%beg, ab_x%end, ab_y%beg, ab_y%end, ab_z%beg, ab_z%end, merge(1, 0, ab_x%beg - buff_size <= -1), merge(1, 0, &
                & ab_x%end + buff_size >= m + 1), merge(1, 0, ab_y%beg - buff_size <= -1), merge(1, 0, &
                & ab_y%end + buff_size >= n + 1), merge(1, 0, ab_z%beg - buff_size <= -1), merge(1, 0, &
                & ab_z%end + buff_size >= p + 1), merge(ab_raw(1), 0, ab_raw_valid), merge(ab_raw(2), -1, ab_raw_valid), &
                & merge(ab_raw(3), 0, ab_raw_valid), merge(ab_raw(4), -1, ab_raw_valid), merge(ab_raw(5), 0, ab_raw_valid), &
                & merge(ab_raw(6), -1, ab_raw_valid)/)

    end function f_active_box_share

    !> Widen the local window across every seam a neighbor's window has reached, so a wave arriving from that neighbor enters cells
    !! the local compute already covers. A periodic axis held by this rank alone is a seam with itself. The band carries the same
    !! growth margin. On the transverse axes it spans the face neighbor's raw box (its window when the raw box is empty), or the
    !! whole width when only a diagonal neighbor reached the seam. Every face whose boundary condition can drive the flow (ab_drive)
    !! is a full-width band too. <share> is identical on every rank, so both ends of the trimmed exchange derive the same windows.
    impure subroutine s_propagate_active_box_seams(share)

        integer, intent(in)   :: share(1:,1:)
        integer               :: j, k, l, r, rf, d, sd, t, edge
        integer               :: lo(1:3), hi(1:3), ext(1:3)
        logical               :: seam(1:3,1:2), full
        type(int_bounds_info) :: wx, wy, wz

        ! A neighbor's edge flag on the shared seam signals that its compute covers the cells across it.

        seam = .false.
#ifdef MFC_MPI
        if (num_procs > 1) then
            do l = nidx(3)%beg, nidx(3)%end
                do k = nidx(2)%beg, nidx(2)%end
                    do j = nidx(1)%beg, nidx(1)%end
                        if (abs(j) + abs(k) + abs(l) == 0) cycle
                        r = neighbor_ranks(j, k, l) + 1
                        if (j == -1 .and. share(8, r) == 1) seam(1, 1) = .true.  ! their x-high edge is our x-low seam
                        if (j == 1 .and. share(7, r) == 1) seam(1, 2) = .true.
                        if (k == -1 .and. share(10, r) == 1) seam(2, 1) = .true.  ! their y-high edge is our y-low seam
                        if (k == 1 .and. share(9, r) == 1) seam(2, 2) = .true.
                        if (l == -1 .and. share(12, r) == 1) seam(3, 1) = .true.  ! their z-high edge is our z-low seam
                        if (l == 1 .and. share(11, r) == 1) seam(3, 2) = .true.
                    end do
                end do
            end do
        end if
#endif
        do d = 1, num_dims
            if (f_ab_self_periodic(d, 1) .and. share(6 + 2*d, proc_rank + 1) == 1) seam(d, 1) = .true.
            if (f_ab_self_periodic(d, 2) .and. share(5 + 2*d, proc_rank + 1) == 1) seam(d, 2) = .true.
        end do
        if (.not. (any(seam) .or. any(ab_drive))) return

        ! The bands as one box, the seam's edge cell on the normal axis and the face neighbor's raw box on the transverse axes, in
        ! local coordinates both ranks share. The union below leaves an already covered seam unchanged.
        ext = (/m, n, p/)
        lo = ext + 1; hi = -1
        do d = 1, num_dims
            do sd = 1, 2
                if (.not. (seam(d, sd) .or. ab_drive(d, sd))) cycle
                edge = merge(0, ext(d), sd == 1)
                lo(d) = min(lo(d), edge); hi(d) = max(hi(d), edge)
                full = ab_drive(d, sd)
                if (.not. full) then
                    rf = f_ab_face_rank(d, sd)
                    ! The face neighbor's flag on the side facing this rank
                    full = share(6 + 2*d - (sd - 1), rf) /= 1
                end if
                do t = 1, 3
                    if (t == d) cycle
                    if (full) then
                        lo(t) = 0; hi(t) = ext(t)
                    else if (share(14, rf) >= share(13, rf)) then
                        lo(t) = min(lo(t), share(12 + 2*t - 1, rf)); hi(t) = max(hi(t), share(12 + 2*t, rf))
                    else
                        lo(t) = min(lo(t), share(2*t - 1, rf)); hi(t) = max(hi(t), share(2*t, rf))
                    end if
                end do
            end do
        end do

        ! Grow the box by the margin, snap it outward to ab_snap, record it as the band and union it into the window.
        call s_ab_grow_snap(lo(1), hi(1), m, ab_band(1), ab_band(2))
        call s_ab_grow_snap(lo(2), hi(2), n, ab_band(3), ab_band(4))
        call s_ab_grow_snap(lo(3), hi(3), p, ab_band(5), ab_band(6))
        ab_seam_band = .true.
        wx%beg = min(ab_x%beg, ab_band(1)); wx%end = max(ab_x%end, ab_band(2))
        wy%beg = min(ab_y%beg, ab_band(3)); wy%end = max(ab_y%end, ab_band(4))
        wz%beg = min(ab_z%beg, ab_band(5)); wz%end = max(ab_z%end, ab_band(6))
        ab_x = wx; ab_y = wy; ab_z = wz

    end subroutine s_propagate_active_box_seams

    !> Whether side <sd> (1 low, 2 high) of axis <d> is a periodic boundary this rank wraps onto itself.
    pure logical function f_ab_self_periodic(d, sd)

        integer, intent(in) :: d, sd

        if (d == 1) then
            f_ab_self_periodic = merge(bc_x%beg, bc_x%end, sd == 1) == BC_PERIODIC
        else if (d == 2) then
            f_ab_self_periodic = merge(bc_y%beg, bc_y%end, sd == 1) == BC_PERIODIC
        else
            f_ab_self_periodic = merge(bc_z%beg, bc_z%end, sd == 1) == BC_PERIODIC
        end if

    end function f_ab_self_periodic

    !> Share column of the rank across side <sd> of axis <d>: this rank on a self-periodic axis, else the face neighbor.
    impure integer function f_ab_face_rank(d, sd)

        integer, intent(in) :: d, sd
        integer             :: dir(1:3)

        f_ab_face_rank = proc_rank + 1
        if (f_ab_self_periodic(d, sd)) return
#ifdef MFC_MPI
        dir = 0
        dir(d) = 2*sd - 3
        f_ab_face_rank = neighbor_ranks(dir(1), dir(2), dir(3)) + 1
#endif

    end function f_ab_face_rank

end module m_active_box
#else
!> Abort when MFC_ACTIVE_BOX is set but the build has no active-box support. Unset or empty continues.
subroutine s_check_active_box_request()

    use m_mpi_common, only: s_mpi_abort

    implicit none

    character(len=64) :: spec
    integer           :: ierr

    call get_environment_variable('MFC_ACTIVE_BOX', value=spec, status=ierr)
    if (ierr == 0 .and. len_trim(spec) > 0) then
        call s_mpi_abort('MFC_ACTIVE_BOX is set but this build has no active-box support.')
    end if

end subroutine s_check_active_box_request
#endif
