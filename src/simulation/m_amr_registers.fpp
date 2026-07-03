!>
!!@file
!!@brief Contains module m_amr_registers

#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine patch boundary (SP4). Depends only on m_derived_types +
!! m_global_parameters so both m_rhs (capture) and m_time_steppers (apply) can use it without cycles. NOTE: m_amr uses m_rhs which
!! uses m_amr_registers, so adding "use m_amr" here would create a compilation cycle. Region info is therefore read from
!! amr_region_lo/hi and amr_isect_lo/hi (m_global_parameters), which s_set_amr_fine_geometry keeps mirroring across regrids. creg
!! uses 0-based transverse indexing relative to the rank's patch INTERSECTION (= the patch at np=1); freg uses 0-based LOCAL fine
!! indexing (aligned: fine children of isect cell t are 2*t and 2*t+1). All arrays are preallocated at max size so regrid requires
!! no reallocation.
module m_amr_registers

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_amr_zero_fine_registers, &
        & s_amr_apply_reflux_state, s_finalize_amr_registers, s_amr_reflux_face_flags, freg

    !> SSP-RK3 effective flux weights: q^{n+1} = q^n + dt*(L(q^n)/6 + L(q^(1))/6 + 2*L(q^(2))/3).
    real(wp), parameter :: rk3_w(3) = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

    !> Registers for the two patch faces normal to one direction: (1:sys_size, transverse-1, transverse-2).
    type t_face_reg
        real(wp), allocatable :: lo(:,:,:)
        real(wp), allocatable :: hi(:,:,:)
    end type t_face_reg

    type(t_face_reg) :: creg(3)  !< coarse flux at the patch boundary faces (relative 0-based transverse)
    type(t_face_reg) :: freg(3)  !< fine flux at the covering fine faces (0-based fine transverse)
    $:GPU_DECLARE(create='[creg, freg]')

contains

    !> Reflux-face participation for THIS rank: own_lo(d)/own_hi(d) = it owns the coarse cell layer just OUTSIDE the patch's
    !! low/high face in dim d (where the coarse capture and both reflux applies run; at an interior face the same rank also holds
    !! the inside cells) - i.e. the outside layer lies in its subdomain in dim d and its patch intersection is nonempty in the
    !! transverse dims. All true at np=1. Also returns sidx/ext (collapsed dims pinned to 0). Reads the COARSE grid state in m/n/p -
    !! never call from inside the fine advance (the swapped fine branch of the capture).
    impure subroutine s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi)

        integer, intent(out) :: sidx(3), ext(3)
        logical, intent(out) :: own_lo(3), own_hi(3)
        logical              :: tv(3), tvd
        integer              :: d, t

        sidx = 0; ext = 0
        sidx(1) = start_idx(1); ext(1) = m
        if (n_glb > 0) then; sidx(2) = start_idx(2); ext(2) = n; end if
        if (p_glb > 0) then; sidx(3) = start_idx(3); ext(3) = p; end if
        tv(1) = amr_isect_lo(1) <= amr_isect_hi(1)
        tv(2) = (n_glb == 0) .or. amr_isect_lo(2) <= amr_isect_hi(2)
        tv(3) = (p_glb == 0) .or. amr_isect_lo(3) <= amr_isect_hi(3)
        own_lo = .false.; own_hi = .false.
        do d = 1, num_dims
            tvd = .true.
            do t = 1, num_dims
                if (t /= d) tvd = tvd .and. tv(t)
            end do
            own_lo(d) = tvd .and. amr_region_lo(d) - 1 >= sidx(d) .and. amr_region_lo(d) - 1 <= sidx(d) + ext(d)
            own_hi(d) = tvd .and. amr_region_hi(d) + 1 >= sidx(d) .and. amr_region_hi(d) + 1 <= sidx(d) + ext(d)
        end do

    end subroutine s_amr_reflux_face_flags

    impure subroutine s_initialize_amr_registers()

        integer :: maxc1, maxc2, maxc3, max_f1, max_f2, max_f3

        if (.not. amr) return
        ! Registers on ALL ranks: regrid moves the patch faces, so any rank can become a participant (fine
        ! cells for freg; outside face layer for creg capture + apply and for RECEIVING freg slices when a
        ! patch face sits on a rank boundary). Participation is re-derived per call from the current box.
        ! max coarse patch cells per dim THIS rank can cover (must match m_amr's preallocation cap). The
        ! transverse extents match the face-neighbor's by construction (cart neighbors share transverse
        ! subdomains), so the whole-array freg sendrecvs in m_mpi_proxy pair up exactly.
        maxc1 = min((m_glb + 1)/2, (m + 1)/2)
        maxc2 = 1; maxc3 = 1
        if (n_glb > 0) maxc2 = min((n_glb + 1)/2, (n + 1)/2)
        if (p_glb > 0) maxc3 = min((p_glb + 1)/2, (p + 1)/2)
        max_f1 = 2*maxc1 - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = 2*maxc2 - 1
        if (p_glb > 0) max_f3 = 2*maxc3 - 1
        ! creg: relative 0-based transverse (0:maxc_t-1); freg: 0-based fine (0:max_f_t).
        ! Device-resident (@:ALLOCATE): capture and both applies run as kernels; no host copies are read.
        @:ALLOCATE(creg(1)%lo(1:sys_size,0:maxc2 - 1,0:maxc3 - 1), creg(1)%hi(1:sys_size,0:maxc2 - 1,0:maxc3 - 1))
        @:ALLOCATE(freg(1)%lo(1:sys_size,0:max_f2,0:max_f3), freg(1)%hi(1:sys_size,0:max_f2,0:max_f3))
        if (n_glb > 0) then
            @:ALLOCATE(creg(2)%lo(1:sys_size,0:maxc1 - 1,0:maxc3 - 1), creg(2)%hi(1:sys_size,0:maxc1 - 1,0:maxc3 - 1))
            @:ALLOCATE(freg(2)%lo(1:sys_size,0:max_f1,0:max_f3), freg(2)%hi(1:sys_size,0:max_f1,0:max_f3))
        end if
        if (p_glb > 0) then
            @:ALLOCATE(creg(3)%lo(1:sys_size,0:maxc1 - 1,0:maxc2 - 1), creg(3)%hi(1:sys_size,0:maxc1 - 1,0:maxc2 - 1))
            @:ALLOCATE(freg(3)%lo(1:sys_size,0:max_f1,0:max_f2), freg(3)%hi(1:sys_size,0:max_f1,0:max_f2))
        end if

    end subroutine s_initialize_amr_registers

    !> Capture the c/f boundary-face fluxes for direction id from the just-finalized flux array. Runs INSIDE s_compute_rhs: coarse
    !! call (amr_in_fine_advance false, coarse globals) fills creg at the patch boundary faces; fine call (flag true, globals
    !! swapped to the fine patch) fills freg at fine faces -1 and m/n/p. creg uses relative 0-based transverse; freg uses 0-based
    !! fine.
    impure subroutine s_amr_capture_boundary_flux(id, flux_dir, stage)

        integer, intent(in)            :: id
        type(vector_field), intent(in) :: flux_dir
        integer, intent(in)            :: stage
        integer                        :: eq, t1, t2, jlo, jhi, t1_hi, t2_hi, o1, o2
        integer                        :: sidx(3), ext(3)
        logical                        :: own_lo(3), own_hi(3), cap_lo, cap_hi
        real(wp)                       :: coef
        logical                        :: accum

        if (.not. amr) return
        if (amr_in_fine_advance .and. .not. amr_rank_owns_patch) return
        ! flux data was just written by device kernels; the face reads below run as device kernels too
        if (amr_subcycle) then
            if (amr_in_fine_advance) then
                coef = 0.5_wp*rk3_w(stage); accum = .true.  ! zeroed by s_amr_zero_fine_registers before substep 1
            else
                coef = rk3_w(stage); accum = (stage > 1)  ! stage 1 overwrites = implicit zero per coarse step
            end if
        else
            coef = 1._wp; accum = .false.  ! overwrite each stage - default behavior, byte-identical
        end if
        if (amr_in_fine_advance) then
            ! fine branch: globals swapped; jlo=-1, jhi=current fine extent in direction id
            select case (id)
            case (1); jlo = -1; jhi = m; t1_hi = n; t2_hi = p
            case (2); jlo = -1; jhi = n; t1_hi = m; t2_hi = p
            case (3); jlo = -1; jhi = p; t1_hi = m; t2_hi = n
            end select
            $:GPU_PARALLEL_LOOP(collapse=3)
            do t2 = 0, t2_hi
                do t1 = 0, t1_hi
                    do eq = 1, sys_size
                        select case (id)
                        case (1)
                            if (accum) then
                                freg(1)%lo(eq, t1, t2) = freg(1)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jlo, t1, t2), wp)
                                freg(1)%hi(eq, t1, t2) = freg(1)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jhi, t1, t2), wp)
                            else
                                freg(1)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jlo, t1, t2), wp)
                                freg(1)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jhi, t1, t2), wp)
                            end if
                        case (2)
                            if (accum) then
                                freg(2)%lo(eq, t1, t2) = freg(2)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(t1, jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2) = freg(2)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(t1, jhi, t2), wp)
                            else
                                freg(2)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(t1, jlo, t2), wp)
                                freg(2)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(t1, jhi, t2), wp)
                            end if
                        case (3)
                            if (accum) then
                                freg(3)%lo(eq, t1, t2) = freg(3)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(t1, t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2) = freg(3)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(t1, t2, jhi), wp)
                            else
                                freg(3)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(t1, t2, jlo), wp)
                                freg(3)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(t1, t2, jhi), wp)
                            end if
                        end select
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else
            ! coarse branch: a face's capture runs on the rank owning the coarse cells just OUTSIDE it (its
            ! flux_n covers that face; at a rank-interior face the same rank also holds the inside cells).
            ! jlo/jhi = LOCAL flux indices of the patch's low/high faces; t1/t2 = 0-based transverse indices
            ! relative to this rank's patch INTERSECTION (o1/o2 = local transverse origins), aligned with the
            ! fine registers: the fine children of isect-relative cell t are faces 2*t and 2*t+1. At np=1 the
            ! intersection is the patch and both flags hold, recovering the single-rank behavior exactly.
            call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi)
            cap_lo = own_lo(id); cap_hi = own_hi(id)
            if (.not. (cap_lo .or. cap_hi)) return
            select case (id)
            case (1); jlo = amr_region_lo(1) - 1 - sidx(1); jhi = amr_region_hi(1) - sidx(1)
                t1_hi = amr_isect_hi(2) - amr_isect_lo(2); o1 = amr_isect_lo(2) - sidx(2)
                t2_hi = amr_isect_hi(3) - amr_isect_lo(3); o2 = amr_isect_lo(3) - sidx(3)
            case (2); jlo = amr_region_lo(2) - 1 - sidx(2); jhi = amr_region_hi(2) - sidx(2)
                t1_hi = amr_isect_hi(1) - amr_isect_lo(1); o1 = amr_isect_lo(1) - sidx(1)
                t2_hi = amr_isect_hi(3) - amr_isect_lo(3); o2 = amr_isect_lo(3) - sidx(3)
            case (3); jlo = amr_region_lo(3) - 1 - sidx(3); jhi = amr_region_hi(3) - sidx(3)
                t1_hi = amr_isect_hi(1) - amr_isect_lo(1); o1 = amr_isect_lo(1) - sidx(1)
                t2_hi = amr_isect_hi(2) - amr_isect_lo(2); o2 = amr_isect_lo(2) - sidx(2)
            end select
            $:GPU_PARALLEL_LOOP(collapse=3)
            do t2 = 0, t2_hi
                do t1 = 0, t1_hi
                    do eq = 1, sys_size
                        select case (id)
                        case (1)
                            if (cap_lo) then
                                if (accum) then
                                    creg(1)%lo(eq, t1, t2) = creg(1)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jlo, o1 + t1, &
                                         & o2 + t2), wp)
                                else
                                    creg(1)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jlo, o1 + t1, o2 + t2), wp)
                                end if
                            end if
                            if (cap_hi) then
                                if (accum) then
                                    creg(1)%hi(eq, t1, t2) = creg(1)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jhi, o1 + t1, &
                                         & o2 + t2), wp)
                                else
                                    creg(1)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jhi, o1 + t1, o2 + t2), wp)
                                end if
                            end if
                        case (2)
                            if (cap_lo) then
                                if (accum) then
                                    creg(2)%lo(eq, t1, t2) = creg(2)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, jlo, &
                                         & o2 + t2), wp)
                                else
                                    creg(2)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, jlo, o2 + t2), wp)
                                end if
                            end if
                            if (cap_hi) then
                                if (accum) then
                                    creg(2)%hi(eq, t1, t2) = creg(2)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, jhi, &
                                         & o2 + t2), wp)
                                else
                                    creg(2)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, jhi, o2 + t2), wp)
                                end if
                            end if
                        case (3)
                            if (cap_lo) then
                                if (accum) then
                                    creg(3)%lo(eq, t1, t2) = creg(3)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, &
                                         & o2 + t2, jlo), wp)
                                else
                                    creg(3)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jlo), wp)
                                end if
                            end if
                            if (cap_hi) then
                                if (accum) then
                                    creg(3)%hi(eq, t1, t2) = creg(3)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(o1 + t1, &
                                         & o2 + t2, jhi), wp)
                                else
                                    creg(3)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(o1 + t1, o2 + t2, jhi), wp)
                                end if
                            end if
                        end select
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Correct the coarse rhs in the first cell OUTSIDE each patch face so the coarse update sees the (child-averaged) fine flux at
    !! every c/f face. Signs follow rhs = (flux_left - flux_right)/dx: low face is the outside cell's RIGHT face => rhs += (F_coarse
    !! - Fbar_fine)/dx; high face is the outside cell's LEFT face => rhs += (Fbar_fine - F_coarse)/dx. Cells INSIDE the patch need
    !! no correction (end-of-step restriction overwrites them). c1/c2 are relative 0-based coarse transverse indices.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch
        integer                                                :: nx_c, ny_c, nz_c, ol1, ol2, ol3, oh1, oh2, oh3, tl1, tl2, tl3
        integer                                                :: dd1_hi, dd2_hi, sidx(3), ext(3)
        logical                                                :: d2, d3, own_lo(3), own_hi(3), has_lo, has_hi
        real(wp)                                               :: fblo, fbhi, mlo, mhi

        if (.not. amr) return
        ! per-face participation: each face's correction runs on the rank owning its OUTSIDE cell layer
        ! (all faces at np=1); freg slices from a rank-boundary face's fine side arrive via
        ! s_mpi_sendrecv_amr_reflux_faces before this is called
        call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi)
        if (.not. (any(own_lo) .or. any(own_hi))) return
        ! device kernels: the coarse rhs stays device-resident for the coarse RK update kernel
        d2 = n_glb > 0; d3 = p_glb > 0
        ! this rank's covered patch cells (isect-relative, 0-based): 0..n{x,y,z}_c; matches creg/freg indexing
        nx_c = amr_isect_hi(1) - amr_isect_lo(1)
        ny_c = 0; nz_c = 0
        if (n_glb > 0) ny_c = amr_isect_hi(2) - amr_isect_lo(2)
        if (p_glb > 0) nz_c = amr_isect_hi(3) - amr_isect_lo(3)
        ! LOCAL cell indices (rhs_vf/dx are rank-local): outside cells ol/oh, isect transverse origins tl
        ol1 = amr_region_lo(1) - 1 - sidx(1); oh1 = amr_region_hi(1) + 1 - sidx(1)
        ol2 = amr_region_lo(2) - 1 - sidx(2); oh2 = amr_region_hi(2) + 1 - sidx(2)
        ol3 = amr_region_lo(3) - 1 - sidx(3); oh3 = amr_region_hi(3) + 1 - sidx(3)
        tl1 = amr_isect_lo(1) - sidx(1); tl2 = amr_isect_lo(2) - sidx(2); tl3 = amr_isect_lo(3) - sidx(3)
        ! x-faces: transverse dims (y, z); children in each active transverse dim
        has_lo = own_lo(1); has_hi = own_hi(1)
        if (has_lo .or. has_hi) then
            nch = 1
            if (n_glb > 0) nch = nch*2
            if (p_glb > 0) nch = nch*2
            dd1_hi = merge(1, 0, n_glb > 0); dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dx(ol1)
            if (has_hi) mhi = dx(oh1)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    do c1 = 0, ny_c
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 0; if (d2) f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, dd1_hi
                                fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(ol1, tl2 + c1, tl3 + c2) = rhs_vf(eq)%sf(ol1, tl2 + c1, &
                            & tl3 + c2) + (creg(1)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(oh1, tl2 + c1, tl3 + c2) = rhs_vf(eq)%sf(oh1, tl2 + c1, &
                            & tl3 + c2) + (fbhi - creg(1)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! y-faces (n_glb > 0): transverse dims (x, z); x is always active (2 children)
        has_lo = own_lo(2); has_hi = own_hi(2)
        if (n_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dy(ol2)
            if (has_hi) mhi = dy(oh2)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    do c1 = 0, nx_c
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(tl1 + c1, ol2, tl3 + c2) = rhs_vf(eq)%sf(tl1 + c1, ol2, &
                            & tl3 + c2) + (creg(2)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(tl1 + c1, oh2, tl3 + c2) = rhs_vf(eq)%sf(tl1 + c1, oh2, &
                            & tl3 + c2) + (fbhi - creg(2)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        ! z-faces (p_glb > 0): transverse dims (x, y); both always active in 3D (4 children)
        has_lo = own_lo(3); has_hi = own_hi(3)
        if (p_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 4
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dz(ol3)
            if (has_hi) mhi = dz(oh3)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, ny_c
                    do c1 = 0, nx_c
                        f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, ol3) = rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, &
                            & ol3) + (creg(3)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, oh3) = rhs_vf(eq)%sf(tl1 + c1, tl2 + c2, &
                            & oh3) + (fbhi - creg(3)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_apply_reflux

    !> Zero the fine registers (called by the subcycle driver before substep 1 - stage-1 overwrite cannot work across two substeps).
    impure subroutine s_amr_zero_fine_registers()

        integer :: d, eq, t1, t2, t1_hi, t2_hi

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        do d = 1, 3
            if (allocated(freg(d)%lo)) then
                t1_hi = ubound(freg(d)%lo, 2); t2_hi = ubound(freg(d)%lo, 3)
                $:GPU_PARALLEL_LOOP(collapse=3)
                do t2 = 0, t2_hi
                    do t1 = 0, t1_hi
                        do eq = 1, sys_size
                            freg(d)%lo(eq, t1, t2) = 0._wp
                            freg(d)%hi(eq, t1, t2) = 0._wp
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end do

    end subroutine s_amr_zero_fine_registers

    !> Berger-Colella state correction (subcycle mode only): after restriction, correct the first coarse cell OUTSIDE each patch
    !! face with the time-accumulated flux mismatch: low face: q += dt*(F_c_eff - Fbar_f_eff)/dx ; high face: q += dt*(Fbar_f_eff -
    !! F_c_eff)/dx. Registers hold EFFECTIVE (rk3_w-weighted, substep-averaged) fluxes in subcycle mode.
    impure subroutine s_amr_apply_reflux_state(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch
        integer                                                :: nx_c, ny_c, nz_c, ol1, ol2, ol3, oh1, oh2, oh3, tl1, tl2, tl3
        integer                                                :: dd1_hi, dd2_hi, sidx(3), ext(3)
        logical                                                :: d2, d3, own_lo(3), own_hi(3), has_lo, has_hi
        real(wp)                                               :: fblo, fbhi, mlo, mhi, dtl

        if (.not. amr) return
        ! per-face participation and index conventions: see s_amr_apply_reflux
        call s_amr_reflux_face_flags(sidx, ext, own_lo, own_hi)
        if (.not. (any(own_lo) .or. any(own_hi))) return
        ! device kernels: the restricted coarse state stays device-resident
        d2 = n_glb > 0; d3 = p_glb > 0
        dtl = dt
        nx_c = amr_isect_hi(1) - amr_isect_lo(1)
        ny_c = 0; nz_c = 0
        if (n_glb > 0) ny_c = amr_isect_hi(2) - amr_isect_lo(2)
        if (p_glb > 0) nz_c = amr_isect_hi(3) - amr_isect_lo(3)
        ol1 = amr_region_lo(1) - 1 - sidx(1); oh1 = amr_region_hi(1) + 1 - sidx(1)
        ol2 = amr_region_lo(2) - 1 - sidx(2); oh2 = amr_region_hi(2) + 1 - sidx(2)
        ol3 = amr_region_lo(3) - 1 - sidx(3); oh3 = amr_region_hi(3) + 1 - sidx(3)
        tl1 = amr_isect_lo(1) - sidx(1); tl2 = amr_isect_lo(2) - sidx(2); tl3 = amr_isect_lo(3) - sidx(3)
        has_lo = own_lo(1); has_hi = own_hi(1)
        if (has_lo .or. has_hi) then
            nch = 1
            if (n_glb > 0) nch = nch*2
            if (p_glb > 0) nch = nch*2
            dd1_hi = merge(1, 0, n_glb > 0); dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dx(ol1)
            if (has_hi) mhi = dx(oh1)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    do c1 = 0, ny_c
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 0; if (d2) f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, dd1_hi
                                fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(ol1, tl2 + c1, tl3 + c2) = q_cons(eq)%sf(ol1, tl2 + c1, &
                            & tl3 + c2) + dtl*(creg(1)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(oh1, tl2 + c1, tl3 + c2) = q_cons(eq)%sf(oh1, tl2 + c1, &
                            & tl3 + c2) + dtl*(fbhi - creg(1)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        has_lo = own_lo(2); has_hi = own_hi(2)
        if (n_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            dd2_hi = merge(1, 0, p_glb > 0)
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dy(ol2)
            if (has_hi) mhi = dy(oh2)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    do c1 = 0, nx_c
                        f20 = 0; if (d3) f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, dd2_hi
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(tl1 + c1, ol2, tl3 + c2) = q_cons(eq)%sf(tl1 + c1, ol2, &
                            & tl3 + c2) + dtl*(creg(2)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(tl1 + c1, oh2, tl3 + c2) = q_cons(eq)%sf(tl1 + c1, oh2, &
                            & tl3 + c2) + dtl*(fbhi - creg(2)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        has_lo = own_lo(3); has_hi = own_hi(3)
        if (p_glb > 0 .and. (has_lo .or. has_hi)) then
            nch = 4
            mlo = 1._wp; mhi = 1._wp
            if (has_lo) mlo = dz(ol3)
            if (has_hi) mhi = dz(oh3)
            $:GPU_PARALLEL_LOOP(collapse=3, private='[f10, f20, dd1, dd2, fblo, fbhi]')
            do eq = 1, sys_size
                do c2 = 0, ny_c
                    do c1 = 0, nx_c
                        f20 = 2*c2
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        if (has_lo) q_cons(eq)%sf(tl1 + c1, tl2 + c2, ol3) = q_cons(eq)%sf(tl1 + c1, tl2 + c2, &
                            & ol3) + dtl*(creg(3)%lo(eq, c1, c2) - fblo)/mlo
                        if (has_hi) q_cons(eq)%sf(tl1 + c1, tl2 + c2, oh3) = q_cons(eq)%sf(tl1 + c1, tl2 + c2, &
                            & oh3) + dtl*(fbhi - creg(3)%hi(eq, c1, c2))/mhi
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_amr_apply_reflux_state

    impure subroutine s_finalize_amr_registers()

        integer :: d

        if (.not. amr) return
        do d = 1, 3
            if (allocated(creg(d)%lo)) then
                @:DEALLOCATE(creg(d)%lo, creg(d)%hi)
            end if
            if (allocated(freg(d)%lo)) then
                @:DEALLOCATE(freg(d)%lo, freg(d)%hi)
            end if
        end do

    end subroutine s_finalize_amr_registers

end module m_amr_registers
