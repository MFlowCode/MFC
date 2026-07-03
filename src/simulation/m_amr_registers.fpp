!>
!!@file
!!@brief Contains module m_amr_registers

#:include 'macros.fpp'

!> @brief AMR flux registers: per-RK-stage refluxing at the coarse/fine patch boundary (SP4). Depends only on m_derived_types +
!! m_global_parameters so both m_rhs (capture) and m_time_steppers (apply) can use it without cycles. NOTE: m_amr uses m_rhs which
!! uses m_amr_registers, so adding "use m_amr" here would create a compilation cycle. Region info is therefore read from
!! amr_region_lo/hi (m_global_parameters), which s_set_amr_fine_geometry keeps mirroring amr_fine%region across regrids. creg uses
!! relative 0-based transverse indexing for regrid readiness; freg uses 0-based fine indexing. All arrays are preallocated at max
!! size so regrid requires no reallocation.
module m_amr_registers

    use m_derived_types
    use m_global_parameters

    implicit none

    private; public :: s_initialize_amr_registers, s_amr_capture_boundary_flux, s_amr_apply_reflux, s_amr_zero_fine_registers, &
        & s_amr_apply_reflux_state, s_finalize_amr_registers

    !> SSP-RK3 effective flux weights: q^{n+1} = q^n + dt*(L(q^n)/6 + L(q^(1))/6 + 2*L(q^(2))/3).
    real(wp), parameter :: rk3_w(3) = [1._wp/6._wp, 1._wp/6._wp, 2._wp/3._wp]

    !> Registers for the two patch faces normal to one direction: (1:sys_size, transverse-1, transverse-2).
    type t_face_reg
        real(wp), allocatable :: lo(:,:,:)
        real(wp), allocatable :: hi(:,:,:)
    end type t_face_reg

    type(t_face_reg) :: creg(3)  !< coarse flux at the patch boundary faces (relative 0-based transverse)
    type(t_face_reg) :: freg(3)  !< fine flux at the covering fine faces (0-based fine transverse)

contains

    impure subroutine s_initialize_amr_registers()

        integer :: maxc1, maxc2, maxc3, max_f1, max_f2, max_f3

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return  ! registers are owner-only (like the fine level itself)
        ! max coarse patch cells per dim (must match m_amr's amr_maxc)
        maxc1 = (m_glb + 1)/2
        maxc2 = 1; maxc3 = 1
        if (n_glb > 0) maxc2 = (n_glb + 1)/2
        if (p_glb > 0) maxc3 = (p_glb + 1)/2
        max_f1 = 2*maxc1 - 1
        max_f2 = 0; max_f3 = 0
        if (n_glb > 0) max_f2 = 2*maxc2 - 1
        if (p_glb > 0) max_f3 = 2*maxc3 - 1
        ! creg: relative 0-based transverse (0:maxc_t-1); freg: 0-based fine (0:max_f_t)
        allocate (creg(1)%lo(1:sys_size,0:maxc2 - 1,0:maxc3 - 1), creg(1)%hi(1:sys_size,0:maxc2 - 1,0:maxc3 - 1))
        allocate (freg(1)%lo(1:sys_size,0:max_f2,0:max_f3), freg(1)%hi(1:sys_size,0:max_f2,0:max_f3))
        if (n_glb > 0) then
            allocate (creg(2)%lo(1:sys_size,0:maxc1 - 1,0:maxc3 - 1), creg(2)%hi(1:sys_size,0:maxc1 - 1,0:maxc3 - 1))
            allocate (freg(2)%lo(1:sys_size,0:max_f1,0:max_f3), freg(2)%hi(1:sys_size,0:max_f1,0:max_f3))
        end if
        if (p_glb > 0) then
            allocate (creg(3)%lo(1:sys_size,0:maxc1 - 1,0:maxc2 - 1), creg(3)%hi(1:sys_size,0:maxc1 - 1,0:maxc2 - 1))
            allocate (freg(3)%lo(1:sys_size,0:max_f1,0:max_f2), freg(3)%hi(1:sys_size,0:max_f1,0:max_f2))
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
        integer                        :: eq, t1, t2, jlo, jhi, t1_hi, t2_hi, al1, al2, al3
        real(wp)                       :: coef
        logical                        :: accum

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        ! flux data was just written by device kernels; the face reads below run on host.
        ! M2: kernel-port removes this transfer (or shrinks it to the boundary faces)
        do eq = 1, sys_size
            $:GPU_UPDATE(host='[flux_dir%vf(eq)%sf]')
        end do
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
        else
            ! coarse branch: jlo/jhi = patch boundary faces; t1/t2 relative 0-based transverse.
            ! amr_region_lo/hi are GLOBAL; the coarse flux array is rank-LOCAL, so index via the
            ! local patch origin al1/al2/al3 (= global origin - start_idx; identical at np=1)
            al1 = amr_region_lo(1) - start_idx(1)
            al2 = amr_region_lo(2); if (n_glb > 0) al2 = al2 - start_idx(2)
            al3 = amr_region_lo(3); if (p_glb > 0) al3 = al3 - start_idx(3)
            select case (id)
            case (1); jlo = al1 - 1; jhi = al1 + amr_region_hi(1) - amr_region_lo(1)
                t1_hi = amr_region_hi(2) - amr_region_lo(2); t2_hi = amr_region_hi(3) - amr_region_lo(3)
            case (2); jlo = al2 - 1; jhi = al2 + amr_region_hi(2) - amr_region_lo(2)
                t1_hi = amr_region_hi(1) - amr_region_lo(1); t2_hi = amr_region_hi(3) - amr_region_lo(3)
            case (3); jlo = al3 - 1; jhi = al3 + amr_region_hi(3) - amr_region_lo(3)
                t1_hi = amr_region_hi(1) - amr_region_lo(1); t2_hi = amr_region_hi(2) - amr_region_lo(2)
            end select
            do t2 = 0, t2_hi
                do t1 = 0, t1_hi
                    do eq = 1, sys_size
                        select case (id)
                        case (1)
                            if (accum) then
                                creg(1)%lo(eq, t1, t2) = creg(1)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jlo, al2 + t1, &
                                     & al3 + t2), wp)
                                creg(1)%hi(eq, t1, t2) = creg(1)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(jhi, al2 + t1, &
                                     & al3 + t2), wp)
                            else
                                creg(1)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jlo, al2 + t1, al3 + t2), wp)
                                creg(1)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(jhi, al2 + t1, al3 + t2), wp)
                            end if
                        case (2)
                            if (accum) then
                                creg(2)%lo(eq, t1, t2) = creg(2)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(al1 + t1, jlo, &
                                     & al3 + t2), wp)
                                creg(2)%hi(eq, t1, t2) = creg(2)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(al1 + t1, jhi, &
                                     & al3 + t2), wp)
                            else
                                creg(2)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(al1 + t1, jlo, al3 + t2), wp)
                                creg(2)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(al1 + t1, jhi, al3 + t2), wp)
                            end if
                        case (3)
                            if (accum) then
                                creg(3)%lo(eq, t1, t2) = creg(3)%lo(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(al1 + t1, &
                                     & al2 + t2, jlo), wp)
                                creg(3)%hi(eq, t1, t2) = creg(3)%hi(eq, t1, t2) + coef*real(flux_dir%vf(eq)%sf(al1 + t1, &
                                     & al2 + t2, jhi), wp)
                            else
                                creg(3)%lo(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(al1 + t1, al2 + t2, jlo), wp)
                                creg(3)%hi(eq, t1, t2) = coef*real(flux_dir%vf(eq)%sf(al1 + t1, al2 + t2, jhi), wp)
                            end if
                        end select
                    end do
                end do
            end do
        end if

    end subroutine s_amr_capture_boundary_flux

    !> Correct the coarse rhs in the first cell OUTSIDE each patch face so the coarse update sees the (child-averaged) fine flux at
    !! every c/f face. Signs follow rhs = (flux_left - flux_right)/dx: low face is the outside cell's RIGHT face => rhs += (F_coarse
    !! - Fbar_fine)/dx; high face is the outside cell's LEFT face => rhs += (Fbar_fine - F_coarse)/dx. Cells INSIDE the patch need
    !! no correction (end-of-step restriction overwrites them). c1/c2 are relative 0-based coarse transverse indices.
    impure subroutine s_amr_apply_reflux(rhs_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch
        integer                                                :: nx_c, ny_c, nz_c, al1, al2, al3, ah1, ah2, ah3
        real(wp)                                               :: fblo, fbhi

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        ! coarse rhs to host for the boundary correction; pushed back below for the coarse RK
        ! update kernel. M2: kernel-port removes this transfer
        do eq = 1, sys_size
            $:GPU_UPDATE(host='[rhs_vf(eq)%sf]')
        end do
        ! current coarse patch extents (relative, 0-based): 0..n{x,y,z}_c
        nx_c = amr_region_hi(1) - amr_region_lo(1)
        ny_c = 0; nz_c = 0
        if (n_glb > 0) ny_c = amr_region_hi(2) - amr_region_lo(2)
        if (p_glb > 0) nz_c = amr_region_hi(3) - amr_region_lo(3)
        ! LOCAL patch bounds (rhs_vf/dx are rank-local; amr_region_lo/hi are global)
        al1 = amr_region_lo(1) - start_idx(1); ah1 = al1 + nx_c
        al2 = amr_region_lo(2); if (n_glb > 0) al2 = al2 - start_idx(2)
        al3 = amr_region_lo(3); if (p_glb > 0) al3 = al3 - start_idx(3)
        ah2 = al2 + ny_c; ah3 = al3 + nz_c
        ! x-faces: transverse dims (y, z); children in each active transverse dim
        nch = 1
        if (n_glb > 0) nch = nch*2
        if (p_glb > 0) nch = nch*2
        do eq = 1, sys_size
            do c2 = 0, nz_c
                f20 = 0; if (p_glb > 0) f20 = 2*c2
                do c1 = 0, ny_c
                    f10 = 0; if (n_glb > 0) f10 = 2*c1
                    fblo = 0._wp; fbhi = 0._wp
                    do dd2 = 0, merge(1, 0, p_glb > 0)
                        do dd1 = 0, merge(1, 0, n_glb > 0)
                            fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2)
                            fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2)
                        end do
                    end do
                    fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                    rhs_vf(eq)%sf(al1 - 1, al2 + c1, al3 + c2) = rhs_vf(eq)%sf(al1 - 1, al2 + c1, al3 + c2) + (creg(1)%lo(eq, c1, &
                           & c2) - fblo)/dx(al1 - 1)
                    rhs_vf(eq)%sf(ah1 + 1, al2 + c1, al3 + c2) = rhs_vf(eq)%sf(ah1 + 1, al2 + c1, &
                           & al3 + c2) + (fbhi - creg(1)%hi(eq, c1, c2))/dx(ah1 + 1)
                end do
            end do
        end do
        ! y-faces (n_glb > 0): transverse dims (x, z); x is always active (2 children)
        if (n_glb > 0) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    f20 = 0; if (p_glb > 0) f20 = 2*c2
                    do c1 = 0, nx_c
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, merge(1, 0, p_glb > 0)
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        rhs_vf(eq)%sf(al1 + c1, al2 - 1, al3 + c2) = rhs_vf(eq)%sf(al1 + c1, al2 - 1, al3 + c2) + (creg(2)%lo(eq, &
                               & c1, c2) - fblo)/dy(al2 - 1)
                        rhs_vf(eq)%sf(al1 + c1, ah2 + 1, al3 + c2) = rhs_vf(eq)%sf(al1 + c1, ah2 + 1, &
                               & al3 + c2) + (fbhi - creg(2)%hi(eq, c1, c2))/dy(ah2 + 1)
                    end do
                end do
            end do
        end if
        ! z-faces (p_glb > 0): transverse dims (x, y); both always active in 3D (4 children)
        if (p_glb > 0) then
            nch = 4
            do eq = 1, sys_size
                do c2 = 0, ny_c
                    f20 = 2*c2
                    do c1 = 0, nx_c
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        rhs_vf(eq)%sf(al1 + c1, al2 + c2, al3 - 1) = rhs_vf(eq)%sf(al1 + c1, al2 + c2, al3 - 1) + (creg(3)%lo(eq, &
                               & c1, c2) - fblo)/dz(al3 - 1)
                        rhs_vf(eq)%sf(al1 + c1, al2 + c2, ah3 + 1) = rhs_vf(eq)%sf(al1 + c1, al2 + c2, &
                               & ah3 + 1) + (fbhi - creg(3)%hi(eq, c1, c2))/dz(ah3 + 1)
                    end do
                end do
            end do
        end if
        ! corrected rhs back to the device for the coarse RK update kernel. M2: kernel-port removes this transfer
        do eq = 1, sys_size
            $:GPU_UPDATE(device='[rhs_vf(eq)%sf]')
        end do

    end subroutine s_amr_apply_reflux

    !> Zero the fine registers (called by the subcycle driver before substep 1 - stage-1 overwrite cannot work across two substeps).
    impure subroutine s_amr_zero_fine_registers()

        integer :: d

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        do d = 1, 3
            if (allocated(freg(d)%lo)) then
                freg(d)%lo = 0._wp; freg(d)%hi = 0._wp
            end if
        end do

    end subroutine s_amr_zero_fine_registers

    !> Berger-Colella state correction (subcycle mode only): after restriction, correct the first coarse cell OUTSIDE each patch
    !! face with the time-accumulated flux mismatch: low face: q += dt*(F_c_eff - Fbar_f_eff)/dx ; high face: q += dt*(Fbar_f_eff -
    !! F_c_eff)/dx. Registers hold EFFECTIVE (rk3_w-weighted, substep-averaged) fluxes in subcycle mode.
    impure subroutine s_amr_apply_reflux_state(q_cons)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons
        integer                                                :: eq, c1, c2, f10, f20, dd1, dd2, nch
        integer                                                :: nx_c, ny_c, nz_c, al1, al2, al3, ah1, ah2, ah3
        real(wp)                                               :: fblo, fbhi

        if (.not. amr) return
        if (.not. amr_rank_owns_patch) return
        nx_c = amr_region_hi(1) - amr_region_lo(1)
        ny_c = 0; nz_c = 0
        if (n_glb > 0) ny_c = amr_region_hi(2) - amr_region_lo(2)
        if (p_glb > 0) nz_c = amr_region_hi(3) - amr_region_lo(3)
        ! LOCAL patch bounds (q_cons/dx are rank-local; amr_region_lo/hi are global)
        al1 = amr_region_lo(1) - start_idx(1); ah1 = al1 + nx_c
        al2 = amr_region_lo(2); if (n_glb > 0) al2 = al2 - start_idx(2)
        al3 = amr_region_lo(3); if (p_glb > 0) al3 = al3 - start_idx(3)
        ah2 = al2 + ny_c; ah3 = al3 + nz_c
        nch = 1
        if (n_glb > 0) nch = nch*2
        if (p_glb > 0) nch = nch*2
        do eq = 1, sys_size
            do c2 = 0, nz_c
                f20 = 0; if (p_glb > 0) f20 = 2*c2
                do c1 = 0, ny_c
                    f10 = 0; if (n_glb > 0) f10 = 2*c1
                    fblo = 0._wp; fbhi = 0._wp
                    do dd2 = 0, merge(1, 0, p_glb > 0)
                        do dd1 = 0, merge(1, 0, n_glb > 0)
                            fblo = fblo + freg(1)%lo(eq, f10 + dd1, f20 + dd2)
                            fbhi = fbhi + freg(1)%hi(eq, f10 + dd1, f20 + dd2)
                        end do
                    end do
                    fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                    q_cons(eq)%sf(al1 - 1, al2 + c1, al3 + c2) = q_cons(eq)%sf(al1 - 1, al2 + c1, al3 + c2) + dt*(creg(1)%lo(eq, &
                           & c1, c2) - fblo)/dx(al1 - 1)
                    q_cons(eq)%sf(ah1 + 1, al2 + c1, al3 + c2) = q_cons(eq)%sf(ah1 + 1, al2 + c1, &
                           & al3 + c2) + dt*(fbhi - creg(1)%hi(eq, c1, c2))/dx(ah1 + 1)
                end do
            end do
        end do
        if (n_glb > 0) then
            nch = 2
            if (p_glb > 0) nch = nch*2
            do eq = 1, sys_size
                do c2 = 0, nz_c
                    f20 = 0; if (p_glb > 0) f20 = 2*c2
                    do c1 = 0, nx_c
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, merge(1, 0, p_glb > 0)
                            do dd1 = 0, 1
                                fblo = fblo + freg(2)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(2)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        q_cons(eq)%sf(al1 + c1, al2 - 1, al3 + c2) = q_cons(eq)%sf(al1 + c1, al2 - 1, &
                               & al3 + c2) + dt*(creg(2)%lo(eq, c1, c2) - fblo)/dy(al2 - 1)
                        q_cons(eq)%sf(al1 + c1, ah2 + 1, al3 + c2) = q_cons(eq)%sf(al1 + c1, ah2 + 1, &
                               & al3 + c2) + dt*(fbhi - creg(2)%hi(eq, c1, c2))/dy(ah2 + 1)
                    end do
                end do
            end do
        end if
        if (p_glb > 0) then
            nch = 4
            do eq = 1, sys_size
                do c2 = 0, ny_c
                    f20 = 2*c2
                    do c1 = 0, nx_c
                        f10 = 2*c1
                        fblo = 0._wp; fbhi = 0._wp
                        do dd2 = 0, 1
                            do dd1 = 0, 1
                                fblo = fblo + freg(3)%lo(eq, f10 + dd1, f20 + dd2)
                                fbhi = fbhi + freg(3)%hi(eq, f10 + dd1, f20 + dd2)
                            end do
                        end do
                        fblo = fblo/real(nch, wp); fbhi = fbhi/real(nch, wp)
                        q_cons(eq)%sf(al1 + c1, al2 + c2, al3 - 1) = q_cons(eq)%sf(al1 + c1, al2 + c2, &
                               & al3 - 1) + dt*(creg(3)%lo(eq, c1, c2) - fblo)/dz(al3 - 1)
                        q_cons(eq)%sf(al1 + c1, al2 + c2, ah3 + 1) = q_cons(eq)%sf(al1 + c1, al2 + c2, &
                               & ah3 + 1) + dt*(fbhi - creg(3)%hi(eq, c1, c2))/dz(ah3 + 1)
                    end do
                end do
            end do
        end if

    end subroutine s_amr_apply_reflux_state

    impure subroutine s_finalize_amr_registers()

        integer :: d

        if (.not. amr) return
        do d = 1, 3
            if (allocated(creg(d)%lo)) deallocate (creg(d)%lo, creg(d)%hi)
            if (allocated(freg(d)%lo)) deallocate (freg(d)%lo, freg(d)%hi)
        end do

    end subroutine s_finalize_amr_registers

end module m_amr_registers
