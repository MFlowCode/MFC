!>
!!@file
!!@brief Contains module m_amr_frame

#! AMD OpenMP lane: assert allocatables present on every kernel here (see OMP_DEFAULT_STR). A kernel naming an unallocated
#! module array aborts, so conditionally allocated ones (sw_jac/jac: igr) are named only under their allocation condition.
#:set MFC_OMP_PRESENT_ALLOCATABLE = True
#:include 'macros.fpp'

!> @brief Block-frame swap (fine grid state in/out of the shared solver).
module m_amr_frame

    use m_derived_types  ! scalar_field, t_box, int_bounds_info
    use m_global_parameters
    use m_mpi_proxy, only: s_mpi_abort  ! @:ASSERT expands to it
    use m_pressure_relaxation, only: s_pressure_relaxation_procedure
    use m_hypoelastic, only: s_hypoelastic_update_fd_coeffs, s_enforce_cont_damage_bounds
    use m_body_forces, only: s_apply_bodyforces
    use m_igr, only: jac, jac_old
    use m_amr_state
    use m_amr_distribution
    use m_amr_store
    use m_amr_exchange

    implicit none

    private
    public :: s_amr_fine_enter, s_amr_fine_leave, s_amr_bodyforces_fine, s_amr_cont_damage_fine, s_amr_pressure_relax_fine, &
        & s_amr_restore_coarse, s_amr_swap_to_fine

contains

    !> Swap to the current fine block and load its solver frame amr_cons_br; paired with s_amr_fine_leave.
    impure subroutine s_amr_fine_enter()

        integer :: loc_cur  !< plain scalar for device routines (NVHPC -Minline)

        call s_amr_swap_to_fine()
        loc_cur = amr_loc_of(amr_cur)
        call s_amr_br_load(loc_cur)

    end subroutine s_amr_fine_enter

    !> Store the solver frame back to the fine block and restore the coarse grid.
    impure subroutine s_amr_fine_leave()

        integer :: loc_cur

        loc_cur = amr_loc_of(amr_cur)
        call s_amr_br_store(loc_cur)
        call s_amr_restore_coarse()

    end subroutine s_amr_fine_leave

    !> 6-equation model: per-stage pressure relaxation on the fine block's interior (cell-local), mirroring the coarse call.
    impure subroutine s_amr_pressure_relax_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_fine_enter()
        call s_pressure_relaxation_procedure(amr_cons_br)
        call s_amr_fine_leave()

    end subroutine s_amr_pressure_relax_fine

    !> Fine-block twin of the coarse post-RK body-force update, applied in the block's own frame; without it the refined region
    !! would feel no body force. bf_spatial_support is rejected with amr (its source arrays are coarse-grid).
    impure subroutine s_amr_bodyforces_fine(ldt)

        real(wp), intent(in) :: ldt  !< the stage's rk_coef(s, 3)*dt/rk_coef(s, 4)

        if (.not. amr_rank_owns_block) return
        call s_amr_fine_enter()
        ! amr_scr_prim is passed for the interface only (the body-force source reads the conserved state and rhoM)
        call s_apply_bodyforces(amr_cons_br, amr_scr_prim, amr_scr_rhs, ldt)
        call s_amr_fine_leave()

    end subroutine s_amr_bodyforces_fine

    !> Fine-block twin of the coarse post-RK continuum-damage clamp (cell-local).
    impure subroutine s_amr_cont_damage_fine()

        if (.not. amr_rank_owns_block) return
        call s_amr_fine_enter()
        call s_enforce_cont_damage_bounds(amr_cons_br)
        call s_amr_fine_leave()

    end subroutine s_amr_cont_damage_fine

    !> Swap the global grid state to the fine block. Must be paired with s_amr_restore_coarse.
    impure subroutine s_amr_swap_to_fine(widen)

        logical, intent(in), optional :: widen  !< also widen idwint to the ghost shells (the RHS conversion range)
        integer                       :: ibm, o, e

        ! Saves are depth-guarded: a nested save would overwrite the sw_* bounce buffers with fine state, and the outermost restore
        ! would install fine extents as the coarse grid. Installs are idempotent, so unguarded.

        amr_swap_depth = amr_swap_depth + 1
        if (amr_swap_depth == 1) then
            sw_m = m; sw_n = n; sw_p = p
            sw_idwint = idwint; sw_idwbuff = idwbuff
        end if
        ! the acoustic source's spatials are coarse-grid cell indices; its support never overlaps a block (checked at startup), so
        ! the fine RHS skips it
        if (amr_swap_depth == 1) sw_acoustic_source = acoustic_source
        acoustic_source = .false.
        m = amr_slots(amr_cur)%m; n = amr_slots(amr_cur)%n; p = amr_slots(amr_cur)%p
        idwint%beg = 0; idwint%end = [m, n, p]
        idwbuff = amr_slots(amr_cur)%idwbuff
        ! save coarse coords, then install the fine coords and extend them into the ghost shell (s_build_level_coords fills 0:m)
        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (amr_dim(${D}$)) then
                if (amr_swap_depth == 1) then; sw_${X}$_cb = ${X}$_cb; sw_${X}$_cc = ${X}$_cc; sw_d${X}$ = d${X}$; end if
                ${X}$_cb(-1:${E}$) = amr_slots(amr_cur)%${X}$_cb(-1:${E}$)
                ${X}$_cc(0:${E}$) = amr_slots(amr_cur)%${X}$_cc(0:${E}$)
                d${X}$(0:${E}$) = amr_slots(amr_cur)%d${X}$(0:${E}$)
                call s_amr_extend_ghost_coords(${D}$, ${X}$_cb, ${X}$_cc, d${X}$, ${E}$, lbound(amr_g${X}$cb, 1), amr_g${X}$cb)
            end if
        #:endfor
        ! batched advance: extend the leader's grid into the slab of amr_bat_n stacked blocks (stride amr_bat_w along amr_bat_sd),
        ! since the flux divergence reads dx/dy/dz at every slab cell. Cell boundaries are not replicated: the batched path does not
        ! read them.
        if (amr_bat_n > 1) then
            do ibm = 2, amr_bat_n
                #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
                    o = amr_bat_off(${D}$, ibm); e = amr_bat_ext(${D}$)
                    if (o > 0) then
                        ${X}$_cc(o - buff_size:o + e + buff_size) = ${X}$_cc(-buff_size:e + buff_size)
                        d${X}$(o - buff_size:o + e + buff_size) = d${X}$(-buff_size:e + buff_size)
                    end if
                #:endfor
            end do
            #:for D, E in [(1, 'm'), (2, 'n'), (3, 'p')]
                if (${D}$ == amr_bat_sd) then
                    e = maxval(amr_bat_off(${D}$,1:amr_bat_n)) + amr_bat_ext(${D}$)
                    ${E}$ = e; idwint(${D}$)%end = e; idwbuff(${D}$)%end = e + buff_size
                end if
            #:endfor
        end if
        if (present(widen)) then
            if (widen) idwint = idwbuff
        end if
        ! sync extents/bounds/coordinates to the GPU_DECLARE'd device copies (stale coarse bounds = OOB kernels)
        call s_amr_sync_grid_state_to_device()
        ! hypoelastic FD coefficients depend on grid spacing: recompute them for the fine grid
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        ! IGR: save the coarse sigma and seed the fine solve with its parent values (best initial guess, and the frozen Dirichlet
        ! ghost data for the block-local Jacobi solve)
        if (igr) call s_amr_igr_swap_sigma()

    end subroutine s_amr_swap_to_fine

    !> Extend one axis of the installed fine grid into the ghost shell with the same parent-cell bisection the interior uses:
    !! boundary jg belongs to parent cell cl = isect_lo + floor(jg/rr), and k = modulo(jg, rr) picks the rr-way split. The parent
    !! boundaries are gcb for a level-1 block; a level>=2 block's parent grid is rebuilt from replicated metadata (the parent's own
    !! arrays exist only on its owner). Blocks stay buff_size inside the domain, so every ghost parent is an in-domain cell.
    impure subroutine s_amr_extend_ghost_coords(d, cb, cc, dd, ext, glb, gcb)

        integer, intent(in)     :: d, ext, glb
        real(wp), intent(inout) :: cb(-1 - buff_size:), cc(-buff_size:), dd(-buff_size:)
        real(wp), intent(in)    :: gcb(glb:)
        real(wp), allocatable   :: pcb(:), tcc(:), tdx(:)
        integer                 :: jg, cl, k, rr, pblk, pnf

        rr = amr_slots(amr_cur)%amr_ref_ratio
        if (amr_block_level(amr_cur) >= 2) then
            pblk = f_amr_parent_block(amr_cur)
            pnf = amr_ref_ratio**amr_block_level(pblk)*(amr_region_hi_all(d, pblk) - amr_region_lo_all(d, pblk) + 1) - 1
            allocate (pcb(-1:pnf), tcc(0:pnf), tdx(0:pnf))
            call s_amr_build_block_coords(pblk, gcb, glb, pcb, tcc, tdx, d)
        else
            allocate (pcb(glb:ubound(gcb, 1))); pcb = gcb
        end if
        do jg = -1 - buff_size, ext + buff_size
            if (jg >= 0 .and. jg <= ext) cycle
            cl = amr_isect_lo(d) + floor(real(jg, wp)/real(rr, wp))
            k = modulo(jg, rr)
            if (k == rr - 1) then
                cb(jg) = pcb(cl)
            else
                cb(jg) = (real(rr - 1 - k, wp)*pcb(cl - 1) + real(k + 1, wp)*pcb(cl))/real(rr, wp)
            end if
        end do
        do jg = -buff_size, ext + buff_size
            if (jg >= 0 .and. jg <= ext) cycle
            dd(jg) = cb(jg) - cb(jg - 1); cc(jg) = 0.5_wp*(cb(jg - 1) + cb(jg))
        end do

    end subroutine s_amr_extend_ghost_coords

    !> Restore the global grid state saved by s_amr_swap_to_fine.
    impure subroutine s_amr_restore_coarse(sync_device)

        !> .false. only from the batched stage when another batch follows at once (its swap re-pushes the grid state)
        logical, intent(in), optional :: sync_device
        logical                       :: sync

        @:ASSERT(amr_swap_depth > 0, "s_amr_restore_coarse without a matching s_amr_swap_to_fine")
        amr_swap_depth = amr_swap_depth - 1
        ! inner restore: the enclosing swap frame still needs its slot installed and the sw_* buffers still hold the coarse state
        if (amr_swap_depth > 0) return
        m = sw_m; n = sw_n; p = sw_p
        idwint = sw_idwint; idwbuff = sw_idwbuff
        acoustic_source = sw_acoustic_source
        ! restore full coarse coords
        #:for D, X, E in [(1, 'x', 'm'), (2, 'y', 'n'), (3, 'z', 'p')]
            if (amr_dim(${D}$)) then; ${X}$_cb = sw_${X}$_cb; ${X}$_cc = sw_${X}$_cc; d${X}$ = sw_d${X}$; end if
        #:endfor
        sync = .true.; if (present(sync_device)) sync = sync_device
        if (sync) call s_amr_sync_grid_state_to_device()
        if (hypoelasticity) call s_hypoelastic_update_fd_coeffs()
        if (igr) call s_amr_igr_restore_sigma()

    end subroutine s_amr_restore_coarse

    !> Save the coarse jac/jac_old and seed the (already swapped-in) fine block's sigma by piecewise-constant parent injection from
    !! the saved coarse sigma: interior = initial guess, ghost shell = frozen Dirichlet coupling data.
    impure subroutine s_amr_igr_swap_sigma()

        integer :: j, k, l, ci, cj, ck, ibm, g, o1, o2, o3, mb1, me1, mb2, me2, mb3, me3, lo1, lo2, lo3, ox, oy, oz
        integer :: cb1, ce1, cb2, ce2, cb3, ce3, lo(3), o(3), mb(3), me(3)

        ! bounds hoisted to scalars: sw_idwbuff etc. are host-only, so naming them in a kernel fails OpenACC's present lookup

        cb1 = sw_idwbuff(1)%beg; ce1 = sw_idwbuff(1)%end
        cb2 = sw_idwbuff(2)%beg; ce2 = sw_idwbuff(2)%end
        cb3 = sw_idwbuff(3)%beg; ce3 = sw_idwbuff(3)%end
        ox = amr_sidx(1); oy = amr_sidx(2); oz = amr_sidx(3)
        ! Save the coarse sigma, outermost swap only (a nested save would overwrite sw_jac with fine state). The seed below is
        ! unguarded: sw_jac still holds the coarse state. jac_old exists only for the Jacobi solver (igr_iter_solver == 1), so it
        ! gets its own kernel: naming an unallocated array in a launched kernel aborts under amdflang.
        if (amr_swap_depth == 1) then
            #:for J in ['jac', 'jac_old']
                #:if J == 'jac_old'
                    if (igr_iter_solver == 1) then
                    #:endif
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = cb3, ce3
                        do k = cb2, ce2
                            do j = cb1, ce1
                                sw_${J}$(j, k, l) = ${J}$(j, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    #:if J == 'jac_old'
                    end if
                #:endif
            #:endfor
        end if
        ! seed each block's buffered range from its own parent; in a batched slab member ibm sits at offset o along amr_bat_sd and
        ! its own extents bound the seed
        do ibm = 1, max(1, amr_bat_n)
            if (amr_bat_n > 1) then
                g = amr_bat_blk(ibm); lo = amr_isect_lo_all(:,g)
                o = amr_bat_off(:,ibm)
                mb = merge(-buff_size, 0, amr_dim); me = merge(amr_bat_mext(:,ibm) + buff_size, 0, amr_dim)
            else
                lo = amr_isect_lo; o = 0; mb = idwbuff%beg; me = idwbuff%end
            end if
            lo1 = lo(1); lo2 = lo(2); lo3 = lo(3); o1 = o(1); o2 = o(2); o3 = o(3)
            mb1 = mb(1); mb2 = mb(2); mb3 = mb(3); me1 = me(1); me2 = me(2); me3 = me(3)
            #:for J in ['jac', 'jac_old']
                #:if J == 'jac_old'
                    if (igr_iter_solver == 1) then
                    #:endif
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, ci, cj, ck]', copyin='[lo1, lo2, lo3, o1, o2, o3, mb1, &
                                        & me1, mb2, me2, mb3, me3]')
                    do l = mb3, me3
                        do k = mb2, me2
                            do j = mb1, me1
                                ci = lo1 + floor(real(j, wp)/real(amr_ref_ratio, wp)) - ox
                                cj = 0; ck = 0
                                if (n_glb > 0) cj = lo2 + floor(real(k, wp)/real(amr_ref_ratio, wp)) - oy
                                if (p_glb > 0) ck = lo3 + floor(real(l, wp)/real(amr_ref_ratio, wp)) - oz
                                ci = min(max(ci, cb1), ce1)
                                cj = min(max(cj, cb2), ce2)
                                ck = min(max(ck, cb3), ce3)
                                ${J}$(j + o1, k + o2, l + o3) = sw_jac(ci, cj, ck)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    #:if J == 'jac_old'
                    end if
                #:endif
            #:endfor
        end do

    end subroutine s_amr_igr_swap_sigma

    !> Restore the coarse jac/jac_old saved by s_amr_igr_swap_sigma (bounds already restored).
    impure subroutine s_amr_igr_restore_sigma()

        integer :: j, k, l, b1, e1, b2, e2, b3, e3

        b1 = idwbuff(1)%beg; e1 = idwbuff(1)%end
        b2 = idwbuff(2)%beg; e2 = idwbuff(2)%end
        b3 = idwbuff(3)%beg; e3 = idwbuff(3)%end
        #:for J in ['jac', 'jac_old']
            #:if J == 'jac_old'
                if (igr_iter_solver == 1) then
                #:endif
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = b3, e3
                    do k = b2, e2
                        do j = b1, e1
                            ${J}$(j, k, l) = sw_${J}$(j, k, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                #:if J == 'jac_old'
                end if
            #:endif
        #:endfor

    end subroutine s_amr_igr_restore_sigma

end module m_amr_frame
