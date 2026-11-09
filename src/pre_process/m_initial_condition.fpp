!>
!! @file
!! @brief Contains module m_initial_condition

!> @brief Assembles initial conditions by layering prioritized patches via constructive solid geometry
module m_initial_condition

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_variables_conversion
    use m_icpp_patches
    use m_assign_variables
    use m_perturbation
    use m_chemistry
    use m_boundary_conditions

    implicit none

    type(ic_context) :: ic  !< Initial-condition state (fields, bc types, patch ids)

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_initial_condition_module

        integer :: i, j, k, l

        allocate (ic%q_prim_vf(1:sys_size))
        allocate (ic%q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (ic%q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
            allocate (ic%q_cons_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do

        if (chemistry) then
            allocate (ic%q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end if

        allocate (ic%patch_id_fp(0:m,0:n,0:p))

        if (qbmm .and. .not. polytropic) then
            allocate (pb%sf(0:m,0:n,0:p,1:nnode,1:nb))
            allocate (mv%sf(0:m,0:n,0:p,1:nnode,1:nb))
        end if

        do i = 1, sys_size
            ic%q_cons_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp) ! TODO :: remove this magic number
            ic%q_prim_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp)
        end do

        allocate (ic%bc_type(1:num_dims,1:2))

        allocate (ic%bc_type(1, 1)%sf(0:0,0:n,0:p))
        allocate (ic%bc_type(1, 2)%sf(0:0,0:n,0:p))

        do l = 0, p
            do k = 0, n
                ic%bc_type(1, 1)%sf(0, k, l) = int(min(bc_x%beg, 0), kind=1)
                ic%bc_type(1, 2)%sf(0, k, l) = int(min(bc_x%end, 0), kind=1)
            end do
        end do

        if (n > 0) then
            allocate (ic%bc_type(2, 1)%sf(-buff_size:m + buff_size,0:0,0:p))
            allocate (ic%bc_type(2, 2)%sf(-buff_size:m + buff_size,0:0,0:p))

            do l = 0, p
                do j = -buff_size, m + buff_size
                    ic%bc_type(2, 1)%sf(j, 0, l) = int(min(bc_y%beg, 0), kind=1)
                    ic%bc_type(2, 2)%sf(j, 0, l) = int(min(bc_y%end, 0), kind=1)
                end do
            end do

            if (p > 0) then
                allocate (ic%bc_type(3, 1)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (ic%bc_type(3, 2)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))

                do k = -buff_size, n + buff_size
                    do j = -buff_size, m + buff_size
                        ic%bc_type(3, 1)%sf(j, k, 0) = int(min(bc_z%beg, 0), kind=1)
                        ic%bc_type(3, 2)%sf(j, k, 0) = int(min(bc_z%end, 0), kind=1)
                    end do
                end do
            end if
        end if

        ! Initial damage state is always zero
        if (cont_damage) then
            ic%q_cons_vf(eqn_idx%damage)%sf = 0._wp
            ic%q_prim_vf(eqn_idx%damage)%sf = 0._wp
        end if

        ! Initial hyper_cleaning state is always zero TODO more general
        if (hyper_cleaning) then
            ic%q_cons_vf(eqn_idx%psi)%sf = 0._wp
            ic%q_prim_vf(eqn_idx%psi)%sf = 0._wp
        end if

        ! Setting default values for patch identities bookkeeping variable. This is necessary to avoid any confusion in the
        ! assessment of the extent of application that the overwrite permissions give a patch when it is being applied in the
        ! domain.
        ic%patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module

    !> Iterate over patches and, depending on the geometry type, call the related subroutine to setup the said geometry on the grid
    !! using the primitive variables included with the patch parameters. The subroutine is complete once the primitive variables are
    !! converted to conservative ones.
    impure subroutine s_generate_initial_condition

        integer :: i

        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(ic%q_cons_vf, ic%q_T_sf, ic%q_prim_vf, idwbuff)
        end if

        call s_apply_icpp_patches(ic%patch_id_fp, ic%q_prim_vf)

        if (num_bc_patches > 0) call s_apply_boundary_patches(ic%q_prim_vf, ic%bc_type)

        if (perturb_flow) call s_perturb_surrounding_flow(ic%q_prim_vf)
        if (perturb_sph) call s_perturb_sphere(ic%q_prim_vf)
        if (mixlayer_perturb) call s_perturb_mixlayer(ic%q_prim_vf)
        if (simplex_perturb) call s_perturb_simplex(ic%q_prim_vf)
        if (chemistry) call s_compute_T_from_primitives(ic%q_T_sf, ic%q_prim_vf, idwint)
        if (chemistry .and. ib) then
            if (any(patch_ib(1:num_ibs)%thermal_layer > 0._wp)) call s_apply_ib_thermal_layers(ic%q_prim_vf, ic%q_T_sf)
        end if

        if (elliptic_smoothing .and. chemistry) then
            call s_elliptic_smoothing(ic%q_prim_vf, ic%bc_type, ic%q_T_sf)
            call s_compute_T_from_primitives(ic%q_T_sf, ic%q_prim_vf, idwint)
        else if (elliptic_smoothing) then
            call s_elliptic_smoothing(ic%q_prim_vf, ic%bc_type)
        end if

        call s_convert_primitive_to_conservative_variables(ic%q_prim_vf, ic%q_cons_vf)

        if (qbmm .and. .not. polytropic) then
            call s_initialize_mv(ic%q_cons_vf, mv%sf)
            call s_initialize_pb(ic%q_cons_vf, mv%sf, pb%sf)
        end if

    end subroutine s_generate_initial_condition

    !> Start the gas around each prescribed-temperature circle/sphere IB in the conduction solution for a wall impulsively brought
    !! to Twall a time t0 ago, T = T_inf + (Twall - T_inf)*(R/r)^((d-1)/2)*erfc((r - R)/delta), delta = 2*sqrt(alpha*t0). Exact for
    !! a sphere and the leading-order term for a cylinder (Carslaw & Jaeger). Pressure and composition are kept, so density scales
    !! as T_inf/T. Where layers overlap, the IB with the largest temperature change sets the cell.
    impure subroutine s_apply_ib_thermal_layers(q_prim_vf, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        type(scalar_field), intent(inout)                      :: q_T_sf
        ! erfc(6) ~ 2e-17: beyond 6 delta a layer changes nothing
        real(wp), parameter :: reach = 6._wp
        real(wp)            :: x(3), c(3), span(3), r, R_ib, dT, dT_i, T_inf
        logical             :: periodic(3)
        integer             :: i, j, k, l, d

        span = [x_domain%end - x_domain%beg, y_domain%end - y_domain%beg, z_domain%end - z_domain%beg]
        periodic = [bc_x%beg == BC_PERIODIC, bc_y%beg == BC_PERIODIC, bc_z%beg == BC_PERIODIC]

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    x = 0._wp
                    x(1) = x_cc(j)
                    if (n > 0) x(2) = y_cc(k)
                    if (p > 0) x(3) = z_cc(l)
                    T_inf = q_T_sf%sf(j, k, l)
                    dT = 0._wp
                    do i = 1, num_ibs
                        if (.not. patch_ib(i)%thermal_layer > 0._wp) cycle
                        c = [patch_ib(i)%x_centroid, patch_ib(i)%y_centroid, patch_ib(i)%z_centroid]
                        do d = 1, num_dims
                            c(d) = x(d) - c(d)
                            if (periodic(d)) c(d) = c(d) - span(d)*nint(c(d)/span(d))
                        end do
                        r = norm2(c(1:num_dims))
                        R_ib = patch_ib(i)%radius
                        if (r - R_ib > reach*patch_ib(i)%thermal_layer) cycle
                        dT_i = (patch_ib(i)%Twall - T_inf)*(R_ib/max(r, R_ib))**(0.5_wp*real(num_dims - 1, &
                                & wp))*erfc(max(r - R_ib, 0._wp)/patch_ib(i)%thermal_layer)
                        if (abs(dT_i) > abs(dT)) dT = dT_i
                    end do
                    if (dT /= 0._wp) then
                        q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l) = q_prim_vf(eqn_idx%cont%beg)%sf(j, k, l)*T_inf/(T_inf + dT)
                        q_T_sf%sf(j, k, l) = T_inf + dT
                    end if
                end do
            end do
        end do

    end subroutine s_apply_ib_thermal_layers

    !> Deallocation procedures for the module
    impure subroutine s_finalize_initial_condition_module

        integer :: i

        do i = 1, sys_size
            deallocate (ic%q_prim_vf(i)%sf)
            deallocate (ic%q_cons_vf(i)%sf)
        end do

        deallocate (ic%q_prim_vf)
        deallocate (ic%q_cons_vf)

        if (chemistry) then
            deallocate (ic%q_T_sf%sf)
        end if

        deallocate (ic%patch_id_fp)

        deallocate (ic%bc_type(1, 1)%sf)
        deallocate (ic%bc_type(1, 2)%sf)

        if (n > 0) then
            deallocate (ic%bc_type(2, 1)%sf)
            deallocate (ic%bc_type(2, 2)%sf)
        end if

        if (p > 0) then
            deallocate (ic%bc_type(3, 1)%sf)
            deallocate (ic%bc_type(3, 2)%sf)
        end if

        deallocate (ic%bc_type)

    end subroutine s_finalize_initial_condition_module

end module m_initial_condition
