!>
!! @file
!! @brief Contains module m_riemann_solvers

!> @brief Approximate and exact Riemann solvers (HLL, HLLC, HLLD, exact) for the multicomponent Navier--Stokes equations

#:include 'case.fpp'
#:include 'macros.fpp'

module m_riemann_solvers

    use m_derived_types
    use m_global_parameters
    use m_riemann_state
    use m_riemann_solver_hllc
    use m_riemann_solver_lf
    use m_riemann_solver_hll
    use m_riemann_solver_hlld
    use m_riemann_solver_hypo_hlld

    implicit none

    private; public :: s_initialize_riemann_solvers_module, s_riemann_solver, s_hll_riemann_solver, s_hllc_riemann_solver, &
        & s_hlld_riemann_solver, s_hypo_hlld_riemann_solver, s_finalize_nc_iface_vel, s_lf_riemann_solver, &
        & s_finalize_riemann_solver_hatR, s_finalize_nc_iface_vel_hatR, s_finalize_riemann_solvers_module

contains

    !> Dispatch to the subroutines that are utilized to compute the Riemann problem solution. For additional information please
    !! reference: 1) s_hll_riemann_solver 2) s_hllc_riemann_solver 3) s_lf_riemann_solver 4) s_hlld_riemann_solver
    subroutine s_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        integer, intent(in)               :: norm_dir
        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! Hypoelasticity enters the Riemann layer in THREE distinct code shapes:
        !   1. HLLC - inline "if (hypoelasticity)" branches inside s_hllc_riemann_solver
        !   2. HLL  - inline "if (hypoelasticity)" branches inside s_hll_riemann_solver
        !   3. HLLD - a separate module (m_riemann_solver_hypo_hlld), called directly from m_rhs (s_compute_directional_rhs)
        !      under hypo_nc_mode_dual_pass
        ! HLLD needs its own path because its anchored dual pass produces BOTH the hat_L and hat_R anchored flux
        ! sets in one fused solve, whose partial RHS are then summed in m_rhs; HLLC/HLL instead add their
        ! non-conservative contribution within a single-pass solve. See
        ! misc/dev_notes/Riemann_and_RHS_source_terms_explanations.md (S5.3).

        #:for NAME, NUM in [('hll', 1), ('hllc', 2), ('hlld', 4), ('lf', 5)]
            if (riemann_solver == ${NUM}$) then
                call s_${NAME}$_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
                                               & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
                                               & q_prim_vf, norm_dir, ix, iy, iz)
            end if
        #:endfor

    end subroutine s_riemann_solver

    !> Initialize the Riemann solvers module
    impure subroutine s_initialize_riemann_solvers_module

        ! Allocating the variables that will be utilized to formulate the left, right, and average states of the Riemann problem, as
        ! well the Riemann problem solution
        integer :: i, j, k, l, src_lo

        @:ALLOCATE(Gs_rs(1:num_fluids))

        do i = 1, num_fluids
            Gs_rs(i) = fluid_pp(i)%G
        end do
        $:GPU_UPDATE(device='[Gs_rs]')

        if (viscous) then
            @:ALLOCATE(Res_gs(1:2, 1:Re_size_max))
        end if

        if (viscous) then
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_gs(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do
            $:GPU_UPDATE(device='[Res_gs, Re_idx, Re_size]')
        end if

        $:GPU_ENTER_DATA(copyin='[is1, is2, is3, isx, isy, isz]')

        is1%beg = -1; is2%beg = 0; is3%beg = 0
        is1%end = m; is2%end = n; is3%end = p

        @:ALLOCATE(flux_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:sys_size))
        @:ALLOCATE(vel_src_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:num_vels))

        ! Size the source-flux buffer to the band that is actually written. These are FULL-DOMAIN arrays, so each unused component
        ! costs (m_alloc+2)(n_alloc+2)(p_alloc+2) reals per rank - 0.37 GB/rank of waste at 400^3 for the five components an
        ! inviscid Cartesian run never touches.
        !   chemistry diffusion : from 1, because m_chemistry lives in src/common, cannot use m_riemann_state, and so takes a flat
        !                         dummy declared `dimension(-1:, -1:, -1:, 1:)` - the lower bounds must agree or every species
        !                         index silently shifts. It is only ever passed this array when diffusion is on.
        !   viscous / surf.tens.: from mom%beg (the viscous stress and work fluxes occupy mom..E)
        !   otherwise           : from adv%beg (the advection source band alone)
        if (chemistry .and. chem_params%diffusion) then
            src_lo = 1
        else if (viscous .or. surface_tension) then
            src_lo = eqn_idx%mom%beg
        else
            src_lo = eqn_idx%adv%beg
        end if
        @:ALLOCATE(flux_src_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, src_lo:sys_size))

        ! The geometric source flux exists only on the cylindrical/axisymmetric paths - every write in the four solvers and both
        ! reads in m_rhs sit under `cyl_coord` or `grid_geometry == 3`, and grid_geometry == 3 implies cyl_coord
        ! (m_global_parameters). A Cartesian run therefore never touches it, so do not pay 6 full-domain arrays for it. It is
        ! zeroed on allocation because the solvers write only the components they touch while m_rhs reads the whole band.
        if (cyl_coord) then
            @:ALLOCATE(flux_gsrc_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:sys_size))
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = -1, p_alloc
                    do k = -1, n_alloc
                        do j = -1, m_alloc
                            flux_gsrc_rsx_vf(j, k, l, i) = 0._wp
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (qbmm) then
            @:ALLOCATE(mom_sp_rsx_vf(-1:m_alloc+1, -1:n_alloc+1, -1:p_alloc+1, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:2))
        end if

        ! _alloc bounds like every rs sibling above: the AMR fine advance swaps m/n/p to fine-block
        ! extents that can exceed the coarse subdomain when amr_max_grid_size pins a larger block
        if (use_nc_iface_vel) then
            @:ALLOCATE(nc_iface_vel_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:num_dims))
        end if

        if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
            @:ALLOCATE(flux_hatR_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:sys_size))
            if (use_nc_iface_vel) then
                @:ALLOCATE(nc_iface_vel_hatR_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:num_dims))
            end if
            if (cyl_coord) then
                @:ALLOCATE(flux_gsrc_hatR_rsx_vf(-1:m_alloc, -1:n_alloc, -1:p_alloc, 1:sys_size))
            end if
        end if

    end subroutine s_initialize_riemann_solvers_module

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_riemann_solvers_module

        if (viscous) then
            @:DEALLOCATE(Re_avg_rsx_vf)
            @:DEALLOCATE(Res_gs)
        end if
        @:DEALLOCATE(vel_src_rsx_vf)
        @:DEALLOCATE(flux_rsx_vf)
        @:DEALLOCATE(flux_src_rsx_vf)
        if (cyl_coord) then
            @:DEALLOCATE(flux_gsrc_rsx_vf)
        end if
        @:DEALLOCATE(Gs_rs)
        if (use_nc_iface_vel) then
            @:DEALLOCATE(nc_iface_vel_rsx_vf)
        end if
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsx_vf)
        end if
        if (hypo_nc_mode == hypo_nc_mode_dual_pass) then
            @:DEALLOCATE(flux_hatR_rsx_vf)
            if (use_nc_iface_vel) then
                @:DEALLOCATE(nc_iface_vel_hatR_rsx_vf)
            end if
            if (cyl_coord) then
                @:DEALLOCATE(flux_gsrc_hatR_rsx_vf)
            end if
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
