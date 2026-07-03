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

    implicit none

    private; public :: s_initialize_riemann_solvers_module, s_riemann_solver, s_hll_riemann_solver, s_hllc_riemann_solver, &
        & s_hlld_riemann_solver, s_lf_riemann_solver, s_finalize_riemann_solvers_module

contains

    !> Dispatch to the subroutines that are utilized to compute the Riemann problem solution. For additional information please
    !! reference: 1) s_hll_riemann_solver 2) s_hllc_riemann_solver 3) s_lf_riemann_solver 4) s_hlld_riemann_solver
    subroutine s_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, qR_prim_rsx_vf, &
                                & dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, q_prim_vf, flux_vf, flux_src_vf, &
                                & flux_gsrc_vf, norm_dir, ix, iy, iz)

        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout) :: qL_prim_rsx_vf, qR_prim_rsx_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), allocatable, dimension(:), intent(inout) :: dqL_prim_dx_vf, dqR_prim_dx_vf, dqL_prim_dy_vf, &
             & dqR_prim_dy_vf, dqL_prim_dz_vf, dqR_prim_dz_vf

        type(scalar_field), dimension(sys_size), intent(inout) :: flux_vf, flux_src_vf, flux_gsrc_vf
        integer, intent(in)                                    :: norm_dir
        type(int_bounds_info), intent(in)                      :: ix, iy, iz

        #:for NAME, NUM in [('hll', 1), ('hllc', 2), ('hlld', 4), ('lf', 5)]
            if (riemann_solver == ${NUM}$) then
                call s_${NAME}$_riemann_solver(qL_prim_rsx_vf, dqL_prim_dx_vf, dqL_prim_dy_vf, dqL_prim_dz_vf, qL_prim_vf, &
                                               & qR_prim_rsx_vf, dqR_prim_dx_vf, dqR_prim_dy_vf, dqR_prim_dz_vf, qR_prim_vf, &
                                               & q_prim_vf, flux_vf, flux_src_vf, flux_gsrc_vf, norm_dir, ix, iy, iz)
            end if
        #:endfor

    end subroutine s_riemann_solver

    !> Initialize the Riemann solvers module
    impure subroutine s_initialize_riemann_solvers_module

        ! Allocating the variables that will be utilized to formulate the left, right, and average states of the Riemann problem, as
        ! well the Riemann problem solution
        integer :: i, j

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

        @:ALLOCATE(flux_rsx_vf(-1:m, -1:n, -1:p, 1:sys_size))
        if (igr) then
            @:ALLOCATE(flux_gsrc_rsx_vf(0:0, 0:0, 0:0, 1:1))
        else
            @:ALLOCATE(flux_gsrc_rsx_vf(-1:m, -1:n, -1:p, 1:sys_size))
        end if
        @:ALLOCATE(flux_src_rsx_vf(-1:m, -1:n, -1:p, eqn_idx%adv%beg:sys_size))
        @:ALLOCATE(vel_src_rsx_vf(-1:m, -1:n, -1:p, 1:num_vels))
        if (qbmm) then
            @:ALLOCATE(mom_sp_rsx_vf(-1:m+1, -1:n+1, -1:p+1, 1:4))
        end if

        if (viscous) then
            @:ALLOCATE(Re_avg_rsx_vf(-1:m, -1:n, -1:p, 1:2))
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
        @:DEALLOCATE(flux_gsrc_rsx_vf)
        @:DEALLOCATE(Gs_rs)
        if (qbmm) then
            @:DEALLOCATE(mom_sp_rsx_vf)
        end if

    end subroutine s_finalize_riemann_solvers_module

end module m_riemann_solvers
