!>
!!@file m_checker_common.f90
!!@brief Contains module m_checker_common

#:include 'macros.fpp'

!> @brief The purpose of the module is to check for compatible input files for.
!!              inputs common to pre-processing, post-processing and simulation
module m_checker_common

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs_common, wp

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_PRE_PROCESS
        call s_check_inputs_time_stepping
        call s_check_inputs_finite_difference
#endif

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif

        ! Run by all three stages
        call s_check_inputs_simulation_domain
        call s_check_inputs_model_eqns_and_num_fluids
        if (igr) then
            call s_check_inputs_igr
        else
#ifndef MFC_POST_PROCESS
            call s_check_inputs_bubbles_euler
            call s_check_inputs_qbmm_and_polydisperse
            call s_check_inputs_adv_n
            call s_check_inputs_hypoelasticity
            call s_check_inputs_phase_change
            call s_check_inputs_ibm
#endif
            if (recon_type == WENO_TYPE) then
                call s_check_inputs_weno
            elseif (recon_type == MUSCL_TYPE) then
                call s_check_inputs_muscl
            end if
            call s_check_inputs_surface_tension
            call s_check_inputs_mhd
        end if
        call s_check_inputs_bc
        call s_check_inputs_stiffened_eos
        call s_check_inputs_moving_bc

    end subroutine s_check_inputs_common

#ifndef MFC_PRE_PROCESS

    !> Checks constraints on the time-stepping parameters.
        !! Called by s_check_inputs_common for simulation and post-processing
    impure subroutine s_check_inputs_time_stepping
        if (cfl_dt) then
            @:PROHIBIT(cfl_target < 0 .or. cfl_target > 1._wp)
            @:PROHIBIT(t_stop <= 0)
            @:PROHIBIT(t_save <= 0)
            @:PROHIBIT(t_save > t_stop)
            @:PROHIBIT(n_start < 0)
        else
            @:PROHIBIT(t_step_start < 0)
            @:PROHIBIT(t_step_stop <= t_step_start)
            @:PROHIBIT(t_step_save > t_step_stop - t_step_start)
        end if
    end subroutine s_check_inputs_time_stepping

    !> Checks constraints on the finite difference parameters.
        !! Called by s_check_inputs_common for simulation and post-processing
    impure subroutine s_check_inputs_finite_difference
        @:PROHIBIT(all(fd_order /= (/dflt_int, 1, 2, 4/)), "fd_order must be 1, 2, or 4")
    end subroutine s_check_inputs_finite_difference

#endif

#ifndef MFC_SIMULATION

    ! Checks constraints on the total number of cells
    impure subroutine s_check_total_cells
        character(len=18) :: numStr !< for int to string conversion
        integer(kind=8) :: min_cells

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(nGlobal < min_cells, &
            "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // &
            "which is currently "//trim(numStr))
    end subroutine s_check_total_cells

#endif

#ifndef MFC_POST_PROCESS

    !> Checks constraints on the bubble parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_bubbles_euler
        @:PROHIBIT(bubbles_euler .and. nb < 1, "The Ensemble-Averaged Bubble Model requires nb >= 1")
        @:PROHIBIT(bubbles_euler .and. polydisperse .and. (nb == 1), "Polydisperse bubble dynamics requires nb > 1")
        @:PROHIBIT(bubbles_euler .and. polydisperse .and. (mod(nb, 2) == 0), "nb must be odd")
        @:PROHIBIT(bubbles_euler .and. (.not. polytropic) .and. f_is_default(R0ref), "R0ref must be set if using bubbles_euler with polytropic = .false.")
        @:PROHIBIT(bubbles_euler .and. nb == dflt_int, "nb must be set if using bubbles_euler")
        @:PROHIBIT(bubbles_euler .and. thermal > 3)
        @:PROHIBIT(bubbles_euler .and. model_eqns == 3, "Bubble models untested with 6-equation model (model_eqns = 3)")
        @:PROHIBIT(bubbles_euler .and. model_eqns == 1, "Bubble models untested with pi-gamma model (model_eqns = 1)")
        @:PROHIBIT(bubbles_euler .and. model_eqns == 4 .and. f_is_default(rhoref), "rhoref must be set if using bubbles_euler with model_eqns = 4")
        @:PROHIBIT(bubbles_euler .and. model_eqns == 4 .and. f_is_default(pref), "pref must be set if using bubbles_euler with model_eqns = 4")
        @:PROHIBIT(bubbles_euler .and. model_eqns == 4 .and. num_fluids /= 1, "4-equation model (model_eqns = 4) is single-component and requires num_fluids = 1")
        @:PROHIBIT(bubbles_euler .and. cyl_coord, "Bubble models untested in cylindrical coordinates")
    end subroutine s_check_inputs_bubbles_euler

    !> Checks constraints on the QBMM and polydisperse bubble parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_qbmm_and_polydisperse
        @:PROHIBIT(polydisperse .and. (.not. bubbles_euler), "Polydisperse bubble modeling requires the bubbles_euler flag to be set")
        @:PROHIBIT(polydisperse .and. f_is_default(poly_sigma), "Polydisperse bubble modeling requires poly_sigma to be set")
        @:PROHIBIT(polydisperse .and. poly_sigma <= 0)
        @:PROHIBIT(qbmm .and. (.not. bubbles_euler), "QBMM requires the bubbles_euler flag to be set")
        @:PROHIBIT(qbmm .and. nnode /= 4)
    end subroutine s_check_inputs_qbmm_and_polydisperse

    !> Checks constraints on the adv_n flag.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_adv_n
        @:PROHIBIT(adv_n .and. (.not. bubbles_euler))
        @:PROHIBIT(adv_n .and. num_fluids /= 1)
        @:PROHIBIT(adv_n .and. qbmm)
    end subroutine s_check_inputs_adv_n

    !> Checks constraints on the hypoelasticity parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_hypoelasticity
        @:PROHIBIT(hypoelasticity .and. model_eqns /= 2)
#ifdef MFC_SIMULATION
        @:PROHIBIT(elasticity .and. fd_order /= 4)
#endif
    end subroutine s_check_inputs_hypoelasticity

    !> Checks constraints on the hyperelasticity parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_hyperelasticity
        @:PROHIBIT(hyperelasticity .and. model_eqns == 1)
        @:PROHIBIT(hyperelasticity .and. model_eqns > 3)
#ifdef MFC_SIMULATION
        @:PROHIBIT(elasticity .and. fd_order /= 4)
#endif
    end subroutine s_check_inputs_hyperelasticity

    !> Checks constraints on the phase change parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_phase_change
        @:PROHIBIT(relax .and. model_eqns /= 3, "phase change requires model_eqns = 3")
        @:PROHIBIT(relax .and. relax_model < 0, "relax_model must be in between 0 and 6")
        @:PROHIBIT(relax .and. relax_model > 6, "relax_model must be in between 0 and 6")
        @:PROHIBIT(relax .and. palpha_eps <= 0._wp, "palpha_eps must be positive")
        @:PROHIBIT(relax .and. palpha_eps >= 1._wp, "palpha_eps must be less than 1")
        @:PROHIBIT(relax .and. ptgalpha_eps <= 0._wp, "ptgalpha_eps must be positive")
        @:PROHIBIT(relax .and. ptgalpha_eps >= 1._wp, "ptgalpha_eps must be less than 1")
        @:PROHIBIT((.not. relax) .and. &
            ((relax_model /= dflt_int) .or. (.not. f_is_default(palpha_eps)) .or. (.not. f_is_default(ptgalpha_eps))), &
            "relax is not set as true, but other phase change parameters have been modified. " // &
            "Either activate phase change or set the values to default")
    end subroutine s_check_inputs_phase_change

    !> Checks constraints on the Immersed Boundaries parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    impure subroutine s_check_inputs_ibm
        @:PROHIBIT(ib .and. n <= 0, "Immersed Boundaries do not work in 1D")
        @:PROHIBIT(ib .and. (num_ibs <= 0 .or. num_ibs > num_patches_max), "num_ibs must be between 1 and num_patches_max")
        @:PROHIBIT((.not. ib) .and. num_ibs > 0, "num_ibs is set, but ib is not enabled")
    end subroutine s_check_inputs_ibm

#endif

    !> Checks constraints on dimensionality and the number of cells for the grid.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_simulation_domain
        @:PROHIBIT(m == dflt_int, "m must be set")
        @:PROHIBIT(n == dflt_int, "n must be set")
        @:PROHIBIT(p == dflt_int, "p must be set")
        @:PROHIBIT(m <= 0)
        @:PROHIBIT(n < 0)
        @:PROHIBIT(p < 0)
        @:PROHIBIT(cyl_coord .and. p > 0 .and. mod(p, 2) /= 1, "p must be odd for cylindrical coordinates")
        @:PROHIBIT(n == 0 .and. p > 0, "p must be 0 if n = 0")
    end subroutine s_check_inputs_simulation_domain

    !> Checks constraints on model equations and number of fluids in the flow.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_model_eqns_and_num_fluids
        @:PROHIBIT(all(model_eqns /= (/1, 2, 3, 4/)), "model_eqns must be 1, 2, 3, or 4")
        @:PROHIBIT(num_fluids /= dflt_int .and. num_fluids < 1, "num_fluids must be positive")
        @:PROHIBIT(model_eqns == 1 .and. num_fluids /= dflt_int, "num_fluids is not supported for model_eqns = 1")
        @:PROHIBIT(model_eqns == 2 .and. num_fluids == dflt_int, "5-equation model (model_eqns = 2) requires num_fluids to be set")
        @:PROHIBIT(model_eqns == 3 .and. num_fluids == dflt_int, "6-equation model (model_eqns = 3) requires num_fluids to be set")
        @:PROHIBIT(model_eqns == 1 .and. mpp_lim)
        @:PROHIBIT(num_fluids == 1 .and. mpp_lim)
        @:PROHIBIT(model_eqns == 3 .and. cyl_coord .and. p /= 0, &
            "6-equation model (model_eqns = 3) does not support cylindrical coordinates (cyl_coord = T and p != 0)")
    end subroutine s_check_inputs_model_eqns_and_num_fluids

    !> Checks constraints regarding IGR order.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_igr
        @:PROHIBIT(all(igr_order /= (/3, 5/)), "igr_order must be 3 or 5")
        @:PROHIBIT(m + 1 < igr_order, "m must be at least igr_order - 1")
        @:PROHIBIT(n > 0 .and. n + 1 < igr_order, "n must be at least igr_order - 1")
        @:PROHIBIT(p > 0 .and. p + 1 < igr_order, "p must be at least igr_order - 1")
    end subroutine s_check_inputs_igr

    !> Checks constraints regarding WENO order.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_weno
        @:PROHIBIT(all(weno_order /= (/1, 3, 5, 7/)), "weno_order must be 1, 3, 5, or 7")
        @:PROHIBIT(m + 1 < weno_order, "m must be at least weno_order - 1")
        @:PROHIBIT(n > 0 .and. n + 1 < weno_order, "n must be at least weno_order - 1")
        @:PROHIBIT(p > 0 .and. p + 1 < weno_order, "p must be at least weno_order - 1")
    end subroutine s_check_inputs_weno

    !> Check constraints regarding MUSCL order
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_muscl
        @:PROHIBIT(all(muscl_order /= (/1, 2/)), "muscl_order must be 1, or 2")
        @:PROHIBIT(m + 1 < muscl_order, "m must be at least muscl_order - 1")
        @:PROHIBIT(n > 0 .and. n + 1 < muscl_order, "n must be at least muscl_order - 1")
        @:PROHIBIT(p > 0 .and. p + 1 < muscl_order, "p must be at least muscl_order - 1")
    end subroutine s_check_inputs_muscl

    !> Checks constraints on the boundary conditions in the x-direction.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_bc
        logical :: skip_check !< Flag to skip the check when iterating over
        !! x, y, and z directions, for special treatment of cylindrical coordinates

        #:for X, VAR in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            #:for BOUND in ['beg', 'end']
                @:PROHIBIT(${VAR}$ == 0 .and. bc_${X}$%${BOUND}$ /= dflt_int, "bc_${X}$%${BOUND}$ is not supported for ${VAR}$ = 0")
                @:PROHIBIT(${VAR}$ > 0 .and. bc_${X}$%${BOUND}$ == dflt_int, "${VAR}$ != 0 but bc_${X}$%${BOUND}$ is not set")
                @:PROHIBIT((bc_${X}$%beg == BC_PERIODIC .and. bc_${X}$%end /= BC_PERIODIC) .or. &
                    (bc_${X}$%end == BC_PERIODIC .and. bc_${X}$%beg /= BC_PERIODIC), &
                    "bc_${X}$%beg and bc_${X}$%end must be both periodic (= -1) or both non-periodic")

                ! For cylindrical coordinates, y and z directions use a different check
                #:if (X == 'y') or (X == 'z')
                    skip_check = cyl_coord
                #:else
                    skip_check = .false.
                #:endif

                if (.not. skip_check) then
                    @:PROHIBIT(bc_${X}$%${BOUND}$ /= dflt_int .and. (bc_${X}$%${BOUND}$ > -1 .or. bc_${X}$%${BOUND}$ < BC_DIRICHLET), &
                        "bc_${X}$%${BOUND}$ must be between -1 and -17")

                    @:PROHIBIT(bc_${X}$%${BOUND}$ /= dflt_int .and. bc_${X}$%${BOUND}$ == BC_AXIS, &
                        "bc_${X}$%${BOUND}$ must not be -14 for non-cylindrical coordinates")
                end if

            #:endfor
        #:endfor

        @:PROHIBIT(any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == BC_NULL), &
            "Boundary condition -13 is not supported")

        ! Check for y and z directions for cylindrical coordinates
        @: PROHIBIT(cyl_coord .and. n == 0, "n must be positive (2D or 3D) for cylindrical coordinates")
        @: PROHIBIT(cyl_coord .and. p == 0 .and. bc_y%beg /= BC_REFLECTIVE, "bc_y%beg must be -2 for 2D cylindrical coordinates (p = 0)")
        @: PROHIBIT(cyl_coord .and. p > 0 .and. bc_y%beg /= BC_AXIS, "bc_y%beg must be -14 for 3D cylindrical coordinates (p > 0)")
        @: PROHIBIT(cyl_coord .and. (bc_y%end > BC_PERIODIC .or. bc_y%end < BC_DIRICHLET), "bc_y%end must be between -1 and -17")
        @: PROHIBIT(cyl_coord .and. bc_y%end == BC_AXIS, "bc_y%end must not be -14")

        ! Check for y and z directions for 3D cylindrical coordinates
        @: PROHIBIT(cyl_coord .and. p > 0 .and. (bc_z%beg /= BC_PERIODIC .and. bc_z%beg /= BC_REFLECTIVE), &
            "bc_z%beg must be -1 or -2 for 3D cylindrical coordinates")

        @: PROHIBIT(cyl_coord .and. p > 0 .and. (bc_z%end /= BC_PERIODIC .and. bc_z%end /= BC_REFLECTIVE), &
            "bc_z%end must be -1 or -2 for 3D cylindrical coordinates")

#ifndef MFC_POST_PROCESS
        if (num_bc_patches > 0) then
            #:for DIR in [('x'), ('y'), ('z')]
                #:for LOC in [('beg'), ('end')]
                    @:PROHIBIT(bc_${DIR}$%${LOC}$ == -1 .or. (bc_${DIR}$%${LOC}$ <= -4 .and. bc_${DIR}$%${LOC}$ >= -14), &
                        "bc_${DIR}$%${LOC}$ is not compatible with num_bc_patches > 0")
                #:endfor
            #:endfor
        end if
#endif

    end subroutine s_check_inputs_bc

    !> Checks constraints on the stiffened equation of state fluids parameters.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_stiffened_eos
        character(len=5) :: iStr !< for int to string conversion
        integer :: bub_fac !< For allowing an extra fluid_pp if there are subgrid bubbles_euler
        integer :: i

        bub_fac = 0
        if (bubbles_euler .and. (num_fluids == 1)) bub_fac = 1

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            @:PROHIBIT(.not. f_is_default(fluid_pp(i)%gamma) .and. fluid_pp(i)%gamma <= 0._wp, &
                "fluid_pp("//trim(iStr)//")%gamma must be positive")

            @:PROHIBIT(model_eqns == 1 .and. (.not. f_is_default(fluid_pp(i)%gamma)), &
                "model_eqns = 1 does not support fluid_pp("//trim(iStr)//")%gamma")

            @:PROHIBIT((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0._wp) .or. &
                (i > num_fluids + bub_fac .and. (.not. f_is_default(fluid_pp(i)%gamma))), &
                "for fluid_pp("//trim(iStr)//")%gamma")

            @:PROHIBIT(.not. f_is_default(fluid_pp(i)%pi_inf) .and. fluid_pp(i)%pi_inf < 0._wp, &
                "fluid_pp("//trim(iStr)//")%pi_inf must be non-negative")

            @:PROHIBIT(model_eqns == 1 .and. (.not. f_is_default(fluid_pp(i)%pi_inf)), &
                "model_eqns = 1 does not support fluid_pp("//trim(iStr)//")%pi_inf")

            @:PROHIBIT((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0._wp) .or. &
                (i > num_fluids + bub_fac .and. (.not. f_is_default(fluid_pp(i)%pi_inf))), &
                "for fluid_pp("//trim(iStr)//")%pi_inf")

            @:PROHIBIT(fluid_pp(i)%cv < 0._wp, &
                "fluid_pp("//trim(iStr)//")%cv must be positive")
        end do
    end subroutine s_check_inputs_stiffened_eos

    !> Checks constraints on the surface tension parameters.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_surface_tension

        integer :: i

        @:PROHIBIT(surface_tension .and. sigma < 0._wp, &
            "sigma must be greater than or equal to zero")

        @:PROHIBIT(surface_tension .and. f_approx_equal(sigma, dflt_real), &
            "sigma must be set if surface_tension is enabled")

        @:PROHIBIT(.not. f_is_default(sigma) .and. .not. surface_tension, &
            "sigma is set but surface_tension is not enabled")

        @:PROHIBIT(surface_tension .and. (model_eqns /= 3 .and. model_eqns /=2), &
            "The surface tension model requires model_eqns=3 or model_eqns=2")

        @:PROHIBIT(surface_tension .and. num_fluids /= 2, &
            "The surface tension model requires num_fluids=2")

#ifdef MFC_PRE_PROCESS
        do i = 1, num_patches
            @:PROHIBIT(surface_tension .and. f_is_default(patch_icpp(i)%cf_val), &
                "patch_icpp(i)%cf_val must be set if surface_tension is enabled")
        end do
#endif MFC_PRE_PROCESS

    end subroutine s_check_inputs_surface_tension

    !> Checks constraints on the inputs for moving boundaries.
        !! Called by s_check_inputs_common for all three stages
    impure subroutine s_check_inputs_moving_bc
        #:for X, VB2, VB3 in [('x', 'vb2', 'vb3'), ('y', 'vb3', 'vb1'), ('z', 'vb1', 'vb2')]
            if (.not. (f_approx_equal(bc_${X}$%vb1, 0._wp) .and. &
                       f_approx_equal(bc_${X}$%vb2, 0._wp) .and. &
                       f_approx_equal(bc_${X}$%vb3, 0._wp))) then
                if (bc_${X}$%beg == BC_SLIP_WALL) then
                    if (.not. (f_approx_equal(bc_${X}$%${VB2}$, 0._wp) .and. &
                               f_approx_equal(bc_${X}$%${VB3}$, 0._wp))) then
                        call s_mpi_abort("bc_${X}$%beg must be -15 if "// &
                                         "bc_${X}$%${VB2}$ or bc_${X}$%${VB3}$ "// &
                                         "is set. Exiting.", CASE_FILE_ERROR_CODE)
                    end if
                elseif (bc_${X}$%beg /= BC_NO_SLIP_WALL) then
                    call s_mpi_abort("bc_${X}$%beg must be -15 or -16 if "// &
                                     "bc_${X}$%vb[1,2,3] is set. Exiting.", CASE_FILE_ERROR_CODE)
                end if
            end if
        #:endfor

        #:for X, VE2, VE3 in [('x', 've2', 've3'), ('y', 've3', 've1'), ('z', 've1', 've2')]
            if (.not. (f_approx_equal(bc_${X}$%ve1, 0._wp) .and. &
                       f_approx_equal(bc_${X}$%ve2, 0._wp) .and. &
                       f_approx_equal(bc_${X}$%ve3, 0._wp))) then
                if (bc_${X}$%end == BC_SLIP_WALL) then
                    if (.not. (f_approx_equal(bc_${X}$%${VE2}$, 0._wp) .and. &
                               f_approx_equal(bc_${X}$%${VE3}$, 0._wp))) then
                        call s_mpi_abort("bc_${X}$%end must be -15 if "// &
                                         "bc_${X}$%${VE2}$ or bc_${X}$%${VE3}$ "// &
                                         "is set. Exiting.", CASE_FILE_ERROR_CODE)
                    end if
                elseif (bc_${X}$%end /= BC_NO_SLIP_WALL) then
                    call s_mpi_abort("bc_${X}$%end must be -15 or -16 if "// &
                                     "bc_${X}$%ve[1,2,3] is set. Exiting.", CASE_FILE_ERROR_CODE)
                end if
            end if
        #:endfor
    end subroutine s_check_inputs_moving_bc

    impure subroutine s_check_inputs_mhd
        @:PROHIBIT(mhd .and. num_fluids /= 1, "MHD is only available for single-component flows")
        @:PROHIBIT(mhd .and. model_eqns /= 2, "MHD is only available for the 5-equation model")

        @:PROHIBIT(relativity .and. .not. mhd)

        @:PROHIBIT(.not. mhd .and. (.not. f_is_default(Bx0)), "Bx0 must not be set if MHD is not enabled")
        @:PROHIBIT(mhd .and. n == 0 .and. f_is_default(Bx0), "Bx0 must be set in 1D MHD simulations")
        @:PROHIBIT(mhd .and. n > 0 .and. (.not. f_is_default(Bx0)), "Bx0 must not be set in 2D/3D MHD simulations")
    end subroutine s_check_inputs_mhd

end module m_checker_common
