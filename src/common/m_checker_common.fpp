!>
!!@file m_checker_common.f90
!!@brief Contains module m_checker_common

!> @brief The purpose of the module is to check for compatible input files for.
!!              inputs common to pre-processing, post-processing and simulation
module m_checker_common

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by all three stages
    subroutine s_check_inputs_common

#ifndef MFC_PRE_PROCESS
        call s_check_inputs_time_stepping
        call s_check_inputs_finite_difference
#endif

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif

#ifndef MFC_POST_PROCESS
        if (bubbles) call s_check_inputs_bubbles
        call s_check_inputs_qbmm_and_polydisperse
        if (adv_n) call s_check_inputs_adv_n
        if (hypoelasticity) call s_check_inputs_hypoelasticity
        call s_check_inputs_phase_change
        call s_check_inputs_ibm
#endif

        ! Run by all three stages
        call s_check_inputs_simulation_domain
        call s_check_inputs_model_eqns_and_num_fluids
        call s_check_inputs_weno
        call s_check_inputs_bc
        call s_check_inputs_stiffened_eos
        call s_check_inputs_surface_tension
        call s_check_inputs_moving_bc

    end subroutine s_check_inputs_common

#ifndef MFC_PRE_PROCESS

    !> Checks constraints on the time-stepping parameters.
        !! Called by s_check_inputs_common for simulation and post-processing
    subroutine s_check_inputs_time_stepping
        if (t_step_start < 0) then
            call s_mpi_abort('t_step_start must be non-negative. Exiting ...')
        elseif (t_step_stop <= t_step_start) then
            call s_mpi_abort('t_step_stop must be greater than t_step_start. '// &
                             'Exiting ...')
        elseif (t_step_save > t_step_stop - t_step_start) then
            call s_mpi_abort('t_step_save must be less or equal to '// &
                             '(t_step_stop - t_step_start). Exiting ...')
        end if
    end subroutine s_check_inputs_time_stepping

    !> Checks constraints on the finite difference parameters.
        !! Called by s_check_inputs_common for simulation and post-processing
    subroutine s_check_inputs_finite_difference
        if (all(fd_order /= (/dflt_int, 1, 2, 4/))) then
            call s_mpi_abort('fd_order must be 1, 2, or 4. Exiting ...')
        end if
    end subroutine s_check_inputs_finite_difference

#endif

#ifndef MFC_SIMULATION

    ! Checks constraints on the total number of cells
    subroutine s_check_total_cells
        character(len=5) :: numStr !< for int to string conversion

        if (nGlobal < 2**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)
            call s_mpi_abort('Total number of cells must be at least '// &
                             '(2^[number of dimensions])*num_procs, which is currently '// &
                             trim(numStr)//'. Exiting ...')
        end if
    end subroutine s_check_total_cells

#endif

#ifndef MFC_POST_PROCESS

    !> Checks constraints on the bubble parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_bubbles
        if (nb < 1) then
            call s_mpi_abort('The Ensemble-Averaged Bubble Model '// &
                             'requires nb >= 1. Exiting ...')
        elseif (polydisperse .and. (nb == 1)) then
            call s_mpi_abort('Polydisperse bubble dynamics requires nb > 1 '// &
                             'Exiting ...')
        elseif (polydisperse .and. (mod(nb, 2) == 0)) then
            call s_mpi_abort('nb must be odd '// &
                             'Exiting ...')
        elseif ((.not. polytropic) .and. f_is_default(R0ref)) then
            call s_mpi_abort('R0ref must be set if using bubbles with '// &
                             'polytropic = .false.. Exiting ...')
        elseif (nb == dflt_int) then
            call s_mpi_abort('nb must be set if using bubbles. Exiting ...')
        elseif (thermal > 3) then
            call s_mpi_abort('thermal must be less than 4 if using bubbles. '// &
                             'Exiting ...')
        elseif (model_eqns == 3) then
            call s_mpi_abort('Bubble models untested with '// &
                             '6-equation model (model_eqns = 3). Exiting ...')
        elseif (model_eqns == 1) then
            call s_mpi_abort('Bubble models untested with '// &
                             'pi-gamma model (model_eqns = 1). Exiting ...')
            !TODO: Comment this out when testing riemann with hll
        elseif (model_eqns == 4 .and. f_is_default(rhoref)) then
            call s_mpi_abort('rhoref must be set if using bubbles with '// &
                             'model_eqns = 4. Exiting ...')
        elseif (model_eqns == 4 .and. f_is_default(pref)) then
            call s_mpi_abort('pref must be set if using bubbles with '// &
                             'model_eqns = 4. Exiting ...')
        elseif (model_eqns == 4 .and. num_fluids > 1) then
            call s_mpi_abort('4-equation model (model_eqns = 4) is '// &
                             'single-component and requires num_fluids = 1. '// &
                             'Exiting ...')
        elseif (cyl_coord) then
            call s_mpi_abort('Bubble models untested in cylindrical coordinates')
        end if
    end subroutine s_check_inputs_bubbles

    !> Checks constraints on the QBMM and polydisperse bubble parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_qbmm_and_polydisperse
        if ((.not. bubbles) .and. polydisperse) then
            call s_mpi_abort('Polydisperse bubble modeling requires the '// &
                             'bubbles flag to be set. Exiting ...')
        elseif (polydisperse .and. f_is_default(poly_sigma)) then
            call s_mpi_abort('Polydisperse bubble modeling requires '// &
                             'poly_sigma > 0. Exiting ...')
        elseif (qbmm .and. (.not. bubbles)) then
            call s_mpi_abort('QBMM is enabled but bubbles are not. Exiting ...')
        elseif (qbmm .and. (nnode /= 4)) then
            call s_mpi_abort('nnode not supported. Exiting ...')
        end if
    end subroutine s_check_inputs_qbmm_and_polydisperse

    !> Checks constraints on the adv_n flag.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_adv_n
        if (.not. bubbles) then
            call s_mpi_abort('adv_n requires bubbles = true.'// &
                             'Exiting ...')
        else if (num_fluids > 1) then
            call s_mpi_abort('adv_n requires num_fluids = 1. '// &
                             'Exiting ...')
        else if (qbmm) then
            call s_mpi_abort('adv_n is incompatible with qbmm.'// &
                             'Exiting ...')
        end if
    end subroutine

    !> Checks constraints on the hypoelasticity parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_hypoelasticity
        if (model_eqns /= 2) then
            call s_mpi_abort('hypoelasticity requires 5-equation model'// &
                             '(model_eqns = 2). Exiting ...')
        end if
    end subroutine s_check_inputs_hypoelasticity

    !> Checks constraints on the phase change parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_phase_change
        if (relax) then
            if (model_eqns /= 3) then
                call s_mpi_abort('phase change requires model_eqns = 3. '// &
                                 'Exiting ...')
            elseif ((relax_model < 0) .or. (relax_model > 6)) then
                call s_mpi_abort('relax_model should be in between 0 and 6. '// &
                                 'Exiting ...')
            elseif ((palpha_eps <= 0d0) .or. (palpha_eps >= 1d0)) then
                call s_mpi_abort('palpha_eps must be in (0,1). Exiting ...')
            elseif ((ptgalpha_eps <= 0d0) .or. (ptgalpha_eps >= 1d0)) then
                call s_mpi_abort('ptgalpha_eps must be in (0,1). Exiting ...')
            end if
        elseif ((relax_model /= dflt_int) .or. (.not. f_is_default(palpha_eps)) &
                .or. (.not. f_is_default(ptgalpha_eps))) then
            call s_mpi_abort('relax is not set as true, but other phase '// &
                             'change parameters have been modified. Either '// &
                             'activate phase change or set the values '// &
                             'to default. Exiting ...')
        end if
    end subroutine s_check_inputs_phase_change

    !> Checks constraints on the Immersed Boundaries parameters.
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_ibm
        if (ib) then
            if (n <= 0) then
                call s_mpi_abort('ib is enabled but n = 0. '// &
                                 'Immersed Boundaries do not work in 1D. '// &
                                 'Exiting ...')
            else if (num_ibs <= 0 .or. num_ibs > num_patches_max) then
                call s_mpi_abort('num_ibs must be between 1 and '// &
                                 'num_patches_max. Exiting ...')
            end if
        end if

        if ((.not. ib) .and. num_ibs > 0) then
            call s_mpi_abort('num_ibs is set, but ib is not enabled. Exiting ...')
        end if
    end subroutine s_check_inputs_ibm

#endif

    !> Checks constraints on dimensionality and the number of cells for the grid.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_simulation_domain
        if (m == dflt_int) then
            call s_mpi_abort('m must be set. Exiting ...')
        elseif (n == dflt_int) then
            call s_mpi_abort('n must be set. Exiting ...')
        elseif (p == dflt_int) then
            call s_mpi_abort('p must be set. Exiting ...')
        elseif (m <= 0) then
            call s_mpi_abort('m must be positive. Exiting ...')
        elseif (n < 0) then
            call s_mpi_abort('n must be non-negative. Exiting ...')
        elseif (p < 0) then
            call s_mpi_abort('p must be non-negative. Exiting ...')
        elseif (cyl_coord .and. p > 0 .and. mod(p, 2) /= 1) then
            call s_mpi_abort('p must be odd for cylindrical coordinates '// &
                             '(cyl_coord = T and p != 0). Exiting ...')
        elseif (n == 0 .and. p > 0) then
            call s_mpi_abort('p must be 0 if n = 0. Exiting ...')
        end if
    end subroutine s_check_inputs_simulation_domain

    !> Checks constraints on model equations and number of fluids in the flow.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_model_eqns_and_num_fluids
        if (all(model_eqns /= (/1, 2, 3, 4/))) then
            call s_mpi_abort('model_eqns must be 1, 2, 3, or 4. Exiting ...')
        elseif (num_fluids /= dflt_int .and. num_fluids < 1) then
            call s_mpi_abort('num_fluids must be positive. Exiting ...')
        elseif (model_eqns == 1 .and. num_fluids /= dflt_int) then
            call s_mpi_abort('num_fluids is not supported for '// &
                             'model_eqns = 1. Exiting ...')
        elseif (model_eqns == 2 .and. num_fluids == dflt_int) then
            call s_mpi_abort('5-equation model (model_eqns = 2) '// &
                             'requires num_fluids to be set. Exiting ...')
        elseif (model_eqns == 3 .and. num_fluids == dflt_int) then
            call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                             'requires num_fluids to be set. Exiting ...')
        elseif (model_eqns == 1 .and. adv_alphan) then
            call s_mpi_abort('adv_alphan is not supported for '// &
                             'model_eqns = 1. Exiting ...')
        elseif (model_eqns == 1 .and. mpp_lim) then
            call s_mpi_abort('mpp_lim is not supported for '// &
                             'model_eqns = 1. Exiting ...')
        elseif (num_fluids == 1 .and. mpp_lim) then
            call s_mpi_abort('mpp_lim is not supported for '// &
                             'num_fluids = 1. Exiting ...')
        elseif (model_eqns == 3 .and. cyl_coord .and. p /= 0) then
            call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                             'does not support cylindrical coordinates '// &
                             '(cyl_coord = T and p != 0). Exiting ...')
        end if
    end subroutine s_check_inputs_model_eqns_and_num_fluids

    !> Checks constraints regarding WENO order.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_weno
        if (all(weno_order /= (/1, 3, 5/))) then
            call s_mpi_abort('weno_order must be 1, 3, or 5. Exiting ...')
        elseif (m + 1 < weno_order) then
            call s_mpi_abort('m must be at least weno_order - 1. Exiting ...')
        elseif (n > 0 .and. n + 1 < weno_order) then
            call s_mpi_abort('n must be at least weno_order - 1. Exiting ...')
        elseif (p > 0 .and. p + 1 < weno_order) then
            call s_mpi_abort('p must be at least weno_order - 1. Exiting ...')
        end if
    end subroutine s_check_inputs_weno

    !> Checks constraints on the boundary conditions in the x-direction.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_bc
        logical :: skip_check !< Flag to skip the check when iterating over
        !! x, y, and z directions, for special treatment of cylindrical coordinates

        #:for X, VAR in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            #:for BOUND in ['beg', 'end']
                if (${VAR}$ == 0 .and. bc_${X}$%${BOUND}$ /= dflt_int) then
                    call s_mpi_abort('bc_${X}$%${BOUND}$ is not '// &
                                     'supported for ${VAR}$ = 0. Exiting ...')
                elseif (${VAR}$ > 0 .and. bc_${X}$%${BOUND}$ == dflt_int) then
                    call s_mpi_abort('${VAR}$ != 0 but bc_${X}$%${BOUND}$ '// &
                                     'is not set. Exiting ...')
                elseif ((bc_${X}$%beg == -1 .and. bc_${X}$%end /= -1) &
                        .or. &
                        (bc_${X}$%end == -1 .and. bc_${X}$%beg /= -1)) then
                    call s_mpi_abort('bc_${X}$%beg and bc_${X}$%end '// &
                                     'must be both periodic (= -1) or both '// &
                                     'non-periodic. Exiting ...')
                end if

                ! For cylindrical coordinates, y and z directions use a different check
                #:if (X == 'y') or (X == 'z')
                    skip_check = cyl_coord
                #:else
                    skip_check = .false.
                #:endif

                if (.not. skip_check) then
                    if (bc_${X}$%${BOUND}$ /= dflt_int) then
                        if (bc_${X}$%${BOUND}$ > -1 .or. bc_${X}$%${BOUND}$ < -16) then
                            call s_mpi_abort('bc_${X}$%${BOUND}$ must be '// &
                                             'between -1 and -16. Exiting ...')
                        elseif (bc_${X}$%${BOUND}$ == -14) then
                            call s_mpi_abort('bc_${X}$%${BOUND}$ must not '// &
                                             'be -14 for non-cylindrical '// &
                                             'coordinates. Exiting ...')
                        end if
                    end if
                end if

            #:endfor
        #:endfor

        if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -13)) then
            call s_mpi_abort('Boundary condition -13 is not supported. Exiting ...')
        end if

        ! Check for y and z directions for cylindrical coordinates
        if (cyl_coord) then
            if (n == 0) then
                call s_mpi_abort('n must be positive for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif (p > 0 .and. bc_y%beg /= -14) then
                call s_mpi_abort('bc_y%beg must be -14 for 3D cylindrical '// &
                                 'coordinates (p > 0). Exiting ...')
            elseif (p == 0 .and. bc_y%beg /= -2) then
                call s_mpi_abort('bc_y%beg must be -2 for 2D cylindrical '// &
                                 'coordinates (p = 0). Exiting ...')
            elseif (bc_y%end > -1 .or. bc_y%end < -16) then
                call s_mpi_abort('bc_y%end must be between -1 and -16. '// &
                                 'Exiting ...')
            elseif (bc_y%end == -14) then
                call s_mpi_abort('bc_y%end must not be -14. Exiting ...')
            end if

            ! 3D cylindrical coordinates
            if (p /= 0) then
                if (bc_z%beg /= -1 .and. bc_z%beg /= -2) then
                    call s_mpi_abort('bc_z%beg must be -1 or -2 for 3D '// &
                                     'cylindrical coordinates. Exiting ...')
                elseif (bc_z%end /= -1 .and. bc_z%end /= -2) then
                    call s_mpi_abort('bc_z%end must be -1 or -2 for 3D '// &
                                     'cylindrical coordinates. Exiting ...')
                end if
            end if
        end if
    end subroutine s_check_inputs_bc

    !> Checks constraints on the stiffened equation of state fluids parameters.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_stiffened_eos
        character(len=5) :: iStr !< for int to string conversion
        integer :: bub_fac !< For allowing an extra fluid_pp if there are subgrid bubbles
        integer :: i

        bub_fac = 0
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (.not. f_is_default(fluid_pp(i)%gamma) &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'gamma must be positive. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    (.not. f_is_default(fluid_pp(i)%gamma))) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. &
                     (.not. f_is_default(fluid_pp(i)%gamma)))) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (.not. f_is_default(fluid_pp(i)%pi_inf) &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf must be non-negative. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    .not. f_is_default(fluid_pp(i)%pi_inf)) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. (.not. f_is_default(fluid_pp(i)%pi_inf)))) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif (fluid_pp(i)%cv < 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'cv must be positive. Exiting ...')
            end if
        end do
    end subroutine s_check_inputs_stiffened_eos

    !> Checks constraints on the surface tension parameters.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_surface_tension
        if (.not. f_is_default(sigma) .and. sigma < 0d0) then
            call s_mpi_abort('The surface tension coefficient must be'// &
                             'greater than or equal to zero. Exiting ...')
        elseif (.not. f_is_default(sigma) .and. model_eqns /= 3) then
            call s_mpi_abort("The surface tension model requires"// &
                             'model_eqns=3. Exiting ...')
        end if
    end subroutine s_check_inputs_surface_tension

    !> Checks constraints on the inputs for moving boundaries.
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_moving_bc
        #:for X, VB2, VB3 in [('x', 'vb2', 'vb3'), ('y', 'vb3', 'vb1'), ('z', 'vb1', 'vb2')]
            if (any((/bc_${X}$%vb1, bc_${X}$%vb2, bc_${X}$%vb3/) /= 0d0)) then
                if (bc_${X}$%beg == -15) then
                    if (any((/bc_${X}$%${VB2}$, bc_${X}$%${VB3}$/) /= 0d0)) then
                        call s_mpi_abort("bc_${X}$%beg must be -15 if "// &
                                         "bc_${X}$%${VB2}$ or bc_${X}$%${VB3}$ "// &
                                         "is set. Exiting ...")
                    end if
                elseif (bc_${X}$%beg /= -16) then
                    call s_mpi_abort("bc_${X}$%beg must be -15 or -16 if "// &
                                     "bc_${X}$%vb[1,2,3] is set. Exiting ...")
                end if
            end if
        #:endfor

        #:for X, VE2, VE3 in [('x', 've2', 've3'), ('y', 've3', 've1'), ('z', 've1', 've2')]
            if (any((/bc_${X}$%ve1, bc_${X}$%ve2, bc_${X}$%ve3/) /= 0d0)) then
                if (bc_${X}$%end == -15) then
                    if (any((/bc_${X}$%${VE2}$, bc_${X}$%${VE3}$/) /= 0d0)) then
                        call s_mpi_abort("bc_${X}$%end must be -15 if "// &
                                         "bc_${X}$%${VE2}$ or bc_${X}$%${VE3}$ "// &
                                         "is set. Exiting ...")
                    end if
                elseif (bc_${X}$%end /= -16) then
                    call s_mpi_abort("bc_${X}$%end must be -15 or -16 if "// &
                                     "bc_${X}$%ve[1,2,3] is set. Exiting ...")
                end if
            end if
        #:endfor
    end subroutine s_check_inputs_moving_bc

end module m_checker_common
