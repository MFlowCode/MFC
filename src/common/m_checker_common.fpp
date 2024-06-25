!>
!!@file m_checker_common.f90
!!@brief Contains module m_checker_common

!> @brief The purpose of the module is to check for compatible input files for
!!              inputs common to pre-processing, post-processing and simulation
module m_checker_common

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file
        !! Called by all three stages
    subroutine s_check_inputs_common

        ! Not run by pre_process
        call s_check_time_stepping
        call s_check_finite_difference

        ! Not run by post_process
        if (bubbles) call s_check_bubbles
        call s_check_qbmm_and_polydisperse
        if (adv_n) call s_check_adv_n
        call s_check_phase_change
        call s_check_inputs_ibm

        ! Run by all three stages
        call s_check_simulation_domain
        call s_check_model_eqn_and_num_fluids
        call s_check_weno
        call s_check_bc
        call s_check_stiffened_eos
        call s_check_inputs_moving_bc

    end subroutine s_check_inputs_common

    !> Checks constraints on the time-stepping parameters
        !! Called by s_check_inputs_common for simulation and post-processing
    subroutine s_check_time_stepping
#ifndef MFC_PRE_PROCESS
        if (t_step_start < 0) then
            call s_mpi_abort('t_step_start must be non-negative. Exiting ...')
        elseif (t_step_stop <= t_step_start) then
            call s_mpi_abort('t_step_stop must be greater than t_step_start. '// &
                             'Exiting ...')
        elseif (t_step_save > t_step_stop - t_step_start) then
            call s_mpi_abort('t_step_save must be less or equal to '// &
                             '(t_step_stop - t_step_start). Exiting ...')
        end if
#endif
    end subroutine s_check_time_stepping

    !> Checks constraints on the finite difference parameters
        !! Called by s_check_inputs_common for simulation and post-processing
    subroutine s_check_finite_difference
#ifndef MFC_PRE_PROCESS
        if (fd_order /= dflt_int .and. all(fd_order /= (/1, 2, 4/))) then
            call s_mpi_abort('fd_order must be 1, 2, or 4. Exiting ...')
        end if
#endif
    end subroutine s_check_finite_difference

    !> Checks constraints on the bubble parameters
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_bubbles
#ifndef MFC_POST_PROCESS
        if (nb < 1) then
            call s_mpi_abort('The Ensemble-Averaged Bubble Model '// &
                             'requires nb >= 1. Exiting ...')
        elseif (polydisperse .and. (nb == 1)) then
            call s_mpi_abort('Polydisperse bubble dynamics requires nb > 1 '// &
                             'Exiting ...')
        elseif (polydisperse .and. (mod(nb, 2) == 0)) then
            call s_mpi_abort('nb must be odd '// &
                             'Exiting ...')
        elseif ((.not. polytropic) .and. R0ref == dflt_real) then
            call s_mpi_abort('R0ref must be set if using bubbles with '// &
                             'polytropic = .false.. Exiting ...')
        elseif (nb == dflt_int) then
            call s_mpi_abort('nb must be set if using bubbles. Exiting ...')
        elseif (thermal > 3) then
            call s_mpi_abort('thermal must be less than 4 if using bubbles. '// &
                             'Exiting ...')
        elseif (model_eqns == 3) then
            call s_mpi_abort('Bubble models untested with '// &
                             '6-equation model. Exiting ...')
        elseif (model_eqns == 1) then
            call s_mpi_abort('Bubble models untested with '// &
                             'pi-gamma model. Exiting ...')
            !TODO: Comment this out when testing riemann with hll
        elseif (model_eqns == 4 .and. rhoref == dflt_real) then
            call s_mpi_abort('rhoref must be set if using bubbles with '// &
                             'model_eqns = 4. Exiting ...')
        elseif (model_eqns == 4 .and. pref == dflt_real) then
            call s_mpi_abort('pref must be set if using bubbles with '// &
                             'model_eqns = 4. Exiting ...')
        elseif (model_eqns == 4 .and. num_fluids > 1) then
            call s_mpi_abort('4-equation model (model_eqns = 4) is '// &
                             'single-component and requires num_fluids = 1. '// &
                             'Exiting ...')
        elseif (cyl_coord) then
            call s_mpi_abort('Bubble models untested in cylindrical coordinates')
        end if
#endif
    end subroutine s_check_bubbles

    !> Checks constraints on the QBMM and polydisperse bubble parameters
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_qbmm_and_polydisperse
#ifndef MFC_POST_PROCESS
        if ((.not. bubbles) .and. polydisperse) then
            call s_mpi_abort('Polydisperse bubble modeling requires the '// &
                             'bubbles flag to be set. Exiting ...')
        elseif (polydisperse .and. f_approx_equal(poly_sigma, dflt_real)) then
            call s_mpi_abort('Polydisperse bubble modeling requires '// &
                             'poly_sigma > 0. Exiting ...')
        elseif (qbmm .and. (.not. bubbles)) then
            call s_mpi_abort('QBMM is enabled but bubbles are not. Exiting ...')
        elseif (qbmm .and. (nnode /= 4)) then
            call s_mpi_abort('nnode not supported. Exiting ...')
        end if
#endif
    end subroutine s_check_qbmm_and_polydisperse

    !> Checks constraints on the adv_n flag
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_adv_n
#ifndef MFC_POST_PROCESS
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
#endif
    end subroutine

    !> Checks constraints on the phase change parameters
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_phase_change
#ifndef MFC_POST_PROCESS
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
        elseif ((relax_model /= dflt_int) .or. (.not. f_approx_equal(palpha_eps, dflt_real)) &
                .or. (.not. f_approx_equal(ptgalpha_eps, dflt_real))) then
            call s_mpi_abort('relax is not set as true, but other phase '// &
                             'change parameters have been modified. Either '// &
                             'activate phase change or set the values '// &
                             'to default. Exiting ...')
        end if
#endif
    end subroutine s_check_phase_change

    !> Checks constraints on the Immersed Boundaries parameters
        !! Called by s_check_inputs_common for pre-processing and simulation
    subroutine s_check_inputs_ibm
#ifndef MFC_POST_PROCESS
        if (num_ibs > 0 .and. .not. ib) then
            call s_mpi_abort('num_ibs is set, but ib is not enabled. Exiting ...')
        end if

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
#endif
    end subroutine s_check_inputs_ibm

    !> Checks constraints on dimensionality and the number of cells for the grid
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_simulation_domain
        if (m <= 0) then
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
    end subroutine s_check_simulation_domain

    !> Checks constraints on model equations and number of fluids in the flow
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_model_eqn_and_num_fluids
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
        end if
    end subroutine s_check_model_eqn_and_num_fluids

    !> Checks constraints regarding WENO order
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_weno
        if (all(weno_order /= (/1, 3, 5/))) then
            call s_mpi_abort('weno_order must be 1, 3, or 5. Exiting ...')
        elseif (m + 1 < weno_order) then
            call s_mpi_abort('m must be at least weno_order - 1. Exiting ...')
        elseif (n > 0 .and. n + 1 < weno_order) then
            call s_mpi_abort('n must be at least weno_order - 1. Exiting ...')
        elseif (p > 0 .and. p + 1 < weno_order) then
            call s_mpi_abort('p must be at least weno_order - 1. Exiting ...')
        end if
    end subroutine s_check_weno

    !> Checks constraints on the boundary conditions in the x-direction
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_bc
        if (bc_x%beg < -16 .or. bc_x%beg > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%beg. Exiting ...')
        elseif (bc_x%end < -16 .or. bc_x%end > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%end. Exiting ...')
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            call s_mpi_abort('bc_x%beg and bc_x%end must be both periodic '// &
                             '(= -1) or both non-periodic. Exiting ...')
        end if

        if (cyl_coord) then ! Cartesian coordinates

            ! Constraints on the boundary conditions in the r-direction
            if (bc_y%beg /= dflt_int &
                .and. &
                ((p > 0 .and. bc_y%beg /= -14) &
                 .or. &
                 (p == 0 .and. bc_y%beg /= -2))) then
                call s_mpi_abort('Unsupported value of bc_y%beg. Exiting ...')
            elseif (bc_y%end /= dflt_int &
                    .and. &
                    (bc_y%end < -16 .or. bc_y%end > -1 .or. bc_y%end == -14)) then
                call s_mpi_abort('Unsupported value of bc_y%end. Exiting ...')
            elseif ((n > 0 .and. bc_y%beg == dflt_int)) then
                call s_mpi_abort('n != 0 but bc_y%beg is not set. Exiting ...')
            elseif ((n > 0 .and. bc_y%end == dflt_int)) then
                call s_mpi_abort('n != 0 but bc_y%end is not set. Exiting ...')

                ! Constraints on the boundary conditions in the theta-direction
            elseif (bc_z%beg /= dflt_int &
                    .and. &
                    (bc_z%beg /= -1 .and. bc_z%beg /= -2)) then
                call s_mpi_abort('Unsupported value of bc_z%beg. Exiting ...')
            elseif (bc_z%end /= dflt_int &
                    .and. &
                    (bc_z%end /= -1 .and. bc_z%end /= -2)) then
                call s_mpi_abort('Unsupported value of bc_z%end. Exiting ...')
            elseif (p == 0 .and. bc_z%beg /= dflt_int) then
                call s_mpi_abort('bc_z%beg is not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. bc_z%beg == dflt_int) then
                call s_mpi_abort('p != 0 but bc_z%beg is not set. Exiting ...')
            elseif (p == 0 .and. bc_z%end /= dflt_int) then
                call s_mpi_abort('bc_z%end is not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. bc_z%end == dflt_int) then
                call s_mpi_abort('p != 0 but bc_z%end is not set. Exiting ...')
            elseif (p > 0 &
                    .and. &
                    ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                     .or. &
                     (bc_z%end == -1 .and. bc_z%beg /= -1))) then
                call s_mpi_abort('bc_z%beg and bc_z%end must be both periodic '// &
                                 '(= -1) or both non-periodic. Exiting ...')
            end if

        else ! not cylindrical

            ! Constraints on the boundary conditions in the y-direction
            if (bc_y%beg /= dflt_int &
                .and. &
                (bc_y%beg < -16 .or. bc_y%beg > -1 .or. bc_y%beg == -14)) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'bc_y%beg. Exiting ...')
            elseif (bc_y%end /= dflt_int &
                    .and. &
                    (bc_y%end < -16 .or. bc_y%end > -1 .or. bc_y%end == -14)) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'bc_y%end. Exiting ...')
            elseif (n == 0 .and. bc_y%beg /= dflt_int) then
                call s_mpi_abort('bc_y%beg is not supported for n = 0. Exiting ...')
            elseif (n > 0 .and. bc_y%beg == dflt_int) then
                call s_mpi_abort('n != 0 but bc_y%beg is not set. Exiting ...')
            elseif (n == 0 .and. bc_y%end /= dflt_int) then
                call s_mpi_abort('bc_y%end is not supported for n = 0. Exiting ...')
            elseif (n > 0 .and. bc_y%end == dflt_int) then
                call s_mpi_abort('n != 0 but bc_y%end is not set. Exiting ...')
            elseif ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                    .or. &
                    (bc_y%end == -1 .and. bc_y%beg /= -1)) then
                call s_mpi_abort('bc_y%beg and bc_y%end must be both periodic '// &
                                 '(= -1) or both non-periodic. Exiting ...')
            end if

            ! Constraints on the boundary conditions in the z-direction
            if (bc_z%beg /= dflt_int &
                .and. &
                (bc_z%beg < -16 .or. bc_z%beg > -1 .or. bc_z%beg == -14)) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'bc_z%beg. Exiting ...')
            elseif (bc_z%end /= dflt_int &
                    .and. &
                    (bc_z%end < -16 .or. bc_z%end > -1 .or. bc_z%end == -14)) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'bc_z%end. Exiting ...')
            elseif (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -13)) then
                call s_mpi_abort('Unsupported choice of boundary condition -13')
            elseif (p == 0 .and. bc_z%beg /= dflt_int) then
                call s_mpi_abort('bc_z%beg is not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. bc_z%beg == dflt_int) then
                call s_mpi_abort('p != 0 but bc_z%beg is not set. Exiting ...')
            elseif (p == 0 .and. bc_z%end /= dflt_int) then
                call s_mpi_abort('bc_z%end is not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. bc_z%end == dflt_int) then
                call s_mpi_abort('p != 0 but bc_z%end is not set. Exiting ...')
            elseif ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                    .or. &
                    (bc_z%end == -1 .and. bc_z%beg /= -1)) then
                call s_mpi_abort('bc_z%beg and bc_z%end must be both periodic '// &
                                 '(= -1) or both non-periodic. Exiting ...')
            end if
        end if
    end subroutine s_check_bc

    !> Checks constraints on the stiffened equation of state fluids parameters
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_stiffened_eos
        character(len=5) :: iStr !< for int to string conversion
        integer :: bub_fac !< For allowing an extra fluid_pp if there are subgrid bubbles
        integer :: i

        bub_fac = 0
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (.not. f_approx_equal(fluid_pp(i)%gamma, dflt_real) &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'gamma must be positive. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    (.not. f_approx_equal(fluid_pp(i)%gamma, dflt_real))) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. &
                     (.not. f_approx_equal(fluid_pp(i)%gamma, dflt_real)))) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (.not. f_approx_equal(fluid_pp(i)%pi_inf, dflt_real) &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf must be non-negative. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    .not. f_approx_equal(fluid_pp(i)%pi_inf, dflt_real)) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. (.not. f_approx_equal(fluid_pp(i)%pi_inf, dflt_real)))) &
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
    end subroutine s_check_stiffened_eos

    !> Checks constraints on the inputs for moving boundaries
        !! Called by s_check_inputs_common for all three stages
    subroutine s_check_inputs_moving_bc
        ! Moving Boundaries Checks: x boundaries
        if (any((/bc_x%vb1, bc_x%vb2, bc_x%vb3/) /= 0d0)) then
            if (bc_x%beg == -15) then
                if (any((/bc_x%vb2, bc_x%vb3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_x%beg and"// &
                                     "bc_x%vb2 or bc_x%vb3. Exiting ...")
                end if
            elseif (bc_x%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_x%beg and"// &
                                 "bc_x%vb1, bc_x%vb2, or bc_x%vb3. Exiting...")
            end if
        end if

        if (any((/bc_x%ve1, bc_x%ve2, bc_x%ve3/) /= 0d0)) then
            if (bc_x%end == -15) then
                if (any((/bc_x%ve2, bc_x%ve3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_x%end and"// &
                                     "bc_x%ve2 or bc_x%ve3. Exiting ...")
                end if
            elseif (bc_x%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_x%end and"// &
                                 "bc_x%ve1, bc_x%ve2, or bc_x%ve3. Exiting...")
            end if
        end if

        ! Moving Boundaries Checks: y boundaries
        if (any((/bc_y%vb1, bc_y%vb2, bc_y%vb3/) /= 0d0)) then
            if (bc_y%beg == -15) then
                if (any((/bc_y%vb1, bc_y%vb3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_y%beg and"// &
                                     "bc_y%vb1 or bc_y%vb3. Exiting ...")
                end if
            elseif (bc_y%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_y%beg and"// &
                                 "bc_y%vb1, bc_y%vb2, or bc_y%vb3. Exiting...")
            end if
        end if

        if (any((/bc_y%ve1, bc_y%ve2, bc_y%ve3/) /= 0d0)) then
            if (bc_y%end == 15) then
                if (any((/bc_y%ve1, bc_y%ve3/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_y%end and"// &
                                     "bc_y%ve1 or bc_y%ve3. Exiting ...")
                end if
            elseif (bc_y%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_y%end and"// &
                                 "bc_y%ve1, bc_y%ve2, or bc_y%ve3. Exiting...")
            end if
        end if

        ! Moving Boundaries Checks: z boundaries
        if (any((/bc_z%vb1, bc_z%vb2, bc_z%vb3/) /= 0d0)) then
            if (bc_z%beg == -15) then
                if (any((/bc_x%vb1, bc_x%vb2/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_z%beg and"// &
                                     "bc_x%vb1 or bc_x%vb1. Exiting ...")
                end if
            elseif (bc_z%beg /= -16) then
                call s_mpi_abort("Unsupported combination of bc_z%beg and"// &
                                 "bc_z%vb1, bc_z%vb2, or bc_z%vb3. Exiting...")
            end if
        end if

        if (any((/bc_z%ve1, bc_z%ve2, bc_z%ve3/) /= 0d0)) then
            if (bc_z%end == -15) then
                if (any((/bc_x%ve1, bc_x%ve2/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                     "bc_z%ve2 or bc_z%ve3. Exiting ...")
                end if
            elseif (bc_z%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                 "bc_z%ve1, bc_z%ve2, or bc_z%ve3. Exiting...")
            end if
        end if

        ! Constraints on the surface tension model
        if (.not. f_approx_equal(sigma, dflt_real) .and. sigma < 0d0) then
            call s_mpi_abort('The surface tension coefficient must be'// &
                             'greater than or equal to zero. Exiting ...')
        elseif (.not. f_approx_equal(sigma, dflt_real) .and. model_eqns /= 3) then
            call s_mpi_abort("The surface tension model requires"// &
                             'model_eqns=3. Exiting ...')
        end if

    end subroutine s_check_inputs_moving_bc

end module m_checker_common
