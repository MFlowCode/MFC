!>
!!@file m_checker.f90
!!@brief Contains module m_checker

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    subroutine s_check_inputs

        integer :: bub_fac !< For allowing an extra fluid_pp if there are subgrid bubbles
        character(len=5) :: iStr, numStr !< for int to string conversion
        integer :: i
        logical :: dir_check !< Logical variable used to test the existence of folders

#ifndef MFC_MPI
        if parallel_io then
        print '(A)', 'MFC built with --no-mpi requires parallel_io=F. '// &
            'Exiting ...'
        call s_mpi_abort()
        end if
#endif

        bub_fac = 0
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1
        ! Startup checks for bubbles and bubble variables
        if (bubbles) then
            if (model_eqns /= 4 .and. model_eqns /= 2) then
                call s_mpi_abort('Bubbles require model_eqns = 4 or 2. '// &
                                 'Exiting ...')
            elseif (nb < 1) then
                call s_mpi_abort('The Ensemble-Averaged Bubble Model requires nb >= 1'// &
                                 'Exiting ...')
            elseif (polydisperse .and. (nb == 1)) then
                call s_mpi_abort('Polydisperse bubble dynamics requires nb > 1 '// &
                                 'Exiting ...')
            elseif (polydisperse .and. (mod(nb, 2) == 0)) then
                call s_mpi_abort('nb must be odd '// &
                                 'Exiting ...')
            elseif (model_eqns == 4 .and. (rhoref == dflt_real)) then
                call s_mpi_abort('rhoref must be set if using bubbles with '// &
                                 'model_eqns = 4. Exiting ...')
            elseif (model_eqns == 4 .and. (pref == dflt_real)) then
                call s_mpi_abort('pref must be set if using bubbles with '// &
                                 'model_eqns = 4. Exiting ...')
            elseif (model_eqns == 4 .and. (num_fluids > 1)) then
                call s_mpi_abort('num_fluids must be 1 if using bubbles with '// &
                                 'model_eqns = 4. Exiting ...')
            elseif ((.not. polytropic) .and. R0ref == dflt_real) then
                call s_mpi_abort('R0ref must be set if using bubbles with '// &
                                 'polytropic = .false.. Exiting ...')
            elseif (nb == dflt_int) then
                call s_mpi_abort('nb must be set if using bubbles. Exiting ...')
            elseif (thermal > 3) then
                call s_mpi_abort('thermal must be less than 4 if using bubbles. '// &
                                 'Exiting ...')
            end if

        end if

        if (adv_n) then
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
        end if

        if (qbmm .and. dist_type == dflt_int) then
            call s_mpi_abort('Dist type must be set if using QBMM. Exiting ...')
        else if (qbmm .and. (dist_type /= 1) .and. rhoRV > 0d0) then
            call s_mpi_abort('rhoRV cannot be used with dist_type \ne 1. Exiting ...')
        else if (polydisperse .and. R0_type == dflt_int) then
            call s_mpi_abort('R0 type must be set if using Polydisperse. Exiting ...')
        end if

        if (hypoelasticity .and. (model_eqns /= 2)) then
            call s_mpi_abort('hypoelasticity requires model_eqns = 2'// &
                             'exiting ...')
        end if
        ! phase change checkers.
        if (relax) then
            if (model_eqns /= 3) then
                call s_mpi_abort('phase change requires model_eqns = 3. '// &
                                 'Exiting ...')
            elseif ((relax_model < 0) .or. (relax_model > 6)) then
                call s_mpi_abort('relax_model should be in between 0 and 6. '// &
                                 'Exiting ...')
            elseif ((palpha_eps <= 0d0) .or. (palpha_eps >= 1d0) .or. &
                    (ptgalpha_eps <= 0d0) .or. (ptgalpha_eps >= 1d0)) then
                call s_mpi_abort('both palpha_eps and ptgalpha_eps must &
 &               be in (0,1). '//'Exiting ...')
            end if

        elseif ((relax_model /= dflt_int) .or. (palpha_eps /= dflt_real) &
                .or. (ptgalpha_eps /= dflt_real)) then
            call s_mpi_abort('relax is not set as true, but other phase change parameters have &
&            been modified. Either activate phase change or set the values to default. '//'Exiting ...')
        end if
        if ((.not. old_grid) .and. old_ic) then
            call s_mpi_abort('old_ic cannot be enabled with old_grid disabled. '// &
                             'Exiting ...')

        elseif ((old_grid .or. old_ic) .and. t_step_old == dflt_int) then
            call s_mpi_abort('old_grid or old_ic enabled, but t_step_old not set. '// &
                             'Exiting ...')

            ! Constraints on dimensionality and the number of cells for the grid
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
        elseif (nGlobal < 2**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)
            call s_mpi_abort('Total number of cells must be at least '// &
                             '(2^[number of dimensions])*num_procs, which is currently '// &
                             trim(numStr)//'. Exiting ...')
            ! Constraints on domain boundaries locations in the x-direction
        elseif ((old_grid .and. x_domain%beg /= dflt_real) &
                .or. &
                ((.not. old_grid) .and. &
                 x_domain%beg == dflt_real)) then
            call s_mpi_abort('old_grid is not enabled but x_domain%beg is '// &
                             'not set. Exiting ...')
        elseif ((old_grid .and. x_domain%end /= dflt_real) &
                .or. &
                ((.not. old_grid) .and. &
                 x_domain%end == dflt_real)) then
            call s_mpi_abort('old_grid is not enabled but x_domain%end is '// &
                             'not set. Exiting ...')
        elseif ((.not. old_grid) &
                .and. &
                x_domain%beg >= x_domain%end) then
            call s_mpi_abort('x_domain%beg must be less than x_domain%end. '// &
                             'Exiting ...')
        end if

        if (cyl_coord) then ! Cartesian coordinates

            ! in case restart of a simulation
            if (old_grid .and. old_ic) then
                ! checking of there is any input to the domains
                if ((x_domain%beg /= dflt_real .or. x_domain%end /= dflt_real) &
                    .or. &
                    (y_domain%beg /= dflt_real .or. y_domain%end /= dflt_real) &
                    .or. &
                    (y_domain%beg /= dflt_real .or. y_domain%end /= dflt_real)) then
                    call s_mpi_abort('x_domain, y_domain, and/or z_domain '// &
                                     'are not supported in restart mode. '// &
                                     '(old_grid = T and old_ic = T). Exiting ...')
                elseif (m == dflt_int .or. n == dflt_int .or. p == dflt_int) then
                    call s_mpi_abort('m, n, and p must be set in restart mode. '// &
                                     '(old_grid = T and old_ic = T). Exiting ...')
                end if
                ! in case it is NOT restart
                ! Constraints on domain boundaries for cylindrical coordinates
            elseif (n == 0) then
                call s_mpi_abort('n must be positive for cylindrical coordinates. '// &
                                 'Exiting ...')
            elseif (y_domain%beg /= 0d0 &
                    .or. &
                    y_domain%end == dflt_real &
                    .or. &
                    y_domain%end < 0d0 &
                    .or. &
                    y_domain%beg >= y_domain%end) then
                call s_mpi_abort('y_domain%beg must be 0 and '// &
                                 'y_domain%end must be positive and '// &
                                 'greater than y_domain%beg for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) &
                    .or. &
                    (p == 0 .and. z_domain%end /= dflt_real)) then
                call s_mpi_abort('z_domain%beg and z_domain%end '// &
                                 'are not supported for p = 0. Exiting ...')
            elseif (p > 0 .and. (z_domain%beg /= 0d0 &
                                 .or. &
                                 z_domain%end /= 2d0*pi)) then
                call s_mpi_abort('z_domain%beg must be 0 and '// &
                                 'z_domain%end must be 2*pi for 3D cylindrical '// &
                                 'coordinates. Exiting ...')
            end if

        else
            ! Constraints on domain boundaries locations in the y-direction
            if ((n == 0 .and. y_domain%beg /= dflt_real) .or. &
                (n > 0 .and. ((old_grid .and. y_domain%beg /= dflt_real) .or. &
                              (.not. old_grid .and. y_domain%beg == dflt_real)))) then
                call s_mpi_abort('y_domain%beg must not be set '// &
                                 'when n = 0 or when n > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif ((n == 0 .and. y_domain%end /= dflt_real) .or. &
                    (n > 0 .and. ((old_grid .and. y_domain%end /= dflt_real) .or. &
                                  (.not. old_grid .and. y_domain%end == dflt_real)))) then
                call s_mpi_abort('y_domain%end must not be set '// &
                                 'when n = 0 or when n > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif (n > 0 .and. .not. old_grid .and. y_domain%beg >= y_domain%end) then
                call s_mpi_abort('y_domain%beg must be less than y_domain%end '// &
                                 'when both are set. '// &
                                 'Exiting ...')

                ! Constraints on domain boundaries locations in the z-direction
            elseif ((p == 0 .and. z_domain%beg /= dflt_real) .or. &
                    (p > 0 .and. ((old_grid .and. z_domain%beg /= dflt_real) .or. &
                                  (.not. old_grid .and. z_domain%beg == dflt_real)))) then
                call s_mpi_abort('z_domain%beg must not be set '// &
                                 'when p = 0 or when p > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif ((p == 0 .and. z_domain%end /= dflt_real) .or. &
                    (p > 0 .and. ((old_grid .and. z_domain%end /= dflt_real) .or. &
                                  (.not. old_grid .and. z_domain%end == dflt_real)))) then
                call s_mpi_abort('z_domain%end must not be set '// &
                                 'when p = 0 or when p > 0 and old_grid = F, and '// &
                                 'must be set otherwise. Exiting ...')
            elseif (p > 0 .and. .not. old_grid .and. z_domain%beg >= z_domain%end) then
                call s_mpi_abort('z_domain%beg must be less than z_domain%end '// &
                                 'when both are set. '// &
                                 'Exiting ...')
            end if
        end if

        if (loops_z < 1) then
            call s_mpi_abort('loops_z must be positive. Exiting ...')
        elseif (loops_y < 1) then
            call s_mpi_abort('loops_y must be positive. Exiting ...')
        end if

        ! Constraints on the grid stretching in the x-direction
        if (stretch_x) then
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_x are incompatible. '// &
                                 'Exiting ...')
            elseif (a_x == dflt_real) then
                call s_mpi_abort('a_x must be set if stretch_x = T. Exiting ...')
            elseif (x_a == dflt_real) then
                call s_mpi_abort('x_a must be set if stretch_x = T. Exiting ...')
            elseif (x_b == dflt_real) then
                call s_mpi_abort('x_b must be set if stretch_x = T. Exiting ...')
            elseif (x_a >= x_b) then
                call s_mpi_abort('x_a must be less than x_b if stretch_x = T. '// &
                                 'Exiting ...')
            elseif ((a_x + log(cosh(a_x*(x_domain%beg - x_a))) &
                     + log(cosh(a_x*(x_domain%beg - x_b))) &
                     - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
                call s_mpi_abort('x_domain%beg is too close to x_a and x_b '// &
                                 'for the given a_x. Exiting ...')
            elseif ((a_x + log(cosh(a_x*(x_domain%end - x_a))) &
                     + log(cosh(a_x*(x_domain%end - x_b))) &
                     - 2d0*log(cosh(0.5d0*a_x*(x_b - x_a))))/a_x <= 0d0) then
                call s_mpi_abort('x_domain%end is too close to x_a and x_b '// &
                                 'for the given a_x. Exiting ...')
            end if
        end if

        if (stretch_y) then
            ! Constraints on the grid stretching in the y-direction
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_y are incompatible. '// &
                                 'Exiting ...')
            elseif (n == 0) then
                call s_mpi_abort('n must be positive if stretch_y = T. Exiting ...')
            elseif (a_y == dflt_real) then
                call s_mpi_abort('a_y must be set if stretch_y = T. Exiting ...')
            elseif (y_a == dflt_real) then
                call s_mpi_abort('y_a must be set if stretch_y = T. Exiting ...')
            elseif (y_b == dflt_real) then
                call s_mpi_abort('y_b must be set if stretch_y = T. Exiting ...')
            elseif (y_a >= y_b) then
                call s_mpi_abort('y_a must be less than y_b if stretch_y = T. '// &
                                 'Exiting ...')
            elseif ((a_y + log(cosh(a_y*(y_domain%beg - y_a))) &
                     + log(cosh(a_y*(y_domain%beg - y_b))) &
                     - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
                call s_mpi_abort('y_domain%beg is too close to y_a and y_b '// &
                                 'for the given a_y. Exiting ...')
            elseif ((a_y + log(cosh(a_y*(y_domain%end - y_a))) &
                     + log(cosh(a_y*(y_domain%end - y_b))) &
                     - 2d0*log(cosh(0.5d0*a_y*(y_b - y_a))))/a_y <= 0d0) then
                call s_mpi_abort('y_domain%end is too close to y_a and y_b '// &
                                 'for the given a_y. Exiting ...')
            end if
        end if

        ! Constraints on the grid stretching in the z-direction
        if (stretch_z) then
            if (old_grid) then
                call s_mpi_abort('old_grid and stretch_z are incompatible. '// &
                                 'Exiting ...')
            elseif (cyl_coord) then
                call s_mpi_abort('stretch_z is not supported for cylindrical '// &
                                 'coordinates. Exiting ...')
            elseif (p == 0) then
                call s_mpi_abort('p must be positive if stretch_z = T. Exiting ...')
            elseif (a_z == dflt_real) then
                call s_mpi_abort('a_z must be set if stretch_z = T. Exiting ...')
            elseif (z_a == dflt_real) then
                call s_mpi_abort('z_a must be set if stretch_z = T. Exiting ...')
            elseif (z_b == dflt_real) then
                call s_mpi_abort('z_b must be set if stretch_z = T. Exiting ...')
            elseif (z_a >= z_b) then
                call s_mpi_abort('z_a must be less than z_b if stretch_z = T. '// &
                                 'Exiting ...')
            elseif ((a_z + log(cosh(a_z*(z_domain%beg - z_a))) &
                     + log(cosh(a_z*(z_domain%beg - z_b))) &
                     - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
                call s_mpi_abort('z_domain%beg is too close to z_a and z_b '// &
                                 'for the given a_z. Exiting ...')
            elseif ((a_z + log(cosh(a_z*(z_domain%end - z_a))) &
                     + log(cosh(a_z*(z_domain%end - z_b))) &
                     - 2d0*log(cosh(0.5d0*a_z*(z_b - z_a))))/a_z <= 0d0) then
                call s_mpi_abort('z_domain%end is too close to z_a and z_b '// &
                                 'for the given a_z. Exiting ...')
            end if
        end if

        ! Constraints on model equations and number of fluids in the flow
        if (all(model_eqns /= (/1, 2, 3, 4/))) then
            call s_mpi_abort('model_eqns must be 1, 2, 3, or 4. Exiting ...')
        elseif (num_fluids /= dflt_int .and. num_fluids < 1) then
            call s_mpi_abort('num_fluids must be positive. Exiting ...')
        elseif (model_eqns == 1 .and. adv_alphan) then
            call s_mpi_abort('adv_alphan is not supported for model_eqns = 1. '// &
                             'Exiting ...')
            ! Constraints on the order of the WENO scheme
        elseif (all(weno_order /= (/1, 3, 5/))) then
            call s_mpi_abort('weno_order must be 1, 3, or 5. Exiting ...')
        elseif (m + 1 < weno_order) then
            call s_mpi_abort('m must be at least weno_order - 1. Exiting ...')
        elseif (n > 0 .and. n + 1 < weno_order) then
            call s_mpi_abort('n must be at least weno_order - 1. Exiting ...')
        elseif (p > 0 .and. p + 1 < weno_order) then
            call s_mpi_abort('p must be at least weno_order - 1. Exiting ...')
        elseif (nGlobal < weno_order**(min(1, m) + min(1, n) + min(1, p))*num_procs) then
            call s_int_to_str(weno_order**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)
            call s_mpi_abort('Total number of cells must be at least '// &
                             '(weno_order^[number of dimensions])*num_procs, '// &
                             'which is currently '//trim(numStr)//'. Exiting ...')

            ! Constraints on the boundary conditions in the x-direction
        elseif (bc_x%beg < -16 .or. bc_x%beg > -1 .or. bc_x%beg == -14) then
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

        else

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
            elseif (n > 0 &
                    .and. &
                    ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                     .or. &
                     (bc_y%end == -1 .and. bc_y%beg /= -1))) then
                call s_mpi_abort('bc_y%beg and bc_y%end must be both periodic '// &
                                 '(= -1) or both non-periodic. Exiting ...')

                ! Constraints on the boundary conditions in the z-direction
            elseif (bc_z%beg /= dflt_int &
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
            elseif (p > 0 &
                    .and. &
                    ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                     .or. &
                     (bc_z%end == -1 .and. bc_z%beg /= -1))) then
                call s_mpi_abort('bc_z%beg and bc_z%end must be both periodic '// &
                                 '(= -1) or both non-periodic. Exiting ...')
            end if

        end if

        ! Constraints on number of patches making up the initial condition
        if (num_patches < 0 .or. &
            (num_patches == 0 .and. t_step_old == dflt_int)) then
            call s_mpi_abort('num_patches must be non-negative for the '// &
                             'non-restart case. Exiting ...')
            ! Constraints on perturbing the initial condition
        elseif (perturb_flow &
                .and. &
                (perturb_flow_fluid == dflt_int .or. perturb_flow_mag == dflt_real)) then
            call s_mpi_abort('perturb_flow_fluid and perturb_flow_mag '// &
                             'must be set with perturb_flow = T. Exiting ...')
        elseif ((.not. perturb_flow) &
                .and. &
                (perturb_flow_fluid /= dflt_int .or. perturb_flow_mag /= dflt_real)) then
            call s_mpi_abort('perturb_flow_fluid and perturb_flow_mag '// &
                             'must not be set with perturb_flow = F. Exiting ...')
        elseif ((perturb_flow_fluid > num_fluids) &
                .or. &
                (perturb_flow_fluid < 0 .and. perturb_flow_fluid /= dflt_int)) then
            call s_mpi_abort('perturb_flow_fluid must be between 0 and '// &
                             'num_fluids. Exiting ...')
        elseif (perturb_sph .and. perturb_sph_fluid == dflt_int) then
            call s_mpi_abort('perturb_sph_fluid must be set with perturb_sph = T. '// &
                             'Exiting ...')
        elseif (.not. perturb_sph .and. perturb_sph_fluid /= dflt_int) then
            call s_mpi_abort('perturb_sph_fluid must not be set with perturb_sph = F. '// &
                             'Exiting ...')
        elseif ((perturb_sph_fluid > num_fluids) &
                .or. &
                (perturb_sph_fluid < 0 .and. perturb_sph_fluid /= dflt_int)) then
            call s_mpi_abort('perturb_sph_fluid must be between 0 and '// &
                             'num_fluids. Exiting ...')
        elseif ((any(fluid_rho /= dflt_real)) .and. (.not. perturb_sph)) then
            call s_mpi_abort('fluid_rho must not be set with perturb_sph = F. '// &
                             'Exiting ...')
        end if

        if (perturb_sph) then
            do i = 1, num_fluids
                call s_int_to_str(i, iStr)
                if (fluid_rho(i) == dflt_real) then
                    call s_mpi_abort('fluid_rho('//trim(iStr)//') must be set '// &
                                     'if perturb_sph = T. Exiting ...')
                end if
            end do
        end if

        ! Constraints on the hypertangent velocity profile
        if (vel_profile .and. (n == 0)) then
            call s_mpi_abort('vel_profile requires n > 0. Exiting ...')
        end if

        ! Constraints on the instability wave
        if (instability_wave .and. (n == 0)) then
            call s_mpi_abort('instability_wave requires n > 0. Exiting ...')
        end if

        ! Constraints on Immersed Boundary Method
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

        if (num_ibs > 0 .and. .not. ib) then
            call s_mpi_abort('num_ibs is set, but ib is not enabled. Exiting ...')
        end if

        ! Constraints on the stiffened equation of state fluids parameters
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'gamma must be positive. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%gamma /= dflt_real) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%gamma <= 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%gamma /= dflt_real)) &
                then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of num_fluids '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (fluid_pp(i)%pi_inf /= dflt_real &
                    .and. &
                    fluid_pp(i)%pi_inf < 0d0) then
                call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf must be non-negative. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%pi_inf /= dflt_real) then
                call s_mpi_abort('model_eqns = 1 does not support '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif ((i <= num_fluids + bub_fac .and. fluid_pp(i)%pi_inf < 0d0) &
                    .or. &
                    (i > num_fluids + bub_fac .and. fluid_pp(i)%pi_inf /= dflt_real)) &
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

        ! Constraints on the surface tension model
        if (sigma /= dflt_real .and. sigma < 0d0) then
            call s_mpi_abort('The surface tension coefficient must be'// &
                             'greater than or equal to zero. Exiting ...')
        elseif (sigma /= dflt_real .and. model_eqns /= 3) then
            call s_mpi_abort("The surface tension model requires"// &
                             'model_eqns=3. Exiting ...')
        end if

        ! Moving Boundaries Checks: x boundaries
        if (any((/bc_x%vb1, bc_x%vb2, bc_x%vb3/) /= 0d0)) then
            if (bc_x%beg == 15) then
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
            if (bc_x%end == 15) then
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
            if (bc_y%beg == 15) then
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
            if (bc_z%beg == 15) then
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
            if (bc_z%end == 15) then
                if (any((/bc_x%ve1, bc_x%ve2/) /= 0d0)) then
                    call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                     "bc_z%ve2 or bc_z%ve3. Exiting ...")
                end if
            elseif (bc_z%end /= -16) then
                call s_mpi_abort("Unsupported combination of bc_z%end and"// &
                                 "bc_z%ve1, bc_z%ve2, or bc_z%ve3. Exiting...")
            end if
        end if

        end subroutine s_check_inputs

    end module m_checker
