!>
!! @file m_start_up.f90
!! @brief Contains module m_checker

#:include 'case.fpp'

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    subroutine s_check_inputs()

        character(len=5) :: iStr, jStr, numStr
        integer :: bub_fac !for allowing an extra fluid_pp if there are bubbles
        integer :: i, j

        bub_fac = 0
        if (bubbles .and. (num_fluids == 1)) bub_fac = 1

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        if (rdma_mpi) then
            call s_mpi_abort('Unsupported value of rdma_mpi. Exiting ...')
        end if
#endif

#ifndef MFC_cuTENSOR
        if (cu_tensor) then
            call s_mpi_abort('Unsupported value of cu_tensor. MFC was not built '// &
                             'with the NVIDIA cuTENSOR library. Exiting ...')
        end if
#endif

        ! Computational Domain Parameters ==================================
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
        elseif (dt <= 0) then
            call s_mpi_abort('dt must be positive. Exiting ...')
        elseif (t_step_start < 0) then
            call s_mpi_abort('t_step_start must be non-negative. Exiting ...')
        elseif (t_step_stop <= t_step_start) then
            call s_mpi_abort('t_step_stop must be greater than t_step_start. '// &
                             'Exiting ...')
        elseif (t_step_save > t_step_stop - t_step_start) then
            call s_mpi_abort('t_step_save must be less or equal to '// &
                             '(t_step_stop - t_step_start). Exiting ...')
        end if
        ! ==================================================================

        ! Simulation Algorithm Parameters ==================================
        if (all(model_eqns /= (/1, 2, 3, 4/))) then
            call s_mpi_abort('model_eqns must be 1, 2, 3, or 4. Exiting ...')
        end if

        if (bubbles) then
            if (model_eqns == 2 .and. bubble_model == 1) then
                call s_mpi_abort('The 5-equation bubbly flow model requires bubble_model = 2 (Keller--Miksis)')
            elseif (nb < 1) then
                call s_mpi_abort('The Ensemble-Averaged Bubble Model requires nb >= 1')
            elseif (bubble_model == 3 .and. (polytropic .neqv. .true.) .and. (.not. qbmm)) then
                call s_mpi_abort('RP bubbles require polytropic compression')
            elseif (cyl_coord) then
                call s_mpi_abort('Bubble models untested in cylindrical coordinates')
            elseif (model_eqns == 3) then
                call s_mpi_abort('Bubble models untested with 6-equation model')
            elseif (model_eqns == 1) then
                call s_mpi_abort('Bubble models untested with pi-gamma model')
                !TODO: Comment this out when testing riemann with hll
            elseif (riemann_solver /= 2) then
                call s_mpi_abort('Bubble modeling requires riemann_solver = 2')
            elseif (avg_state /= 2) then
                call s_mpi_abort('Bubble modeling requires arithmetic average '// &
                                 '(avg_state = 2). Exiting ...')
            end if
        end if

        if (model_eqns == 4 .and. num_fluids /= 1) then
            call s_mpi_abort('4-equation model (model_eqns = 4) is '// &
                             'single-component and requires num_fluids = 1. '// &
                             'Exiting ...')
        end if

        if ((bubbles .neqv. .true.) .and. polydisperse) then
            call s_mpi_abort('Polydisperse bubble modeling requires the bubble switch to be activated')
        elseif (polydisperse .and. f_approx_equal(poly_sigma, dflt_real)) then
            call s_mpi_abort('Polydisperse bubble modeling requires poly_sigma > 0')
        elseif (qbmm .and. (bubbles .neqv. .true.)) then
            call s_mpi_abort('QBMM requires bubbles')
        elseif (qbmm .and. (nnode /= 4)) then
            call s_mpi_abort('nnode not supported')
        end if

        if (model_eqns == 3) then
            if (riemann_solver /= 2) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'requires riemann_solver = 2. Exiting ...')
            elseif (alt_soundspeed) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'does not support alt_sound_speed. Exiting ...')
            elseif (avg_state /= 2) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'requires avg_state = 2. Exiting ...')
            elseif (wave_speeds /= 1) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'requires wave_speeds = 1. Exiting ...')
            elseif (cyl_coord .and. p /= 0) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'does not support cylindrical coordinates '// &
                                 '(cyl_coord = T and p != 0). Exiting ...')
            end if
        end if

        ! phase change checkers.
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

        if (num_fluids /= dflt_int .and. num_fluids < 1) then
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
        elseif (time_stepper < 1 .or. time_stepper > 5) then
            if (time_stepper /= 23) then
                call s_mpi_abort('time_stepper must be between 1 and 5. '// &
                                 'Exiting ...')
            end if
        elseif (all(weno_order /= (/1, 3, 5/))) then
            call s_mpi_abort('weno_order must be 1, 3, or 5. Exiting ...')
        elseif (m + 1 < num_stcls_min*weno_order) then
            call s_int_to_str(num_stcls_min*weno_order, numStr)
            call s_mpi_abort('m must be greater than or equal to '// &
                             '(num_stcls_min*weno_order - 1), whose value is '// &
                             trim(numStr)//'. Exiting ...')
        elseif (n + 1 < min(1, n)*num_stcls_min*weno_order) then
            call s_mpi_abort('For 2D simulation, n must be greater than or '// &
                             'equal to (num_stcls_min*weno_order - 1), '// &
                             'whose value is '//trim(numStr)//'. Exiting ...')
        elseif (p + 1 < min(1, p)*num_stcls_min*weno_order) then
            call s_mpi_abort('For 3D simulation, p must be greater than or '// &
                             'equal to (num_stcls_min*weno_order - 1), '// &
                             'whose value is '//trim(numStr)//'. Exiting ...')
        elseif (weno_order /= 1 .and. f_approx_equal(weno_eps, dflt_real)) then
            call s_mpi_abort('weno_order != 1, but weno_eps is not set. '// &
                             'A typical value of weno_eps is 1e-6. '// &
                             'Exiting ...')
        elseif (weno_eps <= 0d0) then
            call s_mpi_abort('weno_eps must be positive. '// &
                             'A typical value of weno_eps is 1e-6. '// &
                             'Exiting ...')
        elseif (teno .and. f_approx_equal(teno_CT, dflt_real)) then
            call s_mpi_abort('teno is used, but teno_CT is not set. '// &
                             'A typical value of teno_CT is 1e-6. '// &
                             'Exiting ...')
        elseif (teno .and. teno_CT <= 0d0) then
            call s_mpi_abort('teno_CT must be positive. '// &
                             'A typical value of teno_CT is 1e-6. '// &
                             'Exiting ...')
        elseif (count([mapped_weno, wenoz, teno]) >= 2) then
            call s_mpi_abort('Only one of mapped_weno, wenoz, or teno'// &
                             'can be set to true. Exiting ...')
        elseif (weno_order == 1 .and. mapped_weno) then
            call s_mpi_abort('mapped_weno is not supported for '// &
                             'weno_order = 1. Exiting ...')
        elseif (weno_order == 1 .and. wenoz) then
            call s_mpi_abort('wenoz is not supported for '// &
                             'weno_order = 1. Exiting ...')
        elseif (weno_order /= 5 .and. teno) then
            call s_mpi_abort('teno is only supported for '// &
                             'weno_order = 5. Exiting ...')
        elseif (weno_order /= 5 .and. mp_weno) then
            call s_mpi_abort('mp_weno is only supported for '// &
                             'weno_order = 5. Exiting ...')
        elseif (model_eqns == 1 .and. weno_avg) then
            call s_mpi_abort('weno_avg is not supported for '// &
                             'model_eqns = 1. Exiting ...')
        elseif (riemann_solver < 1 .or. riemann_solver > 3) then
            call s_mpi_abort('riemann_solver must be 1, 2, or 3. Exiting ...')
        elseif (all(wave_speeds /= (/dflt_int, 1, 2/))) then
            call s_mpi_abort('wave_speeds must be 1 or 2. Exiting ...')
        elseif (riemann_solver == 3 .and. wave_speeds /= dflt_int) then
            call s_mpi_abort('Exact Riemann (riemann_solver = 3) '// &
                             'does not support wave_speeds. Exiting ...')
            call s_mpi_abort('Unsupported value of avg_state. Exiting ...')
        elseif (riemann_solver /= 3 .and. wave_speeds == dflt_int) then
            call s_mpi_abort('wave_speeds must be set if '// &
                             'riemann_solver != 3. Exiting ...')
        elseif (riemann_solver /= 3 .and. avg_state == dflt_int) then
            call s_mpi_abort('avg_state must be set if '// &
                             'riemann_solver != 3. Exiting ...')
        elseif (bc_x%beg < -16 .or. bc_x%beg > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%beg. Exiting ...')
        elseif (bc_x%end < -16 .or. bc_x%end > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%end. Exiting ...')
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            call s_mpi_abort('bc_x%beg and bc_x%end must be both periodic '// &
                             '(= -1) or both non-periodic. Exiting ...')
        elseif (bc_y%beg /= dflt_int &
                .and. &
                (((cyl_coord .neqv. .true.) .and. (bc_y%beg < -16 .or. bc_y%beg > -1 .or. bc_y%beg == -14)) &
                 .or. &
                 (cyl_coord .and. p == 0 .and. bc_y%beg /= -2) &
                 .or. &
                 (cyl_coord .and. p > 0 .and. bc_y%beg /= -14))) then
            call s_mpi_abort('Unsupported value of bc_y%beg. Exiting ...')
        elseif (bc_y%end /= dflt_int &
                .and. &
                (bc_y%end < -16 .or. bc_y%end > -1 .or. bc_y%end == -14)) then
            call s_mpi_abort('Unsupported value of bc_y%end. Exiting ...')
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
        elseif (bc_z%beg /= dflt_int &
                .and. &
                (bc_z%beg < -16 .or. bc_z%beg > -1 .or. bc_z%beg == -14)) then
            call s_mpi_abort('Unsupported value of bc_z%beg. Exiting ...')
        elseif (bc_z%end /= dflt_int &
                .and. &
                (bc_z%end < -16 .or. bc_z%end > -1 .or. bc_z%end == -14)) then
            call s_mpi_abort('Unsupported value of bc_z%end. Exiting ...')
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
        elseif (model_eqns == 1 .and. alt_soundspeed) then
            call s_mpi_abort('model_eqns = 1 does not support alt_soundspeed. '// &
                             'Exiting ...')
        elseif (model_eqns == 4 .and. alt_soundspeed) then
            call s_mpi_abort('4-equation model (model_eqns = 4) does not '// &
                             'support alt_soundspeed. Exiting ...')
        elseif ((num_fluids /= 2 .and. num_fluids /= 3) .and. alt_soundspeed) then
            call s_mpi_abort('alt_soundspeed requires num_fluids = 2 or 3. '// &
                             'Exiting ...')
        elseif (riemann_solver /= 2 .and. alt_soundspeed) then
            call s_mpi_abort('alt_soundspeed requires HLLC Riemann solver '// &
                             '(riemann_solver = 2). Exiting ...')
        elseif (hypoelasticity .and. (riemann_solver /= 1)) then
            call s_mpi_abort('hypoelasticity requires HLL Riemann solver '// &
                             '(riemann_solver = 1). Exiting ...')
        end if

        if (adap_dt) then
            if (time_stepper /= 3) then
                call s_mpi_abort('adapt_dt requires Runge-Kutta 3 '// &
                                 '(time_stepper = 3). Exiting ...')
            else if (qbmm) then
                call s_mpi_abort('adapt_dt is not supported with QBMM. Exiting ...')
            else if (.not. polytropic) then
                call s_mpi_abort('adapt_dt is enabled, but polytropic is not. '// &
                                 'Exiting ...')
            else if (.not. adv_n) then
                call s_mpi_abort('adapt_dt is enabled, but adv_n is not. '// &
                                 'Exiting ...')
            end if
        end if
        ! END: Simulation Algorithm Parameters =============================

        ! Finite Difference Parameters =====================================
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            call s_mpi_abort('fd_order must be 1, 2, or 4. Exiting ...')
        elseif (probe_wrt .and. fd_order == dflt_int) then
            call s_mpi_abort('probe_wrt is enabled, but fd_order is not set. '// &
                             'Exiting ...')
        elseif (integral_wrt .and. (.not. bubbles)) then
            call s_mpi_abort('integral_wrt is enabled, but bubbles is not. '// &
                             'Exiting ...')
        end if
        ! END: Finite Difference Parameters ================================

        ! Fluids Physical Parameters =======================================
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

            do j = 1, 2
                call s_int_to_str(j, jStr)
                if ((.not. f_approx_equal(fluid_pp(i)%Re(j), dflt_real)) &
                    .and. &
                    fluid_pp(i)%Re(j) <= 0d0) then
                    call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//') must be positive. '// &
                                     'Exiting ...')
                end if

                if (model_eqns == 1 &
                    .and. &
                    (.not. f_approx_equal(fluid_pp(i)%Re(j), dflt_real))) then
                    call s_mpi_abort('model_eqns = 1 does not support '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

                if (i > num_fluids &
                    .and. &
                    (.not. f_approx_equal(fluid_pp(i)%Re(j), dflt_real))) then
                    call s_mpi_abort('First index ('//trim(iStr)//') of '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//') exceeds '// &
                                     'num_fluids. Exiting ...')
                end if

                if (weno_order == 1 &
                    .and. &
                    (weno_avg .neqv. .true.) &
                    .and. &
                    (.not. f_approx_equal(fluid_pp(i)%Re(j), dflt_real))) then
                    call s_mpi_abort('weno_order = 1 without weno_avg '// &
                                     'does not support '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

            end do

        end do
        ! END: Fluids Physical Parameters ==================================

        ! Constraints on the surface tension model
        if (.not. f_approx_equal(sigma, dflt_real) .and. sigma < 0d0) then
            call s_mpi_abort('The surface tension coefficient must be'// &
                             'greater than or equal to zero. Exiting ...')
        elseif (.not. f_approx_equal(sigma, dflt_real) .and. model_eqns /= 3) then
            call s_mpi_abort("The surface tension model requires"// &
                             'model_eqns=3. Exiting ...')
        end if

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

        ! Check IB parameters
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

        #:for DIR in ['x', 'y', 'z']
            if (bf_${DIR}$ .and. f_approx_equal(k_${DIR}$, dflt_real)) then
                call s_mpi_abort('k_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_approx_equal(w_${DIR}$, dflt_real)) then
                call s_mpi_abort('w_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_approx_equal(p_${DIR}$, dflt_real)) then
                call s_mpi_abort('p_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_approx_equal(g_${DIR}$, dflt_real)) then
                call s_mpi_abort('g_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            end if
        #:endfor

    end subroutine s_check_inputs

end module m_checker
