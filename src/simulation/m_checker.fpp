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
        integer :: bub_fac !< For allowing an extra fluid_pp if there are subgrid bubbles
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
        if (dt <= 0) then
            call s_mpi_abort('dt must be positive. Exiting ...')
        end if
        ! END: Computational Domain Parameters =============================

        ! Simulation Algorithm Parameters ==================================
        ! Constraints on WENO scheme parameters
        if (all(weno_order /= (/1, 3, 5/))) then
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
        end if

        ! Constraints on Riemann solver parameters
        if (riemann_solver < 1 .or. riemann_solver > 3) then
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
        end if

        ! Constraints on time stepping parameters
        if (time_stepper < 1 .or. time_stepper > 5) then
            if (time_stepper /= 23) then
                call s_mpi_abort('time_stepper must be between 1 and 5. '// &
                                 'Exiting ...')
            end if
        end if

        ! Constraints on pairing parameters with 6-equation model
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

        ! Constraints on hypoelasticity parameters
        if (hypoelasticity .and. (model_eqns /= 2)) then
            call s_mpi_abort('hypoelasticity requires 5-equation model'// &
                             '(model_eqns = 2). Exiting ...')
        elseif (hypoelasticity .and. (riemann_solver /= 1)) then
            call s_mpi_abort('hypoelasticity requires HLL Riemann solver '// &
                             '(riemann_solver = 1). Exiting ...')
        end if

        ! Constraints on bubble parameters
        if (bubbles) then
            if (riemann_solver /= 2) then
                call s_mpi_abort('Bubble modeling requires riemann_solver = 2')
            elseif (avg_state /= 2) then
                call s_mpi_abort('Bubble modeling requires arithmetic average '// &
                                 '(avg_state = 2). Exiting ...')
            elseif (model_eqns == 2 .and. bubble_model == 1) then
                call s_mpi_abort('The 5-equation bubbly flow model requires '// &
                                 'bubble_model = 2 (Keller--Miksis). Exiting ...')
            end if
        end if

        ! Constraints on adaptive time stepping
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

        ! Constraints on alt_soundspeed
        if (model_eqns == 1 .and. alt_soundspeed) then
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
        end if
        ! END: Simulation Algorithm Parameters =============================

        ! Fluids Physical Parameters =======================================
        do i = 1, num_fluids
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
        ! END: Fluids Physical Parameters ===================================

        ! Probe Parameters =================================================
        if (probe_wrt .and. fd_order == dflt_int) then
            call s_mpi_abort('probe_wrt is enabled, but fd_order is not set. '// &
                             'Exiting ...')
        elseif (integral_wrt .and. (.not. bubbles)) then
            call s_mpi_abort('integral_wrt is enabled, but bubbles is not. '// &
                             'Exiting ...')
        end if
        ! END: Probe Parameters ============================================

        ! Constraints on body force parameters
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
