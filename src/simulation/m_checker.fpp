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

    !> Checks compatibility of parameters in the input file.
        !! Used by the simulation stage
    subroutine s_check_inputs

        call s_check_inputs_compilers

        call s_check_inputs_weno
        call s_check_inputs_riemann_solver
        call s_check_inputs_time_stepping
        call s_check_inputs_model_eqns
        if (monopole) call s_check_inputs_monopole
        if (hypoelasticity) call s_check_inputs_hypoelasticity
        if (bubbles) call s_check_inputs_bubbles
        if (adap_dt) call s_check_inputs_adapt_dt
        if (alt_soundspeed) call s_check_inputs_alt_soundspeed
        call s_check_inputs_stiffened_eos_viscosity
        call s_check_inputs_body_forces
        call s_check_inputs_misc

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    subroutine s_check_inputs_compilers
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
    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    subroutine s_check_inputs_weno
        character(len=5) :: numStr !< for int to string conversion

        if (m + 1 < num_stcls_min*weno_order) then
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
        elseif (weno_order /= 1 .and. f_is_default(weno_eps)) then
            call s_mpi_abort('weno_order != 1, but weno_eps is not set. '// &
                             'A typical value of weno_eps is 1e-6. '// &
                             'Exiting ...')
        elseif (weno_eps <= 0d0) then
            call s_mpi_abort('weno_eps must be positive. '// &
                             'A typical value of weno_eps is 1e-6. '// &
                             'Exiting ...')
        elseif (teno .and. f_is_default(teno_CT)) then
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
    end subroutine s_check_inputs_weno

    !> Checks constraints on Riemann solver parameters
    subroutine s_check_inputs_riemann_solver
        if (riemann_solver /= 2 .and. model_eqns == 3) then
            call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                             'requires riemann_solver = 2. Exiting ...')
        elseif (riemann_solver < 1 .or. riemann_solver > 3) then
            call s_mpi_abort('riemann_solver must be 1, 2, or 3. Exiting ...')
        elseif (all(wave_speeds /= (/dflt_int, 1, 2/))) then
            call s_mpi_abort('wave_speeds must be 1 or 2. Exiting ...')
        elseif (riemann_solver == 3 .and. wave_speeds /= dflt_int) then
            call s_mpi_abort('Exact Riemann (riemann_solver = 3) '// &
                             'does not support wave_speeds. Exiting ...')
        elseif (all(avg_state /= (/dflt_int, 1, 2/))) then
            call s_mpi_abort('Unsupported value of avg_state. Exiting ...')
        elseif (riemann_solver /= 3 .and. wave_speeds == dflt_int) then
            call s_mpi_abort('wave_speeds must be set if '// &
                             'riemann_solver != 3. Exiting ...')
        elseif (riemann_solver /= 3 .and. avg_state == dflt_int) then
            call s_mpi_abort('avg_state must be set if '// &
                             'riemann_solver != 3. Exiting ...')
        end if
    end subroutine s_check_inputs_riemann_solver

    !> Checks constraints on time stepping parameters
    subroutine s_check_inputs_time_stepping
        if (dt <= 0) then
            call s_mpi_abort('dt must be positive. Exiting ...')
        end if

        if (time_stepper < 1 .or. time_stepper > 5) then
            if (time_stepper /= 23) then
                call s_mpi_abort('time_stepper must be between 1 and 5. '// &
                                 'Exiting ...')
            end if
        end if
    end subroutine s_check_inputs_time_stepping

    !> Checks constraints on parameters related to 6-equation model
    subroutine s_check_inputs_model_eqns
        if (model_eqns == 3) then
            if (avg_state /= 2) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'requires avg_state = 2. Exiting ...')
            elseif (wave_speeds /= 1) then
                call s_mpi_abort('6-equation model (model_eqns = 3) '// &
                                 'requires wave_speeds = 1. Exiting ...')
            end if
        end if
    end subroutine s_check_inputs_model_eqns

    !> Checks constraints on monopole parameters
    subroutine s_check_inputs_monopole
        integer :: j
        character(len=5) :: jStr

        if (num_mono == dflt_int) then
            call s_mpi_abort('num_mono must be specified for monopole. Exiting ...')
        elseif (num_mono < 0) then
            call s_mpi_abort('num_mono must be non-negative. Exiting ...')
        end if

        do j = 1, num_mono
            call s_int_to_str(j, jStr)
            if (mono(j)%support == dflt_int) then
                call s_mpi_abort('mono('//trim(jStr)//')%support must be '// &
                                 'specified. Exiting ...')
            elseif (f_is_default(mono(j)%mag)) then
                call s_mpi_abort('mono('//trim(jStr)//')%mag must be '// &
                                 'specified. Exiting ...')
            elseif (mono(j)%mag <= 0d0) then
                call s_mpi_abort('mono('//trim(jStr)//')%mag must be '// &
                                 'positive. Exiting ...')
            end if

            if (n == 0) then ! 1D
                if (.not. any(mono(j)%support == (/0, 1/))) then ! undocumented support 0
                    call s_mpi_abort('Only Mono(i)support = 1 is allowed for '// &
                                     '1D simulations. Exiting ...')
                end if
                if (mono(j)%support == 1 .and. f_is_default(mono(j)%loc(1))) then
                    call s_mpi_abort('mono_loc(1) must be specified for '// &
                                     'Mono(i)support = 1. Exiting ...')
                end if
            elseif (p == 0) then ! 2D
                if (.not. any(mono(j)%support == (/1, 2, 3, 4, 5/))) then
                    call s_mpi_abort('Only Mono(i)support = 1, 2, 3, 4, or 5 is '// &
                                     'allowed for 2D simulations. Exiting ...')
                end if
                if (any(mono(j)%support == (/1, 2, 3, 5/)) .and. &
                    (f_is_default(mono(j)%loc(1)) .or. &
                     f_is_default(mono(j)%loc(2)))) then
                    call s_mpi_abort('mono('//trim(jStr)//')%loc(1:2) must be '// &
                                     'specified for Mono(i)support = 1, 3, or 5. '// &
                                     'Exiting ...')
                elseif (mono(j)%support == 4 .and. f_is_default(mono(j)%loc(1))) then
                    call s_mpi_abort('mono('//trim(jStr)//')%loc(1) must be '// &
                                     'specified for Mono(i)support = 4. Exiting ...')
                end if
            else ! 3D
                if (.not. any(mono(j)%support == (/3, 4, 5, 6/))) then
                    call s_mpi_abort('Only Mono(i)support = 3, 4, 5, or 6 is '// &
                                     'allowed for 3D simulations. Exiting ...')
                elseif (mono(j)%support == 6 .and. (.not. cyl_coord)) then
                    call s_mpi_abort('Mono(i)support = 6 requires cyl_coord = true. '// &
                                     'Exiting ...')
                elseif (cyl_coord .and. mono(j)%support /= 6) then
                    call s_mpi_abort('cyl_coord = true requires Mono(i)support = 6. '// &
                                     'Exiting ...')
                end if
                if (any(mono(j)%support == (/3, 5, 6/)) .and. &
                    (f_is_default(mono(j)%loc(1)) .or. &
                     f_is_default(mono(j)%loc(2)) .or. &
                     f_is_default(mono(j)%loc(3)))) then
                    call s_mpi_abort('mono('//trim(jStr)//')%loc(1:3) must be '// &
                                     'specified for Mono(i)support = 3, 5, or 6. '// &
                                     'Exiting ...')
                elseif (mono(j)%support == 4 .and. &
                        (f_is_default(mono(j)%loc(3)))) then
                    call s_mpi_abort('mono('//trim(jStr)//')%loc(3) must be '// &
                                     'specified for Mono(i)support = 4. Exiting ...')
                end if
            end if
        end do

    end subroutine s_check_inputs_monopole

    !> Checks constraints on hypoelasticity parameters
    subroutine s_check_inputs_hypoelasticity
        if (riemann_solver /= 1) then
            call s_mpi_abort('hypoelasticity requires HLL Riemann solver '// &
                             '(riemann_solver = 1). Exiting ...')
        end if
    end subroutine

    !> Checks constraints on bubble parameters
    subroutine s_check_inputs_bubbles
        if (riemann_solver /= 2) then
            call s_mpi_abort('Bubble modeling requires riemann_solver = 2')
        elseif (avg_state /= 2) then
            call s_mpi_abort('Bubble modeling requires arithmetic average '// &
                             '(avg_state = 2). Exiting ...')
        elseif (model_eqns == 2 .and. bubble_model == 1) then
            call s_mpi_abort('The 5-equation bubbly flow model requires '// &
                             'bubble_model = 2 (Keller--Miksis). Exiting ...')
        end if
    end subroutine s_check_inputs_bubbles

    !> Checks constraints on adaptive time stepping parameters (adap_dt)
    subroutine s_check_inputs_adapt_dt
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
    end subroutine s_check_inputs_adapt_dt

    !> Checks constraints on alternative sound speed parameters (alt_soundspeed)
    subroutine s_check_inputs_alt_soundspeed
        if (model_eqns /= 2) then
            call s_mpi_abort('5-equation model (model_eqns = 2) '// &
                             'is required for alt_soundspeed. Exiting ...')
        elseif (num_fluids /= 2 .and. num_fluids /= 3) then
            call s_mpi_abort('alt_soundspeed requires num_fluids = 2 or 3. '// &
                             'Exiting ...')
        elseif (riemann_solver /= 2) then
            call s_mpi_abort('alt_soundspeed requires HLLC Riemann solver '// &
                             '(riemann_solver = 2). Exiting ...')
        end if
    end subroutine s_check_inputs_alt_soundspeed

    !> Checks constraints on viscosity parameters (fluid_pp(i)%Re(1:2))
        !! of the stiffened gas equation of state
    subroutine s_check_inputs_stiffened_eos_viscosity
        character(len=5) :: iStr, jStr
        integer :: i, j

        do i = 1, num_fluids
            do j = 1, 2
                call s_int_to_str(j, jStr)
                if ((.not. f_is_default(fluid_pp(i)%Re(j))) &
                    .and. &
                    fluid_pp(i)%Re(j) <= 0d0) then
                    call s_mpi_abort('fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//') must be positive. '// &
                                     'Exiting ...')
                else if (model_eqns == 1 &
                         .and. &
                         (.not. f_is_default(fluid_pp(i)%Re(j)))) then
                    call s_mpi_abort('model_eqns = 1 does not support '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                else if (i > num_fluids &
                         .and. &
                         (.not. f_is_default(fluid_pp(i)%Re(j)))) then
                    call s_mpi_abort('First index ('//trim(iStr)//') of '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//') exceeds '// &
                                     'num_fluids. Exiting ...')
                else if (weno_order == 1 .and. (.not. weno_avg) &
                         .and. &
                         (.not. f_is_default(fluid_pp(i)%Re(j)))) then
                    call s_mpi_abort('weno_order = 1 without weno_avg '// &
                                     'does not support '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if
            end do
        end do
    end subroutine s_check_inputs_stiffened_eos_viscosity

    !> Checks constraints on body forces parameters (bf_x[y,z], etc.)
    subroutine s_check_inputs_body_forces
        #:for DIR in ['x', 'y', 'z']
            if (bf_${DIR}$ .and. f_is_default(k_${DIR}$)) then
                call s_mpi_abort('k_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_is_default(w_${DIR}$)) then
                call s_mpi_abort('w_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_is_default(p_${DIR}$)) then
                call s_mpi_abort('p_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. f_is_default(g_${DIR}$)) then
                call s_mpi_abort('g_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            end if
        #:endfor
    end subroutine s_check_inputs_body_forces

    !> Checks miscellaneous constraints,
        !! including constraints on probe_wrt and integral_wrt
    subroutine s_check_inputs_misc
        ! Write probe data
        if (probe_wrt .and. fd_order == dflt_int) then
            call s_mpi_abort('probe_wrt is enabled, but fd_order is not set. '// &
                             'Exiting ...')
            ! Write integral data for bubbles
        elseif (integral_wrt .and. (.not. bubbles)) then
            call s_mpi_abort('integral_wrt is enabled, but bubbles is not. '// &
                             'Exiting ...')
        end if
    end subroutine s_check_inputs_misc

end module m_checker
