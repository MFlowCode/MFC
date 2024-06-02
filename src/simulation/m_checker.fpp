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

        character(len=5) :: iStr
        character(len=5) :: jStr
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
            call s_mpi_abort('Unsupported value of m. Exiting ...')
        elseif (n < 0) then
            call s_mpi_abort('Unsupported value of n. Exiting ...')
        elseif (p < 0) then
            call s_mpi_abort('Unsupported value of p. Exiting ...')
        elseif (cyl_coord .and. p > 0 .and. mod(p, 2) /= 1) then
            call s_mpi_abort('Unsupported value of p. Exiting ...')
        elseif (n == 0 .and. p > 0) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'n and p. Exiting ...')
        elseif (dt <= 0) then
            call s_mpi_abort('Unsupported value of dt. Exiting ...')
        elseif (t_step_start < 0) then
            call s_mpi_abort('Unsupported value of t_step_start. Exiting ...')
        elseif (t_step_stop <= t_step_start) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             't_step_start and t_step_stop. '// &
                             'Exiting ...')
        elseif (t_step_save > t_step_stop - t_step_start) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             't_step_start, t_step_stop and '// &
                             't_step_save. Exiting ...')
        end if
        ! ==================================================================

        ! Simulation Algorithm Parameters ==================================
        if (all(model_eqns /= (/1, 2, 3, 4/))) then
            call s_mpi_abort('Unsupported value of model_eqns. Exiting ...')
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
            elseif (avg_state == 1) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'bubbles and Roe average (please use avg_state = 2). '// &
                                 'Exiting ...')
            end if
        end if

        if (model_eqns == 4 .and. num_fluids /= 1) then
            call s_mpi_abort('The 4-equation model implementation is not a multi-component and requires num_fluids = 1')
        end if

        if ((bubbles .neqv. .true.) .and. polydisperse) then
            call s_mpi_abort('Polydisperse bubble modeling requires the bubble switch to be activated')
        elseif (polydisperse .and. (poly_sigma == dflt_real)) then
            call s_mpi_abort('Polydisperse bubble modeling requires poly_sigma > 0')
        elseif (qbmm .and. (bubbles .neqv. .true.)) then
            call s_mpi_abort('QBMM requires bubbles')
        elseif (qbmm .and. (nnode /= 4)) then
            call s_mpi_abort('nnode not supported')
        end if

        if (model_eqns == 3) then
            if (riemann_solver /= 2) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'model_eqns (6-eq) and riemann_solver (please use riemann_solver = 2). '// &
                                 'Exiting ...')
            elseif (alt_soundspeed) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'model_eqns (6-eq) and alt_soundspeed. '// &
                                 'Exiting ...')
            elseif (avg_state == 1) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'model_eqns (6-eq) and Roe average (please use avg_state = 2). '// &
                                 'Exiting ...')
            elseif (wave_speeds == 2) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'model_eqns (6-eq) and wave_speeds (please use wave_speeds = 1). '// &
                                 'Exiting ...')
            elseif (cyl_coord .and. p /= 0) then
                call s_mpi_abort('Unsupported combination of values of '// &
                                 'model_eqns (6-eq) and cylindrical coordinates. '// &
                                 'Exiting ...')
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
            elseif ((palpha_eps <= 0d0) .or. (palpha_eps >= 1d0) .or. &
                    (ptgalpha_eps <= 0d0) .or. (ptgalpha_eps >= 1d0)) then
                call s_mpi_abort('both palpha_eps and ptgalpha_eps must &
 &               be in (0,1). '//'Exiting ...')
            end if
        elseif ((relax_model /= dflt_int) .or. (palpha_eps /= dflt_real) &
                .or. (ptgalpha_eps /= dflt_real)) then
            call s_mpi_abort('relax is not set as true, but other phase change parameters have &
&               been modified. Either activate phase change or set the values to default. '//'Exiting ...')
        end if

        if (num_fluids /= dflt_int &
            .and. &
            (num_fluids < 1 .or. num_fluids > num_fluids)) then
            call s_mpi_abort('Unsupported value of num_fluids. Exiting ...')
        elseif ((model_eqns == 1 .and. num_fluids /= dflt_int) &
                .or. &
                (model_eqns == 2 .and. num_fluids == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and num_fluids. '// &
                             'Exiting ...')
        elseif (model_eqns == 1 .and. adv_alphan) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and adv_alphan. '// &
                             'Exiting ...')
        elseif (model_eqns == 1 .and. mpp_lim) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and mpp_lim. Exiting ...')
        elseif (num_fluids == 1 .and. mpp_lim) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'num_fluids and mpp_lim. Exiting ...')
        elseif (time_stepper < 1 .or. time_stepper > 5) then
            if (time_stepper /= 23) then
                call s_mpi_abort('Unsupported value of time_stepper. Exiting ...')
            end if
        elseif (all(weno_order /= (/1, 3, 5/))) then
            call s_mpi_abort('Unsupported value of weno_order. Exiting ...')
        elseif (m + 1 < num_stcls_min*weno_order) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'm and weno_order. Exiting ...')
        elseif (n + 1 < min(1, n)*num_stcls_min*weno_order) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'n and weno_order. Exiting ...')
        elseif (p + 1 < min(1, p)*num_stcls_min*weno_order) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'p and weno_order. Exiting ...')
        elseif (weno_eps <= 0d0 .or. weno_eps > 1d-6) then
            call s_mpi_abort('Unsupported value of weno_eps. Exiting ...')
        elseif (weno_order == 1 .and. mapped_weno) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'weno_order and mapped_weno. '// &
                             'Exiting ...')
        elseif (weno_order /= 5 .and. mp_weno) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'weno_order and mp_weno. Exiting ...')
        elseif (model_eqns == 1 .and. weno_avg) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'model_eqns and weno_avg. Exiting ...')
        elseif (riemann_solver < 1 .or. riemann_solver > 3) then
            call s_mpi_abort('Unsupported value of riemann_solver. Exiting ...')
        elseif (all(wave_speeds /= (/dflt_int, 1, 2/))) then
            call s_mpi_abort('Unsupported value of wave_speeds. Exiting ...')
        elseif ((riemann_solver /= 3 .and. wave_speeds == dflt_int) &
                .or. &
                (riemann_solver == 3 .and. wave_speeds /= dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'riemann_solver and wave_speeds. '// &
                             'Exiting ...')
        elseif (all(avg_state /= (/dflt_int, 1, 2/))) then
            call s_mpi_abort('Unsupported value of avg_state. Exiting ...')
        elseif (riemann_solver /= 3 .and. avg_state == dflt_int) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'riemann_solver and avg_state. '// &
                             'Exiting ...')
        elseif (bc_x%beg < -16 .or. bc_x%beg > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%beg. Exiting ...')
        elseif (bc_x%end < -16 .or. bc_x%end > -1 .or. bc_x%beg == -14) then
            call s_mpi_abort('Unsupported value of bc_x%end. Exiting ...')
        elseif ((bc_x%beg == -1 .and. bc_x%end /= -1) &
                .or. &
                (bc_x%end == -1 .and. bc_x%beg /= -1)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'bc_x%beg and bc_x%end. Exiting ...')
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
        elseif ((n == 0 .and. bc_y%beg /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%beg == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'n and bc_y%beg. Exiting ...')
        elseif ((n == 0 .and. bc_y%end /= dflt_int) &
                .or. &
                (n > 0 .and. bc_y%end == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'n and bc_y%end. Exiting ...')
        elseif ((bc_y%beg == -1 .and. bc_y%end /= -1) &
                .or. &
                (bc_y%end == -1 .and. bc_y%beg /= -1)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'bc_y%beg and bc_y%end. Exiting ...')
        elseif (bc_z%beg /= dflt_int &
                .and. &
                (bc_z%beg < -16 .or. bc_z%beg > -1 .or. bc_z%beg == -14)) then
            call s_mpi_abort('Unsupported value of bc_z%beg. Exiting ...')
        elseif (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -13)) then
            call s_mpi_abort('Unsupported choice of boundary condition -13')
        elseif (bc_z%end /= dflt_int &
                .and. &
                (bc_z%end < -16 .or. bc_z%end > -1 .or. bc_z%end == -14)) then
            call s_mpi_abort('Unsupported value of bc_z%end. Exiting ...')
        elseif ((p == 0 .and. bc_z%beg /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%beg == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'p and bc_z%beg. Exiting ...')
        elseif ((p == 0 .and. bc_z%end /= dflt_int) &
                .or. &
                (p > 0 .and. bc_z%end == dflt_int)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'p and bc_z%end. Exiting ...')
        elseif ((bc_z%beg == -1 .and. bc_z%end /= -1) &
                .or. &
                (bc_z%end == -1 .and. bc_z%beg /= -1)) then
            call s_mpi_abort('Unsupported combination of values of '// &
                             'bc_z%beg and bc_z%end. Exiting ...')
        elseif (model_eqns == 1 .and. alt_soundspeed) then
            call s_mpi_abort('Unsupported combination of model_eqns '// &
                             'and alt_soundspeed. Exiting ...')
        elseif (model_eqns == 4 .and. alt_soundspeed) then
            call s_mpi_abort('Unsupported combination of model_eqns '// &
                             'and alt_soundspeed. Exiting ...')
        elseif ((num_fluids /= 2 .and. num_fluids /= 3) .and. alt_soundspeed) then
            call s_mpi_abort('Unsupported combination of num_fluids '// &
                             'and alt_soundspeed. Exiting ...')
        elseif (riemann_solver /= 2 .and. alt_soundspeed) then
            call s_mpi_abort('Unsupported combination of riemann_solver '// &
                             'and alt_soundspeed. Exiting ...')
        elseif (hypoelasticity .and. (riemann_solver /= 1)) then
            call s_mpi_abort('hypoelasticity requires riemann_solver = 1'// &
                             'Exiting ...')
        end if

        if (adap_dt) then
            if (time_stepper /= 3) then
                call s_mpi_abort('Unsupported combination of adap_dt '// &
                                 'and time_stepper. Exiting ...')
            else if (qbmm) then
                call s_mpi_abort('Unsupported combination of adap_dt '// &
                                 'and qbmm. Exiting ...')
            else if (.not. polytropic) then
                call s_mpi_abort('Unsupported combination of adap_dt '// &
                                 'and polytropic. Exiting ...')
            else if (.not. adv_n) then
                call s_mpi_abort('Unsupported combination of adap_dt '// &
                                 'and adv_n. Exiting ...')
            end if
        end if
        ! END: Simulation Algorithm Parameters =============================

        ! Finite Difference Parameters =====================================
        if (fd_order /= dflt_int &
            .and. &
            fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'fd_order. Exiting ...')
        elseif (probe_wrt .and. fd_order == dflt_int) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for probe_wrt, and fd_order. '// &
                             'Exiting ...')
        elseif (integral_wrt .and. (bubbles .neqv. .true.)) then
            call s_mpi_abort('Unsupported choice of the combination of '// &
                             'values for integral_wrt, and bubbles. '// &
                             'Exiting ...')
        end if
        ! END: Finite Difference Parameters ================================

        ! Fluids Physical Parameters =======================================
        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            if (fluid_pp(i)%gamma /= dflt_real &
                .and. &
                fluid_pp(i)%gamma <= 0d0) then
                call s_mpi_abort('Unsupported value of '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'gamma. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%gamma /= dflt_real) then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of model_eqns '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
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
                call s_mpi_abort('Unsupported value of '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'pi_inf. Exiting ...')
            elseif (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%pi_inf /= dflt_real) then
                call s_mpi_abort('Unsupported combination '// &
                                 'of values of model_eqns '// &
                                 'and fluid_pp('//trim(iStr)//')%'// &
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
                call s_mpi_abort('Unsupported value of '// &
                                 'fluid_pp('//trim(iStr)//')%'// &
                                 'cv. Make sure cv is positive. Exiting ...')
            end if

            do j = 1, 2
                call s_int_to_str(j, jStr)
                if (fluid_pp(i)%Re(j) /= dflt_real &
                    .and. &
                    fluid_pp(i)%Re(j) <= 0d0) then
                    call s_mpi_abort('Unsupported value of '// &
                                     'fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

                if (model_eqns == 1 &
                    .and. &
                    fluid_pp(i)%Re(j) /= dflt_real) then
                    call s_mpi_abort('Unsupported combination '// &
                                     'of values of model_eqns '// &
                                     'and fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

                if (i > num_fluids &
                    .and. &
                    fluid_pp(i)%Re(j) /= dflt_real) then
                    call s_mpi_abort('Unsupported combination '// &
                                     'of values of num_fluids '// &
                                     'and fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

                if (weno_order == 1 &
                    .and. &
                    (weno_avg .neqv. .true.) &
                    .and. &
                    fluid_pp(i)%Re(j) /= dflt_real) then
                    call s_mpi_abort('Unsupported combination '// &
                                     'of values of weno_order, '// &
                                     'weno_avg and fluid_pp('//trim(iStr)//')%'// &
                                     'Re('//trim(jStr)//'). Exiting ...')
                end if

            end do

        end do
        ! END: Fluids Physical Parameters ==================================

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
                call s_mpi_abort('Unsupported choices of the combination of values for '// &
                                 'ib and n. Immersed Boundaries do not work in 1D. Exiting ...')
            else if (num_ibs <= 0 .or. num_ibs > num_patches_max) then
                call s_mpi_abort('Unsupported choice for the value of '// &
                                 'num_ibs. Exiting ...')
            end if
        end if

        if (num_ibs > 0 .and. .not. ib) then
            call s_mpi_abort('Unsupported choices of the combination of values for '// &
                             'num_ibs and ib. Exiting ...')
        end if

        #:for DIR in ['x', 'y', 'z']
            if (bf_${DIR}$ .and. k_${DIR}$ == dflt_real) then
                call s_mpi_abort('k_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. w_${DIR}$ == dflt_real) then
                call s_mpi_abort('w_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. p_${DIR}$ == dflt_real) then
                call s_mpi_abort('p_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            elseif (bf_${DIR}$ .and. g_${DIR}$ == dflt_real) then
                call s_mpi_abort('g_${DIR}$ must be specified if bf_${DIR}$ is true '// &
                                 'Exiting ...')
            end if
        #:endfor

    end subroutine s_check_inputs

end module m_checker
