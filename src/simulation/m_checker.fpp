!>
!!@file m_start_up.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'
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
        call s_check_inputs_acoustic_src
        call s_check_inputs_hypoelasticity
        call s_check_inputs_bubbles
        call s_check_inputs_adapt_dt
        call s_check_inputs_alt_soundspeed
        call s_check_inputs_stiffened_eos_viscosity
        call s_check_inputs_body_forces
        call s_check_inputs_misc

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    subroutine s_check_inputs_compilers
#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

#ifndef MFC_cuTENSOR
        @:PROHIBIT(cu_tensor, "MFC was not built with the NVIDIA cuTENSOR library")
#endif
    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    subroutine s_check_inputs_weno
        character(len=5) :: numStr !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
            "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
            "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
            "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(weno_order /= 1 .and. f_is_default(weno_eps), &
            "weno_order != 1, but weno_eps is not set. A typical value of weno_eps is 1e-6")
        @:PROHIBIT(weno_eps <= 0d0, "weno_eps must be positive. A typical value of weno_eps is 1e-6")
        @:PROHIBIT(teno .and. f_is_default(teno_CT), "teno is used, but teno_CT is not set. A typical value of teno_CT is 1e-6")
        @:PROHIBIT(teno .and. teno_CT <= 0d0, "teno_CT must be positive. A typical value of teno_CT is 1e-6")
        @:PROHIBIT(count([mapped_weno, wenoz, teno]) >= 2, "Only one of mapped_weno, wenoz, or teno can be set to true")
        @:PROHIBIT(weno_order == 1 .and. mapped_weno)
        @:PROHIBIT(weno_order == 1 .and. wenoz)
        @:PROHIBIT(weno_order /= 5 .and. teno)
        @:PROHIBIT(weno_order /= 5 .and. mp_weno)
        @:PROHIBIT(model_eqns == 1 .and. weno_avg)
    end subroutine s_check_inputs_weno

    !> Checks constraints on Riemann solver parameters
    subroutine s_check_inputs_riemann_solver
        @:PROHIBIT(riemann_solver /= 2 .and. model_eqns == 3, "6-equation model (model_eqns = 3) requires riemann_solver = 2")
        @:PROHIBIT(riemann_solver < 1 .or. riemann_solver > 3, "riemann_solver must be 1, 2, or 3")
        @:PROHIBIT(all(wave_speeds /= (/dflt_int, 1, 2/)), "wave_speeds must be 1 or 2")
        @:PROHIBIT(riemann_solver == 3 .and. wave_speeds /= dflt_int, "Exact Riemann (riemann_solver = 3) does not support wave_speeds")
        @:PROHIBIT(all(avg_state /= (/dflt_int, 1, 2/)), "Unsupported value of avg_state")
        @:PROHIBIT(riemann_solver /= 3 .and. wave_speeds == dflt_int, "wave_speeds must be set if riemann_solver != 3")
        @:PROHIBIT(riemann_solver /= 3 .and. avg_state == dflt_int, "avg_state must be set if riemann_solver != 3")
        @:PROHIBIT(all(low_Mach /= (/0, 1, 2/)), "low_Mach must be 0, 1 or 2")
        @:PROHIBIT(riemann_solver /= 2 .and. low_Mach /= 0, "low_Mach = 1 or 2 requires riemann_solver = 2")
        @:PROHIBIT(low_Mach /= 0 .and. model_eqns /= 2, "low_Mach = 1 or 2 requires model_eqns = 2")
    end subroutine s_check_inputs_riemann_solver

    !> Checks constraints on time stepping parameters
    subroutine s_check_inputs_time_stepping
        @:PROHIBIT(dt <= 0)
        @:PROHIBIT(time_stepper < 1 .or. time_stepper > 5)
    end subroutine s_check_inputs_time_stepping

    !> Checks constraints on parameters related to 6-equation model
    subroutine s_check_inputs_model_eqns
        @:PROHIBIT(model_eqns == 3 .and. avg_state /= 2, "6-equation model (model_eqns = 3) requires avg_state = 2")
        @:PROHIBIT(model_eqns == 3 .and. wave_speeds /= 1, "6-equation model (model_eqns = 3) requires wave_speeds = 1")
    end subroutine s_check_inputs_model_eqns

    !> Checks constraints on acoustic_source parameters
    subroutine s_check_inputs_acoustic_src

        integer :: j, dim
        character(len=5) :: jStr

        !! When it's obvious that the checks are only relevant if acoustic_source is enabled,
        !! `acoustic_source .and.` is removed from the conditions for clarity.
        !! `if (.not. acoustic_source) return` ensures equivalent behavior
        if (.not. acoustic_source) return

        if (n == 0) then
            dim = 1
        else if (p == 0) then
            dim = 2
        else
            dim = 3
        end if

        @:PROHIBIT(acoustic_source .and. num_source == dflt_int, "num_source must be specified for acoustic_source")
        @:PROHIBIT(acoustic_source .and. num_source < 0, "num_source must be non-negative")

        do j = 1, num_source
            call s_int_to_str(j, jStr)

            @:PROHIBIT(acoustic_source .and. acoustic(j)%support == dflt_int, &
                "acoustic("//trim(jStr)//")%support must be specified for acoustic_source")

            @:PROHIBIT(dim == 1 .and. acoustic(j)%support /= 1, &
                "Only acoustic("//trim(jStr)//")%support = 1 is allowed for 1D simulations")
            @:PROHIBIT(dim == 1 .and. acoustic(j)%support == 1 .and. f_is_default(acoustic(j)%loc(1)), &
                "acoustic("//trim(jStr)//")%loc(1) must be specified for acoustic("//trim(jStr)//")%support = 1")
            @:PROHIBIT(dim == 2 .and. (.not. any(acoustic(j)%support == (/2, 5, 6, 9, 10/))), &
                "Only acoustic("//trim(jStr)//")%support = 2, 5, 6, 9, or 10 is allowed for 2D simulations")
            @:PROHIBIT(dim == 2 .and. (.not. any(acoustic(j)%support == (/6, 10/))) .and. cyl_coord, &
                "Only acoustic("//trim(jStr)//")%support = 6 or 10 is allowed for 2D axisymmetric simulations")
            @:PROHIBIT(dim == 2 .and. any(acoustic(j)%support == (/2, 5, 6, 9, 10/)) .and. &
                (f_is_default(acoustic(j)%loc(1)) .or. f_is_default(acoustic(j)%loc(2))), &
                "acoustic("//trim(jStr)//")%loc(1:2) must be specified for acoustic("//trim(jStr)//")%support = 2")
            @:PROHIBIT(dim == 3 .and. (.not. any(acoustic(j)%support == (/3, 7, 11/))), &
                "Only acoustic("//trim(jStr)//")%support = 3, 7, or 11 is allowed for 3D simulations")
            @:PROHIBIT(dim == 3 .and. cyl_coord, &
                "Acoustic source is not supported in 3D cylindrical simulations")
            @:PROHIBIT(dim == 3 .and. acoustic(j)%support == 3 .and. &
                (f_is_default(acoustic(j)%loc(1)) .or. f_is_default(acoustic(j)%loc(2))), &
                "acoustic("//trim(jStr)//")%loc(1:2) must be specified for acoustic("//trim(jStr)//")%support = 3")
            @:PROHIBIT(dim == 3 .and. any(acoustic(j)%support == (/7, 11/)) .and. &
                (f_is_default(acoustic(j)%loc(1)) .or. &
                 f_is_default(acoustic(j)%loc(2)) .or. &
                 f_is_default(acoustic(j)%loc(3))), &
                "acoustic("//trim(jStr)//")%loc(1:3) must be specified for acoustic("//trim(jStr)//")%support = 7 or 11")

            @:PROHIBIT(f_is_default(acoustic(j)%mag), &
                "acoustic("//trim(jStr)//")%mag must be specified")
            @:PROHIBIT(acoustic(j)%pulse == dflt_int, &
                "acoustic("//trim(jStr)//")%pulse must be specified")
            @:PROHIBIT(.not. any(acoustic(j)%pulse == (/1, 2, 3/)), &
                "Only acoustic("//trim(jStr)//")%pulse = 1, 2, or 3 is allowed")

            @:PROHIBIT(any(acoustic(j)%pulse == (/1, 3/)) .and. &
                (f_is_default(acoustic(j)%frequency) .eqv. f_is_default(acoustic(j)%wavelength)), &
                "One and only one of acoustic("//trim(jStr)//")%frequency "// &
                "or acoustic("//trim(jStr)//")%wavelength must be specified for pulse = 1 or 3")
            @:PROHIBIT(acoustic(j)%pulse == 2 .and. &
                (f_is_default(acoustic(j)%gauss_sigma_time) .eqv. f_is_default(acoustic(j)%gauss_sigma_dist)), &
                "One and only one of acoustic("//trim(jStr)//")%gauss_sigma_time "// &
                "or acoustic("//trim(jStr)//")%gauss_sigma_dist must be specified for pulse = 2")

            @:PROHIBIT(f_is_default(acoustic(j)%npulse), &
                "acoustic("//trim(jStr)//")%npulse must be specified")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. (.not. f_is_integer(acoustic(j)%npulse)), &
                "acoustic("//trim(jStr)//")%npulse must be an integer for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%npulse >= 5 .and. acoustic(j)%dipole, &
                "acoustic("//trim(jStr)//")%dipole is not supported for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support < 5 .and. f_is_default(acoustic(j)%dir), &
                "acoustic("//trim(jStr)//")%dir must be specified for support < 5 (planer support)")
            @:PROHIBIT(acoustic(j)%support == 1 .and. f_approx_equal(acoustic(j)%dir, 0d0), &
                "acoustic("//trim(jStr)//")dir must be non-zero for support = 1")
            @:PROHIBIT(acoustic(j)%pulse == 2 .and. f_is_default(acoustic(j)%delay), &
                "acoustic("//trim(jStr)//")%delay must be specified for pulse = 2 (Gaussian)")
            @:PROHIBIT(acoustic(j)%pulse == 3 .and. acoustic(j)%support >= 5, &
                "acoustic("//trim(jStr)//")%support >= 5 (Cylindrical or Spherical support) is not allowed for pulse = 3 (square wave)")

            @:PROHIBIT((acoustic(j)%support == 2 .or. acoustic(j)%support == 3) .and. f_is_default(acoustic(j)%length), &
                "acoustic("//trim(jStr)//")%length must be specified for support = 2 or 3")
            @:PROHIBIT((acoustic(j)%support == 2 .or. acoustic(j)%support == 3) .and. acoustic(j)%length <= 0d0, &
                "acoustic("//trim(jStr)//")%length must be positive for support = 2 or 3")
            @:PROHIBIT(acoustic(j)%support == 3 .and. f_is_default(acoustic(j)%height), &
                "acoustic("//trim(jStr)//")%height must be specified for support = 3")
            @:PROHIBIT(acoustic(j)%support == 3 .and. acoustic(j)%height <= 0d0, &
                "acoustic("//trim(jStr)//")%height must be positive for support = 3")

            @:PROHIBIT(acoustic(j)%support >= 5 .and. f_is_default(acoustic(j)%foc_length), &
                "acoustic("//trim(jStr)//")%foc_length must be specified for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. acoustic(j)%foc_length <= 0d0, &
                "acoustic("//trim(jStr)//")%foc_length must be positive for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. f_is_default(acoustic(j)%aperture), &
                "acoustic("//trim(jStr)//")%aperture must be specified for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. acoustic(j)%aperture <= 0d0, &
                "acoustic("//trim(jStr)//")%aperture must be positive for support >= 5 (non-planar supports)")

            @:PROHIBIT(any(acoustic(j)%support == (/9, 10, 11/)) .and. acoustic(j)%num_elements == dflt_int, &
                "acoustic("//trim(jStr)//")%num_elements must be specified for support = 9, 10, or 11 (transducer array)")
            @:PROHIBIT(any(acoustic(j)%support == (/9, 10, 11/)) .and. acoustic(j)%num_elements <= 0, &
                "acoustic("//trim(jStr)//")%num_elements must be positive for support = 9, 10, or 11 (transducer array)")
            @:PROHIBIT(acoustic(j)%element_on /= dflt_int .and. acoustic(j)%element_on < 0, &
                "acoustic("//trim(jStr)//")%element_on must be non-negative for support = 9, 10, or 11 (transducer array)")
            @:PROHIBIT(acoustic(j)%element_on /= dflt_int .and. acoustic(j)%element_on > acoustic(j)%num_elements, &
                "acoustic("//trim(jStr)//")%element_on must be less than or equal to num_elements for support = 9, 10, or 11 (transducer array)")
            @:PROHIBIT(any(acoustic(j)%support == (/9, 10/)) .and. f_is_default(acoustic(j)%element_spacing_angle), &
                "acoustic("//trim(jStr)//")%element_spacing_angle must be specified for support = 9 or 10 (2D transducer array)")
            @:PROHIBIT(any(acoustic(j)%support == (/9, 10/)) .and. acoustic(j)%element_spacing_angle < 0d0, &
                "acoustic("//trim(jStr)//")%element_spacing_angle must be non-negative for support = 9 or 10 (2D transducer array)")
            @:PROHIBIT(acoustic(j)%support == 11 .and. f_is_default(acoustic(j)%element_polygon_ratio), &
                "acoustic("//trim(jStr)//")%element_polygon_ratio must be specified for support = 11 (3D transducer array)")
            @:PROHIBIT(acoustic(j)%support == 11 .and. acoustic(j)%element_polygon_ratio <= 0d0, &
                "acoustic("//trim(jStr)//")%element_polygon_ratio must be positive for support = 11 (3D transducer array)")
        end do

    end subroutine s_check_inputs_acoustic_src

    !> Checks constraints on hypoelasticity parameters
    subroutine s_check_inputs_hypoelasticity
        @:PROHIBIT(hypoelasticity .and. riemann_solver /= 1, "hypoelasticity requires HLL Riemann solver (riemann_solver = 1)")
    end subroutine

    !> Checks constraints on bubble parameters
    subroutine s_check_inputs_bubbles
        @:PROHIBIT(bubbles .and. riemann_solver /= 2, "Bubble modeling requires HLLC Riemann solver (riemann_solver = 2)")
        @:PROHIBIT(bubbles .and. avg_state /= 2, "Bubble modeling requires arithmetic average (avg_state = 2)")
        @:PROHIBIT(bubbles .and. model_eqns == 2 .and. bubble_model == 1, &
            "The 5-equation bubbly flow model does not support bubble_model = 1 (Gilmore)")
    end subroutine s_check_inputs_bubbles

    !> Checks constraints on adaptive time stepping parameters (adap_dt)
    subroutine s_check_inputs_adapt_dt
        @:PROHIBIT(adap_dt .and. time_stepper /= 3, "adapt_dt requires Runge-Kutta 3 (time_stepper = 3)")
        @:PROHIBIT(adap_dt .and. qbmm)
        @:PROHIBIT(adap_dt .and. (.not. polytropic))
        @:PROHIBIT(adap_dt .and. (.not. adv_n))
    end subroutine s_check_inputs_adapt_dt

    !> Checks constraints on alternative sound speed parameters (alt_soundspeed)
    subroutine s_check_inputs_alt_soundspeed
        @:PROHIBIT(alt_soundspeed .and. model_eqns /= 2, "5-equation model (model_eqns = 2) is required for alt_soundspeed")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver /= 2, "alt_soundspeed requires HLLC Riemann solver (riemann_solver = 2)")
        @:PROHIBIT(alt_soundspeed .and. num_fluids /= 2 .and. num_fluids /= 3)
    end subroutine s_check_inputs_alt_soundspeed

    !> Checks constraints on viscosity parameters (fluid_pp(i)%Re(1:2))
        !! of the stiffened gas equation of state
    subroutine s_check_inputs_stiffened_eos_viscosity
        character(len=5) :: iStr, jStr
        integer :: i, j

        do i = 1, num_fluids
            do j = 1, 2
                call s_int_to_str(j, jStr)
                @:PROHIBIT((.not. f_is_default(fluid_pp(i)%Re(j))) .and. fluid_pp(i)%Re(j) <= 0d0, &
                    "fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//") must be positive.")
                @:PROHIBIT(model_eqns == 1 .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                    "model_eqns = 1 does not support fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//")")
                @:PROHIBIT(i > num_fluids .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                    "First index ("//trim(iStr)//") of fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//") exceeds num_fluids")
                @:PROHIBIT(weno_order == 1 .and. (.not. weno_avg) .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                    "weno_order = 1 without weno_avg does not support fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//")")
            end do
        end do
    end subroutine s_check_inputs_stiffened_eos_viscosity

    !> Checks constraints on body forces parameters (bf_x[y,z], etc.)
    subroutine s_check_inputs_body_forces
        #:for DIR in ['x', 'y', 'z']
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(k_${DIR}$), "k_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(w_${DIR}$), "w_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(p_${DIR}$), "p_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(g_${DIR}$), "g_${DIR}$ must be specified if bf_${DIR}$ is true")
        #:endfor
    end subroutine s_check_inputs_body_forces

    !> Checks miscellaneous constraints,
        !! including constraints on probe_wrt and integral_wrt
    subroutine s_check_inputs_misc
        @:PROHIBIT(probe_wrt .and. fd_order == dflt_int, "fd_order must be specified for probe_wrt")
        @:PROHIBIT(integral_wrt .and. (.not. bubbles))
    end subroutine s_check_inputs_misc

end module m_checker
