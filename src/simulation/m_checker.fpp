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

    use m_helper_basic         !< Functions to compare floating point numbers

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_igr
        else
            if (recon_type == WENO_TYPE) then
                call s_check_inputs_weno
            else if (recon_type == MUSCL_TYPE) then
                call s_check_inputs_muscl
            end if
            call s_check_inputs_riemann_solver
            call s_check_inputs_model_eqns
            call s_check_inputs_acoustic_src
            call s_check_inputs_hypoelasticity
            call s_check_inputs_bubbles_euler
            call s_check_inputs_bubbles_lagrange
            call s_check_inputs_adapt_dt
            call s_check_inputs_alt_soundspeed
            call s_check_inputs_grcbc
            call s_check_inputs_geometry_precision
            call s_check_inputs_mhd
            call s_check_inputs_continuum_damage
        end if

        call s_check_inputs_time_stepping
        call s_check_inputs_stiffened_eos_viscosity
        call s_check_inputs_body_forces
        call s_check_inputs_misc

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers
#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif
    end subroutine s_check_inputs_compilers

    impure subroutine s_check_inputs_igr
        @:PROHIBIT(num_igr_iters < 0, "num_igr_iters must be greater than 0")
        @:PROHIBIT(num_igr_warm_start_iters < 0, "num_igr_warm_start_iters must be greater than 0")
        @:PROHIBIT((igr_iter_solver /= 1 .and. igr_iter_solver /= 2), &
            "igr_iter_solver must be 1 or 2")
        @:PROHIBIT(alf_factor < 0, "alf factor must be positive")
        @:PROHIBIT(model_eqns /= 2, "IGR only supports model_eqns = 2")
        @:PROHIBIT(ib, "IGR does not support the immersed boundary method")
        @:PROHIBIT(bubbles_euler, "IGR does not support Euler-Euler bubble models")
        @:PROHIBIT(bubbles_lagrange, "IGR does not support Euler-Lagrange bubbles models")
        @:PROHIBIT(alt_soundspeed, "IGR does not support alt_soundspeed = T")
        @:PROHIBIT(surface_tension, "IGR does not support surface tension")
        @:PROHIBIT(hypoelasticity, "IGR does not support hypoelasticity")
        @:PROHIBIT(acoustic_source, "IGR does not support acoustic sources")
        @:PROHIBIT(relax, "IGR does not support phase change")
        @:PROHIBIT(mhd, "IGR does not support magnetohydrodynamics")
        @:PROHIBIT(hyperelasticity, "IGR does not support hyperelasticity")
        @:PROHIBIT(cyl_coord, "IGR does not support cylindrical or axisymmetric coordinates")
        @:PROHIBIT(probe_wrt, "IGR does not support probe writes")

        #:for DIR in [('x'), ('y'), ('z')]
            #:for LOC in [('beg'), ('end')]
                @:PROHIBIT((bc_${DIR}$%${LOC}$ <= -4 .and. bc_${DIR}$%${LOC}$ >= -14), &
                    "Characteristic boundary condition bc_${DIR}$%${LOC}$ is not compatible with IGR")
            #:endfor
        #:endfor

    end subroutine s_check_inputs_igr

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno
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

        @:PROHIBIT(weno_eps <= 0._wp, "weno_eps must be positive. A typical value of weno_eps is 1e-6")
        @:PROHIBIT(wenoz .and. weno_order == 7 .and. f_is_default(real(wenoz_q, wp)), &
            "wenoz is used at 7th order, but wenoz_q is not set. It should be either 2, 3, or 4")
        @:PROHIBIT(wenoz .and. weno_order == 7 .and. .not. (f_approx_equal(real(wenoz_q, wp), real(2, wp)) .or. &
            f_approx_equal(real(wenoz_q, wp), real(3, wp)) .or. f_approx_equal(real(wenoz_q, wp), real(4, wp))), &
            "wenoz_q must be either 2, 3, or 4")
        @:PROHIBIT(teno .and. f_is_default(teno_CT), "teno is used, but teno_CT is not set. A typical value of teno_CT is 1e-6")
        @:PROHIBIT(teno .and. teno_CT <= 0._wp, "teno_CT must be positive. A typical value of teno_CT is 1e-6")
        @:PROHIBIT(count([mapped_weno, wenoz, teno]) >= 2, "Only one of mapped_weno, wenoz, or teno can be set to true")
        @:PROHIBIT(weno_order == 1 .and. mapped_weno)
        @:PROHIBIT(weno_order == 1 .and. wenoz)
        @:PROHIBIT((weno_order == 1 .or. weno_order == 3) .and. teno)
        @:PROHIBIT(weno_order /= 5 .and. mp_weno)
        @:PROHIBIT(model_eqns == 1 .and. weno_avg)
    end subroutine s_check_inputs_weno

    impure subroutine s_check_inputs_muscl
        character(len=5) :: numStr !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
            "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
            "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
            "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is "//trim(numStr))
        @:PROHIBIT(muscl_order == 2 .and. muscl_lim == dflt_int, "muscl_lim must be defined if using muscl_order = 2")
        @:PROHIBIT(muscl_order /= 1 .and. muscl_lim < 1 .or. muscl_lim > 5, "muscl_lim must be 1,2,3,4, or 5")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on Riemann solver parameters
    impure subroutine s_check_inputs_riemann_solver
        @:PROHIBIT(riemann_solver /= 2 .and. model_eqns == 3, "6-equation model (model_eqns = 3) requires riemann_solver = 2")
        @:PROHIBIT(riemann_solver < 1 .or. riemann_solver > 4, "riemann_solver must be 1, 2, 3, or 4")
        @:PROHIBIT(all(wave_speeds /= (/dflt_int, 1, 2/)), "wave_speeds must be 1 or 2")
        @:PROHIBIT(riemann_solver == 3 .and. wave_speeds /= dflt_int, "Exact Riemann (riemann_solver = 3) does not support wave_speeds")
        @:PROHIBIT(all(avg_state /= (/dflt_int, 1, 2/)), "Unsupported value of avg_state")
        @:PROHIBIT(riemann_solver /= 3 .and. wave_speeds == dflt_int, "wave_speeds must be set if riemann_solver != 3")
        @:PROHIBIT(riemann_solver /= 3 .and. avg_state == dflt_int, "avg_state must be set if riemann_solver != 3")
        @:PROHIBIT(all(low_Mach /= (/0, 1, 2/)), "low_Mach must be 0, 1 or 2")
        @:PROHIBIT(riemann_solver /= 2 .and. low_Mach == 2, "low_Mach = 2 requires riemann_solver = 2")
        @:PROHIBIT(low_Mach /= 0 .and. all(model_eqns /= (/2, 3/)), "low_Mach = 1 or 2 requires model_eqns = 2 or 3")
    end subroutine s_check_inputs_riemann_solver

    !> Checks constraints on geometry and precision
    impure subroutine s_check_inputs_geometry_precision
        ! Prevent spherical geometry in single precision
#ifdef MFC_SINGLE_PRECISION
        @:PROHIBIT(.not. (cyl_coord .neqv. .true. .or. (cyl_coord .and. p == 0)), "Fully 3D cylindrical grid (geometry = 3) is not supported in single precision.")
#endif
    end subroutine s_check_inputs_geometry_precision

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping
        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if
        @:PROHIBIT(time_stepper < 1 .or. time_stepper > 3)
    end subroutine s_check_inputs_time_stepping

    !> Checks constraints on parameters related to 6-equation model
    impure subroutine s_check_inputs_model_eqns
        @:PROHIBIT(model_eqns == 3 .and. avg_state /= 2, "6-equation model (model_eqns = 3) requires avg_state = 2")
        @:PROHIBIT(model_eqns == 3 .and. wave_speeds /= 1, "6-equation model (model_eqns = 3) requires wave_speeds = 1")
    end subroutine s_check_inputs_model_eqns

    !> Checks constraints for GRCBC
    impure subroutine s_check_inputs_grcbc
        #:for DIR in ['x', 'y', 'z']
            @:PROHIBIT(bc_${DIR}$%grcbc_in .and. (bc_${DIR}$%beg /= -7 .and. bc_${DIR}$%end /= -7), "Subsonic Inflow requires bc = -7")
            @:PROHIBIT(bc_${DIR}$%grcbc_out .and. (bc_${DIR}$%beg /= -8 .and. bc_${DIR}$%end /= -8), "Subsonic Outflow requires bc = -8")
            @:PROHIBIT(bc_${DIR}$%grcbc_vel_out .and. (bc_${DIR}$%beg /= -8 .and. bc_${DIR}$%end /= -8), "Subsonic Outflow requires bc = -8")
        #:endfor
    end subroutine s_check_inputs_grcbc

    !> Checks constraints on acoustic_source parameters
    impure subroutine s_check_inputs_acoustic_src

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
            @:PROHIBIT((dim == 2  .and. .not. cyl_coord) .and. (.not. any(acoustic(j)%support == (/2, 5, 9/))), &
                "Only acoustic("//trim(jStr)//")%support = 2, 5, 6, 9, or 10 is allowed for 2D simulations")
            @:PROHIBIT((dim == 2  .and. cyl_coord) .and. (.not. any(acoustic(j)%support == (/2, 6, 10/))), &
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
            @:PROHIBIT(.not. any(acoustic(j)%pulse == (/1, 2, 3, 4/)), &
                "Only acoustic("//trim(jStr)//")%pulse = 1, 2, 3 or 4 is allowed")

            @:PROHIBIT(any(acoustic(j)%pulse == (/1, 3/)) .and. &
                (f_is_default(acoustic(j)%frequency) .eqv. f_is_default(acoustic(j)%wavelength)), &
                "One and only one of acoustic("//trim(jStr)//")%frequency "// &
                "or acoustic("//trim(jStr)//")%wavelength must be specified for pulse = 1 or 3")
            @:PROHIBIT(acoustic(j)%pulse == 2 .and. &
                (f_is_default(acoustic(j)%gauss_sigma_time) .eqv. f_is_default(acoustic(j)%gauss_sigma_dist)), &
                "One and only one of acoustic("//trim(jStr)//")%gauss_sigma_time "// &
                "or acoustic("//trim(jStr)//")%gauss_sigma_dist must be specified for pulse = 2")
            @:PROHIBIT(acoustic(j)%pulse == 4 .and. acoustic(j)%bb_num_freq == dflt_int, &
                "The number of broadband frequencies acoustic("//trim(jStr)//")%bb_num_freq must be specified for pulse = 4")
            @:PROHIBIT(acoustic(j)%pulse == 4 .and. f_is_default(acoustic(j)%bb_bandwidth), &
                "The broadband wave band width acoustic("//trim(jStr)//")%bb_bandwidth must be specified for pulse = 4")
            @:PROHIBIT(acoustic(j)%pulse == 4 .and. f_is_default(acoustic(j)%bb_lowest_freq), &
                "The broadband wave lower frequency bound acoustic("//trim(jStr)//")%bb_lowest_freq must be specified for pulse = 4")

            @:PROHIBIT(f_is_default(acoustic(j)%npulse), &
                "acoustic("//trim(jStr)//")%npulse must be specified")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. (.not. f_is_integer(acoustic(j)%npulse)), &
                "acoustic("//trim(jStr)//")%npulse must be an integer for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%npulse >= 5 .and. acoustic(j)%dipole, &
                "acoustic("//trim(jStr)//")%dipole is not supported for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support < 5 .and. f_is_default(acoustic(j)%dir), &
                "acoustic("//trim(jStr)//")%dir must be specified for support < 5 (planer support)")
            @:PROHIBIT(acoustic(j)%support == 1 .and. f_approx_equal(acoustic(j)%dir, 0._wp), &
                "acoustic("//trim(jStr)//")dir must be non-zero for support = 1")
            @:PROHIBIT(acoustic(j)%pulse == 2 .and. f_is_default(acoustic(j)%delay), &
                "acoustic("//trim(jStr)//")%delay must be specified for pulse = 2 (Gaussian)")
            @:PROHIBIT(acoustic(j)%pulse == 3 .and. acoustic(j)%support >= 5, &
                "acoustic("//trim(jStr)//")%support >= 5 (Cylindrical or Spherical support) is not allowed for pulse = 3 (square wave)")

            @:PROHIBIT((acoustic(j)%support == 2 .or. acoustic(j)%support == 3) .and. f_is_default(acoustic(j)%length), &
                "acoustic("//trim(jStr)//")%length must be specified for support = 2 or 3")
            @:PROHIBIT((acoustic(j)%support == 2 .or. acoustic(j)%support == 3) .and. acoustic(j)%length <= 0._wp, &
                "acoustic("//trim(jStr)//")%length must be positive for support = 2 or 3")
            @:PROHIBIT(acoustic(j)%support == 3 .and. f_is_default(acoustic(j)%height), &
                "acoustic("//trim(jStr)//")%height must be specified for support = 3")
            @:PROHIBIT(acoustic(j)%support == 3 .and. acoustic(j)%height <= 0._wp, &
                "acoustic("//trim(jStr)//")%height must be positive for support = 3")

            @:PROHIBIT(acoustic(j)%support >= 5 .and. f_is_default(acoustic(j)%foc_length), &
                "acoustic("//trim(jStr)//")%foc_length must be specified for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. acoustic(j)%foc_length <= 0._wp, &
                "acoustic("//trim(jStr)//")%foc_length must be positive for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. f_is_default(acoustic(j)%aperture), &
                "acoustic("//trim(jStr)//")%aperture must be specified for support >= 5 (non-planar supports)")
            @:PROHIBIT(acoustic(j)%support >= 5 .and. acoustic(j)%aperture <= 0._wp, &
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
            @:PROHIBIT(any(acoustic(j)%support == (/9, 10/)) .and. acoustic(j)%element_spacing_angle < 0._wp, &
                "acoustic("//trim(jStr)//")%element_spacing_angle must be non-negative for support = 9 or 10 (2D transducer array)")
            @:PROHIBIT(acoustic(j)%support == 11 .and. f_is_default(acoustic(j)%element_polygon_ratio), &
                "acoustic("//trim(jStr)//")%element_polygon_ratio must be specified for support = 11 (3D transducer array)")
            @:PROHIBIT(acoustic(j)%support == 11 .and. acoustic(j)%element_polygon_ratio <= 0._wp, &
                "acoustic("//trim(jStr)//")%element_polygon_ratio must be positive for support = 11 (3D transducer array)")
        end do

    end subroutine s_check_inputs_acoustic_src

    !> Checks constraints on hypoelasticity parameters
    impure subroutine s_check_inputs_hypoelasticity
        @:PROHIBIT(hypoelasticity .and. riemann_solver /= 1, "hypoelasticity requires HLL Riemann solver (riemann_solver = 1)")
    end subroutine

    !> Checks constraints on bubble parameters
    impure subroutine s_check_inputs_bubbles_euler
        @:PROHIBIT(bubbles_euler .and. bubbles_lagrange, "Activate only one of the bubble subgrid models")
        @:PROHIBIT(bubbles_euler .and. riemann_solver /= 2, "Bubble modeling requires HLLC Riemann solver (riemann_solver = 2)")
        @:PROHIBIT(bubbles_euler .and. avg_state /= 2, "Bubble modeling requires arithmetic average (avg_state = 2)")
        @:PROHIBIT(bubbles_euler .and. model_eqns == 2 .and. bubble_model == 1, &
            "The 5-equation bubbly flow model does not support bubble_model = 1 (Gilmore)")
    end subroutine s_check_inputs_bubbles_euler

    !> Checks constraints on adaptive time stepping parameters (adap_dt)
    impure subroutine s_check_inputs_adapt_dt
        @:PROHIBIT(adap_dt .and. time_stepper /= 3, "adapt_dt requires Runge-Kutta 3 (time_stepper = 3)")
        @:PROHIBIT(adap_dt .and. qbmm)
        @:PROHIBIT(adap_dt .and. (.not. polytropic) .and. (.not. bubbles_lagrange))
        @:PROHIBIT(adap_dt .and. (.not. adv_n) .and. (.not. bubbles_lagrange))
    end subroutine s_check_inputs_adapt_dt

    !> Checks constraints on alternative sound speed parameters (alt_soundspeed)
    impure subroutine s_check_inputs_alt_soundspeed
        @:PROHIBIT(alt_soundspeed .and. model_eqns /= 2, "5-equation model (model_eqns = 2) is required for alt_soundspeed")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver /= 2, "alt_soundspeed requires HLLC Riemann solver (riemann_solver = 2)")
        @:PROHIBIT(alt_soundspeed .and. num_fluids /= 2 .and. num_fluids /= 3)
    end subroutine s_check_inputs_alt_soundspeed

    !> Checks constraints on viscosity parameters (fluid_pp(i)%Re(1:2))
        !! of the stiffened gas equation of state
    impure subroutine s_check_inputs_stiffened_eos_viscosity
        character(len=5) :: iStr, jStr
        integer :: i, j

        do i = 1, num_fluids
            do j = 1, 2
                call s_int_to_str(j, jStr)
                @:PROHIBIT((.not. f_is_default(fluid_pp(i)%Re(j))) .and. fluid_pp(i)%Re(j) <= 0._wp, &
                    "fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//") must be positive.")
                @:PROHIBIT(model_eqns == 1 .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                    "model_eqns = 1 does not support fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//")")
                @:PROHIBIT(i > num_fluids .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                    "First index ("//trim(iStr)//") of fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//") exceeds num_fluids")
                if (.not. igr) then
                    @:PROHIBIT(weno_order == 1 .and. (.not. weno_avg) .and. (.not. f_is_default(fluid_pp(i)%Re(j))), &
                        "weno_order = 1 without weno_avg does not support fluid_pp("//trim(iStr)//")%"// "Re("//trim(jStr)//")")
                end if
            end do
            @:PROHIBIT(.not. f_is_default(fluid_pp(i)%Re(1)) .and. .not. viscous, &
                "Re(1) is specified, but viscous is not set to true")
            @:PROHIBIT(.not. f_is_default(fluid_pp(i)%Re(2)) .and. .not. viscous, &
                "Re(2) is specified, but viscous is not set to true")
            @:PROHIBIT(f_is_default(fluid_pp(i)%Re(1)) .and. viscous, &
                "Re(1) is not specified, but viscous is set to true")
        end do

    end subroutine s_check_inputs_stiffened_eos_viscosity

    !> Checks constraints on body forces parameters (bf_x[y,z], etc.)
    impure subroutine s_check_inputs_body_forces
        #:for DIR in ['x', 'y', 'z']
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(k_${DIR}$), "k_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(w_${DIR}$), "w_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(p_${DIR}$), "p_${DIR}$ must be specified if bf_${DIR}$ is true")
            @:PROHIBIT(bf_${DIR}$ .and. f_is_default(g_${DIR}$), "g_${DIR}$ must be specified if bf_${DIR}$ is true")
        #:endfor
    end subroutine s_check_inputs_body_forces

    !> Checks constraints on lagrangian bubble parameters
    impure subroutine s_check_inputs_bubbles_lagrange
        @:PROHIBIT(bubbles_lagrange .and. file_per_process, "file_per_process must be false for bubbles_lagrange")
        @:PROHIBIT(bubbles_lagrange .and. n==0, "bubbles_lagrange accepts 2D and 3D simulations only")
        @:PROHIBIT(bubbles_lagrange .and. model_eqns==3, "The 6-equation flow model does not support bubbles_lagrange")
        @:PROHIBIT(bubbles_lagrange .and. lag_params%cluster_type>=2 .and. lag_params%smooth_type/=1, "cluster_type=2 requires smooth_type=1")
    end subroutine s_check_inputs_bubbles_lagrange

    !> Checks constraints on continuum damage model parameters
    impure subroutine s_check_inputs_continuum_damage
        @:PROHIBIT(cont_damage .and. f_is_default(tau_star))
        @:PROHIBIT(cont_damage .and. f_is_default(cont_damage_s))
        @:PROHIBIT(cont_damage .and. f_is_default(alpha_bar))
    end subroutine s_check_inputs_continuum_damage

    !> Checks miscellaneous constraints,
        !! including constraints on probe_wrt and integral_wrt
    impure subroutine s_check_inputs_misc
        @:PROHIBIT(probe_wrt .and. fd_order == dflt_int, "fd_order must be specified for probe_wrt")
        @:PROHIBIT(integral_wrt .and. (.not. bubbles_euler))
    end subroutine s_check_inputs_misc

    impure subroutine s_check_inputs_mhd
        @:PROHIBIT(mhd .and. (riemann_solver /= 1 .and. riemann_solver /= 4), &
            "MHD simulations require riemann_solver = 1 (HLL) or riemann_solver = 4 (HLLD)")
        @:PROHIBIT(riemann_solver == 4 .and. .not. mhd, "HLLD is only available for MHD simulations")
        @:PROHIBIT(riemann_solver == 4 .and. relativity, "HLLD is not available for RMHD")
        @:PROHIBIT(powell .and. .not. mhd)
        @:PROHIBIT(powell .and. n == 0, "Powell's method is not supported for 1D simulations")
        @:PROHIBIT(powell .and. fd_order == dflt_int, "fd_order must be set if Powell's method is enabled")
    end subroutine s_check_inputs_mhd

end module m_checker
