!>
!!@file m_checker.f90
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief The purpose of the module is to check for compatible input files
module m_checker

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file.
        !! Used by the pre_process stage
    impure subroutine s_check_inputs

        call s_check_parallel_io
        call s_check_inputs_restart
        call s_check_inputs_grid_stretching
        call s_check_inputs_qbmm_and_polydisperse
        call s_check_inputs_perturb_density
        call s_check_inputs_chemistry
        call s_check_inputs_misc
        call s_check_bc

    end subroutine s_check_inputs

    !> Checks if mpi is enabled with parallel_io
    impure subroutine s_check_parallel_io
#ifndef MFC_MPI
        @:PROHIBIT(parallel_io, "MFC built with --no-mpi requires parallel_io=F")
#endif
        @:PROHIBIT(down_sample .and. (.not. parallel_io), "down sample with parallel_io = T")
        @:PROHIBIT(down_sample .and. (.not. igr), "down sample with igr = T")
        @:PROHIBIT(down_sample .and. (p == 0), "down sample with 3D only")
        @:PROHIBIT(down_sample .and. (.not. file_per_process), "down sample with file_per_process = T")
        @:PROHIBIT(down_sample .and. (MOD(m+1,3) > 0), "down sample with m divisible by 3")
        @:PROHIBIT(down_sample .and. (MOD(n+1,3) > 0), "down sample with n divisible by 3")
        @:PROHIBIT(down_sample .and. (MOD(p+1,3) > 0), "down sample with p divisible by 3")

    end subroutine s_check_parallel_io

    !> Checks constraints on the restart parameters
        !! (old_grid, old_ic, etc.)
    impure subroutine s_check_inputs_restart
        logical :: skip_check !< Flag to skip the check when iterating over
        !! x, y, and z directions, for special treatment of cylindrical coordinates
        integer :: i

        @:PROHIBIT((.not. old_grid) .and. old_ic, &
            "old_ic can only be enabled with old_grid enabled")

        @:PROHIBIT(old_grid .and. t_step_old == dflt_int, &
            "old_grid requires t_step_old to be set")
        @:PROHIBIT(old_grid .and. ((.not. f_is_default(x_domain%beg)) .or. (.not. f_is_default(x_domain%end))), &
            "x_domain is not supported with old_grid enabled")
        @:PROHIBIT(old_grid .and. ((.not. f_is_default(y_domain%beg)) .or. (.not. f_is_default(y_domain%end))), &
            "y_domain is not supported with old_grid enabled")
        @:PROHIBIT(old_grid .and. ((.not. f_is_default(z_domain%beg)) .or. (.not. f_is_default(z_domain%end))), &
            "z_domain is not supported with old_grid enabled")

        #:for DIR, VAR in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            ! For cylindrical coordinates, the y and z directions use a different check
            #:if (DIR == 'y') or (DIR == 'z')
                skip_check = cyl_coord
            #:else
                skip_check = .false.
            #:endif

            if (.not. skip_check) then
                #:for BOUND in ['beg', 'end']
                    @:PROHIBIT(${VAR}$ == 0 .and. (.not. f_is_default((${DIR}$_domain%${BOUND}$))), &
                        "${DIR}$_domain%${BOUND}$ must not be set when ${VAR}$ = 0")
                    @:PROHIBIT(${VAR}$ > 0 .and. old_grid .and. (.not. f_is_default(${DIR}$_domain%${BOUND}$)), &
                        "${DIR}$_domain%${BOUND}$ must not be set when ${VAR}$ > 0 and old_grid = T")
                    @:PROHIBIT(${VAR}$ > 0 .and. (.not. old_grid) .and. f_is_default(${DIR}$_domain%${BOUND}$), &
                        "${DIR}$_domain%${BOUND}$ must be set when ${VAR}$ > 0 and old_grid = F")
                    @:PROHIBIT(${VAR}$ > 0 .and. (.not. old_grid) .and. ${DIR}$_domain%beg >= ${DIR}$_domain%end, &
                        "${DIR}$_domain%beg must be less than ${DIR}$_domain%end when both are set")
                #:endfor
            end if

        #:endfor

        @:PROHIBIT(cyl_coord .and. n == 0, &
            "n must be positive (2D or 3D) for cylindrical coordinates")
        @:PROHIBIT(cyl_coord .and. (f_is_default(y_domain%beg) .or. f_is_default(y_domain%end)), &
            "y_domain%beg and y_domain%end must be set for n = 0 (2D cylindrical coordinates)")
        @:PROHIBIT(cyl_coord .and. ((.not. f_approx_equal(y_domain%beg , 0._wp)) .or. y_domain%end <= 0._wp), &
            "y_domain%beg must be 0 and y_domain%end must be positive for cylindrical coordinates")
        @:PROHIBIT(cyl_coord .and. p == 0 .and. ((.not. f_is_default(z_domain%beg)) .or. (.not. f_is_default(z_domain%end))), &
            "z_domain%beg and z_domain%end are not supported for p = 0 (2D cylindrical coordinates)")
        @:PROHIBIT(cyl_coord .and. p > 0 .and. (.not. (f_approx_equal(z_domain%beg, 0._wp) .and. f_approx_equal(z_domain%end, 2._wp*pi))), &
            "z_domain%beg must be 0 and z_domain%end must be 2*pi for 3D cylindrical coordinates")

        @:PROHIBIT(num_patches < 0)
        @:PROHIBIT(num_patches == 0 .and. t_step_old == dflt_int, &
            "num_patches must be positive for the non-restart case")

    end subroutine s_check_inputs_restart

    !> Checks constraints on grid stretching parameters
        !! (loops_x[y,z], stretch_x[y,z], etc.)
    impure subroutine s_check_inputs_grid_stretching
        ! Constraints on loops for grid stretching
        @:PROHIBIT(loops_x < 1)
        @:PROHIBIT(loops_y < 1)

        ! Constraints specific to stretch_y
        @:PROHIBIT(stretch_y .and. n == 0)

        ! Constraints specific to stretch_z
        @:PROHIBIT(stretch_z .and. p == 0)
        @:PROHIBIT(stretch_z .and. cyl_coord)

        ! Common checks for all directions (stretch_x, stretch_y, and stretch_z)
        #:for X in ['x', 'y', 'z']
            @:PROHIBIT(stretch_${X}$ .and. old_grid, "old_grid and stretch_${X}$ are incompatible")
            @:PROHIBIT(stretch_${X}$ .and. f_is_default(a_${X}$), "a_${X}$ must be set with stretch_${X}$ enabled")
            @:PROHIBIT(stretch_${X}$ .and. f_is_default(${X}$_a), "${X}$_a must be set with stretch_${X}$ enabled")
            @:PROHIBIT(stretch_${X}$ .and. f_is_default(${X}$_b), "${X}$_b must be set with stretch_${X}$ enabled")
            @:PROHIBIT(stretch_${X}$ .and. ${X}$_a >= ${X}$_b, "${X}$_a must be less than ${X}$_b with stretch_${X}$ enabled")
            !&< Deactivate prettify
            @:PROHIBIT(stretch_${X}$ .and. (a_${X}$ + log(cosh(a_${X}$*(${X}$_domain%beg - ${X}$_a))) &
                                                    + log(cosh(a_${X}$*(${X}$_domain%beg - ${X}$_b))) &
                                                    - 2._wp*log(cosh(0.5_wp*a_${X}$*(${X}$_b - ${X}$_a)))) / a_${X}$ <= 0._wp, &
                "${X}$_domain%beg is too close to ${X}$_a and ${X}$_b for the given a_${X}$")
            @:PROHIBIT(stretch_${X}$ .and. (a_${X}$ + log(cosh(a_${X}$*(${X}$_domain%end - ${X}$_a))) &
                                                    + log(cosh(a_${X}$*(${X}$_domain%end - ${X}$_b))) &
                                                    - 2._wp*log(cosh(0.5_wp*a_${X}$*(${X}$_b - ${X}$_a)))) / a_${X}$ <= 0._wp, &
                "${X}$_domain%end is too close to ${X}$_a and ${X}$_b for the given a_${X}$")
            !&>
        #:endfor
    end subroutine s_check_inputs_grid_stretching

    !> Checks constraints on the QBMM and polydisperse bubble parameters
        !! (qbmm, polydisperse, dist_type and rhoRV)
    impure subroutine s_check_inputs_qbmm_and_polydisperse
        @:PROHIBIT(qbmm .and. dist_type == dflt_int, "dist_type must be set if using QBMM")
        @:PROHIBIT(qbmm .and. dist_type /= 1 .and. rhoRV > 0._wp, "rhoRV cannot be used with dist_type != 1")
    end subroutine s_check_inputs_qbmm_and_polydisperse

    !> Checks constraints on initial partial density perturbation
        !! (perturb_flow, perturb_flow_fluid, perturb_flow_mag, perturb_sph,
        !! perturb_sph_fluid, and fluid_rho)
    impure subroutine s_check_inputs_perturb_density
        character(len=5) :: iStr !< for int to string conversion
        integer :: i

        @:PROHIBIT(perturb_flow .and. (perturb_flow_fluid == dflt_int .or. f_is_default(perturb_flow_mag)), &
            "perturb_flow_fluid and perturb_flow_mag must be set with perturb_flow = T")
        @:PROHIBIT((.not. perturb_flow) .and. (perturb_flow_fluid /= dflt_int .or. (.not. f_is_default(perturb_flow_mag))), &
            "perturb_flow_fluid and perturb_flow_mag must not be set with perturb_flow = F")
        @:PROHIBIT(perturb_flow_fluid > num_fluids .or. (perturb_flow_fluid < 0 .and. perturb_flow_fluid /= dflt_int), &
            "perturb_flow_fluid must be between 0 and num_fluids")
        @:PROHIBIT(perturb_sph .and. perturb_sph_fluid == dflt_int, &
            "perturb_sph_fluid must be set with perturb_sph = T")
        @:PROHIBIT((.not. perturb_sph) .and. perturb_sph_fluid /= dflt_int, &
            "perturb_sph_fluid must not be set with perturb_sph = F")
        @:PROHIBIT(perturb_sph_fluid > num_fluids .or. (perturb_sph_fluid < 0 .and. perturb_sph_fluid /= dflt_int), &
            "perturb_sph_fluid must be between 0 and num_fluids")
        @:PROHIBIT((.not. perturb_sph) .and. (.not. f_all_default(fluid_rho)), &
            "fluid_rho must not be set with perturb_sph = F")

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)
            @:PROHIBIT(perturb_sph .and. f_is_default(fluid_rho(i)), &
                "fluid_rho("//trim(iStr)//") must be set if perturb_sph = T")
        end do
    end subroutine s_check_inputs_perturb_density

    impure subroutine s_check_inputs_chemistry

        if (chemistry) then
            @:ASSERT(num_species > 0)
        end if

    end subroutine s_check_inputs_chemistry

    !> Checks miscellaneous constraints
        !! (mixlayer_vel_profile and mixlayer_perturb)
    impure subroutine s_check_inputs_misc
        ! Hypertangent velocity profile
        @:PROHIBIT(mixlayer_vel_profile .and. (n == 0), &
            "mixlayer_vel_profile requires n > 0")

        ! Instability wave
        @:PROHIBIT(mixlayer_perturb .and. p == 0, "mixlayer_perturb requires p > 0")
        @:PROHIBIT(elliptic_smoothing .and. elliptic_smoothing_iters < 1, &
            "elliptic_smoothing_iters must be positive")

    end subroutine s_check_inputs_misc

    impure subroutine s_check_bc

        integer :: i
        character(len=5) :: iStr !< for int to string conversion

        call s_int_to_str(num_bc_patches_max, iStr)
        @:PROHIBIT(num_bc_patches > num_bc_patches_max, "num_bc_patches must be <= "//trim(iStr))

        do i = 1, num_bc_patches
            call s_int_to_str(i, iStr)

            ! Line Segment BC
            if (patch_bc(i)%geometry == 1) then
                @:PROHIBIT( .not. f_is_default(patch_bc(i)%radius), "Line Segment Patch can't have radius defined")
                #:for DIR in [('1'), ('2')]
                    @:PROHIBIT(patch_bc(i)%dir == ${DIR}$ .and. (.not. f_is_default(patch_bc(i)%centroid(${DIR}$)) .or. .not. f_is_default(patch_bc(i)%centroid(3))), &
                        "Line Segment Patch of Dir ${DIR}$ can't have a centroid in Dir ${DIR}$ or 3" )
                    @:PROHIBIT(patch_bc(i)%dir == ${DIR}$ .and. (.not. f_is_default(patch_bc(i)%length(${DIR}$)) .or. .not. f_is_default(patch_bc(i)%length(3))), &
                        "Line Segment Patch of Dir ${DIR}$ can't have a length in Dir ${DIR}$ or 3" )
                #:endfor
            end if

            ! Circle BC
            if (patch_bc(i)%geometry == 2) then
                @:PROHIBIT(f_is_default(patch_bc(i)%radius), "Circle Patch must have radius defined")
                @:PROHIBIT(.not. f_is_default(patch_bc(i)%length(1)) .or. .not. f_is_default(patch_bc(i)%length(2)) .or. .not. f_is_default(patch_bc(i)%length(3)), &
                    "Circle Patch can't have lengths defined")

                #:for DIR in [('1'), ('2'), ('3')]
                    @:PROHIBIT(patch_bc(i)%dir == ${DIR}$ .and. .not. f_is_default(patch_bc(i)%centroid(${DIR}$)), &
                        "Circle Patch of Dir ${DIR}$ can't have a centroid in Dir ${DIR}$")
                #:endfor
            end if

            ! Rectangle BC
            if (patch_bc(i)%geometry == 3) then
                @:PROHIBIT( .not. f_is_default(patch_bc(i)%radius), "Rectangle Patch can't have radius defined")

                #:for DIR, VAR in [('1', '2', '3'), ('2', '3', '1'), ('3', '1', '2')]
                    @:PROHIBIT(patch_bc(i)%dir == ${DIR}$ .and. .not. f_is_default(patch_bc(i)%centroid(${DIR}$)), &
                        "Rectangle Patch of Dir ${DIR}$ can't have a centroid in Dir ${DIR}$")
                    @:PROHIBIT(patch_bc(i)%dir == ${DIR}$ .and. .not. f_is_default(patch_bc(i)%length(${DIR}$)), &
                        "Rectangle Patch of Dir ${DIR}$ can't have a length in Dir ${DIR}$")
                #:endfor
            end if

            ! Incompatible BC check
            @:PROHIBIT(((patch_bc(i)%type >= BC_AXIS .and. patch_bc(i)%type <= BC_RIEMANN_EXTRAP) .or. &
                (patch_bc(i)%type == BC_PERIODIC) .or. patch_bc(i)%type < BC_DIRICHLET), &
                "Incompatible BC type for boundary condition patch "//trim(iStr))
        end do

    end subroutine s_check_bc

    impure subroutine s_check_moving_IBM

    end subroutine s_check_moving_IBM

end module m_checker
