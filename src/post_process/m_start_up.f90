!>
!! @file m_start_up.f90
!! @brief  Contains module m_start_up

!> @brief This module contains the subroutines that read in and check the
!!              consistency of the user provided inputs. This module also allocates, initializes and
!!              deallocates the relevant variables and sets up the time stepping,
!!              MPI decomposition and I/O procedures
module m_start_up

    ! Dependencies

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_mpi_common            !< Common MPI subroutines

    use m_boundary_common       !< Common boundary conditions subroutines

    use m_variables_conversion  !< Subroutines to change the state variables from
                                !! one form to another

    use m_data_input            !< Procedures reading raw simulation data to fill
                                !! the conservative, primitive and grid variables

    use m_data_output           !< Procedures that write the grid and chosen flow
                                !! variable(s) to the formatted database file(s)

    use m_derived_variables     !< Procedures used to compute quantities derived
                                !! from the conservative and primitive variables
    use m_helper

    use m_compile_specific

    use m_checker_common

    use m_checker

    use m_thermochem, only: num_species, species_names

    use m_finite_differences

    use m_chemistry

    implicit none

contains

    !>  Reads the configuration file post_process.inp, in order
        !!      to populate parameters in module m_global_parameters.f90
        !!      with the user provided inputs
    impure subroutine s_read_input_file

        character(LEN=name_len) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: file_check !<
            !! Generic logical used for the purpose of asserting whether a file
            !! is or is not present in the designated location

        integer :: iostatus
            !! Integer to check iostat of file read

        character(len=1000) :: line

        ! Namelist for all of the parameters to be inputted by the user
        namelist /user_inputs/ case_dir, m, n, p, t_step_start, &
            t_step_stop, t_step_save, model_eqns, &
            num_fluids, mpp_lim, &
            weno_order, bc_x, &
            bc_y, bc_z, fluid_pp, format, precision, &
            output_partial_domain, x_output, y_output, z_output, &
            hypoelasticity, G, mhd, &
            chem_wrt_Y, chem_wrt_T, avg_state, &
            alpha_rho_wrt, rho_wrt, mom_wrt, vel_wrt, &
            E_wrt, pres_wrt, alpha_wrt, gamma_wrt, &
            heat_ratio_wrt, pi_inf_wrt, pres_inf_wrt, &
            cons_vars_wrt, prim_vars_wrt, c_wrt, &
            omega_wrt, qm_wrt, schlieren_wrt, schlieren_alpha, &
            fd_order, mixture_err, alt_soundspeed, &
            flux_lim, flux_wrt, cyl_coord, &
            parallel_io, rhoref, pref, bubbles_euler, qbmm, sigR, &
            R0ref, nb, polytropic, thermal, Ca, Web, Re_inv, &
            polydisperse, poly_sigma, file_per_process, relax, &
            relax_model, cf_wrt, sigma, adv_n, ib, num_ibs, &
            cfl_adap_dt, cfl_const_dt, t_save, t_stop, n_start, &
            cfl_target, surface_tension, bubbles_lagrange, &
            sim_data, hyperelasticity, Bx0, relativity, cont_damage, &
            num_bc_patches, igr, igr_order, recon_type, muscl_order

        ! Inquiring the status of the post_process.inp file
        file_loc = 'post_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        ! Checking whether the input file is there. If it is, the input file
        ! is read. If not, the program is terminated.
        if (file_check) then
            open (1, FILE=trim(file_loc), FORM='formatted', &
                  STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: '//trim(line)
                call s_mpi_abort('Invalid line in post_process.inp. It is '// &
                                 'likely due to a datatype mismatch. Exiting.')
            end if

            close (1)

            call s_update_cell_bounds(cells_bounds, m, n, p)

            ! Store m,n,p into global m,n,p
            m_glb = m
            n_glb = n
            p_glb = p

            nGlobal = (m_glb + 1)*(n_glb + 1)*(p_glb + 1)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -17) .or. &
                num_bc_patches > 0) then
                bc_io = .true.
            end if

        else
            call s_mpi_abort('File post_process.inp is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !>  Checking that the user inputs make sense, i.e. that the
        !!      individual choices are compatible with the code's options
        !!      and that the combination of these choices results into a
        !!      valid configuration for the post-process
    impure subroutine s_check_input_file

        character(LEN=len_trim(case_dir)) :: file_loc !<
            !! Generic string used to store the address of a particular file

        logical :: dir_check !<
            !! Logical variable used to test the existence of folders

        ! Checking the existence of the case folder
        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir)//'/.'

        call my_inquire(file_loc, dir_check)

        ! Constraint on the location of the case directory
        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Unsupported choice for the value of '// &
                             'case_dir. Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

    end subroutine s_check_input_file

    impure subroutine s_perform_time_step(t_step)

        integer, intent(inout) :: t_step
        if (proc_rank == 0) then
            if (cfl_dt) then
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, "")', &
                    int(ceiling(100._wp*(real(t_step - n_start)/(n_save)))), &
                    t_step, n_save
            else
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, " @ t_step = ", I0, "")', &
                    int(ceiling(100._wp*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                    (t_step - t_step_start)/t_step_save + 1, &
                    (t_step_stop - t_step_start)/t_step_save + 1, &
                    t_step
            end if
        end if

        ! Populating the grid and conservative variables
        call s_read_data_files(t_step)

        ! Populating the buffer regions of the grid and conservative variables
        if (buff_size > 0) then
            call s_populate_grid_variables_buffers()
            call s_populate_variables_buffers(bc_type, q_cons_vf)
        end if

        ! Initialize the Temperature cache.
        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwbuff)

        ! Converting the conservative variables to the primitive ones
        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwbuff)

    end subroutine s_perform_time_step

    impure subroutine s_save_data(t_step, varname, pres, c, H)

        integer, intent(inout) :: t_step
        character(LEN=name_len), intent(inout) :: varname
        real(wp), intent(inout) :: pres, c, H

        integer :: i, j, k, l

        integer :: x_beg, x_end, y_beg, y_end, z_beg, z_end

        if (output_partial_domain) then
            call s_define_output_region
            x_beg = -offset_x%beg + x_output_idx%beg
            x_end = offset_x%end + x_output_idx%end
            y_beg = -offset_y%beg + y_output_idx%beg
            y_end = offset_y%end + y_output_idx%end
            z_beg = -offset_z%beg + z_output_idx%beg
            z_end = offset_z%end + z_output_idx%end
        else
            x_beg = -offset_x%beg
            x_end = offset_x%end + m
            y_beg = -offset_y%beg
            y_end = offset_y%end + n
            z_beg = -offset_z%beg
            z_end = offset_z%end + p
        end if

        ! Opening a new formatted database file
        call s_open_formatted_database_file(t_step)

        if (sim_data .and. proc_rank == 0) then
            call s_open_intf_data_file()
            call s_open_energy_data_file()
        end if

        if (sim_data) then
            call s_write_intf_data_file(q_prim_vf)
            call s_write_energy_data_file(q_prim_vf, q_cons_vf)
        end if

        ! Adding the grid to the formatted database file
        call s_write_grid_to_formatted_database_file(t_step)

        ! Computing centered finite-difference coefficients in x-direction
        if (omega_wrt(2) .or. omega_wrt(3) .or. qm_wrt .or. schlieren_wrt) then
            call s_compute_finite_difference_coefficients(m, x_cc, &
                                                          fd_coeff_x, buff_size, &
                                                          fd_number, fd_order, offset_x)
        end if

        ! Computing centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. omega_wrt(3) .or. qm_wrt .or. (n > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(n, y_cc, &
                                                          fd_coeff_y, buff_size, &
                                                          fd_number, fd_order, offset_y)
        end if

        ! Computing centered finite-difference coefficients in z-direction
        if (omega_wrt(1) .or. omega_wrt(2) .or. qm_wrt .or. (p > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(p, z_cc, &
                                                          fd_coeff_z, buff_size, &
                                                          fd_number, fd_order, offset_z)
        end if

        ! Adding the partial densities to the formatted database file
        if ((model_eqns == 2) .or. (model_eqns == 3) .or. (model_eqns == 4)) then
            do i = 1, num_fluids
                if (alpha_rho_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:, :, :) = q_cons_vf(i)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    if (model_eqns /= 4) then
                        write (varname, '(A,I0)') 'alpha_rho', i
                    else
                        write (varname, '(A,I0)') 'rho', i
                    end if
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '

                end if
            end do
        end if

        ! Adding the density to the formatted database file
        if ((rho_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) .and. (.not. relativity)) then
            q_sf(:, :, :) = rho_sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (relativity .and. (rho_wrt .or. prim_vars_wrt)) then
            q_sf(:, :, :) = q_prim_vf(1)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (relativity .and. (rho_wrt .or. cons_vars_wrt)) then
            ! For relativistic flow, conservative and primitive densities are different
            ! Hard-coded single-component for now
            q_sf(:, :, :) = q_cons_vf(1)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'D'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        ! Adding the momentum to the formatted database file
        do i = 1, E_idx - mom_idx%beg
            if (mom_wrt(i) .or. cons_vars_wrt) then
                q_sf(:, :, :) = q_cons_vf(i + cont_idx%end)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A,I0)') 'mom', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if
        end do

        ! Adding the velocity to the formatted database file
        do i = 1, E_idx - mom_idx%beg
            if (vel_wrt(i) .or. prim_vars_wrt) then
                q_sf(:, :, :) = q_prim_vf(i + cont_idx%end)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A,I0)') 'vel', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if
        end do

        ! Adding the species' concentrations to the formatted database file
        if (chemistry) then
            do i = 1, num_species
                if (chem_wrt_Y(i) .or. prim_vars_wrt) then
                    q_sf(:, :, :) = q_prim_vf(chemxb + i - 1)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,A)') 'Y_', trim(species_names(i))
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '

                end if
            end do

            if (chem_wrt_T) then
                q_sf(:, :, :) = q_T_sf%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A)') 'T'
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end if

        ! Adding the flux limiter function to the formatted database file
        do i = 1, E_idx - mom_idx%beg
            if (flux_wrt(i)) then

                call s_derive_flux_limiter(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'flux', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        ! Adding the energy to the formatted database file
        if (E_wrt .or. cons_vars_wrt) then
            q_sf(:, :, :) = q_cons_vf(E_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'E'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the magnetic field to the formatted database file
        if (mhd .and. prim_vars_wrt) then
            do i = B_idx%beg, B_idx%end
                q_sf(:, :, :) = q_prim_vf(i)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)

                ! 1D: output By, Bz
                if (n == 0) then
                    if (i == B_idx%beg) then
                        write (varname, '(A)') 'By'
                    else
                        write (varname, '(A)') 'Bz'
                    end if
                    ! 2D/3D: output Bx, By, Bz
                else
                    if (i == B_idx%beg) then
                        write (varname, '(A)') 'Bx'
                    elseif (i == B_idx%beg + 1) then
                        write (varname, '(A)') 'By'
                    else
                        write (varname, '(A)') 'Bz'
                    end if
                end if

                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
        end if

        ! Adding the elastic shear stresses to the formatted database file
        if (elasticity) then
            do i = 1, stress_idx%end - stress_idx%beg + 1
                if (prim_vars_wrt) then
                    q_sf(:, :, :) = q_prim_vf(i - 1 + stress_idx%beg)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,I0)') 'tau', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                end if
                varname(:) = ' '
            end do
        end if

        if (hyperelasticity) then
            do i = 1, xiend - xibeg + 1
                if (prim_vars_wrt) then
                    q_sf(:, :, :) = q_prim_vf(i - 1 + xibeg)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,I0)') 'xi', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                end if
                varname(:) = ' '
            end do
        end if

        if (cont_damage) then
            q_sf(:, :, :) = q_cons_vf(damage_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'damage_state'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        ! Adding the pressure to the formatted database file
        if (pres_wrt .or. prim_vars_wrt) then
            q_sf(:, :, :) = q_prim_vf(E_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'pres'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the volume fraction(s) to the formatted database file
        if (((model_eqns == 2) .and. (bubbles_euler .neqv. .true.)) &
            .or. (model_eqns == 3) &
            ) then

            do i = 1, num_fluids - 1
                if (alpha_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:, :, :) = q_cons_vf(i + E_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,I0)') 'alpha', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '

                end if
            end do

            if (alpha_wrt(num_fluids) &
                .or. &
                (cons_vars_wrt .or. prim_vars_wrt)) then
                if (igr) then
                    do k = z_beg, z_end
                        do j = y_beg, y_end
                            do i = x_beg, x_end
                                q_sf(i, j, k) = 1._wp
                                do l = 1, num_fluids - 1
                                    q_sf(i, j, k) = q_sf(i, j, k) - q_cons_vf(E_idx + l)%sf(i, j, k)
                                end do
                            end do
                        end do
                    end do
                else
                    q_sf(:, :, :) = q_cons_vf(adv_idx%end)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                end if
                write (varname, '(A,I0)') 'alpha', num_fluids
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '

            end if

        end if

        ! Adding specific heat ratio function to formatted database file
        if (gamma_wrt &
            .or. &
            (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
            q_sf(:, :, :) = gamma_sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'gamma'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the specific heat ratio to the formatted database file
        if (heat_ratio_wrt) then

            call s_derive_specific_heat_ratio(q_sf)

            write (varname, '(A)') 'heat_ratio'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding liquid stiffness function to formatted database file
        if (pi_inf_wrt &
            .or. &
            (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
            q_sf(:, :, :) = pi_inf_sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A)') 'pi_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the liquid stiffness to the formatted database file
        if (pres_inf_wrt) then

            call s_derive_liquid_stiffness(q_sf)

            write (varname, '(A)') 'pres_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the sound speed to the formatted database file
        if (c_wrt) then
            do k = -offset_z%beg, p + offset_z%end
                do j = -offset_y%beg, n + offset_y%end
                    do i = -offset_x%beg, m + offset_x%end
                        do l = 1, adv_idx%end - E_idx
                            adv(l) = q_prim_vf(E_idx + l)%sf(i, j, k)
                        end do

                        pres = q_prim_vf(E_idx)%sf(i, j, k)

                        H = ((gamma_sf(i, j, k) + 1._wp)*pres + &
                             pi_inf_sf(i, j, k))/rho_sf(i, j, k)

                        call s_compute_speed_of_sound(pres, rho_sf(i, j, k), &
                                                      gamma_sf(i, j, k), pi_inf_sf(i, j, k), &
                                                      H, adv, 0._wp, 0._wp, c)

                        q_sf(i, j, k) = c
                    end do
                end do
            end do

            write (varname, '(A)') 'c'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the vorticity to the formatted database file
        do i = 1, 3
            if (omega_wrt(i)) then

                call s_derive_vorticity_component(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'omega', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (ib) then
            q_sf(:, :, :) = real(ib_markers%sf(-offset_x%beg:m + offset_x%end, -offset_y%beg:n + offset_y%end, -offset_z%beg:p + offset_z%end))
            varname = 'ib_markers'
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! Adding Q_M to the formatted database file
        if (p > 0 .and. qm_wrt) then
            call s_derive_qm(q_prim_vf, q_sf)

            write (varname, '(A)') 'qm'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        ! Adding numerical Schlieren function to formatted database file
        if (schlieren_wrt) then

            call s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

            write (varname, '(A)') 'schlieren'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

        end if

        ! Adding the color function to formatted database file
        if (cf_wrt) then
            q_sf(:, :, :) = q_cons_vf(c_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
            write (varname, '(A,I0)') 'color_function'
            call s_write_variable_to_formatted_database_file(varname, t_step)
            varname(:) = ' '

        end if

        ! Adding the volume fraction(s) to the formatted database file
        if (bubbles_euler) then
            do i = adv_idx%beg, adv_idx%end
                q_sf(:, :, :) = q_cons_vf(i)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A,I0)') 'alpha', i - E_idx
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
        end if

        ! Adding the bubble variables  to the formatted database file
        if (bubbles_euler) then
            !nR
            do i = 1, nb
                q_sf(:, :, :) = q_cons_vf(bub_idx%rs(i))%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A,I3.3)') 'nR', i
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do

            !nRdot
            do i = 1, nb
                q_sf(:, :, :) = q_cons_vf(bub_idx%vs(i))%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A,I3.3)') 'nV', i
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
            if ((polytropic .neqv. .true.) .and. (.not. qbmm)) then
                !nP
                do i = 1, nb
                    q_sf(:, :, :) = q_cons_vf(bub_idx%ps(i))%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,I3.3)') 'nP', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                    varname(:) = ' '
                end do

                !nM
                do i = 1, nb
                    q_sf(:, :, :) = q_cons_vf(bub_idx%ms(i))%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                    write (varname, '(A,I3.3)') 'nM', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                    varname(:) = ' '
                end do
            end if

            ! number density
            if (adv_n) then
                q_sf(:, :, :) = q_cons_vf(n_idx)%sf(x_beg:x_end, y_beg:y_end, z_beg:z_end)
                write (varname, '(A)') 'n'
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end if
        end if

        ! Adding the lagrangian subgrid variables  to the formatted database file
        if (bubbles_lagrange) then
            !! Void fraction field
            q_sf(:, :, :) = 1._wp - q_cons_vf(beta_idx)%sf( &
                            -offset_x%beg:m + offset_x%end, &
                            -offset_y%beg:n + offset_y%end, &
                            -offset_z%beg:p + offset_z%end)
            write (varname, '(A)') 'voidFraction'
            call s_write_variable_to_formatted_database_file(varname, t_step)
            varname(:) = ' '

            call s_write_lag_bubbles_results(t_step) !! Individual bubble evolution
        end if

        if (sim_data .and. proc_rank == 0) then
            call s_close_intf_data_file()
            call s_close_energy_data_file()
        end if

        ! Closing the formatted database file
        call s_close_formatted_database_file()

    end subroutine s_save_data

    impure subroutine s_initialize_modules
        ! Computation of parameters, allocation procedures, and/or any other tasks
        ! needed to properly setup the modules
        call s_initialize_global_parameters_module()
        if (bubbles_euler .and. nb > 1) then
            call s_simpson(weight, R0)
        end if
        if (bubbles_euler .and. .not. polytropic) then
            call s_initialize_nonpoly()
        end if
        if (num_procs > 1) then
            call s_initialize_mpi_proxy_module()
            call s_initialize_mpi_common_module()
        end if
        call s_initialize_boundary_common_module()
        call s_initialize_variables_conversion_module()
        call s_initialize_data_input_module()
        call s_initialize_derived_variables_module()
        call s_initialize_data_output_module()

        ! Associate pointers for serial or parallel I/O
        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if
    end subroutine s_initialize_modules

    impure subroutine s_initialize_mpi_domain
        ! Initialization of the MPI environment
        call s_mpi_initialize()

        ! Processor with rank 0 assigns default user input values prior to reading
        ! those in from the input file. Next, the user inputs are read in and their
        ! consistency is checked. The detection of any inconsistencies automatically
        ! leads to the termination of the post-process.
        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()

            print '(" Post-processing a ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s)")', m, n, p, num_procs
        end if

        ! Broadcasting the user inputs to all of the processors and performing the
        ! parallel computational domain decomposition. Neither procedure has to be
        ! carried out if the simulation is in fact not truly executed in parallel.
        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    impure subroutine s_finalize_modules
        ! Disassociate pointers for serial and parallel I/O
        s_read_data_files => null()

!        if (sim_data .and. proc_rank == 0) then
!            call s_close_intf_data_file()
!            call s_close_energy_data_file()
!        end if

        ! Deallocation procedures for the modules
        call s_finalize_data_output_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_input_module()
        call s_finalize_variables_conversion_module()
        if (num_procs > 1) then
            call s_finalize_mpi_proxy_module()
            call s_finalize_mpi_common_module()
        end if
        call s_finalize_global_parameters_module()

        ! Finalizing the MPI environment
        call s_mpi_finalize()
    end subroutine s_finalize_modules

end module m_start_up
