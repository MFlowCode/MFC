!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive ones and vice versa. In
!!              addition, the module also contains the subroutines used to obtain
!!              the mixture variables.
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code
    ! ==========================================================================

    implicit none

    !> Abstract interface to two subroutines designed for the transfer/conversion
    !! of the mixture/species variables to the mixture variables
    abstract interface

        !> Structure of the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines
        !!  @param q_vf Conservative or primitive variables
        !!  @param i First-coordinate cell index
        !!  @param j First-coordinate cell index
        !!  @param k First-coordinate cell index
        !!  @param rho Density
        !!  @param gamma Specific heat ratio function
        !!  @param pi_inf Liquid stiffness function
        subroutine s_convert_xxxxx_to_mixture_variables(q_vf, i, j, k, &
                                                        rho, gamma, pi_inf, G)

            ! Importing the derived type scalar_field from m_derived_types.f90
            ! and global variable sys_size, from m_global_variables.f90, as
            ! the abstract interface does not inherently have access to them
            import :: scalar_field, sys_size

            type(scalar_field), dimension(sys_size), intent(IN) :: q_vf

            integer, intent(IN) :: i, j, k

            real(kind(0d0)), intent(OUT) :: rho
            real(kind(0d0)), intent(OUT) :: gamma
            real(kind(0d0)), intent(OUT) :: pi_inf

            real(kind(0d0)), optional, intent(OUT) :: G

        end subroutine s_convert_xxxxx_to_mixture_variables

    end interface

    ! NOTE: These abstract interfaces allow for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the mixture quantities are needed.

    procedure(s_convert_xxxxx_to_mixture_variables), &
        pointer :: s_convert_to_mixture_variables => null() !<
    !! Pointer referencing the subroutine s_convert_mixture_to_mixture_variables
    !! or s_convert_species_to_mixture_variables, based on model equations choice

contains

    !>  This subroutine is designed for the gamma/pi_inf model
        !!      and provided a set of either conservative or primitive
        !!      variables, transfers the density, specific heat ratio
        !!      function and the liquid stiffness function from q_vf to
        !!      rho, gamma and pi_inf.
        !! @param q_vf conservative or primitive variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        !! @param rho density
        !! @param gamma  specific heat ratio function
        !! @param pi_inf liquid stiffness
    subroutine s_convert_mixture_to_mixture_variables(q_vf, i, j, k, &
                                                      rho, gamma, pi_inf, G)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_vf
        integer, intent(IN) :: i, j, k

        real(kind(0d0)), intent(OUT) :: rho
        real(kind(0d0)), intent(OUT) :: gamma
        real(kind(0d0)), intent(OUT) :: pi_inf

        real(kind(0d0)), optional, intent(OUT) :: G

        ! Transfering the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        rho = q_vf(1)%sf(i, j, k)
        gamma = q_vf(gamma_idx)%sf(i, j, k)
        pi_inf = q_vf(pi_inf_idx)%sf(i, j, k)

    end subroutine s_convert_mixture_to_mixture_variables ! ----------------

    !>  This procedure is used alongside with the gamma/pi_inf
        !!      model to transfer the density, the specific heat ratio
        !!      function and liquid stiffness function from the vector
        !!      of conservative or primitive variables to their scalar
        !!      counterparts. Specifially designed for when subgrid bubbles
        !!      must be included.
        !! @param qK_vf primitive variables
        !! @param rho_K density
        !! @param gamma_K specific heat ratio
        !! @param pi_inf_K liquid stiffness
        !! @param j Cell index
        !! @param k Cell index
        !! @param l Cell index
    subroutine s_convert_species_to_mixture_variables_bubbles(qK_vf, &
                                                              j, k, l, &
                                                              rho_K, gamma_K, pi_inf_K, &
                                                              G)

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K !<
            !! Partial densities and volume fractions

        real(kind(0d0)), optional, intent(OUT) :: G

        integer, intent(IN) :: j, k, l
        integer :: i

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        ! alpha_rho_K(1) = qK_vf(i)%sf(i,j,k)
        ! alpha_K(1)     = qK_vf(E_idx+i)%sf(i,j,k)

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively

        if (model_eqns == 4) then
            rho_K = qK_vf(1)%sf(j, k, l)
            gamma_K = fluid_pp(1)%gamma    !qK_vf(gamma_idx)%sf(i,j,k)
            pi_inf_K = fluid_pp(1)%pi_inf   !qK_vf(pi_inf_idx)%sf(i,j,k)
        else if ((model_eqns == 2) .and. bubbles) then
            rho_k = 0d0; gamma_k = 0d0; pi_inf_k = 0d0

            if (mpp_lim .and. (num_fluids > 2)) then
                do i = 1, num_fluids
                    rho_k = rho_k + qK_vf(i)%sf(j, k, l)
                    gamma_k = gamma_k + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf_k = pi_inf_k + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                end do
            else if (num_fluids == 2) then
                rho_K = qK_vf(1)%sf(j, k, l)
                gamma_K = fluid_pp(1)%gamma
                pi_inf_K = fluid_pp(1)%pi_inf
            else if (num_fluids > 2) then
                !TODO: This may need fixing for hypo + bubbles
                do i = 1, num_fluids - 1 !leave out bubble part of mixture
                    rho_k = rho_k + qK_vf(i)%sf(j, k, l)
                    gamma_k = gamma_k + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf_k = pi_inf_k + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                end do
                !rho_K    = qK_vf(1)%sf(j,k,l)
                !gamma_K  = fluid_pp(1)%gamma
                !pi_inf_K = fluid_pp(1)%pi_inf
            else
                rho_K = qK_vf(1)%sf(j, k, l)
                gamma_K = fluid_pp(1)%gamma
                pi_inf_K = fluid_pp(1)%pi_inf
            end if
        end if

    end subroutine s_convert_species_to_mixture_variables_bubbles ! ----------------

    !>  This subroutine is designed for the volume fraction model
        !!              and provided a set of either conservative or primitive
        !!              variables, computes the density, the specific heat ratio
        !!              function and the liquid stiffness function from q_vf and
        !!              stores the results into rho, gamma and pi_inf.
        !! @param q_vf primitive variables
        !! @param rho density
        !! @param gamma specific heat ratio
        !! @param pi_inf liquid stiffness
        !! @param j Cell index
        !! @param k Cell index
        !! @param l Cell index
    subroutine s_convert_species_to_mixture_variables(q_vf, j, k, l, &
                                                      rho, gamma, pi_inf, G)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_vf

        integer, intent(IN) :: j, k, l !<
        ! Indices of the cell for which to compute the mixture variables

        real(kind(0d0)), intent(OUT) :: rho
        real(kind(0d0)), intent(OUT) :: gamma
        real(kind(0d0)), intent(OUT) :: pi_inf

        real(kind(0d0)), optional, intent(OUT) :: G

        integer :: i !< Generic loop iterator

        ! Computing the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        rho = 0d0; gamma = 0d0; pi_inf = 0d0

        do i = 1, num_fluids
            rho = rho + q_vf(i)%sf(j, k, l)
            gamma = gamma + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
            pi_inf = pi_inf + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
        end do

        if (present(G)) then
            G = 0d0
            do i = 1, num_fluids
                G = G + q_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%G
            end do
            G = max(0d0, G)
        end if

    end subroutine s_convert_species_to_mixture_variables ! ----------------

    !> Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_variables_conversion_module() ! -------------------

        ! Depending on the model selection for the equations of motion, the
        ! appropriate procedure for the conversion to the mixture variables
        ! is targeted by the procedure pointer

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables

        else if (bubbles) then
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else
            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

    end subroutine s_initialize_variables_conversion_module ! -----------------

    !> Converts the conservative variables to the primitive ones
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                             q_prim_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf

        ! Density, specific heat ratio function, liquid stiffness function
        ! and dynamic pressure, as defined in the incompressible flow sense,
        ! respectively
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: nbub, nR3, vftmp
        real(kind(0d0)), dimension(nb) :: nRtmp

        real(kind(0d0)) :: G

        ! Generic loop iterators
        integer :: i, j, k, l, q

        ! Converting the conservative variables to the primitive variables
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Obtaining the density, specific heat ratio function
                    ! and the liquid stiffness function, respectively
                    if (model_eqns /= 4) then
                        if (hypoelasticity) then
                            call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                rho, gamma, pi_inf, G)
                        else
                            call s_convert_to_mixture_variables(q_cons_vf, j, k, l, &
                                                                rho, gamma, pi_inf)
                        end if
                    end if

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, cont_idx%end
                        q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed
                    ! iteratively by cycling through the momentum equations
                    dyn_pres = 0d0

                    ! Computing velocity and dynamic pressure from momenta
                    do i = mom_idx%beg, mom_idx%end
                        if (model_eqns /= 4) then
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/rho
                            dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)* &
                                       q_prim_vf(i)%sf(j, k, l)/2d0
                        else
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) &
                                                       /q_cons_vf(1)%sf(j, k, l)
                        end if
                    end do

                    if ((model_eqns /= 4) .and. (bubbles .neqv. .true.)) then
                        ! Computing the pressure from the energy
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            (q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres - pi_inf)/gamma
                    else if ((model_eqns /= 4) .and. bubbles) then
                        print *, 'getting model_eqns 2 with bubbles. cons to prim'
                        ! p = ( E/(1-alf) - 0.5 rho u u/(1-alf) - pi_inf_k )/gamma_k
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            ((q_cons_vf(E_idx)%sf(j, k, l) - dyn_pres)/(1.d0 - q_cons_vf(alf_idx)%sf(j, k, l)) &
                             - pi_inf)/gamma
                    else
                        ! Tait EOS
                        ! p = (pl0 + pi_infty)(rho/(rho_l0(1-alf)))^gamma - pi_infty
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            (pref + fluid_pp(1)%pi_inf)* &
                            ( &
                            q_prim_vf(1)%sf(j, k, l)/ &
                            (rhoref*(1 - q_prim_vf(alf_idx)%sf(j, k, l))) &
                            )**(1/fluid_pp(1)%gamma + 1) - fluid_pp(1)%pi_inf
                    end if

                    ! Set partial pressures to mixture pressure
                    if (model_eqns == 3) then
                        do i = internalEnergies_idx%beg, internalEnergies_idx%end
                            q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
                        end do
                    end if

                    ! Transfer the advection equation(s) variable(s)
                    do i = adv_idx%beg, adv_idx%end
                        q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do

                    if (bubbles) then
                        ! From cons: ntmp = DSQRT( (4.d0*pi/3.d0)*nR3/vftmp )
                        do i = 1, nb
                            nRtmp(i) = q_cons_vf(bub_idx%rs(i))%sf(j, k, l)
                        end do
                        vftmp = q_cons_vf(alf_idx)%sf(j, k, l)
                        nR3 = 0d0
                        do q = 1, nb
                            nR3 = nR3 + weight(q)*(nRtmp(q)**3d0)
                        end do
                        nbub = DSQRT((4.d0*pi/3.d0)*nR3/vftmp)
                        do i = bub_idx%beg, bub_idx%end
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/nbub
                        end do
                    end if

                    if (hypoelasticity) then
                        do i = stress_idx%beg, stress_idx%end
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/rho
                            ! subtracting elastic contribution for pressure calculation
                            if (G > 1000) then !TODO: Change to >0 if stable
                                q_prim_vf(E_idx)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l) - &
                                                               ((q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G))/gamma
                                ! extra terms in 2 and 3D
                                if ((i == stress_idx%beg + 1) .or. &
                                    (i == stress_idx%beg + 3) .or. &
                                    (i == stress_idx%beg + 4)) then
                                    q_prim_vf(E_idx)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l) - &
                                                                   ((q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G))/gamma
                                end if
                            end if
                        end do
                    end if

                end do
            end do
        end do

    end subroutine s_convert_conservative_to_primitive_variables ! ---------

    !> Converts the primitive variables to the conservative ones.
        !!  Used when initializing patches.
        !! @param q_cons_vf Conservative variables
        !! @param q_prim_vf Primitive variables
    subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, &
                                                             q_cons_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function
        ! and dynamic pressure, as defined in the incompressible flow sense,
        ! respectively
        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: dyn_pres
        real(kind(0d0)) :: nbub, R3, vftmp
        real(kind(0d0)), dimension(nb) :: Rtmp

        real(kind(0d0)) :: G

        integer :: i, j, k, l, q !< Generic loop iterators

        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    ! Obtaining the density, specific heat ratio function
                    ! and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                                                        rho, gamma, pi_inf)

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, cont_idx%end
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed
                    ! iteratively by cycling through the velocity equations
                    dyn_pres = 0d0

                    ! Computing momenta and dynamic pressure from velocity
                    do i = mom_idx%beg, mom_idx%end
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)* &
                                   q_prim_vf(i)%sf(j, k, l)/2d0
                    end do

                    ! Computing the energy from the pressure
                    if ((model_eqns /= 4) .and. (bubbles .neqv. .true.)) then
                        ! E = Gamma*P + \rho u u /2 + \pi_inf
                        q_cons_vf(E_idx)%sf(j, k, l) = &
                            gamma*q_prim_vf(E_idx)%sf(j, k, l) + dyn_pres + pi_inf
                    else if ((model_eqns /= 4) .and. (bubbles)) then
                        ! \tilde{E} = dyn_pres + (1-\alf)(\Gamma p_l + \Pi_inf)
                        q_cons_vf(E_idx)%sf(j, k, l) = dyn_pres + &
                                                       (1.d0 - q_prim_vf(alf_idx)%sf(j, k, l))* &
                                                       (gamma*q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)
                    else
                        !Tait EOS, no conserved energy variable
                        q_cons_vf(E_idx)%sf(j, k, l) = 0.
                    end if

                    ! Computing the internal energies from the pressure and continuities
                    if (model_eqns == 3) then
                        do i = internalEnergies_idx%beg, internalEnergies_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i - adv_idx%end)%sf(j, k, l)* &
                                                       fluid_pp(i - adv_idx%end)%gamma* &
                                                       q_prim_vf(E_idx)%sf(j, k, l) + &
                                                       fluid_pp(i - adv_idx%end)%pi_inf
                        end do
                    end if

                    ! Transferring the advection equation(s) variable(s)
                    do i = adv_idx%beg, adv_idx%end
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    if (bubbles) then
                        ! From prim: Compute nbub = (3/4pi) * \alpha / \bar{R^3}
                        do i = 1, nb
                            Rtmp(i) = q_prim_vf(bub_idx%rs(i))%sf(j, k, l)
                        end do
                        !call s_comp_n_from_prim_cpu(q_prim_vf(alf_idx)%sf(j, k, l), Rtmp, nbub)
                        vftmp = q_prim_vf(alf_idx)%sf(j, k, l)
                        R3 = 0d0
                        do q = 1, nb
                            R3 = R3 + weight(q)*(Rtmp(q)**3d0)
                        end do
                        nbub = (3.d0/(4.d0*pi))*vftmp/R3
                        if (j == 0 .and. k == 0 .and. l == 0) print *, 'In convert, nbub:', nbub
                        do i = bub_idx%beg, bub_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)*nbub
                            ! IF( j==0 .and. k==0 .and. l==0) THEN
                            !     PRINT*, 'nmom', i, q_cons_vf(i)%sf(j,k,l)
                            ! END IF
                        end do
                    end if

                    if (hypoelasticity) then
                        do i = stress_idx%beg, stress_idx%end
                            q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                            ! adding elastic contribution
                            if (G > 1000) then
                                q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                               (q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G)
                                ! extra terms in 2 and 3D
                                if ((i == stress_idx%beg + 1) .or. &
                                    (i == stress_idx%beg + 3) .or. &
                                    (i == stress_idx%beg + 4)) then
                                    q_cons_vf(E_idx)%sf(j, k, l) = q_cons_vf(E_idx)%sf(j, k, l) + &
                                                                   (q_prim_vf(i)%sf(j, k, l)**2d0)/(4d0*G)
                                end if
                            end if
                        end do
                    end if
                end do
            end do
        end do

    end subroutine s_convert_primitive_to_conservative_variables ! ---------

    !> Deallocation procedures for the module
    subroutine s_finalize_variables_conversion_module() ! ----------------

        ! Nullifying the procedure pointer to the subroutine transfering/
        ! computing the mixture/species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module ! --------------

end module m_variables_conversion
