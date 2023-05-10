!>
!! @file m_assign_variables.f90
!! @brief Contains module m_assign_variables
module m_assign_variables

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another
    ! ==========================================================================

    implicit none


    type(scalar_field) :: alf_sum

    procedure(s_assign_patch_xxxxx_primitive_variables), &
    pointer :: s_assign_patch_primitive_variables => null() !<
    !! Depending on the multicomponent flow model, this variable is a pointer to
    !! either the subroutine s_assign_patch_mixture_primitive_variables, or the
    !! subroutine s_assign_patch_species_primitive_variables

    !> Abstract interface to the two subroutines that assign the patch primitive
    !! variables, either mixture or species, depending on the subroutine, to a
    !! particular cell in the computational domain
    abstract interface

        !> Skeleton of s_assign_patch_mixture_primitive_variables
        !!      and s_assign_patch_species_primitive_variables
        !! @param patch_id is the patch identifier
        !! @param j (x) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param k (y,th) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param l (z) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        subroutine s_assign_patch_xxxxx_primitive_variables(patch_id, j, k, l, &
                                                eta, q_prim_vf, patch_id_fp)

            import :: scalar_field, sys_size, n, m, p, wp

            integer, intent(IN) :: patch_id
            integer, intent(IN) :: j, k, l
            integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
            type(scalar_field), dimension(1:sys_size) :: q_prim_vf
            real(wp) :: eta !<

        end subroutine s_assign_patch_xxxxx_primitive_variables

    end interface

    private; public :: s_initialize_assign_variables_module, &
        s_assign_patch_primitive_variables, &
        s_assign_patch_mixture_primitive_variables, &
        s_assign_patch_species_primitive_variables, &
        s_assign_patch_species_primitive_variables_bubbles, &
        s_finialize_assign_variables_module

contains

    subroutine s_initialize_assign_variables_module()

        allocate (alf_sum%sf(0:m, 0:n, 0:p))

        ! Depending on multicomponent flow model, the appropriate procedure
        ! for assignment of the patch mixture or species primitive variables
        ! to a cell in the domain is targeted by the procedure pointer

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_assign_patch_primitive_variables => &
                s_assign_patch_mixture_primitive_variables
        else if (bubbles) then
            s_assign_patch_primitive_variables => &
                s_assign_patch_species_primitive_variables_bubbles
        else ! Volume fraction model
            s_assign_patch_primitive_variables => &
                s_assign_patch_species_primitive_variables
        end if
    
    end subroutine s_initialize_assign_variables_module

    !>  This subroutine assigns the mixture primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_mixture_primitive_variables(patch_id, j, k, l, &
                                         eta, q_prim_vf, patch_id_fp)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: eta !<

        integer, intent(IN) :: j, k, l

        real(wp) :: rho    !< density
        real(wp), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(wp) :: pres   !< pressure
        real(wp) :: gamma  !< specific heat ratio function
        real(wp) :: x_centroid, y_centroid
        real(wp) :: epsilon, beta

        integer :: smooth_patch_id
        integer :: i !< generic loop operator

        ! Assigning the mixture primitive variables of a uniform state patch
        if (patch_icpp(patch_id)%geometry /= 6) then

            ! Transferring the identity of the smoothing patch
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

            ! Density
            q_prim_vf(1)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%rho &
                + (1._wp - eta)*patch_icpp(smooth_patch_id)%rho

            ! Velocity
            do i = 1, E_idx - mom_idx%beg
                q_prim_vf(i + 1)%sf(j, k, l) = &
                    1._wp/q_prim_vf(1)%sf(j, k, l)* &
                    (eta*patch_icpp(patch_id)%rho &
                     *patch_icpp(patch_id)%vel(i) &
                     + (1._wp - eta)*patch_icpp(smooth_patch_id)%rho &
                     *patch_icpp(smooth_patch_id)%vel(i))
            end do

            ! Specific heat ratio function
            q_prim_vf(gamma_idx)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%gamma &
                + (1._wp - eta)*patch_icpp(smooth_patch_id)%gamma

            ! Pressure
            q_prim_vf(E_idx)%sf(j, k, l) = &
                1._wp/q_prim_vf(gamma_idx)%sf(j, k, l)* &
                (eta*patch_icpp(patch_id)%gamma &
                 *patch_icpp(patch_id)%pres &
                 + (1._wp - eta)*patch_icpp(smooth_patch_id)%gamma &
                 *patch_icpp(smooth_patch_id)%pres)

            ! Liquid stiffness function
            q_prim_vf(pi_inf_idx)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%pi_inf &
                + (1._wp - eta)*patch_icpp(smooth_patch_id)%pi_inf

            ! Assigning mixture primitive variables of isentropic vortex patch
        else

            ! Transferring the centroid of the isentropic vortex patch, the
            ! amplitude of its disturbance and also its domain of influence
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            epsilon = patch_icpp(patch_id)%epsilon
            beta = patch_icpp(patch_id)%beta

            ! Reference density, velocity, pressure and specific heat ratio
            ! function of the isentropic vortex patch
            rho = patch_icpp(patch_id)%rho
            vel = patch_icpp(patch_id)%vel
            pres = patch_icpp(patch_id)%pres
            gamma = patch_icpp(patch_id)%gamma

            ! Density
            q_prim_vf(1)%sf(j, k, 0) = &
                rho*(1._wp - (rho/pres)*(epsilon/(2._wp*pi))* &
                     (epsilon/(8._wp*beta*(gamma + 1._wp)*pi))* &
                     exp(2._wp*beta*(1._wp - (x_cc(j) - x_centroid)**2 &
                                   - (y_cc(k) - y_centroid)**2)) &
                     )**gamma

            ! Velocity
            q_prim_vf(2)%sf(j, k, 0) = &
                vel(1) - (y_cc(k) - y_centroid)*(epsilon/(2._wp*pi))* &
                exp(beta*(1._wp - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))
            q_prim_vf(3)%sf(j, k, 0) = &
                vel(2) + (x_cc(j) - x_centroid)*(epsilon/(2._wp*pi))* &
                exp(beta*(1._wp - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))

            ! Pressure
            q_prim_vf(4)%sf(j, k, 0) = &
                pres*(1._wp - (rho/pres)*(epsilon/(2._wp*pi))* &
                      (epsilon/(8._wp*beta*(gamma + 1._wp)*pi))* &
                      exp(2._wp*beta*(1._wp - (x_cc(j) - x_centroid)**2 &
                                    - (y_cc(k) - y_centroid)**2)) &
                      )**(gamma + 1._wp)

            ! Specific heat ratio function
            q_prim_vf(5)%sf(j, k, 0) = gamma

            ! Liquid stiffness function
            q_prim_vf(6)%sf(j, k, 0) = 0._wp

        end if

        ! Updating the patch identities bookkeeping variable
        if (1._wp - eta < (1._wp * (10._wp ** -(16)))) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_mixture_primitive_variables ! ------------

    !>  This subroutine assigns the species primitive variables. This follows
        !!  s_assign_patch_species_primitive_variables with adaptation for
        !!  ensemble-averaged bubble modeling
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_species_primitive_variables_bubbles(patch_id, j, k, l, &
                                                eta, q_prim_vf, patch_id_fp)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: eta !<
        integer, intent(IN) :: j, k, l

        ! Density, the specific heat ratio function and the liquid stiffness
        ! function, respectively, obtained from the combination of primitive
        ! variables of the current and smoothing patches
        real(wp) :: rho          !< density
        real(wp) :: gamma
        real(wp) :: lit_gamma    !< specific heat ratio
        real(wp) :: pi_inf       !< stiffness from SEOS
        real(wp) :: orig_rho
        real(wp) :: orig_gamma
        real(wp) :: orig_pi_inf
        real(wp) :: muR, muV

        real(wp), dimension(sys_size) :: orig_prim_vf !<
        real(wp), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(wp) :: pres   !< pressure
        real(wp) :: x_centroid, y_centroid
        real(wp) :: epsilon, beta

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
            !! Vector to hold original values of cell for smoothing purposes

        integer :: i  !< Generic loop iterator
        integer :: smooth_patch_id

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! Transferring original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0._wp
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1._wp - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Computing Mixture Variables from Original Primitive Variables
        ! call s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            orig_rho, &
            orig_gamma, &
            orig_pi_inf)

        ! Computing Mixture Variables of Current Patch =====================

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha(i - E_idx)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0._wp
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1._wp - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Partial densities
        if (model_eqns /= 4) then
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha_rho(i)
            end do
        end if

        ! Bubbles variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(patch_id)%r0 ! = R0(i)
                muV = patch_icpp(patch_id)%v0 ! = 0

                if (qbmm) then
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = exp((sigR**2)*2)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = exp((sigR**2)/2)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if

                    if (j == 0 .and. k == 0 .and. l == 0) then
                        print *, 'moments @ (0,0,0): ', &
                            q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l), &
                            q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l), &
                            q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l), &
                            q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l), &
                            q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l), &
                            q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l)
                    end if
                else
                    q_prim_vf(bub_idx%rs(i))%sf(j, k, l) = muR
                    q_prim_vf(bub_idx%vs(i))%sf(j, k, l) = muV
                    if (.not. polytropic) then
                        q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = patch_icpp(patch_id)%p0
                        q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = patch_icpp(patch_id)%m0
                    end if
                end if
            end do
        end if
  
        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf)

        ! ==================================================================

        ! Computing Mixture Variables of Smoothing Patch ===================

        if (model_eqns /= 4) then
            ! Partial densities
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha_rho(i)
            end do
        end if

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha(i - E_idx)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0._wp
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1._wp - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Bubbles variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(smooth_patch_id)%r0 ! = R0(i)
                muV = V0(i)*patch_icpp(smooth_patch_id)%v0 ! = 0
                if (qbmm) then
                    ! Initialize the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = exp((sigR**2)*2._wp)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if
                else
                    q_prim_vf(bub_idx%rs(i))%sf(j, k, l) = muR
                    q_prim_vf(bub_idx%vs(i))%sf(j, k, l) = muV
                    if (.not. polytropic) then
                        q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = patch_icpp(patch_id)%p0
                        q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = patch_icpp(patch_id)%m0
                    end if
                end if
            end do
        end if

        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf)

        ! ==================================================================

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1._wp - eta)*orig_prim_vf(E_idx))

        ! Volume fractions \alpha
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1._wp - eta)*orig_prim_vf(i)
        end do

        ! Elastic Shear Stress
        if (hypoelasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1._wp - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0._wp
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1._wp - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Partial densities \alpha \rho
        if (model_eqns /= 4) then
            !mixture density is an input
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%alpha_rho(i) &
                    + (1._wp - eta)*orig_prim_vf(i)
            end do
        else
            !get mixture density from pressure via Tait EOS
            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1._wp + gamma)/gamma

            ! \rho = (( p_l + pi_inf)/( p_ref + pi_inf))**(1/little_gam) * rhoref(1-alf)
            q_prim_vf(1)%sf(j, k, l) = &
                (((q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)/(pref + pi_inf))**(1/lit_gamma))* &
                rhoref*(1 - q_prim_vf(alf_idx)%sf(j, k, l))
        end if

        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1._wp - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Set streamwise velocity to hypertangent function of y
        if (vel_profile) then
            q_prim_vf(1 + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(1)*tanh(y_cc(k)) &
                + (1._wp - eta)*orig_prim_vf(1 + cont_idx%end))
        end if

        ! Smoothed bubble variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(patch_id)%r0 ! = 1*R0(i)
                muV = V0(i)*patch_icpp(patch_id)%v0 ! = 1*V0(i)
                if (qbmm) then
                    ! Initialize the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = exp((sigR**2)*2._wp)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if
                else
                    ! q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = &
                    !     (eta * R0(i)*patch_icpp(patch_id)%r0 &
                    !     + (1._wp-eta)*orig_prim_vf(bub_idx%rs(i)))
                    ! q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = &
                    !     (eta * V0(i)*patch_icpp(patch_id)%v0 &
                    !     + (1._wp-eta)*orig_prim_vf(bub_idx%vs(i)))
                    q_prim_vf(bub_idx%rs(i))%sf(j, k, l) = muR
                    q_prim_vf(bub_idx%vs(i))%sf(j, k, l) = muV

                    if (.not. polytropic) then
                        q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = patch_icpp(patch_id)%p0
                        q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = patch_icpp(patch_id)%m0
                    end if

                end if
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0._wp
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1._wp - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        if (bubbles .and. (.not. polytropic)) then
            do i = 1, nb
                if (q_prim_vf(bub_idx%ps(i))%sf(j, k, l) == dflt_real) then
                    q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = pb0(i)
                    print *, 'setting to pb0'
                end if
                if (q_prim_vf(bub_idx%ms(i))%sf(j, k, l) == dflt_real) then
                    q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = mass_v0(i)
                end if
            end do
        end if
      
        if (patch_icpp(patch_id)%geometry == 6) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            epsilon = patch_icpp(patch_id)%epsilon
            beta = patch_icpp(patch_id)%beta

            ! Reference density, velocity, pressure and specific heat ratio
            ! function of the isentropic vortex patch
            rho = patch_icpp(patch_id)%rho
            vel = patch_icpp(patch_id)%vel
            pres = patch_icpp(patch_id)%pres
            gamma = patch_icpp(patch_id)%gamma

            ! Density
            q_prim_vf(1)%sf(j, k, 0) = &
                rho*(1d0 - (rho/pres)*(epsilon/(2d0*pi))* &
                     (epsilon/(8d0*beta*(gamma + 1d0)*pi))* &
                     exp(2d0*beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                                   - (y_cc(k) - y_centroid)**2)) &
                     )**gamma

            ! Velocity
            q_prim_vf(2)%sf(j, k, 0) = &
                vel(1) - (y_cc(k) - y_centroid)*(epsilon/(2d0*pi))* &
                exp(beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))
            q_prim_vf(3)%sf(j, k, 0) = &
                vel(2) + (x_cc(j) - x_centroid)*(epsilon/(2d0*pi))* &
                exp(beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))

            ! Pressure
            q_prim_vf(4)%sf(j, k, 0) = &
                pres*(1d0 - (rho/pres)*(epsilon/(2d0*pi))* &
                      (epsilon/(8d0*beta*(gamma + 1d0)*pi))* &
                      exp(2d0*beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                                    - (y_cc(k) - y_centroid)**2)) &
                      )**(gamma + 1d0)

        end if

        ! Updating the patch identities bookkeeping variable
        if (1._wp - eta < (1._wp * (10._wp ** -(16)))) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_species_primitive_variables_bubbles ! ------------

    !>      This subroutine assigns the species primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_species_primitive_variables(patch_id, j, k, l, &
                                            eta, q_prim_vf, patch_id_fp)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(wp) :: eta !<
        integer, intent(IN) :: j, k, l

        real(wp) :: rho
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: orig_rho
        real(wp) :: orig_gamma
        real(wp) :: orig_pi_inf !<
            !! Density, the specific heat ratio function and the liquid stiffness
            !! function, respectively, obtained from the combination of primitive
            !! variables of the current and smoothing patches

        real(wp), dimension(sys_size) :: orig_prim_vf !<
        real(wp), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(wp) :: pres   !< pressure
        real(wp) :: x_centroid, y_centroid
        real(wp) :: epsilon, beta

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
        ! Vector to hold original values of cell for smoothing purposes

        integer :: smooth_patch_id

        integer :: i !< generic loop iterator

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! Transferring original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        ! Computing Mixture Variables from Original Primitive Variables
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            orig_rho, &
            orig_gamma, &
            orig_pi_inf)

        ! Computing Mixture Variables of Current Patch =====================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha_rho(i)
        end do

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha(i - E_idx)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf)

        ! ==================================================================

        ! Computing Mixture Variables of Smoothing Patch ===================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                patch_icpp(smooth_patch_id)%alpha_rho(i)
        end do

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                patch_icpp(smooth_patch_id)%alpha(i - E_idx)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf)

        ! ==================================================================

        ! Partial densities
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha_rho(i) &
                + (1._wp - eta)*orig_prim_vf(i)
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1._wp - eta)*orig_prim_vf(i)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1._wp - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Set streamwise velocity to hypertangent function of y
        if (vel_profile) then
            q_prim_vf(1 + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(1)*tanh(y_cc(k)) &
                + (1._wp - eta)*orig_prim_vf(1 + cont_idx%end))
        end if

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1._wp - eta)*orig_prim_vf(E_idx))

        ! Elastic Shear Stress
        if (hypoelasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1._wp - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        ! Set partial pressures to mixture pressure
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
            end do
        end if
        
        if (patch_icpp(patch_id)%geometry == 6) then
            x_centroid = patch_icpp(patch_id)%x_centroid
            y_centroid = patch_icpp(patch_id)%y_centroid
            epsilon = patch_icpp(patch_id)%epsilon
            beta = patch_icpp(patch_id)%beta

            ! Reference density, velocity, pressure and specific heat ratio
            ! function of the isentropic vortex patch
            rho = patch_icpp(patch_id)%rho
            vel = patch_icpp(patch_id)%vel
            pres = patch_icpp(patch_id)%pres
            gamma = patch_icpp(patch_id)%gamma

            ! Density
            q_prim_vf(1)%sf(j, k, 0) = &
                rho*(1d0 - (rho/pres)*(epsilon/(2d0*pi))* &
                     (epsilon/(8d0*beta*(gamma + 1d0)*pi))* &
                     exp(2d0*beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                                   - (y_cc(k) - y_centroid)**2)) &
                     )**gamma

            ! Velocity
            q_prim_vf(2)%sf(j, k, 0) = &
                vel(1) - (y_cc(k) - y_centroid)*(epsilon/(2d0*pi))* &
                exp(beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))
            q_prim_vf(3)%sf(j, k, 0) = &
                vel(2) + (x_cc(j) - x_centroid)*(epsilon/(2d0*pi))* &
                exp(beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                          - (y_cc(k) - y_centroid)**2))

            ! Pressure
            q_prim_vf(4)%sf(j, k, 0) = &
                pres*(1d0 - (rho/pres)*(epsilon/(2d0*pi))* &
                      (epsilon/(8d0*beta*(gamma + 1d0)*pi))* &
                      exp(2d0*beta*(1d0 - (x_cc(j) - x_centroid)**2 &
                                    - (y_cc(k) - y_centroid)**2)) &
                      )**(gamma + 1d0)

        end if

        ! Updating the patch identities bookkeeping variable
        if (1._wp - eta < (1._wp * (10._wp ** -(16)))) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_species_primitive_variables ! ------------

    subroutine s_finialize_assign_variables_module

        ! Nullifying procedure pointer to the subroutine assigning either
        ! the patch mixture or species primitive variables to a cell in the
        ! computational domain
        s_assign_patch_primitive_variables => null()

    end subroutine s_finialize_assign_variables_module

end module
