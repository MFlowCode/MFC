module m_assign_patches

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another

    ! ==========================================================================

    implicit none

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
        subroutine s_assign_patch_xxxxx_primitive_variables(patch_id, j, k, l)

            integer, intent(IN) :: patch_id
            integer, intent(IN) :: j, k, l

        end subroutine s_assign_patch_xxxxx_primitive_variables

    end interface

    ! NOTE: The abstract interface allows for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the patch primitive variables are to be assigned in
    ! a cell in the computational domain.
    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !< primitive variables
    type(scalar_field), allocatable, dimension(:) :: q_cons_vf !< conservative variables
    type(scalar_field) :: alf_sum

    real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
    real(kind(0d0)) :: epsilon, beta
    integer :: smooth_patch_id

    real(kind(0d0)) :: eta !<
    !! In the case that smoothing of patch boundaries is enabled and the boundary
    !! between two adjacent patches is to be smeared out, this variable's purpose
    !! is to act as a pseudo volume fraction to indicate the contribution of each
    !! patch toward the composition of a cell's fluid state.

    integer, allocatable, dimension(:, :, :) :: patch_id_fp !<
    !! Bookkepping variable used to track the patch identities (id) associated
    !! with each of the cells in the computational domain. Note that only one
    !! patch identity may be associated with any one cell.

    procedure(s_assign_patch_xxxxx_primitive_variables), &
        pointer :: s_assign_patch_primitive_variables => null() !<
    !! Depending on the multicomponent flow model, this variable is a pointer to
    !! either the subroutine s_assign_patch_mixture_primitive_variables, or the
    !! subroutine s_assign_patch_species_primitive_variables
    
contains

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
    subroutine s_assign_patch_mixture_primitive_variables(patch_id, j, k, l)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: rho    !< density
        real(kind(0d0)), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(kind(0d0)) :: pres   !< pressure
        real(kind(0d0)) :: gamma  !< specific heat ratio function

        integer :: i !< generic loop operator

        ! Assigning the mixture primitive variables of a uniform state patch
        if (patch_icpp(patch_id)%geometry /= 6) then

            ! Transferring the identity of the smoothing patch
            smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

            ! Density
            q_prim_vf(1)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%rho &
                + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho

            ! Velocity
            do i = 1, E_idx - mom_idx%beg
                q_prim_vf(i + 1)%sf(j, k, l) = &
                    1d0/q_prim_vf(1)%sf(j, k, l)* &
                    (eta*patch_icpp(patch_id)%rho &
                     *patch_icpp(patch_id)%vel(i) &
                     + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho &
                     *patch_icpp(smooth_patch_id)%vel(i))
            end do

            ! Specific heat ratio function
            q_prim_vf(gamma_idx)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%gamma &
                + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma

            ! Pressure
            q_prim_vf(E_idx)%sf(j, k, l) = &
                1d0/q_prim_vf(gamma_idx)%sf(j, k, l)* &
                (eta*patch_icpp(patch_id)%gamma &
                 *patch_icpp(patch_id)%pres &
                 + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma &
                 *patch_icpp(smooth_patch_id)%pres)

            ! Liquid stiffness function
            q_prim_vf(pi_inf_idx)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%pi_inf &
                + (1d0 - eta)*patch_icpp(smooth_patch_id)%pi_inf

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

            ! Specific heat ratio function
            q_prim_vf(5)%sf(j, k, 0) = gamma

            ! Liquid stiffness function
            q_prim_vf(6)%sf(j, k, 0) = 0d0

        end if

        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_mixture_primitive_variables ! ------------

    !>  This subroutine assigns the species primitive variables. This follows
        !!  s_assign_patch_species_primitive_variables with adaptation for
        !!  ensemble-averaged bubble modeling
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_species_primitive_variables_bubbles(patch_id, j, k, l)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(IN) :: j, k, l

        ! Density, the specific heat ratio function and the liquid stiffness
        ! function, respectively, obtained from the combination of primitive
        ! variables of the current and smoothing patches
        real(kind(0d0)) :: rho          !< density
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: lit_gamma    !< specific heat ratio
        real(kind(0d0)) :: pi_inf       !< stiffness from SEOS
        real(kind(0d0)) :: orig_rho
        real(kind(0d0)) :: orig_gamma
        real(kind(0d0)) :: orig_pi_inf
        real(kind(0d0)) :: muR, muV

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
            !! Vector to hold original values of cell for smoothing purposes

        integer :: i  !< Generic loop iterator

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! Transferring original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Computing Mixture Variables from Original Primitive Variables
        call s_convert_species_to_mixture_variables( &
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
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
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
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = dexp((sigR**2)*2)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = dexp((sigR**2)/2)*muR*muV
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
        call s_convert_species_to_mixture_variables( &
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
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
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
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = dexp((sigR**2)*2d0)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR*muV
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
        call s_convert_species_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf)

        ! ==================================================================

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1d0 - eta)*orig_prim_vf(E_idx))

        ! Volume fractions \alpha
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1d0 - eta)*orig_prim_vf(i)
        end do

        ! Elastic Shear Stress
        if (hypoelasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1d0 - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Partial densities \alpha \rho
        if (model_eqns /= 4) then
            !mixture density is an input
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%alpha_rho(i) &
                    + (1d0 - eta)*orig_prim_vf(i)
            end do
        else
            !get mixture density from pressure via Tait EOS
            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1.d0 + gamma)/gamma

            ! \rho = (( p_l + pi_inf)/( p_ref + pi_inf))**(1/little_gam) * rhoref(1-alf)
            q_prim_vf(1)%sf(j, k, l) = &
                (((q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)/(pref + pi_inf))**(1/lit_gamma))* &
                rhoref*(1 - q_prim_vf(alf_idx)%sf(j, k, l))
        end if

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1d0 - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Smoothed bubble variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(patch_id)%r0 ! = 1*R0(i)
                muV = V0(i)*patch_icpp(patch_id)%v0 ! = 1*V0(i)
                if (qbmm) then
                    ! Initialize the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = dexp((sigR**2)*2d0)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if
                else
                    ! q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = &
                    !     (eta * R0(i)*patch_icpp(patch_id)%r0 &
                    !     + (1d0-eta)*orig_prim_vf(bub_idx%rs(i)))
                    ! q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = &
                    !     (eta * V0(i)*patch_icpp(patch_id)%v0 &
                    !     + (1d0-eta)*orig_prim_vf(bub_idx%vs(i)))
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
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
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

        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k, l) = patch_id

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
    subroutine s_assign_patch_species_primitive_variables(patch_id, j, k, l)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: rho
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: pi_inf
        real(kind(0d0)) :: orig_rho
        real(kind(0d0)) :: orig_gamma
        real(kind(0d0)) :: orig_pi_inf !<
            !! Density, the specific heat ratio function and the liquid stiffness
            !! function, respectively, obtained from the combination of primitive
            !! variables of the current and smoothing patches

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
        ! Vector to hold original values of cell for smoothing purposes

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
                + (1d0 - eta)*orig_prim_vf(i)
        end do
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1d0 - eta)*orig_prim_vf(i)
        end do

        ! Density and the specific heat ratio and liquid stiffness functions
        call s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1d0 - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1d0 - eta)*orig_prim_vf(E_idx))

        ! Elastic Shear Stress
        if (hypoelasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1d0 - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        ! Set partial pressures to mixture pressure
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
            end do
        end if

        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_species_primitive_variables ! ------------

    !> Computation of parameters, allocation procedures, and/or
        !!              any other tasks needed to properly setup the module
    subroutine s_initialize_assign_patches_module() ! -------------------

        integer :: i !< generic loop iterator

        ! Allocating the primitive and conservative variables
        allocate (q_prim_vf(1:sys_size))
        allocate (q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (q_prim_vf(i)%sf(0:m, 0:n, 0:p))
            allocate (q_cons_vf(i)%sf(0:m, 0:n, 0:p))
        end do
        allocate (alf_sum%sf(0:m, 0:n, 0:p))

        ! Allocating the patch identities bookkeeping variable
        allocate (patch_id_fp(0:m, 0:n, 0:p))

        ! Setting default values for conservative and primitive variables so
        ! that in the case that the initial condition is wrongly laid out on
        ! the grid the simulation component will catch the problem on start-
        ! up. The conservative variables do not need to be similarly treated
        ! since they are computed directly from the primitive variables.
        do i = 1, sys_size
            q_cons_vf(i)%sf = dflt_real
            q_prim_vf(i)%sf = dflt_real
        end do

        ! Setting default values for patch identities bookkeeping variable.
        ! This is necessary to avoid any confusion in the assessment of the
        ! extent of application that the overwrite permissions give a patch
        ! when it is being applied in the domain.
        patch_id_fp = 0

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


    end subroutine s_initialize_assign_patches_module ! -----------------


    !>  Deallocation procedures for the module
    subroutine s_finalize_assign_patches_module() ! ---------------------

        integer :: i !< Generic loop iterator

        ! Dellocating the primitive and conservative variables
        do i = 1, sys_size
            deallocate (q_prim_vf(i)%sf)
            deallocate (q_cons_vf(i)%sf)
        end do

        deallocate (q_prim_vf)
        deallocate (q_cons_vf)

        ! Deallocating the patch identities bookkeeping variable
        deallocate (patch_id_fp)

        ! Nullifying procedure pointer to the subroutine assigning either
        ! the patch mixture or species primitive variables to a cell in the
        ! computational domain
        s_assign_patch_primitive_variables => null()

    end subroutine s_finalize_assign_patches_module ! -------------------

end module m_assign_patches
