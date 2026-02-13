!>
!!
!! module m_assign_variables

#:include 'case.fpp'
#:include 'macros.fpp'

module m_assign_variables

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from

    use m_helper_basic          !< Functions to compare floating point numbers

    use m_thermochem, only: num_species, gas_constant, get_mixture_molecular_weight

    implicit none

    public :: s_perturb_primitive

    type(scalar_field) :: alf_sum

    procedure(s_assign_patch_xxxxx_primitive_variables), &
        pointer :: s_assign_patch_primitive_variables => null() !<
    !! the multicomponent flow model, this variable is a pointer to
    !! subroutine s_assign_patch_mixture_primitive_variables, or the
    !!

    !> Abstract interface to the two subroutines that assign the patch primitive
    !! mixture or species, depending on the subroutine, to a
    !! in the computational domain
    abstract interface

        !> Skeleton of s_assign_patch_mixture_primitive_variables
        !!
        !! is the patch identifier
        !! (x) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! (y,th) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! (z) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! pseudo volume fraction
        !! Primitive variables
        !! Array to track patch ids
        subroutine s_assign_patch_xxxxx_primitive_variables(patch_id, j, k, l, &
                                                            eta, q_prim_vf, patch_id_fp)

            import :: scalar_field, sys_size, n, m, p, wp

            integer, intent(in) :: patch_id
            integer, intent(in) :: j, k, l
            real(wp), intent(in) :: eta
            type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_MIXED_PRECISION
            integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
            integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif

        end subroutine s_assign_patch_xxxxx_primitive_variables

    end interface

    private; 
    public :: s_initialize_assign_variables_module, &
              s_assign_patch_primitive_variables, &
              s_assign_patch_mixture_primitive_variables, &
              s_assign_patch_species_primitive_variables, &
              s_finalize_assign_variables_module

contains

    impure subroutine s_initialize_assign_variables_module

        if (.not. igr) then
            allocate (alf_sum%sf(0:m, 0:n, 0:p))
        end if

        ! on multicomponent flow model, the appropriate procedure
        ! assignment of the patch mixture or species primitive variables
        ! a cell in the domain is targeted by the procedure pointer

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_assign_patch_primitive_variables => &
                s_assign_patch_mixture_primitive_variables
        else ! Volume fraction model
            s_assign_patch_primitive_variables => &
                s_assign_patch_species_primitive_variables
        end if

    end subroutine s_initialize_assign_variables_module

    !>  This subroutine assigns the mixture primitive variables
        !! patch designated by the patch_id, to the cell that
        !! by the indexes (j,k,l). In addition, the
        !! the patch identities in the entire
        !! updated with the new assignment. Note that if
        !! of the patch's boundaries is employed, the
        !! variables in the cell will be a type of
        !! the current patch's primitive variables
        !! of the smoothing patch. The specific details
        !! combination may be found in Shyue's work (1998).
        !! the patch identifier
        !! the x-dir node index
        !! the y-dir node index
        !! the z-dir node index
        !! pseudo volume fraction
        !! Primitive variables
        !! Array to track patch ids
    subroutine s_assign_patch_mixture_primitive_variables(patch_id, j, k, l, &
                                                          eta, q_prim_vf, patch_id_fp)
        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: patch_id
        integer, intent(in) :: j, k, l
        real(wp), intent(in) :: eta
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif

        real(wp) :: Ys(1:num_species)

        integer :: smooth_patch_id
        integer :: i !< generic loop operator

        ! the mixture primitive variables of a uniform state patch

        ! the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        !
        q_prim_vf(1)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%rho &
            + (1._wp - eta)*patch_icpp(smooth_patch_id)%rho

        !
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + 1)%sf(j, k, l) = &
                1._wp/q_prim_vf(1)%sf(j, k, l)* &
                (eta*patch_icpp(patch_id)%rho &
                 *patch_icpp(patch_id)%vel(i) &
                 + (1._wp - eta)*patch_icpp(smooth_patch_id)%rho &
                 *patch_icpp(smooth_patch_id)%vel(i))
        end do

        ! heat ratio function
        q_prim_vf(gamma_idx)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%gamma &
            + (1._wp - eta)*patch_icpp(smooth_patch_id)%gamma

        !
        q_prim_vf(E_idx)%sf(j, k, l) = &
            1._wp/q_prim_vf(gamma_idx)%sf(j, k, l)* &
            (eta*patch_icpp(patch_id)%gamma &
             *patch_icpp(patch_id)%pres &
             + (1._wp - eta)*patch_icpp(smooth_patch_id)%gamma &
             *patch_icpp(smooth_patch_id)%pres)

        ! stiffness function
        q_prim_vf(pi_inf_idx)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%pi_inf &
            + (1._wp - eta)*patch_icpp(smooth_patch_id)%pi_inf

        ! Concentrations
        if (chemistry) then
            block
                real(wp) :: sum, term

                ! the species concentrations
                sum = 0._wp
                do i = 1, num_species
                    term = &
                        eta*patch_icpp(patch_id)%Y(i) &
                        + (1._wp - eta)*patch_icpp(smooth_patch_id)%Y(i)
                    q_prim_vf(chemxb + i - 1)%sf(j, k, l) = term
                    sum = sum + term
                end do

                sum = max(sum, verysmall)

                ! the species concentrations
                do i = 1, num_species
                    q_prim_vf(chemxb + i - 1)%sf(j, k, l) = &
                        q_prim_vf(chemxb + i - 1)%sf(j, k, l)/sum
                    Ys(i) = q_prim_vf(chemxb + i - 1)%sf(j, k, l)
                end do
            end block
        end if

        ! the patch identities bookkeeping variable
        if (1._wp - eta < 1.e-16_wp) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_mixture_primitive_variables

    !Stable perturbation in pressure (Ando)
    !! the x-dir node index
    !! the y-dir node index
    !! the z-dir node index
    !! Primitive variables
    subroutine s_perturb_primitive(j, k, l, q_prim_vf)

        integer, intent(in) :: j, k, l
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        integer :: i
        real(wp) :: pres_mag, loc, n_tait, B_tait, p0
        real(wp) :: R3bar, n0, ratio, nH, vfH, velH, rhoH, deno

        p0 = 101325._wp
        pres_mag = 1.e-1_wp
        loc = x_cc(177)
        n_tait = gs_min(1)
        B_tait = ps_inf(1)

        if (j < 177) then
            q_prim_vf(E_idx)%sf(j, k, l) = 0.5_wp*q_prim_vf(E_idx)%sf(j, k, l)
        end if

        if (qbmm) then
            do i = 1, nb
                q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l) = q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l)*((p0 - bub_pp%pv)/(q_prim_vf(E_idx)%sf(j, k, l)*p0 - bub_pp%pv))**(1._wp/3._wp)
            end do
        end if

        R3bar = 0._wp

        if (qbmm) then
            do i = 1, nb
                R3bar = R3bar + weight(i)*0.5_wp*(q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l))**3._wp
                R3bar = R3bar + weight(i)*0.5_wp*(q_prim_vf(bubxb + 1 + (i - 1)*nmom)%sf(j, k, l))**3._wp
            end do
        else
            do i = 1, nb
                if (polytropic) then
                    R3bar = R3bar + weight(i)*(q_prim_vf(bubxb + (i - 1)*2)%sf(j, k, l))**3._wp
                else
                    R3bar = R3bar + weight(i)*(q_prim_vf(bubxb + (i - 1)*4)%sf(j, k, l))**3._wp
                end if
            end do
        end if

        n0 = 3._wp*q_prim_vf(alf_idx)%sf(j, k, l)/(4._wp*pi*R3bar)

        ratio = ((1._wp + B_tait)/(q_prim_vf(E_idx)%sf(j, k, l) + B_tait))**(1._wp/n_tait)

        nH = n0/((1._wp - q_prim_vf(alf_idx)%sf(j, k, l))*ratio + (4._wp*pi/3._wp)*n0*R3bar)
        vfH = (4._wp*pi/3._wp)*nH*R3bar
        rhoH = (1._wp - vfH)/ratio
        deno = 1._wp - (1._wp - q_prim_vf(alf_idx)%sf(j, k, l))/rhoH

        if (f_approx_equal(deno, 0._wp)) then
            velH = 0._wp
        else
            velH = (q_prim_vf(E_idx)%sf(j, k, l) - 1._wp)/(1._wp - q_prim_vf(alf_idx)%sf(j, k, l))/deno
            velH = sqrt(velH)
            velH = velH*deno
        end if

        do i = cont_idx%beg, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = rhoH
        end do

        do i = mom_idx%beg, mom_idx%end
            q_prim_vf(i)%sf(j, k, l) = velH
        end do

        q_prim_vf(alf_idx)%sf(j, k, l) = vfH

    end subroutine s_perturb_primitive

    !>  This subroutine assigns the species primitive variables. This follows
        !! adaptation for
        !! modeling
        !! the patch identifier
        !! the x-dir node index
        !! the y-dir node index
        !! the z-dir node index
        !! pseudo volume fraction
        !! Primitive variables
        !! Array to track patch ids
    impure subroutine s_assign_patch_species_primitive_variables(patch_id, j, k, l, &
                                                                 eta, q_prim_vf, patch_id_fp)
        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: patch_id
        integer, intent(in) :: j, k, l
        real(wp), intent(in) :: eta
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! the specific heat ratio function and the liquid stiffness
        ! respectively, obtained from the combination of primitive
        ! of the current and smoothing patches
        real(wp) :: rho         !< density
        real(wp) :: gamma
        real(wp) :: lit_gamma   !< specific heat ratio
        real(wp) :: pi_inf      !< stiffness from SEOS
        real(wp) :: qv          !< reference energy from SEOS
        real(wp) :: orig_rho
        real(wp) :: orig_gamma
        real(wp) :: orig_pi_inf
        real(wp) :: orig_qv
        real(wp) :: muR, muV
        real(wp) :: R3bar
        real(wp) :: rcoord, theta, phi, xi_sph
        real(wp), dimension(3) :: xi_cart

        real(wp) :: Ys(1:num_species)

        real(stp), dimension(sys_size) :: orig_prim_vf !<
            !! hold original values of cell for smoothing purposes

        integer :: i  !< Generic loop iterator
        integer :: smooth_patch_id

        ! the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        if (mpp_lim .and. bubbles_euler) then
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

        ! Mixture Variables from Original Primitive Variables
        ! s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            orig_rho, &
            orig_gamma, &
            orig_pi_inf, orig_qv)

        ! Mixture Variables of Current Patch

        if (.not. igr .or. num_fluids > 1) then
            ! fraction(s)
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha(i - E_idx)
            end do
        end if

        if (mpp_lim .and. bubbles_euler) then
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

        ! densities
        if (model_eqns /= 4) then
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha_rho(i)
            end do
        end if

        ! and the specific heat ratio and liquid stiffness functions
        ! s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf, &
            patch_icpp(patch_id)%qv)

        ! Mixture Variables of Smoothing Patch

        if (model_eqns /= 4) then
            ! densities
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha_rho(i)
            end do
        end if

        if (.not. igr .or. num_fluids > 1) then
            ! fraction(s)
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha(i - E_idx)
            end do
        end if

        if (mpp_lim .and. bubbles_euler) then
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

        ! euler variables
        if (bubbles_euler) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(smooth_patch_id)%r0/R0ref
                muV = patch_icpp(smooth_patch_id)%v0
                if (qbmm) then
                    ! the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + (sigR*R0ref)**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*(sigR*R0ref)*(sigV*sqrt(p0ref/rho0ref))
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + (sigV*sqrt(p0ref/rho0ref))**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = exp((sigR**2)*2._wp)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + (sigV*sqrt(p0ref/rho0ref))**2
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

            if (adv_n) then
                ! number density
                R3bar = 0._wp
                do i = 1, nb
                    R3bar = R3bar + weight(i)*(q_prim_vf(bub_idx%rs(i))%sf(j, k, l))**3._wp
                end do
                q_prim_vf(n_idx)%sf(j, k, l) = 3*q_prim_vf(alf_idx)%sf(j, k, l)/(4*pi*R3bar)
            end if
        end if

        ! and the specific heat ratio and liquid stiffness functions
        ! s_convert_species_to_mixture_variables( &
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf, &
            patch_icpp(smooth_patch_id)%qv)

        !
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1._wp - eta)*orig_prim_vf(E_idx))

        if (.not. igr .or. num_fluids > 1) then
            ! fractions \alpha
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                    + (1._wp - eta)*orig_prim_vf(i)
            end do
        end if

        if (mhd) then
            if (n == 0) then ! 1D: By, Bz
                q_prim_vf(B_idx%beg)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%By &
                    + (1._wp - eta)*orig_prim_vf(B_idx%beg)
                q_prim_vf(B_idx%beg + 1)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%Bz &
                    + (1._wp - eta)*orig_prim_vf(B_idx%beg + 1)
            else ! 2D/3D: Bx, By, Bz
                q_prim_vf(B_idx%beg)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%Bx &
                    + (1._wp - eta)*orig_prim_vf(B_idx%beg)
                q_prim_vf(B_idx%beg + 1)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%By &
                    + (1._wp - eta)*orig_prim_vf(B_idx%beg + 1)
                q_prim_vf(B_idx%beg + 2)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%Bz &
                    + (1._wp - eta)*orig_prim_vf(B_idx%beg + 2)
            end if
        end if

        ! Shear Stress
        if (elasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1._wp - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        ! Shear Stress
        if (hyperelasticity) then

            if (pre_stress) then ! pre stressed initial condition in spatial domain
                rcoord = sqrt((x_cc(j)**2 + y_cc(k)**2 + z_cc(l)**2))
                theta = atan2(y_cc(k), x_cc(j))
                phi = atan2(sqrt(x_cc(j)**2 + y_cc(k)**2), z_cc(l))
                !spherical coord, assuming Rmax=1
                xi_sph = (rcoord**3 - R0ref**3 + 1._wp)**(1._wp/3._wp)
                xi_cart(1) = xi_sph*sin(phi)*cos(theta)
                xi_cart(2) = xi_sph*sin(phi)*sin(theta)
                xi_cart(3) = xi_sph*cos(phi)
            else
                xi_cart(1) = x_cc(j)
                xi_cart(2) = y_cc(k)
                xi_cart(3) = z_cc(l)
            end if

            ! the reference map to the q_prim vector field
            do i = 1, num_dims
                q_prim_vf(i + xibeg - 1)%sf(j, k, l) = eta*xi_cart(i) + &
                                                       (1._wp - eta)*orig_prim_vf(i + xibeg - 1)
            end do
        end if

        if (mpp_lim .and. bubbles_euler) then
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

        ! densities \alpha \rho
        if (model_eqns /= 4) then
            !mixture density is an input
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%alpha_rho(i) &
                    + (1._wp - eta)*orig_prim_vf(i)
            end do
        else
            !get mixture density from pressure via Tait EOS
            pi_inf = pi_infs(1)
            gamma = gammas(1)
            lit_gamma = gs_min(1)

            ! = (( p_l + pi_inf)/( p_ref + pi_inf))**(1/little_gam) * rhoref(1-alf)
            q_prim_vf(1)%sf(j, k, l) = &
                (((q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)/(pref + pi_inf))**(1/lit_gamma))* &
                rhoref*(1 - q_prim_vf(alf_idx)%sf(j, k, l))
        end if

        ! and the specific heat ratio and liquid stiffness functions
        ! s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                                            rho, gamma, pi_inf, qv)

        !
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1._wp - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Concentrations
        if (chemistry) then
            block
                real(wp) :: sum, term

                ! the species concentrations
                sum = 0._wp
                do i = 1, num_species
                    term = &
                        eta*patch_icpp(patch_id)%Y(i) &
                        + (1._wp - eta)*patch_icpp(smooth_patch_id)%Y(i)
                    q_prim_vf(chemxb + i - 1)%sf(j, k, l) = term
                    sum = sum + term
                end do

                if (sum < verysmall) then
                    sum = 1._wp
                end if

                ! the species concentrations
                do i = 1, num_species
                    q_prim_vf(chemxb + i - 1)%sf(j, k, l) = &
                        q_prim_vf(chemxb + i - 1)%sf(j, k, l)/sum
                    Ys(i) = q_prim_vf(chemxb + i - 1)%sf(j, k, l)
                end do
            end block
        end if

        ! streamwise velocity to hyperbolic tangent function of y
        if (mixlayer_vel_profile) then
            q_prim_vf(1 + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(1)*tanh(y_cc(k)*mixlayer_vel_coef) &
                 + (1._wp - eta)*orig_prim_vf(1 + cont_idx%end))
        end if

        ! partial pressures to mixture pressure for the 6-eqn model
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
            end do
        end if

        ! bubble variables
        if (bubbles_euler) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(patch_id)%r0/R0ref
                muV = patch_icpp(patch_id)%v0
                if (qbmm) then
                    ! the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + (sigR*R0ref)**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*(sigR*R0ref)*(sigV*sqrt(p0ref/rho0ref))
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + (sigV*sqrt(p0ref/rho0ref))**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1._wp
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = exp((sigR**2)*2._wp)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = exp((sigR**2)/2._wp)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + (sigV*sqrt(p0ref/rho0ref))**2
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

            if (adv_n) then
                ! number density
                R3bar = 0._wp
                do i = 1, nb
                    R3bar = R3bar + weight(i)*(q_prim_vf(bub_idx%rs(i))%sf(j, k, l))**3._wp
                end do
                q_prim_vf(n_idx)%sf(j, k, l) = 3*q_prim_vf(alf_idx)%sf(j, k, l)/(4*pi*R3bar)
            end if
        end if

        if (mpp_lim .and. bubbles_euler) then
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

        if (bubbles_euler .and. (.not. polytropic) .and. (.not. qbmm)) then
            do i = 1, nb
                if (f_is_default(real(q_prim_vf(bub_idx%ps(i))%sf(j, k, l), kind=wp))) then
                    q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = pb0(i)
                    ! *, 'setting to pb0'
                end if
                if (f_is_default(real(q_prim_vf(bub_idx%ms(i))%sf(j, k, l), kind=wp))) then
                    q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = mass_v0(i)
                end if
            end do
        end if

        if (surface_tension) then
            q_prim_vf(c_idx)%sf(j, k, l) = eta*patch_icpp(patch_id)%cf_val + &
                                           (1._wp - eta)*orig_prim_vf(c_idx)
        end if

        ! the patch identities bookkeeping variable
        if (1._wp - eta < 1.e-16_wp) patch_id_fp(j, k, l) = patch_id

        ! (j == 1) then
        !     *, (q_prim_vf(bub_idx%rs(i))%sf(j, k, l), i = 1, nb)
        !     *, (q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l), i = 1, nb)
        !     *, (R0(i), i = 1, nb)
        !     *, patch_icpp(patch_id)%r0
        !     *, (bub_idx%rs(i), i = 1, nb)
        !     *, (bub_idx%fullmom(i, 1, 0), i = 1, nb)
        ! if

    end subroutine s_assign_patch_species_primitive_variables

    impure subroutine s_finalize_assign_variables_module

        ! procedure pointer to the subroutine assigning either
        ! patch mixture or species primitive variables to a cell in the
        ! domain
        s_assign_patch_primitive_variables => null()

    end subroutine s_finalize_assign_variables_module

end module m_assign_variables
