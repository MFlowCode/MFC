!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features a database of subroutines that allow for the
!!              conversion of state variables from one type into another. At this
!!              time, the state variables type conversions below are available:
!!                             1) Mixture        => Mixture
!!                             2) Species        => Mixture
!!                             3) Conservative   => Primitive
!!                             5) Conservative   => Flux
!!                             6) Primitive      => Conservative
!!                             8) Primitive      => Flux
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use openacc

    use nvtx
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_variables_conversion_module, &
         s_convert_to_mixture_variables, &
         s_convert_mixture_to_mixture_variables, &
         s_convert_species_to_mixture_variables_bubbles, &
         s_convert_species_to_mixture_variables_bubbles_acc, &
         s_convert_species_to_mixture_variables, &
         s_convert_species_to_mixture_variables_acc, &
         s_convert_conservative_to_primitive_variables, &
         s_convert_primitive_to_conservative_variables, &
         s_convert_primitive_to_flux_variables, &
         s_finalize_variables_conversion_module

    abstract interface ! =======================================================

        !> The abstract interface to the procedures that are utilized to convert
        !! the mixture and the species variables into the mixture variables. For
        !! more information, refer to:
        !!               1) s_convert_mixture_to_mixture_variables
        !!               2) s_convert_species_to_mixture_variables
        !! @param qK_vf Conservative/primitive variables
        !! @param rho_K Mixture density
        !! @param gamma_K Mixture sp. heat ratio
        !! @param pi_inf_K Mixture stiffness function
        !! @param Re_K Reynolds number
        !! @param i Cell location first index
        !! @param j Cell location second index
        !! @param k Cell location third index
        subroutine s_convert_abstract_to_mixture_variables(qK_vf, rho_K, &
                                                           gamma_K, pi_inf_K, &
                                                           Re_K, i, j, k)

            import :: scalar_field, sys_size, num_fluids

            type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

            real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

            real(kind(0d0)), dimension(2), intent(OUT) :: Re_K

            integer, intent(IN) :: i, j, k

        end subroutine s_convert_abstract_to_mixture_variables

        !> The abstract interface to the procedures that are used to compute the
        !! Roe and the arithmetic average states. For additional information see:
        !!                 1) s_compute_roe_average_state
        !!                 2) s_compute_arithmetic_average_state
        !! @param i Cell location first index
        !! @param j Cell location second index
        !! @param k Cell location third index
        subroutine s_compute_abstract_average_state(i, j, k)

            integer, intent(IN) :: i, j, k

        end subroutine s_compute_abstract_average_state

    end interface ! ============================================================

    !> @name  Left/right states
    !> @{
    real(kind(0d0))                              ::    rho_L, rho_R      !< left/right states density
    real(kind(0d0)), allocatable, dimension(:)   ::    vel_L, vel_R      !< left/right states velocity
    real(kind(0d0))                              ::   pres_L, pres_R      !< left/right states pressure
    real(kind(0d0))                              ::      E_L, E_R      !< left/right states total energy
    real(kind(0d0))                              ::      H_L, H_R      !< left/right states enthalpy
    real(kind(0d0)), allocatable, dimension(:)   ::     mf_L, mf_R      !< left/right states mass fraction
    real(kind(0d0))                              ::  gamma_L, gamma_R      !< left/right states specific heat ratio
    real(kind(0d0))                              :: pi_inf_L, pi_inf_R      !< left/right states liquid stiffness
    real(kind(0d0)), dimension(2)                ::     Re_L, Re_R      !< left/right states Reynolds number
    real(kind(0d0))                              ::   alpha_L, alpha_R    !< left/right states void fraction
    !> @}

    !> @name Averaged states
    !> @{
    real(kind(0d0)), allocatable, dimension(:, :, :) :: rho_avg_sf !< averaged (Roe/arithmetic) density
    real(kind(0d0)), allocatable, dimension(:)     :: vel_avg    !< averaged (Roe/arithmetic) velocity
    real(kind(0d0))                                   :: H_avg      !< averaged (Roe/arithmetic) enthalpy
    type(scalar_field), allocatable, dimension(:)     :: mf_avg_vf  !< averaged (Roe/arithmetic) mass fraction
    real(kind(0d0))                                   :: gamma_avg  !< averaged (Roe/arithmetic) specific heat ratio
    real(kind(0d0)), allocatable, dimension(:, :, :) :: c_avg_sf   !< averaged (Roe/arithmetic) speed of sound

    real(kind(0d0))                                   :: alpha_avg !< averaging for bubbly mixture speed of sound
    real(kind(0d0))                                   :: pres_avg  !< averaging for bubble mixture speed of sound
    !> @}

    real(kind(0d0)) :: ixb,ixe,iyb,iye,izb,ize
    real(kind(0d0)) :: momxb, momxe
    real(kind(0d0)) :: contxb, contxe
    real(kind(0d0)) :: bubxb, bubxe
    real(kind(0d0)) :: advxb, advxe
    real(kind(0d0)),allocatable, dimension(:) :: gammas, pi_infs, bubrs
!$acc declare create(ixb, ixe, iyb, iye, iye, izb, ize, momxb, momxe, bubxb, bubxe, contxb, contxe, advxb, advxe, gammas, pi_infs, bubrs)


    integer :: is1b, is2b, is3b, is1e, is2e, is3e
!$acc declare create(is1b, is2b, is3b, is1e, is2e, is3e)

    ! Density, dynamic pressure, surface energy, specific heat ratio
    ! function, liquid stiffness function, shear and volume Reynolds
    ! numbers and the Weber numbers




    procedure(s_convert_abstract_to_mixture_variables), &
        pointer :: s_convert_to_mixture_variables => null() !<
    !! Pointer to the procedure utilized to convert either the mixture or the
    !! species variables into the mixture variables, based on model equations

    procedure(s_compute_abstract_average_state), &
        pointer :: s_compute_average_state => null() !<
    !! Pointer to the subroutine utilized to calculate either the Roe or the
    !! arithmetic average state variables, based on the chosen average state

contains

    !> This procedure is used alongside with the gamma/pi_inf
        !!      model to transfer the density, the specific heat ratio
        !!      function and liquid stiffness function from the vector
        !!      of conservative or primitive variables to their scalar
        !!      counterparts.
        !! @param qK_vf conservative or primitive variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        !! @param rho_K density
        !! @param gamma_K  specific heat ratio function
        !! @param pi_inf_K liquid stiffness
        !! @param Re_k Reynolds number
    subroutine s_convert_mixture_to_mixture_variables(qK_vf, rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      Re_K, i, j, k)
!$acc routine seq

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(2), intent(OUT) :: Re_K

        integer, intent(IN) :: i, j, k

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively
        rho_K = qK_vf(1)%sf(i, j, k)
        gamma_K = qK_vf(gamma_idx)%sf(i, j, k)
        pi_inf_K = qK_vf(pi_inf_idx)%sf(i, j, k)

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
        !! @param Re_K mixture Reynolds number
        !! @param i Cell index
        !! @param j Cell index
        !! @param k Cell index
    subroutine s_convert_species_to_mixture_variables_bubbles(qK_vf, rho_K, &
                                                              gamma_K, pi_inf_K, &
                                                              Re_K, i, j, k)
        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(2), intent(OUT) :: Re_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K  !<
            !! Partial densities and volume fractions

        integer, intent(IN) :: i, j, k
        integer :: l

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        ! alpha_rho_K(1) = qK_vf(i)%sf(i,j,k)
        ! alpha_K(1)     = qK_vf(E_idx+i)%sf(i,j,k)

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively

        if (model_eqns == 4) then
            rho_K = qK_vf(1)%sf(i, j, k)
            gamma_K = gammas(1)!qK_vf(gamma_idx)%sf(i,j,k)
            pi_inf_K = pi_infs(1) !qK_vf(pi_inf_idx)%sf(i,j,k)
        else if ((model_eqns == 2) .and. bubbles) then
            rho_k = 0d0; gamma_k = 0d0; pi_inf_k = 0d0

            if (mpp_lim .and. (num_fluids > 2)) then
                do l = 1, num_fluids
                    rho_k = rho_k + qK_vf(l)%sf(i, j, k)
                    gamma_k = gamma_k + qK_vf(l + E_idx)%sf(i, j, k)*gammas(l)
                    pi_inf_k = pi_inf_k + qK_vf(l + E_idx)%sf(i, j, k)*pi_infs(l)
                end do
            else if (num_fluids == 2) then
                rho_K = qK_vf(1)%sf(i, j, k)
                gamma_K = gammas(1)
                pi_inf_K = pi_infs(1)
            else if (num_fluids > 2) then
                do l = 1, num_fluids - 1 !leave out bubble part of mixture
                    rho_k = rho_k + qK_vf(l)%sf(i, j, k)
                    gamma_k = gamma_k + qK_vf(l + E_idx)%sf(i, j, k)*gammas(l)
                    pi_inf_k = pi_inf_k + qK_vf(l + E_idx)%sf(i, j, k)*pi_infs(l)
                end do
            else
                rho_K = qK_vf(1)%sf(i, j, k)
                gamma_K = gammas(1)
                pi_inf_K = pi_infs(1)
            end if
        end if

    end subroutine s_convert_species_to_mixture_variables_bubbles ! ----------------

    !>  This subroutine is designed for the volume fraction model
        !!              and provided a set of either conservative or primitive
        !!              variables, computes the density, the specific heat ratio
        !!              function and the liquid stiffness function from q_vf and
        !!              stores the results into rho, gamma and pi_inf.
        !! @param qK_vf primitive variables
        !! @param rho_K density
        !! @param gamma_K specific heat ratio
        !! @param pi_inf_K liquid stiffness
        !! @param Re_K mixture Reynolds number
        !! @param k Cell index
        !! @param l Cell index
        !! @param r Cell index
    subroutine s_convert_species_to_mixture_variables(qK_vf, rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      Re_K,  k, l, r)

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(2), intent(OUT) :: Re_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K !<
            !! Partial densities and volume fractions

        integer, intent(IN) :: k, l, r

        integer :: i, j !< Generic loop iterators

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture

        do i = 1, num_fluids
            alpha_rho_K(i) = qK_vf(i)%sf(k, l, r)
            alpha_K(i) = qK_vf(advxb + i - 1)%sf(k, l, r)
        end do

        if (mpp_lim) then

            do i = 1, num_fluids
                alpha_rho_K(i) = max(0d0, alpha_rho_K(i))
                alpha_K(i) = min(max(0d0, alpha_K(i)), 1d0)
            end do

            alpha_K = alpha_K/max(sum(alpha_K),1d-16)

        end if

        ! Calculating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively, from the species analogs
        rho_K = 0d0; gamma_K = 0d0; pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do

        ! Computing the shear and bulk Reynolds numbers from species analogs
        do i = 1, 2

            Re_K(i) = -1d6; !if (Re_size(i) > 0) Re_K(i) = 0d0

            !do j = 1, Re_size(i)
            !    Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) &
            !              + Re_K(i)
            !end do

            Re_K(i) = 1d0/max(Re_K(i), 1d-16)

        end do

    end subroutine s_convert_species_to_mixture_variables ! ----------------


    subroutine s_convert_species_to_mixture_variables_acc( rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                       alpha_K, alpha_rho_K,  k, l, r)
!$acc routine seq

        real(kind(0d0)), intent(INOUT) :: rho_K, gamma_K, pi_inf_K

        ! TODO: Ensure these were meant to be "INOUT" and not just "IN"
        real(kind(0d0)), dimension(:), intent(INOUT) :: alpha_rho_K, alpha_K !<
            !! Partial densities and volume fractions

        integer, intent(IN) :: k, l, r

        integer :: i, j !< Generic loop iterators
        real(kind(0d0)) :: alpha_K_sum

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        rho_K = 0d0
        gamma_K = 0d0
        pi_inf_K = 0d0

        alpha_K_sum = 0d0

        if (mpp_lim) then
!$acc loop seq
            do i = 1, num_fluids
                alpha_rho_K(i) = max(0d0, alpha_rho_K(i))
                alpha_K(i) = min(max(0d0, alpha_K(i)), 1d0)
                alpha_K_sum = alpha_K_sum + alpha_K(i)
            end do

            alpha_K = alpha_K/max(alpha_K_sum,sgm_eps)

        end if

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do


    end subroutine s_convert_species_to_mixture_variables_acc ! ----------------

    subroutine s_convert_species_to_mixture_variables_bubbles_acc( rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                       alpha_K, alpha_rho_K,  k, l, r)
!$acc routine seq

        real(kind(0d0)), intent(INOUT) :: rho_K, gamma_K, pi_inf_K


        real(kind(0d0)), dimension(:), intent(IN) :: alpha_rho_K, alpha_K !<
            !! Partial densities and volume fractions
        integer, intent(IN) :: k, l, r
        integer :: i, j !< Generic loop iterators

        rho_K = 0d0
        gamma_K = 0d0
        pi_inf_K = 0d0

        if(mpp_lim .and. (model_eqns == 2) .and. (num_fluids > 2)) then
            do i = 1, num_fluids
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            end do
        else if((model_eqns == 2) .and. (num_fluids > 2)) then
            do i = 1, num_fluids - 1
                rho_K = rho_K + alpha_rho_K(i)
                gamma_K = gamma_K + alpha_K(i)*gammas(i)
                pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            end do 
        else           
            rho_K = alpha_rho_K(1)
            gamma_K = gammas(1)
            pi_inf_K = pi_infs(1)
        end if

    end subroutine s_convert_species_to_mixture_variables_bubbles_acc


    !> The goal of this subroutine is to calculate the Roe
        !!      average density, velocity, enthalpy, mass fractions,
        !!      specific heat ratio function and the speed of sound,
        !!      at the cell-boundaries, from the left and the right
        !!      cell-average variables.
        !!  @param j Cell index
        !!  @param k Cell index
        !!  @param l Cell index
    subroutine s_compute_roe_average_state(j, k, l) ! ------------------------

        integer, intent(IN) :: j, k, l
        integer :: i

        rho_avg_sf(j, k, l) = sqrt(rho_L*rho_R)

        vel_avg = (sqrt(rho_L)*vel_L + sqrt(rho_R)*vel_R)/ &
                  (sqrt(rho_L) + sqrt(rho_R))

        H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/ &
                (sqrt(rho_L) + sqrt(rho_R))

        do i = 1, cont_idx%end
            mf_avg_vf(i)%sf(j, k, l) = (sqrt(rho_L)*mf_L(i) &
                                        + sqrt(rho_R)*mf_R(i)) &
                                       /(sqrt(rho_L) &
                                         + sqrt(rho_R))
        end do

        gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/ &
                    (sqrt(rho_L) + sqrt(rho_R))

        c_avg_sf(j, k, l) = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)

    end subroutine s_compute_roe_average_state ! ---------------------------

    !>  The goal of this subroutine is to compute the arithmetic
        !!      average density, velocity, enthalpy, mass fractions, the
        !!      specific heat ratio function and the sound speed, at the
        !!      cell-boundaries, from the left and right cell-averages.
        !!  @param j Cell index
        !!  @param k Cell index
        !!  @param l Cell index
    subroutine s_compute_arithmetic_average_state(j, k, l) ! -----------------

        integer, intent(IN) :: j, k, l
        integer :: i

        rho_avg_sf(j, k, l) = 5d-1*(rho_L + rho_R)
        vel_avg = 5d-1*(vel_L + vel_R)

        do i = 1, cont_idx%end
            mf_avg_vf(i)%sf(j, k, l) = 5d-1*(mf_L(i) + mf_R(i))
        end do

        H_avg = 5d-1*(H_L + H_R)

        gamma_avg = 5d-1*(gamma_L + gamma_R)

        alpha_avg = 5d-1*(alpha_L + alpha_R)
        pres_avg = 5d-1*(pres_L + pres_R)

        if (model_eqns .ne. 4) then
            c_avg_sf(j, k, l) = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
        else
            ! For Tait EOS
            c_avg_sf(j, k, l) = sqrt( &
                                (1d0/fluid_pp(1)%gamma + 1d0)* &
                                (pres_avg + fluid_pp(1)%pi_inf)/ &
                                (rho_avg_sf(j, k, l)*(1d0 - alpha_avg)) &
                                )
        end if

    end subroutine s_compute_arithmetic_average_state ! --------------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_variables_conversion_module() ! ----------------


        integer :: i, j
        
        momxb = mom_idx%beg; momxe = mom_idx%end
        bubxb = bub_idx%beg; bubxe = bub_idx%end
        advxb = adv_idx%beg; advxe = adv_idx%end
        contxb = cont_idx%beg; contxe = cont_idx%end
!$acc update device(momxb, momxe, bubxb, bubxe, advxb, advxe, contxb, contxe)

        ixb = -buff_size
        ixe = m - ixb

        iyb = 0; iye = 0; izb = 0; ize = 0;

        if(n > 0) then
            iyb = -buff_size
            iye = n - iyb
        end if

        if(p > 0) then
            izb = -buff_size
            ize = p - izb
        end if

!$acc update device(ixb, ixe, iyb, iye, izb, ize)
        
        allocate(gammas(1:num_fluids))
        allocate(pi_infs(1:num_fluids))

        allocate(bubrs(1:nb))

        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do
!$acc update device(gammas, pi_infs)

        if (bubbles) then

            do i = 1, nb
                bubrs(i) = bub_idx%rs(i)
            end do

        end if

!$acc update device(bubrs)




        
!$acc update device(small_alf, dflt_real, dflt_int)
!$acc update device(pi, dt, sys_size, pref, rhoref, gamma_idx, pi_inf_idx, E_idx, alf_idx, mpp_lim, bubbles, alt_soundspeed, avg_state, num_fluids, model_eqns, num_dims, mixture_err, nb, weight, grid_geometry, cyl_coord, mapped_weno, mp_weno, weno_eps)
!$acc update device(nb, R0ref, Ca, Web, Re_inv, weight, R0, V0, bubbles, polytropic, polydisperse, qbmm, nmom, nnode, nmomsp, nmomtot, R0_type, ptil, bubble_model, thermal, poly_sigma, sgm_eps)


!$acc update device(R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v, k_n, k_v, pb0, mass_n0, mass_v0, Pe_T, Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN , mul0, ss, gamma_v, mu_v, gamma_m, gamma_n, mu_n, gam)


!$acc update device(monopole, num_mono, mono)
        do i = 1, num_mono
!$acc update device(mono(i)%mag)
!$acc update device(mono(i)%length)
!$acc update device(mono(i)%npulse)
!$acc update device(mono(i)%dir)
!$acc update device(mono(i)%delay)

        end do

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        if (model_eqns == 1) then        ! gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        elseif (bubbles) then           ! Volume fraction for bubbles
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if



    end subroutine s_initialize_variables_conversion_module ! --------------

    !> The following procedure handles the conversion between
        !!      the conservative variables and the primitive variables.
        !! @param qK_cons_vf Conservative variables
        !! @param qK_prim_vf Primitive variables
        !! @param gm_alphaK_vf Gradient magnitude of the volume fraction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, &
                                                             qK_prim_vf, &
                                                             gm_alphaK_vf, &
                                                             ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qK_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qK_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(IN) :: gm_alphaK_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)),   dimension(2) :: alpha_K, alpha_rho_K
        real(kind(0d0)) :: rho_K, gamma_K, pi_inf_K, dyn_pres_K
        real(kind(0d0)), dimension(nb) :: nRtmp
        real(kind(0d0)) :: vftmp, nR3, nbub_sc

        integer :: i, j, k, l !< Generic loop iterators


        if((model_eqns .ne. 4) .and. (bubbles .neqv. .true.)) then 
!$acc parallel loop collapse(3) gang vector default(present) private( alpha_K, alpha_rho_K)
            do l = izb, ize
                do k = iyb, iye
                    do j = ixb, ixe
                        dyn_pres_K = 0d0  
                   
!$acc loop seq
                        do i = 1, num_fluids
                            alpha_rho_K(i) = qK_cons_vf(i)%sf(j, k, l)
                            alpha_K(i) = qK_cons_vf(advxb + i - 1)%sf(j, k, l)
                        end do

                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, alpha_K, alpha_rho_K, j, k, l)


!$acc loop seq
                        do i = momxb, momxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                            /max( rho_K, sgm_eps)
                                dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf(i)%sf(j, k, l) &
                                             *qK_prim_vf(i)%sf(j, k, l)
                        end do 


                        qK_prim_vf(E_idx)%sf(j, k, l) = (qK_cons_vf(E_idx)%sf(j, k, l) &
                                 - dyn_pres_K - pi_inf_K )/gamma_K
                    end do
                end do
            end do
!$acc end parallel loop 
        elseif((model_eqns .ne. 4)) then 
!$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_K, alpha_K, nRtmp)
            do l = izb, ize
                do k = iyb, iye
                    do j = ixb, ixe
                        dyn_pres_K = 0d0
!$acc loop seq 
                        do i = 1, num_fluids
                            alpha_rho_K(i) = qK_cons_vf(i)%sf(j, k, l)
                            alpha_K(i) = qK_cons_vf(advxb + i - 1)%sf(j, k, l)
                        end do

                        call s_convert_species_to_mixture_variables_bubbles_acc(rho_K, gamma_K, pi_inf_K,  alpha_K, alpha_rho_K, j, k, l)


!$acc loop seq
                        do i = momxb, momxe
                                qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                            /max(rho_K, sgm_eps)
                                dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf(i)%sf(j, k, l) &
                                             *qK_prim_vf(i)%sf(j, k, l)
                        end do 

                       qK_prim_vf(E_idx)%sf(j, k, l) = (((qK_cons_vf(E_idx)%sf(j, k, l) - dyn_pres_K)/(1.d0 - qK_cons_vf(alf_idx)%sf(j, k, l)))  - pi_inf_K  )/gamma_K

!$acc loop seq 
                        do i = 1, nb
                            nRtmp(i) = qK_cons_vf(bubrs(i))%sf(j, k, l)
                        end do

                        vftmp = qK_cons_vf(alf_idx)%sf(j, k, l)

                        call s_comp_n_from_cons(vftmp, nRtmp, nbub_sc)

!$acc loop seq 
                        do i = bubxb, bubxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/nbub_sc
                        end do

                    end do
                end do
            end do
!$acc end parallel loop 
        else 
!$acc parallel loop collapse(3) gang vector default(present) private(rho_K, gamma_K, pi_inf_K,  dyn_pres_K)
            do l = izb, ize
                do k = iyb, iye
                    do j = ixb, ixe
                        dyn_pres_K = 0d0
!$acc loop 
                        do i = momxb, momxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                        /qK_cons_vf(1)%sf(j, k, l)
                        end do 
!$acc end loop
                        qK_prim_vf(E_idx)%sf(j, k, l) = (pref + pi_infs(1))* &
                            (( qK_cons_vf(1)%sf(j, k, l)/ (rhoref*(1.d0 - qK_cons_vf(E_idx + 1)%sf(j, k, l))) &
                             )**(1.d0/gammas(1) + 1.d0)) - pi_infs(1)
                    end do
                end do
            end do
!$acc end parallel loop 
        end if


    end subroutine s_convert_conservative_to_primitive_variables ! ---------



    !>  The following procedure handles the conversion between
        !!      the primitive variables and the conservative variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Conservative variables
        !!  @param gm_alphaK_vf Gradient magnitude of the volume fractions
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_conservative_variables(qK_prim_vf, &
                                                             qK_cons_vf, &
                                                             gm_alphaK_vf, &
                                                             ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qK_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qK_cons_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(IN) :: gm_alphaK_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ! Calculating the momentum and energy from the velocity and pressure
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    if (proc_rank == 0) then
                        print '(A)', 'Conversion from primitive to '// &
                            'conservative variables not '// &
                            'implemented. Exiting ...'
                        call s_mpi_abort()
                    end if

                end do
            end do
        end do

    end subroutine s_convert_primitive_to_conservative_variables ! ---------




    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, & ! ------
                                                     FK_vf, &
                                                     FK_src_vf, &
                                                     is1, is2, is3, s2b, s3b)

        integer :: s2b, s3b
        real(kind(0d0)), dimension(0:, s2b:, s3b:, 1:), intent(IN) :: qK_prim_vf
        real(kind(0d0)), dimension(0:, s2b:, s3b:, 1:), intent(INOUT) :: FK_vf 
        real(kind(0d0)), dimension(0:, s2b:, s3b:, INT(advxb):), intent(INOUT) :: FK_src_vf
 
        type(bounds_info), intent(IN) :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(kind(0d0)), dimension(INT(contxe))           :: alpha_rho_K
        real(kind(0d0)), dimension(num_fluids)            :: alpha_K  
        real(kind(0d0))                                   ::       rho_K
        real(kind(0d0)), dimension(num_dims)              ::       vel_K
        real(kind(0d0))                                   :: vel_K_sum  
        real(kind(0d0))                                   ::      pres_K
        real(kind(0d0))                                   ::         E_K
        real(kind(0d0))                                   ::     gamma_K
        real(kind(0d0))                                   ::    pi_inf_K
        real(kind(0d0)), dimension(2)                     ::        Re_K

        integer :: i, j, k, l !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        !$acc update device(is1b, is2b, is3b, is1e, is2e, is3e)

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity

!$acc parallel loop collapse(3) gang vector default(present) private(alpha_rho_K, vel_K, alpha_K)
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e

!$acc loop seq
                    do i = 1, contxe
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

!$acc loop seq
                    do i = advxb, advxe
                        alpha_K(i - E_idx) = qK_prim_vf(j, k, l, i)
                    end do
!$acc loop seq
                    do i = 1, num_dims
                        vel_K(i) = qK_prim_vf(j, k, l, contxe + i)
                    end do

                    vel_K_sum = 0d0
!$acc loop seq
                    do i = 1, num_dims
                        vel_K_sum = vel_K_sum + vel_K(i)**2d0
                    end do

                    pres_K = qK_prim_vf(j, k, l, E_idx)

                    if(bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho_K, gamma_K, pi_inf_K, alpha_K, alpha_rho_K, j, k, l)

                    else
                        call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, alpha_K, alpha_rho_K, j, k, l)
                    end if


                    ! Computing the energy from the pressure
                    E_K = gamma_K*pres_K + pi_inf_K &
                          + 5d-1*rho_K*vel_K_sum


                    ! mass flux, this should be \alpha_i \rho_i u_i
!$acc loop seq
                    do i = 1, contxe
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

!$acc loop seq
                    do i = 1, num_dims
                        FK_vf(j, k, l, contxe + dir_idx(i)) = &
                            rho_K*vel_K(dir_idx(1)) &
                            *vel_K(dir_idx(i)) &
                            + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, E_idx) = vel_K(dir_idx(1))*(E_K + pres_K)

                    ! have been using == 2
                    if (riemann_solver == 1) then
!$acc loop seq
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = 0d0
                            FK_src_vf(j, k, l, i) = alpha_K(i - E_idx)
                        end do

                    else
                        ! Could be bubbles!
!$acc loop seq
                        do i = advxb, advxe
                            FK_vf(j, k, l, i) = vel_K(dir_idx(1))*alpha_K(i - E_idx)
                        end do

!$acc loop seq
                        do i = advxb, advxe
                            FK_src_vf(j, k, l, i) = vel_K(dir_idx(1))
                        end do
                        
                    end if
                end do
            end do
        end do

    end subroutine s_convert_primitive_to_flux_variables ! -----------------




    subroutine s_finalize_variables_conversion_module() ! ------------------

        deallocate(gammas, pi_infs)
        deallocate(bubrs)

        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module ! ----------------

end module m_variables_conversion
