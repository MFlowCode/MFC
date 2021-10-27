!>
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module consists of subroutines used in the conversion of the
!!              conservative variables into the primitive variables. In addition,
!!              the module also contains subroutines used to compute the mixture
!!              variables. For a specific time-step undergoing the post-process,
!!              the mixture variables are stored everywhere on the grid so that
!!              they may possibly be outputted or used to compute other derived
!!              flow quantities.
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_variables_conversion_module, &
 s_convert_to_mixture_variables, &
 s_convert_mixture_to_mixture_variables, &
 s_convert_species_to_mixture_variables_bubbles, &
 s_convert_species_to_mixture_variables, &
 s_convert_conservative_to_primitive_variables, &
 s_finalize_variables_conversion_module

    !> Abstract interface to two subroutines designed for the transfer/conversion
    !! of the mixture/species variables to the mixture variables
    abstract interface

        !>  Structure of the s_convert_mixture_to_mixture_variables
        !!      and s_convert_species_to_mixture_variables subroutines
        !! @param q_cons_vf Conservative variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        subroutine s_convert_xxxxx_to_mixture_variables(q_cons_vf, i, j, k)

            ! Importing the derived type scalar_field from m_derived_types.f90
            ! and global variable sys_size, from m_global_variables.f90, as
            ! the abstract interface does not inherently have access to them
            import :: scalar_field, sys_size

            type(scalar_field), &
                dimension(sys_size), &
                intent(IN) :: q_cons_vf

            integer, intent(IN) :: i, j, k

        end subroutine s_convert_xxxxx_to_mixture_variables

    end interface

    ! NOTE: These abstract interfaces allow for the declaration of a pointer to
    ! a procedure such that the choice of the model equations does not have to
    ! be queried every time the mixture quantities are needed.

    real(kind(0d0)), allocatable, dimension(:, :, :), public :: rho_sf !< Scalar density function
    real(kind(0d0)), allocatable, dimension(:, :, :), public :: gamma_sf !< Scalar sp. heat ratio function
    real(kind(0d0)), allocatable, dimension(:, :, :), public :: pi_inf_sf !< Scalar liquid stiffness function

    procedure(s_convert_xxxxx_to_mixture_variables), &
        pointer :: s_convert_to_mixture_variables => null() !<
    !! Pointer referencing the subroutine s_convert_mixture_to_mixture_variables
    !! or s_convert_species_to_mixture_variables, based on model equations choice

    integer, private :: flg !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.

contains

    !>  This subroutine is constructed for the gamma/pi_inf model
        !!      and provided a set of conservative variables, transfers
        !!      the density, specific heat ratio function and the liquid
        !!      stiffness function from q_cons_vf to rho_sf, gamma_sf and
        !!      pi_inf_sf.
        !! @param q_cons_vf Conservative variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
    subroutine s_convert_mixture_to_mixture_variables(q_cons_vf, i, j, k) ! --

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        integer, intent(IN) :: i, j, k

        ! Transfering the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        rho_sf(i, j, k) = q_cons_vf(1)%sf(i, j, k)
        gamma_sf(i, j, k) = q_cons_vf(gamma_idx)%sf(i, j, k)
        pi_inf_sf(i, j, k) = q_cons_vf(pi_inf_idx)%sf(i, j, k)

    end subroutine s_convert_mixture_to_mixture_variables ! ----------------

    !>  This subroutine is designed for the volume fraction model
        !!      with sub-grid ensemble averaged bubbles (Bryngelson 2019)
        !!      and provided a set of conservative variables, calculates
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function from q_cons_vf and stores the results
        !!      into rho_sf, gamma_sf and pi_inf_sf.
        !!  @param qK_vf Conservative variables
        !!  @param j cell index to transfer mixture variables
        !!  @param k cell index to transfer mixture variables
        !!  @param l cell index to transfer mixture variables
    subroutine s_convert_species_to_mixture_variables_bubbles(qK_vf, j, k, l)

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf
        integer, intent(IN) :: j, k, l

        integer :: i !< Generic loop iterator

        if (model_eqns == 4) then
            rho_sf = qK_vf(1)%sf(j, k, l)
            gamma_sf = fluid_pp(1)%gamma
            pi_inf_sf = fluid_pp(1)%pi_inf
        else if (model_eqns == 2 .and. bubbles .and. adv_alphan) then
            rho_sf = 0d0; gamma_sf = 0d0; pi_inf_sf = 0d0

            if (mpp_lim .and. num_fluids > 1) then
                do i = 1, num_fluids
                    rho_sf = rho_sf + qK_vf(i + E_idx)%sf(j, k, l)*qK_vf(i)%sf(j, k, l)
                    gamma_sf = gamma_sf + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf_sf = pi_inf_sf + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                end do
            else if (num_fluids > 1) then
                do i = 1, num_fluids - 1 !leave out bubble part of mixture
                    rho_sf = rho_sf + qK_vf(i + E_idx)%sf(j, k, l)*qK_vf(i)%sf(j, k, l)
                    gamma_sf = gamma_sf + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%gamma
                    pi_inf_sf = pi_inf_sf + qK_vf(i + E_idx)%sf(j, k, l)*fluid_pp(i)%pi_inf
                end do
            else
                rho_sf = qK_vf(1)%sf(j, k, l)
                gamma_sf = fluid_pp(1)%gamma
                pi_inf_sf = fluid_pp(1)%pi_inf
            end if
        end if

    end subroutine s_convert_species_to_mixture_variables_bubbles ! ----------------

    !>  This subroutine is designed for the volume fraction model
        !!      and provided a set of conservative variables, calculates
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function from q_cons_vf and stores the results
        !!      into rho_sf, gamma_sf and pi_inf_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param j cell index to transfer mixture variables
        !!  @param k cell index to transfer mixture variables
        !!  @param l cell index to transfer mixture variables
    subroutine s_convert_species_to_mixture_variables(q_cons_vf, j, k, l) ! --

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        integer, intent(IN) :: j, k, l

        integer :: i !< Generic loop iterator

        ! Computing the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively
        if (adv_alphan) then
            if (bubbles .neqv. .true.) then
                rho_sf(j, k, l) = 0d0
                gamma_sf(j, k, l) = 0d0
                pi_inf_sf(j, k, l) = 0d0

                do i = 1, num_fluids
                    rho_sf(j, k, l) = rho_sf(j, k, l) &
                                      + q_cons_vf(i)%sf(j, k, l)
                    gamma_sf(j, k, l) = gamma_sf(j, k, l) &
                                        + q_cons_vf(i + E_idx)%sf(j, k, l) &
                                        *fluid_pp(i)%gamma
                    pi_inf_sf(j, k, l) = pi_inf_sf(j, k, l) &
                                         + q_cons_vf(i + E_idx)%sf(j, k, l) &
                                         *fluid_pp(i)%pi_inf
                end do
            else
                rho_sf(j, k, l) = q_cons_vf(1)%sf(j, k, l)
                gamma_sf(j, k, l) = fluid_pp(1)%gamma
                pi_inf_sf(j, k, l) = fluid_pp(1)%pi_inf
            end if
        else
            rho_sf(j, k, l) = q_cons_vf(num_fluids)%sf(j, k, l)
            gamma_sf(j, k, l) = fluid_pp(num_fluids)%gamma
            pi_inf_sf(j, k, l) = fluid_pp(num_fluids)%pi_inf

            do i = 1, num_fluids - 1
                rho_sf(j, k, l) = rho_sf(j, k, l) &
                                  + q_cons_vf(i)%sf(j, k, l)
                gamma_sf(j, k, l) = gamma_sf(j, k, l) &
                                    + q_cons_vf(i + E_idx)%sf(j, k, l) &
                                    *(fluid_pp(i)%gamma &
                                      - fluid_pp(num_fluids)%gamma)
                pi_inf_sf(j, k, l) = pi_inf_sf(j, k, l) &
                                     + q_cons_vf(i + E_idx)%sf(j, k, l) &
                                     *(fluid_pp(i)%pi_inf &
                                       - fluid_pp(num_fluids)%pi_inf)
            end do

        end if

    end subroutine s_convert_species_to_mixture_variables ! ----------------

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    subroutine s_initialize_variables_conversion_module() ! -------------------

        ! Allocating the density, the specific heat ratio function and the
        ! liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then

            ! Simulation is 3D
            if (p > 0) then

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 -buff_size:p + buff_size))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   -buff_size:p + buff_size))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    -buff_size:p + buff_size))

                ! Simulation is 2D
            else

                allocate (rho_sf(-buff_size:m + buff_size, &
                                 -buff_size:n + buff_size, &
                                 0:0))
                allocate (gamma_sf(-buff_size:m + buff_size, &
                                   -buff_size:n + buff_size, &
                                   0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                    -buff_size:n + buff_size, &
                                    0:0))

            end if

            ! Simulation is 1D
        else

            allocate (rho_sf(-buff_size:m + buff_size, &
                             0:0, &
                             0:0))
            allocate (gamma_sf(-buff_size:m + buff_size, &
                               0:0, &
                               0:0))
            allocate (pi_inf_sf(-buff_size:m + buff_size, &
                                0:0, &
                                0:0))

        end if

        ! Depending on the model selection for the equations of motion, the
        ! appropriate procedure for the conversion to the mixture variables
        ! is targeted by the procedure pointer

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        elseif (bubbles) then           ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables_bubbles
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

        ! Annotating the dimensionality of the dataset undergoing the post-
        ! process. A flag value of 1 indicates that the dataset is 3D, while
        ! a flag value of 0 indicates that it is not.
        if (p > 0) then
            flg = 1
        else
            flg = 0
        end if

    end subroutine s_initialize_variables_conversion_module ! -----------------

    !> Converts the conservative variables to the primitive ones
        !!  @param q_cons_vf Conservative variabels
        !!  @param q_prim_vf Primitive variables
    subroutine s_convert_conservative_to_primitive_variables(q_cons_vf, &
                                                             q_prim_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: q_prim_vf

        ! Dynamic pressure, as defined in the incompressible flow sense
        real(kind(0d0)) :: dyn_pres

        ! Bubble parameters
        real(kind(0d0)) :: nbub
        real(kind(0d0)), dimension(:), allocatable :: nRtmp

        integer :: i, j, k, l !< Generic loop iterators

        allocate (nRtmp(nb))

        ! Converting the conservative variables to the primitive variables
        do l = -buff_size*flg, (p + buff_size)*flg
            do k = -buff_size, n + buff_size
                do j = -buff_size, m + buff_size

                    ! Obtaining the density, specific heat ratio function
                    ! and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_cons_vf, j, k, l)

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, cont_idx%end
                        q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed
                    ! iteratively by cycling through the momentum equations
                    dyn_pres = 0d0

                    ! Computing velocity and dynamic pressure from momenta
                    do i = mom_idx%beg, mom_idx%end
                        if (model_eqns == 4) then
                            !u = \rho u/ \rho
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/ &
                                                       q_cons_vf(1)%sf(j, k, l)
                        else
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/ &
                                                       rho_sf(j, k, l)
                            dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)* &
                                       q_prim_vf(i)%sf(j, k, l)/2d0
                        end if
                    end do

                    if (model_eqns == 4) then
                        ! Computing the pressure from the energy
                        ! Tait EOS
                        ! p_l = (pl0 + pi_infty)(rho/(rho_l0(1-alf)))^gamma
                        ! - pi_infty
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            (pref + fluid_pp(1)%pi_inf)* &
                            (( &
                             q_cons_vf(1)%sf(j, k, l)/ &
                             (rhoref*(1.d0 - q_cons_vf(alf_idx)%sf(j, k, l))) &
                             )**(1.d0/fluid_pp(1)%gamma + 1.d0)) - fluid_pp(1)%pi_inf
                    else if ((model_eqns .ne. 4) .and. bubbles .neqv. .true.) then
                        ! Computing the pressure from the energy
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            (q_cons_vf(E_idx)%sf(j, k, l) &
                             - dyn_pres - pi_inf_sf(j, k, l))/gamma_sf(j, k, l)
                    else
                        ! p = ( E/(1-alf) - 0.5 rho u u/(1-alf) - pi_inf_k )/gamma_k
                        q_prim_vf(E_idx)%sf(j, k, l) = &
                            ((q_cons_vf(E_idx)%sf(j, k, l) &
                              - dyn_pres)/(1.d0 - q_cons_vf(alf_idx)%sf(j, k, l)) &
                             - pi_inf_sf(j, k, l) &
                             )/gamma_sf(j, k, l)
                    end if

                    ! Set partial pressures to mixture pressure
                    if (model_eqns == 3) then
                        do i = internalEnergies_idx%beg, internalEnergies_idx%end
                            q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
                        end do
                    end if

                    ! Transferring the advection equation(s) variable(s)
                    do i = adv_idx%beg, adv_idx%end
                        q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do

                    ! \phi = (n\phi)/n  (n = nbub)
                    if (bubbles) then
                        ! n = sqrt( 4pi/(3 alpha) * (nR)**3 )
                        do i = 1, nb
                            nRtmp(i) = q_cons_vf(bub_idx%rs(i))%sf(j, k, l)
                        end do
                        call s_comp_n_from_cons(q_cons_vf(alf_idx)%sf(j, k, l), nRtmp, nbub)
                        do i = bub_idx%beg, sys_size
                            q_prim_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)/nbub
                        end do
                    end if
                end do
            end do
        end do

        deallocate (nRtmp)

    end subroutine s_convert_conservative_to_primitive_variables ! ---------

    !> Deallocation procedures for the module
    subroutine s_finalize_variables_conversion_module() ! ----------------

        ! Deallocating the density, the specific heat ratio function and the
        ! liquid stiffness function
        deallocate (rho_sf)
        deallocate (gamma_sf)
        deallocate (pi_inf_sf)

        ! Nullifying the procedure pointer to the subroutine transfering/
        ! computing the mixture/species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module ! --------------

end module m_variables_conversion
