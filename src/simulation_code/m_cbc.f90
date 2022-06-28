!>
!! @file m_cbc.f90
!! @brief Contains module m_cbc
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The module features a large database of characteristic boundary
!!              conditions (CBC) for the Euler system of equations. This system
!!              is augmented by the appropriate advection equations utilized to
!!              capture the material interfaces. The closure is achieved by the
!!              stiffened equation of state and mixture relations. At this time,
!!              the following CBC are available:
!!                           1) Slip Wall
!!                           2) Nonreflecting Subsonic Buffer
!!                           3) Nonreflecting Subsonic Inflow
!!                           4) Nonreflecting Subsonic Outflow
!!                           5) Force-Free Subsonic Outflow
!!                           6) Constant Pressure Subsonic Outflow
!!                           7) Supersonic Inflow
!!                           8) Supersonic Outflow
!!              Please refer to Thompson (1987, 1990) for detailed descriptions.
module m_cbc

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_cbc_module, s_cbc, s_finalize_cbc_module

    abstract interface ! =======================================================

        !> Abstract interface to the procedures that are utilized to calculate
        !! the L variables. For additional information refer to the following:
        !!            1) s_compute_slip_wall_L
        !!            2) s_compute_nonreflecting_subsonic_buffer_L
        !!            3) s_compute_nonreflecting_subsonic_inflow_L
        !!            4) s_compute_nonreflecting_subsonic_outflow_L
        !!            5) s_compute_force_free_subsonic_outflow_L
        !!            6) s_compute_constant_pressure_subsonic_outflow_L
        !!            7) s_compute_supersonic_inflow_L
        !!            8) s_compute_supersonic_outflow_L
        !! @param dflt_int Default null integer
        subroutine s_compute_abstract_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------

        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        end subroutine s_compute_abstract_L

    end interface ! ============================================================

    !! The cell-average primitive variables. They are obtained by reshaping (RS)
    !! q_prim_vf in the coordinate direction normal to the domain boundary along
    !! which the CBC is applied.

    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsx_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: q_prim_rsz_vf


    type(scalar_field), allocatable, dimension(:) :: F_rs_vf, F_src_rs_vf !<
    !! Cell-average fluxes (src - source). These are directly determined from the
    !! cell-average primitive variables, q_prims_rs_vf, and not a Riemann solver.

    real(kind(0d0)), allocatable, dimension(:, :, :, :)  :: F_rsx_vf, F_src_rsx_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :)  :: F_rsy_vf, F_src_rsy_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :)  :: F_rsz_vf, F_src_rsz_vf !<

    real(kind(0d0)), allocatable, dimension(:, :, :, :)  :: flux_rsx_vf, flux_src_rsx_vf !<
    real(kind(0d0)), allocatable, dimension(:, :, :, :)  :: flux_rsy_vf, flux_src_rsy_vf
    real(kind(0d0)), allocatable, dimension(:, :, :, :) :: flux_rsz_vf, flux_src_rsz_vf
    !! The cell-boundary-average of the fluxes. They are initially determined by
    !! reshaping flux_vf and flux_src_vf in a coordinate direction normal to the
    !! domain boundary along which CBC is applied. flux_rs_vf and flux_src_rs_vf
    !! are subsequently modified based on the selected CBC.
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho   !< Cell averaged partial densiy
    real(kind(0d0))                              :: rho         !< Cell averaged density
    real(kind(0d0)), allocatable, dimension(:)   :: vel         !< Cell averaged velocity
    real(kind(0d0))                              :: pres        !< Cell averaged pressure
    real(kind(0d0))                              :: E           !< Cell averaged energy
    real(kind(0d0))                              :: H           !< Cell averaged enthalpy
    real(kind(0d0)), allocatable, dimension(:)   :: adv         !< Cell averaged advected variables
    real(kind(0d0)), allocatable, dimension(:)   :: mf          !< Cell averaged mass fraction
    real(kind(0d0))                              :: gamma       !< Cell averaged specific heat ratio
    real(kind(0d0))                              :: pi_inf      !< Cell averaged liquid stiffness
    real(kind(0d0))                              :: c           !< Cell averaged speed of sound
    real(kind(0d0)), dimension(2)   :: Re          !< Cell averaged Reynolds numbers

    real(kind(0d0)), allocatable, dimension(:) :: dalpha_rho_ds !< Spatial derivatives in s-dir of partial density
    real(kind(0d0)), allocatable, dimension(:) ::       dvel_ds !< Spatial derivatives in s-dir of velocity
    real(kind(0d0))                            ::      dpres_ds !< Spatial derivatives in s-dir of pressure
    real(kind(0d0)), allocatable, dimension(:) ::       dadv_ds !< Spatial derivatives in s-dir of advection variables

    !! Note that these are only obtained in those cells on the domain boundary along which the
    !! CBC is applied by employing finite differences (FD) on the cell-average primitive variables, q_prim_rs_vf.

    real(kind(0d0)), dimension(3) :: lambda !< Eigenvalues (see Thompson 1987,1990)
    real(kind(0d0)), allocatable, dimension(:) :: L      !< L matrix (see Thompson 1987,1990)

    real(kind(0d0)), allocatable, dimension(:) :: ds !< Cell-width distribution in the s-direction

    ! CBC Coefficients =========================================================
    real(kind(0d0)),  allocatable, dimension(:, :) :: fd_coef_x !< Finite diff. coefficients x-dir
    real(kind(0d0)),  allocatable, dimension(:, :) :: fd_coef_y !< Finite diff. coefficients y-dir
    real(kind(0d0)),  allocatable, dimension(:, :) :: fd_coef_z !< Finite diff. coefficients z-dir
    !! The first dimension identifies the location of a coefficient in the FD
    !! formula, while the last dimension denotes the location of the CBC.

! Bug with NVHPC when using nullified pointers in a declare create
!    real(kind(0d0)), pointer, dimension(:, :) :: fd_coef => null()

    real(kind(0d0)),  allocatable, dimension(:, :, :) :: pi_coef_x !< Polynominal interpolant coefficients in x-dir
    real(kind(0d0)),  allocatable, dimension(:, :, :) :: pi_coef_y !< Polynominal interpolant coefficients in y-dir
    real(kind(0d0)),  allocatable, dimension(:, :, :) :: pi_coef_z !< Polynominal interpolant coefficients in z-dir
    !! The first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the location of the CBC.

!    real(kind(0d0)), pointer, dimension(:, :, :) :: pi_coef => null()
    ! ==========================================================================

   ! procedure(s_compute_abstract_L), pointer :: s_compute_L => null()  !<
!    procedure(s_compute_abstract_L), pointer :: s_compute_L => null() !<
    !! Pointer to procedure used to calculate L variables, based on choice of CBC

    type(bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions

    integer :: dj

    integer :: momxb, momxe, advxb, advxe, contxb, contxe, bubxb, bubxe
    integer :: bcxb, bcxe, bcyb, bcye, bczb, bcze
    real(kind(0d0)), allocatable, dimension(:) :: gammas, pi_infs

!$acc declare create(q_prim_rsx_vf, q_prim_rsy_vf, q_prim_rsz_vf,  F_rsx_vf, F_src_rsx_vf,flux_rsx_vf, flux_src_rsx_vf, &
!$acc                 F_rsy_vf, F_src_rsy_vf,flux_rsy_vf, flux_src_rsy_vf, F_rsz_vf, F_src_rsz_vf,flux_rsz_vf, flux_src_rsz_vf,alpha_rho,vel,adv,mf,Re, &
!$acc                dalpha_rho_ds,dvel_ds,dadv_ds,lambda,L,ds,fd_coef_x,fd_coef_y,fd_coef_z,      &
!$acc                pi_coef_x,pi_coef_y,pi_coef_z, momxb, momxe, advxb, advxe, contxb, contxe, bubxb, bubxe, gammas, pi_infs, bcxb, bcxe, bcyb, bcye, bczb, bcze, is1, is2, is3, dj)

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_cbc_module() ! ---------------------------------

        integer :: i

        if (all((/bc_x%beg, bc_x%end/) > -5) &
            .and. &
            (n > 0 .and. all((/bc_y%beg, bc_y%end/) > -5)) &
            .and. &
            (p > 0 .and. all((/bc_z%beg, bc_z%end/) > -5))) return

        ! Allocating the cell-average primitive variables
        !allocate (q_prim_rs_vf(1:sys_size))



        ! Allocating the cell-average and cell-boundary-average fluxes
        !allocate (F_rs_vf(1:sys_size), F_src_rs_vf(1:sys_size))
        !allocate (flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))
        if(n == 0) then
            is2%beg = 0
        
        else
            is2%beg = -buff_size
        end if
        
        is2%end = n - is2%beg



        if(p == 0) then
            is3%beg = 0
        
        else
            is3%beg = -buff_size
        end if
        is3%end = p - is3%beg

        allocate (q_prim_rsx_vf(0:buff_size, &
                                         is2%beg:is2%end, &
                                         is3%beg:is3%end, 1:sys_size))

        if (weno_order > 1) then

            allocate (F_rsx_vf(0:buff_size, &
                                        is2%beg:is2%end, &
                                        is3%beg:is3%end, 1:adv_idx%end))


            allocate (F_src_rsx_vf(0:buff_size, &
                                                  is2%beg:is2%end, &
                                                  is3%beg:is3%end, adv_idx%beg:adv_idx%end))


        end if

        allocate (flux_rsx_vf(-1:buff_size, &
                                       is2%beg:is2%end, &
                                       is3%beg:is3%end, 1:adv_idx%end))

        allocate (flux_src_rsx_vf(-1:buff_size, &
                                                 is2%beg:is2%end, &
                                                 is3%beg:is3%end, adv_idx%beg:adv_idx%end))


        if(n > 0) then
            
            if(m == 0) then
                is2%beg = 0
            
            else
                is2%beg = -buff_size
            end if
            
            is2%end = m - is2%beg



            if(p == 0) then
                is3%beg = 0
            
            else
                is3%beg = -buff_size
            end if
            is3%end = p - is3%beg

            allocate (q_prim_rsy_vf(0:buff_size, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1) then

                allocate (F_rsy_vf(0:buff_size, &
                                            is2%beg:is2%end, &
                                            is3%beg:is3%end, 1:adv_idx%end))


                allocate (F_src_rsy_vf(0:buff_size, &
                                                      is2%beg:is2%end, &
                                                      is3%beg:is3%end, adv_idx%beg:adv_idx%end))


            end if

            allocate (flux_rsy_vf(-1:buff_size, &
                                           is2%beg:is2%end, &
                                           is3%beg:is3%end, 1:adv_idx%end))

            allocate (flux_src_rsy_vf(-1:buff_size, &
                                                     is2%beg:is2%end, &
                                                     is3%beg:is3%end, adv_idx%beg:adv_idx%end))


        end if

        if(p > 0) then

            if(n == 0) then
                is2%beg = 0
            
            else
                is2%beg = -buff_size
            end if
            
            is2%end = n - is2%beg



            if(m == 0) then
                is3%beg = 0
            
            else
                is3%beg = -buff_size
            end if
            is3%end = m - is3%beg

            allocate (q_prim_rsz_vf(0:buff_size, &
                                             is2%beg:is2%end, &
                                             is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1) then

                allocate (F_rsz_vf(0:buff_size, &
                                            is2%beg:is2%end, &
                                            is3%beg:is3%end, 1:adv_idx%end))


                allocate (F_src_rsz_vf(0:buff_size, &
                                                      is2%beg:is2%end, &
                                                      is3%beg:is3%end, adv_idx%beg:adv_idx%end))


            end if

            allocate (flux_rsz_vf(-1:buff_size, &
                                           is2%beg:is2%end, &
                                           is3%beg:is3%end, 1:adv_idx%end))

            allocate (flux_src_rsz_vf(-1:buff_size, &
                                                     is2%beg:is2%end, &
                                                     is3%beg:is3%end, adv_idx%beg:adv_idx%end))

         end if       



        ! Allocating the cell-average partial densities, the velocity, the
        ! advected variables, the mass fractions, as well as Weber numbers
        allocate (alpha_rho(1:cont_idx%end))
        allocate (vel(1:num_dims))
        allocate (adv(1:adv_idx%end - E_idx))
        allocate (mf(1:cont_idx%end))

        ! Allocating the first-order spatial derivatives in the s-direction
        ! of the partial densities, the velocity and the advected variables
        allocate (dalpha_rho_ds(1:cont_idx%end))
        allocate (dvel_ds(1:num_dims))
        allocate (dadv_ds(1:adv_idx%end - E_idx))

        ! Allocating L, see Thompson (1987, 1990)
        allocate (L(1:adv_idx%end))

        ! Allocating the cell-width distribution in the s-direction
        allocate (ds(0:buff_size))

        allocate(gammas(1:num_fluids), pi_infs(1:num_fluids))


        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do

        !$acc update device(gammas, pi_infs)



        ! Allocating/Computing CBC Coefficients in x-direction =============
        if (all((/bc_x%beg, bc_x%end/) <= -5)) then

            allocate (fd_coef_x(0:buff_size, -1:1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
            end if

            call s_compute_cbc_coefficients(1, -1)
            call s_compute_cbc_coefficients(1, 1)

        elseif (bc_x%beg <= -5) then

            allocate (fd_coef_x(0:buff_size, -1:-1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
            end if

            call s_compute_cbc_coefficients(1, -1)

        elseif (bc_x%end <= -5) then

            allocate (fd_coef_x(0:buff_size, 1:1))

            if (weno_order > 1) then
                allocate (pi_coef_x(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
            end if

            call s_compute_cbc_coefficients(1, 1)

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in y-direction =============
        if (n > 0) then

            if (all((/bc_y%beg, bc_y%end/) <= -5)) then

                allocate (fd_coef_y(0:buff_size, -1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(2, -1)
                call s_compute_cbc_coefficients(2, 1)

            elseif (bc_y%beg <= -5) then

                allocate (fd_coef_y(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(2, -1)

            elseif (bc_y%end <= -5) then

                allocate (fd_coef_y(0:buff_size, 1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_y(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(2, 1)

            end if

        end if
        ! ==================================================================

        ! Allocating/Computing CBC Coefficients in z-direction =============
        if (p > 0) then

            if (all((/bc_z%beg, bc_z%end/) <= -5)) then

                allocate (fd_coef_z(0:buff_size, -1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:1))
                end if

                call s_compute_cbc_coefficients(3, -1)
                call s_compute_cbc_coefficients(3, 1)

            elseif (bc_z%beg <= -5) then

                allocate (fd_coef_z(0:buff_size, -1:-1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, -1:-1))
                end if

                call s_compute_cbc_coefficients(3, -1)

            elseif (bc_z%end <= -5) then

                allocate (fd_coef_z(0:buff_size, 1:1))

                if (weno_order > 1) then
                    allocate (pi_coef_z(0:weno_polyn - 1, 0:weno_order - 3, 1:1))
                end if

                call s_compute_cbc_coefficients(3, 1)

            end if

        end if
        ! ==================================================================

        !$acc update device(fd_coef_x, fd_coef_y, fd_coef_z, pi_coef_x, pi_coef_y, pi_coef_z)

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end
        bubxb = bub_idx%beg
        bubxe = bub_idx%end


        !$acc update device(momxb, momxe, advxb, advxe, contxb, contxe, bubxb, bubxe)

        bcxb = bc_x%beg
        bcxe = bc_x%end

        !$acc update device(bcxb, bcxe)

        if(n > 0) then
            bcyb = bc_y%beg
            bcye = bc_y%end

            !$acc update device(bcyb, bcye)
        end if

        if(p > 0) then
            bczb = bc_z%beg
            bcze = bc_z%end

            !$acc update device(bczb, bcze)
        end if

    end subroutine s_initialize_cbc_module ! -------------------------------

    subroutine s_compute_cbc_coefficients(cbc_dir, cbc_loc) ! --------------
        ! Description: The purpose of this subroutine is to compute the grid
        !              dependent FD and PI coefficients, or CBC coefficients,
        !              provided the CBC coordinate direction and location.

        ! CBC coordinate direction and location
        integer, intent(IN) :: cbc_dir, cbc_loc

        ! Cell-boundary locations in the s-direction
        real(kind(0d0)), dimension(0:buff_size + 1) :: s_cb

        ! Generic loop iterator
        integer :: i

        ! Associating CBC coefficients pointers
        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)

        ! Determining the cell-boundary locations in the s-direction
        s_cb(0) = 0d0

        do i = 0, buff_size
            s_cb(i + 1) = s_cb(i) + ds(i)
        end do

        ! Computing CBC1 Coefficients ======================================
        if(cbc_dir == 1) then
            if (weno_order == 1) then

                fd_coef_x(:, cbc_loc) = 0d0
                fd_coef_x(0, cbc_loc) = -2d0/(ds(0) + ds(1))
                fd_coef_x(1, cbc_loc) = -fd_coef_x(0, cbc_loc)

                ! ==================================================================

                ! Computing CBC2 Coefficients ======================================
            elseif (weno_order == 3) then

                fd_coef_x(:, cbc_loc) = 0d0
                fd_coef_x(0, cbc_loc) = -6d0/(3d0*ds(0) + 2d0*ds(1) - ds(2))
                fd_coef_x(1, cbc_loc) = -4d0*fd_coef_x(0, cbc_loc)/3d0
                fd_coef_x(2, cbc_loc) = fd_coef_x(0, cbc_loc)/3d0

                pi_coef_x(0, 0, cbc_loc) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

                ! ==================================================================

                ! Computing CBC4 Coefficients ======================================
            else

                fd_coef_x(:, cbc_loc) = 0d0
                fd_coef_x(0, cbc_loc) = -50d0/(25d0*ds(0) + 2d0*ds(1) &
                                             - 1d1*ds(2) + 1d1*ds(3) &
                                             - 3d0*ds(4))
                fd_coef_x(1, cbc_loc) = -48d0*fd_coef_x(0, cbc_loc)/25d0
                fd_coef_x(2, cbc_loc) = 36d0*fd_coef_x(0, cbc_loc)/25d0
                fd_coef_x(3, cbc_loc) = -16d0*fd_coef_x(0, cbc_loc)/25d0
                fd_coef_x(4, cbc_loc) = 3d0*fd_coef_x(0, cbc_loc)/25d0

                pi_coef_x(0, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                     (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
                pi_coef_x(0, 1, cbc_loc) = &
                    ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                      (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                           (s_cb(4) - s_cb(1)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_x(0, 2, cbc_loc) = &
                    (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                         (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                                (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                    ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                     (s_cb(0) - s_cb(4)))
                pi_coef_x(1, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                     (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
                pi_coef_x(1, 1, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                      (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                           (s_cb(2) - s_cb(4)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_x(1, 2, cbc_loc) = &
                    ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                     (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                           (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

            end if
        elseif(cbc_dir == 2) then
            if (weno_order == 1) then

                fd_coef_y(:, cbc_loc) = 0d0
                fd_coef_y(0, cbc_loc) = -2d0/(ds(0) + ds(1))
                fd_coef_y(1, cbc_loc) = -fd_coef_y(0, cbc_loc)

                ! ==================================================================

                ! Computing CBC2 Coefficients ======================================
            elseif (weno_order == 3) then

                fd_coef_y(:, cbc_loc) = 0d0
                fd_coef_y(0, cbc_loc) = -6d0/(3d0*ds(0) + 2d0*ds(1) - ds(2))
                fd_coef_y(1, cbc_loc) = -4d0*fd_coef_y(0, cbc_loc)/3d0
                fd_coef_y(2, cbc_loc) = fd_coef_y(0, cbc_loc)/3d0

                pi_coef_y(0, 0, cbc_loc) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

                ! ==================================================================

                ! Computing CBC4 Coefficients ======================================
            else

                fd_coef_y(:, cbc_loc) = 0d0
                fd_coef_y(0, cbc_loc) = -50d0/(25d0*ds(0) + 2d0*ds(1) &
                                             - 1d1*ds(2) + 1d1*ds(3) &
                                             - 3d0*ds(4))
                fd_coef_y(1, cbc_loc) = -48d0*fd_coef_y(0, cbc_loc)/25d0
                fd_coef_y(2, cbc_loc) = 36d0*fd_coef_y(0, cbc_loc)/25d0
                fd_coef_y(3, cbc_loc) = -16d0*fd_coef_y(0, cbc_loc)/25d0
                fd_coef_y(4, cbc_loc) = 3d0*fd_coef_y(0, cbc_loc)/25d0

                pi_coef_y(0, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                     (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
                pi_coef_y(0, 1, cbc_loc) = &
                    ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                      (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                           (s_cb(4) - s_cb(1)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_y(0, 2, cbc_loc) = &
                    (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                         (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                                (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                    ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                     (s_cb(0) - s_cb(4)))
                pi_coef_y(1, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                     (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
                pi_coef_y(1, 1, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                      (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                           (s_cb(2) - s_cb(4)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_y(1, 2, cbc_loc) = &
                    ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                     (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                           (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

            end if
        else
            if (weno_order == 1) then

                fd_coef_z(:, cbc_loc) = 0d0
                fd_coef_z(0, cbc_loc) = -2d0/(ds(0) + ds(1))
                fd_coef_z(1, cbc_loc) = -fd_coef_z(0, cbc_loc)

                ! ==================================================================

                ! Computing CBC2 Coefficients ======================================
            elseif (weno_order == 3) then

                fd_coef_z(:, cbc_loc) = 0d0
                fd_coef_z(0, cbc_loc) = -6d0/(3d0*ds(0) + 2d0*ds(1) - ds(2))
                fd_coef_z(1, cbc_loc) = -4d0*fd_coef_z(0, cbc_loc)/3d0
                fd_coef_z(2, cbc_loc) = fd_coef_z(0, cbc_loc)/3d0

                pi_coef_z(0, 0, cbc_loc) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

                ! ==================================================================

                ! Computing CBC4 Coefficients ======================================
            else

                fd_coef_z(:, cbc_loc) = 0d0
                fd_coef_z(0, cbc_loc) = -50d0/(25d0*ds(0) + 2d0*ds(1) &
                                             - 1d1*ds(2) + 1d1*ds(3) &
                                             - 3d0*ds(4))
                fd_coef_z(1, cbc_loc) = -48d0*fd_coef_z(0, cbc_loc)/25d0
                fd_coef_z(2, cbc_loc) = 36d0*fd_coef_z(0, cbc_loc)/25d0
                fd_coef_z(3, cbc_loc) = -16d0*fd_coef_z(0, cbc_loc)/25d0
                fd_coef_z(4, cbc_loc) = 3d0*fd_coef_z(0, cbc_loc)/25d0

                pi_coef_z(0, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                     (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
                pi_coef_z(0, 1, cbc_loc) = &
                    ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                      (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                           (s_cb(4) - s_cb(1)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_z(0, 2, cbc_loc) = &
                    (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                         (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                                (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                    ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                     (s_cb(0) - s_cb(4)))
                pi_coef_z(1, 0, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                     (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                           (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
                pi_coef_z(1, 1, cbc_loc) = &
                    ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                     ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                      (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                           (s_cb(2) - s_cb(4)))))/ &
                    ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                     (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                pi_coef_z(1, 2, cbc_loc) = &
                    ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                     (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                           (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

            end if
        end if
        ! END: Computing CBC4 Coefficients =================================

        ! Nullifying CBC coefficients


    end subroutine s_compute_cbc_coefficients ! ----------------------------

        !!  The goal of the procedure is to associate the FD and PI
        !!      coefficients, or CBC coefficients, with the appropriate
        !!      targets, based on the coordinate direction and location
        !!      of the CBC.
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
    subroutine s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc) ! ---

        integer, intent(IN) :: cbc_dir, cbc_loc

        integer :: i !< Generic loop iterator

        ! Associating CBC Coefficients in x-direction ======================
        if (cbc_dir == 1) then

            !fd_coef => fd_coef_x; if (weno_order > 1) pi_coef => pi_coef_x

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dx(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dx(m - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in y-direction ======================
        elseif (cbc_dir == 2) then

            !fd_coef => fd_coef_y; if (weno_order > 1) pi_coef => pi_coef_y

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dy(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dy(n - i)
                end do
            end if
            ! ==================================================================

            ! Associating CBC Coefficients in z-direction ======================
        else

            !fd_coef => fd_coef_z; if (weno_order > 1) pi_coef => pi_coef_z

            if (cbc_loc == -1) then
                do i = 0, buff_size
                    ds(i) = dz(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dz(p - i)
                end do
            end if

        end if

        !$acc update device(ds)

        ! ==================================================================

    end subroutine s_associate_cbc_coefficients_pointers ! -----------------

    !>  The following is the implementation of the CBC based on
        !!      the work of Thompson (1987, 1990) on hyperbolic systems.
        !!      The CBC is indirectly applied in the computation of the
        !!      right-hand-side (RHS) near the relevant domain boundary
        !!      through the modification of the fluxes.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_cbc(q_prim_vf, flux_vf, flux_src_vf, & ! -----------------
                     cbc_dir, cbc_loc, &
                     ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! First-order time derivatives of the partial densities, density,
        ! velocity, pressure, advection variables, and the specific heat
        ! ratio and liquid stiffness functions
        real(kind(0d0)), dimension(contxe)   :: dalpha_rho_dt
        real(kind(0d0))                            ::       drho_dt
        real(kind(0d0)), dimension(num_dims)       ::       dvel_dt
        real(kind(0d0))                            ::      dpres_dt
        real(kind(0d0)), dimension(advxe - E_idx) ::       dadv_dt
        real(kind(0d0))                            ::     dgamma_dt
        real(kind(0d0))                            ::    dpi_inf_dt
        real(kind(0d0)), dimension(contxe)   :: alpha_rho, dalpha_rho_ds, mf
        real(kind(0d0)), dimension(2) :: Re_cbc
        real(kind(0d0)), dimension(1:num_dims)     :: vel, dvel_ds
        real(kind(0d0)), dimension(1:advxe-E_idx) :: adv, dadv_ds
        real(kind(0d0)), dimension(1:advxe)       :: L
        real(kind(0d0)), dimension(3)             ::lambda  


        real(kind(0d0))                              :: rho         !< Cell averaged density
        real(kind(0d0))                              :: pres        !< Cell averaged pressure
        real(kind(0d0))                              :: E           !< Cell averaged energy
        real(kind(0d0))                              :: H           !< Cell averaged enthalpy
        real(kind(0d0))                              :: gamma       !< Cell averaged specific heat ratio
        real(kind(0d0))                              :: pi_inf      !< Cell averaged liquid stiffness
        real(kind(0d0))                              :: c

        real(kind(0d0)) :: vel_K_sum, vel_dv_dt_sum

        integer :: i, j, k, r !< Generic loop iterators

        real(kind(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for Wood mixture sound speed

        ! Reshaping of inputted data and association of the FD and PI
        ! coefficients, or CBC coefficients, respectively, hinging on
        ! selected CBC coordinate direction


        ! Allocating L, see Thompson (1987, 1990)



        call s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                              cbc_dir, cbc_loc, &
                              ix, iy, iz) 

        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)       

        if(cbc_dir == 1) then

        ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
            if (weno_order == 3) then

                call s_convert_primitive_to_flux_variables(q_prim_rsx_vf, &
                                                           F_rsx_vf, &
                                                           F_src_rsx_vf, &
                                                           is1, is2, is3, starty, startz)

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            flux_rsx_vf(0, k, r, i) = F_rsx_vf(0, k, r, i) &
                                                        + pi_coef_x(0, 0, cbc_loc)* &
                                                        (F_rsx_vf(1, k, r, i) - &
                                                         F_rsx_vf(0, k, r, i))
                        end do
                    end do
                end do

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    flux_src_rsx_vf(0, k, r, i) = F_src_rsx_vf(0, k, r, i) + &
                                                    (F_src_rsx_vf(1, k, r, i) - &
                                                     F_src_rsx_vf(0, k, r, i)) &
                                                    *pi_coef_x(0, 0, cbc_loc)
                        end do
                    end do
                end do
                ! ==================================================================

                ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
            elseif (weno_order == 5) then

                call s_convert_primitive_to_flux_variables(q_prim_rsx_vf, &
                                                           F_rsx_vf, &
                                                           F_src_rsx_vf, &
                                                           is1, is2, is3, starty, startz)


    !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_rsx_vf(j, k, r, i) = F_rsx_vf(j, k, r, i) &
                                                    + pi_coef_x(j, 0, cbc_loc)* &
                                                    (F_rsx_vf(3, k, r, i) - &
                                                     F_rsx_vf(2, k, r, i)) &
                                                    + pi_coef_x(j, 1, cbc_loc)* &
                                                    (F_rsx_vf(2, k, r, i) - &
                                                     F_rsx_vf(1, k, r, i)) &
                                                    + pi_coef_x(j, 2, cbc_loc)* &
                                                    (F_rsx_vf(1, k, r, i) - &
                                                     F_rsx_vf(0, k, r, i))
                            end do
                        end do
                    end do
                end do

    !$acc parallel loop collapse(4) gang vector default(present) 
                do i = advxb, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_src_rsx_vf(j, k, r, i) = F_src_rsx_vf(j, k, r, i) + &
                                                        (F_src_rsx_vf(3, k, r, i) - &
                                                         F_src_rsx_vf(2, k, r, i)) &
                                                        *pi_coef_x(j, 0, cbc_loc) + &
                                                        (F_src_rsx_vf(2, k, r, i) - &
                                                         F_src_rsx_vf(1, k, r, i)) &
                                                        *pi_coef_x(j, 1, cbc_loc) + &
                                                        (F_src_rsx_vf(1, k, r, i) - &
                                                         F_src_rsx_vf(0, k, r, i)) &
                                                        *pi_coef_x(j, 2, cbc_loc)
                            end do
                        end do
                    end do
                end do

            end if
            ! ==================================================================

            ! FD2 or FD4 of RHS at j = 0 =======================================
    !$acc parallel loop collapse(2) gang vector default(present) private(alpha_rho, vel, adv, mf, dvel_ds, dadv_ds, Re_cbc, dalpha_rho_ds,dvel_dt, dadv_dt, dalpha_rho_dt,L, lambda)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end

                    ! Transferring the Primitive Variables =======================
    !$acc loop seq
                    do i = 1, contxe
                        alpha_rho(i) = q_prim_rsx_vf(0, k, r, i)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        vel(i) = q_prim_rsx_vf(0, k, r, contxe + i)
                    end do

                    vel_K_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_K_sum = vel_K_sum + vel(i)**2d0
                    end do

                    pres = q_prim_rsx_vf(0, k, r, E_idx)

    !$acc loop seq
                    do i = 1, advxe - E_idx
                        adv(i) = q_prim_rsx_vf(0, k, r, E_idx + i)
                    end do

                    if(bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, adv, alpha_rho, 0, k, r)

                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, adv, alpha_rho, Re_cbc, 0, k, r)
                    end if

                    E = gamma*pres + pi_inf + 5d-1*rho*vel_K_sum

                    H = (E + pres)/rho

    !$acc loop seq
                    do i = 1, contxe
                        mf(i)  = alpha_rho(i)/rho
                    end do

                    ! Compute mixture sound speed
                    if (alt_soundspeed) then
                        blkmod1 = ((gammas(1) + 1d0)*pres + &
                                   pi_infs(1))/gammas(1)
                        blkmod2 = ((gammas(2) + 1d0)*pres + &
                                   pi_infs(2))/gammas(2)
                        c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
                    elseif (model_eqns == 3) then
                        c = 0d0
    !$acc loop seq
                        do i = 1, num_fluids
                            c = c + q_prim_rsx_vf(0, k, r, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                (pres + pi_infs(i)/(gammas(i) + 1d0))
                        end do
                        c = c/rho
                    else
                        c = ((H - 5d-1*vel_K_sum)/gamma)
                    end if

                    c = sqrt(c)

    !                  IF (mixture_err .AND. c < 0d0) THEN
    !                    c = sgm_eps
    !                  ELSE
    !                    c = SQRT(c)
    !                  END IF

                    ! ============================================================

                    ! First-Order Spatial Derivatives of Primitive Variables =====

    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = 0d0
                        end do

    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = 0d0
                        end do

                        dpres_ds = 0d0
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = 0d0
                        end do

    !$acc loop seq
                    do j = 0, buff_size

    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = q_prim_rsx_vf(j, k, r, i)* &
                                               fd_coef_x(j, cbc_loc) + &
                                               dalpha_rho_ds(i)
                        end do
    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = q_prim_rsx_vf(j, k, r, contxe + i)* &
                                         fd_coef_x(j, cbc_loc) + &
                                         dvel_ds(i)
                        end do

                        dpres_ds = q_prim_rsx_vf(j, k, r, E_idx)* &
                                   fd_coef_x(j, cbc_loc) + &
                                   dpres_ds
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = q_prim_rsx_vf(j, k, r, E_idx + i)* &
                                         fd_coef_x(j, cbc_loc) + &
                                         dadv_ds(i)
                        end do
                    end do
                    ! ============================================================

                    ! First-Order Temporal Derivatives of Primitive Variables ====
                    lambda(1) = vel(dir_idx(1)) - c
                    lambda(2) = vel(dir_idx(1))
                    lambda(3) = vel(dir_idx(1)) + c

                   ! call s_compute_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------

                    if((cbc_loc == -1 .and. bcxb == -5) .or. (cbc_loc == 1 .and. bcxe == -5)) then
                        call s_compute_slip_wall_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -6) .or. (cbc_loc == 1 .and. bcxe == -6)) then
                        call s_compute_nonreflecting_subsonic_buffer_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -7) .or. (cbc_loc == 1 .and. bcxe == -7)) then
                        call s_compute_nonreflecting_subsonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -8) .or. (cbc_loc == 1 .and. bcxe == -8)) then
                        call s_compute_nonreflecting_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -9) .or. (cbc_loc == 1 .and. bcxe == -9)) then
                        call s_compute_force_free_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -10) .or. (cbc_loc == 1 .and. bcxe == -10)) then
                        call s_compute_constant_pressure_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcxb == -11) .or. (cbc_loc == 1 .and. bcxe == -11)) then  
                        call s_compute_supersonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else  
                        call s_compute_supersonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    end if


                    ! Be careful about the cylindrical coordinate!
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                        dpres_dt = -5d-1*(L(advxe) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                   /y_cc(n)
                    else
                        dpres_dt = -5d-1*(L(advxe) + L(1))
                    end if

    !$acc loop seq
                    do i = 1, contxe
                        dalpha_rho_dt(i) = &
                            -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                              (L(1) - L(advxe))/(2d0*rho*c) + &
                                              (dir_flg(dir_idx(i)) - 1d0)* &
                                              L(momxb + i)
                    end do

                    vel_dv_dt_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_dv_dt_sum = vel_dv_dt_sum +vel(i)*dvel_dt(i)
                    end do 

                    ! The treatment of void fraction source is unclear
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                        end do
                    else
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i)
                        end do
                    end if

                    drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0

                    if (model_eqns == 1) then
                        drho_dt = dalpha_rho_dt(1)
                        dgamma_dt = dadv_dt(1)
                        dpi_inf_dt = dadv_dt(2)
                    else
    !$acc loop seq
                        do i = 1, num_fluids
                            drho_dt = drho_dt + dalpha_rho_dt(i)
                            dgamma_dt = dgamma_dt + dadv_dt(i)*gammas(i)
                            dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_infs(i)
                        end do
                    end if
                    ! ============================================================

                    ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
    !$acc loop seq
                    do i = 1, contxe
                        flux_rsx_vf(-1, k, r, i) = flux_rsx_vf(0, k, r, i) &
                                                     + ds(0)*dalpha_rho_dt(i)
                    end do

    !$acc loop seq
                    do i = momxb, momxe
                        flux_rsx_vf(-1, k, r, i) = flux_rsx_vf(0, k, r, i) &
                                                     + ds(0)*(vel(i - contxe)*drho_dt &
                                                              + rho*dvel_dt(i - contxe))
                    end do

                    flux_rsx_vf(-1, k, r, E_idx) = flux_rsx_vf(0, k, r, E_idx) &
                                                     + ds(0)*(pres*dgamma_dt &
                                                              + gamma*dpres_dt &
                                                              + dpi_inf_dt &
                                                              + rho*vel_dv_dt_sum &
                                                              + 5d-1*drho_dt*vel_K_sum)

                    if (riemann_solver == 1) then
    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsx_vf(-1, k, r, i) = 0d0
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsx_vf(-1, k, r, i) = &
                                1d0/max(abs(vel(dir_idx(1))), sgm_eps) &
                                *sign(1d0, vel(dir_idx(1))) &
                                *(flux_rsx_vf(0, k, r, i) &
                                  + vel(dir_idx(1)) &
                                  *flux_src_rsx_vf(0, k, r, i) &
                                  + ds(0)*dadv_dt(i - E_idx))
                        end do

                    else

    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsx_vf(-1, k, r, i) = flux_rsx_vf(0, k, r, i) - &
                                                         adv(i - E_idx)*flux_src_rsx_vf(0, k, r, i) + &
                                                         ds(0)*dadv_dt(i - E_idx)
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsx_vf(-1, k, r, i) = 0d0
                        end do

                    end if
                    ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============

                end do
            end do

        else if(cbc_dir == 2) then
        ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
            if (weno_order == 3) then

                call s_convert_primitive_to_flux_variables(q_prim_rsy_vf, &
                                                           F_rsy_vf, &
                                                           F_src_rsy_vf, &
                                                           is1, is2, is3, startx, startz)

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            flux_rsy_vf(0, k, r, i) = F_rsy_vf(0, k, r, i) &
                                                        + pi_coef_y(0, 0, cbc_loc)* &
                                                        (F_rsy_vf(1, k, r, i) - &
                                                         F_rsy_vf(0, k, r, i))
                        end do
                    end do
                end do

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    flux_src_rsy_vf(0, k, r, i) = F_src_rsy_vf(0, k, r, i) + &
                                                    (F_src_rsy_vf(1, k, r, i) - &
                                                     F_src_rsy_vf(0, k, r, i)) &
                                                    *pi_coef_y(0, 0, cbc_loc)
                        end do
                    end do
                end do
                ! ==================================================================

                ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
            elseif (weno_order == 5) then

                call s_convert_primitive_to_flux_variables(q_prim_rsy_vf, &
                                                           F_rsy_vf, &
                                                           F_src_rsy_vf, &
                                                           is1, is2, is3, startx, startz)


    !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_rsy_vf(j, k, r, i) = F_rsy_vf(j, k, r, i) &
                                                    + pi_coef_y(j, 0, cbc_loc)* &
                                                    (F_rsy_vf(3, k, r, i) - &
                                                     F_rsy_vf(2, k, r, i)) &
                                                    + pi_coef_y(j, 1, cbc_loc)* &
                                                    (F_rsy_vf(2, k, r, i) - &
                                                     F_rsy_vf(1, k, r, i)) &
                                                    + pi_coef_y(j, 2, cbc_loc)* &
                                                    (F_rsy_vf(1, k, r, i) - &
                                                     F_rsy_vf(0, k, r, i))
                            end do
                        end do
                    end do
                end do

    !$acc parallel loop collapse(4) gang vector default(present) 
                do i = advxb, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_src_rsy_vf(j, k, r, i) = F_src_rsy_vf(j, k, r, i) + &
                                                        (F_src_rsy_vf(3, k, r, i) - &
                                                         F_src_rsy_vf(2, k, r, i)) &
                                                        *pi_coef_y(j, 0, cbc_loc) + &
                                                        (F_src_rsy_vf(2, k, r, i) - &
                                                         F_src_rsy_vf(1, k, r, i)) &
                                                        *pi_coef_y(j, 1, cbc_loc) + &
                                                        (F_src_rsy_vf(1, k, r, i) - &
                                                         F_src_rsy_vf(0, k, r, i)) &
                                                        *pi_coef_y(j, 2, cbc_loc)
                            end do
                        end do
                    end do
                end do

            end if
            ! ==================================================================

            ! FD2 or FD4 of RHS at j = 0 =======================================
    !$acc parallel loop collapse(2) gang vector default(present) private(alpha_rho, vel, adv, mf, dvel_ds, dadv_ds, Re_cbc, dalpha_rho_ds,dvel_dt, dadv_dt, dalpha_rho_dt,L, lambda)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end

                    ! Transferring the Primitive Variables =======================
    !$acc loop seq
                    do i = 1, contxe
                        alpha_rho(i) = q_prim_rsy_vf(0, k, r, i)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        vel(i) = q_prim_rsy_vf(0, k, r, contxe + i)
                    end do

                    vel_K_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_K_sum = vel_K_sum + vel(i)**2d0
                    end do

                    pres = q_prim_rsy_vf(0, k, r, E_idx)

    !$acc loop seq
                    do i = 1, advxe - E_idx
                        adv(i) = q_prim_rsy_vf(0, k, r, E_idx + i)
                    end do

                    if(bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, adv, alpha_rho, 0, k, r)

                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, adv, alpha_rho, Re_cbc, 0, k, r)
                    end if

                    E = gamma*pres + pi_inf + 5d-1*rho*vel_K_sum

                    H = (E + pres)/rho



    !$acc loop seq
                    do i = 1, contxe
                        mf(i)  = alpha_rho(i)/rho
                    end do

                    ! Compute mixture sound speed
                    if (alt_soundspeed) then
                        blkmod1 = ((gammas(1) + 1d0)*pres + &
                                   pi_infs(1))/gammas(1)
                        blkmod2 = ((gammas(2) + 1d0)*pres + &
                                   pi_infs(2))/gammas(2)
                        c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
                    elseif (model_eqns == 3) then
                        c = 0d0
    !$acc loop seq
                        do i = 1, num_fluids
                            c = c + q_prim_rsy_vf(0, k, r, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                (pres + pi_infs(i)/(gammas(i) + 1d0))
                        end do
                        c = c/rho
                    else
                        c = ((H - 5d-1*vel_K_sum)/gamma)
                    end if

                    c = sqrt(c)

    !                  IF (mixture_err .AND. c < 0d0) THEN
    !                    c = sgm_eps
    !                  ELSE
    !                    c = SQRT(c)
    !                  END IF

                    ! ============================================================

                    ! First-Order Spatial Derivatives of Primitive Variables =====

                    
    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = 0d0
                        end do

    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = 0d0
                        end do

                        dpres_ds = 0d0
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = 0d0
                        end do

    !$acc loop seq
                    do j = 0, buff_size

    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = q_prim_rsy_vf(j, k, r, i)* &
                                               fd_coef_y(j, cbc_loc) + &
                                               dalpha_rho_ds(i)
                        end do
    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = q_prim_rsy_vf(j, k, r, contxe + i)* &
                                         fd_coef_y(j, cbc_loc) + &
                                         dvel_ds(i)
                        end do

                        dpres_ds = q_prim_rsy_vf(j, k, r, E_idx)* &
                                   fd_coef_y(j, cbc_loc) + &
                                   dpres_ds
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = q_prim_rsy_vf(j, k, r, E_idx + i)* &
                                         fd_coef_y(j, cbc_loc) + &
                                         dadv_ds(i)
                        end do

                    end do
                    ! ============================================================

                    ! First-Order Temporal Derivatives of Primitive Variables ====
                    lambda(1) = vel(dir_idx(1)) - c
                    lambda(2) = vel(dir_idx(1))
                    lambda(3) = vel(dir_idx(1)) + c

                 !   call s_compute_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------

                    if((cbc_loc == -1 .and. bcyb == -5) .or. (cbc_loc == 1 .and. bcye == -5)) then
                        call s_compute_slip_wall_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -6) .or. (cbc_loc == 1 .and. bcye == -6)) then
                        call s_compute_nonreflecting_subsonic_buffer_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -7) .or. (cbc_loc == 1 .and. bcye == -7)) then
                        call s_compute_nonreflecting_subsonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -8) .or. (cbc_loc == 1 .and. bcye == -8)) then
                        call s_compute_nonreflecting_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -9) .or. (cbc_loc == 1 .and. bcye == -9)) then
                        call s_compute_force_free_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -10) .or. (cbc_loc == 1 .and. bcye == -10)) then
                        call s_compute_constant_pressure_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bcyb == -11) .or. (cbc_loc == 1 .and. bcye == -11)) then  
                        call s_compute_supersonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else  
                        call s_compute_supersonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    end if


                    ! Be careful about the cylindrical coordinate!
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                        dpres_dt = -5d-1*(L(advxe) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                   /y_cc(n)
                    else
                        dpres_dt = -5d-1*(L(advxe) + L(1))
                    end if

    !$acc loop seq
                    do i = 1, contxe
                        dalpha_rho_dt(i) = &
                            -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                              (L(1) - L(advxe))/(2d0*rho*c) + &
                                              (dir_flg(dir_idx(i)) - 1d0)* &
                                              L(momxb + i)
                    end do

                    vel_dv_dt_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_dv_dt_sum = vel_dv_dt_sum +vel(i)*dvel_dt(i)
                    end do 

                    ! The treatment of void fraction source is unclear
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                        end do
                    else
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i)
                        end do
                    end if

                    drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0

                    if (model_eqns == 1) then
                        drho_dt = dalpha_rho_dt(1)
                        dgamma_dt = dadv_dt(1)
                        dpi_inf_dt = dadv_dt(2)
                    else
    !$acc loop seq
                        do i = 1, num_fluids
                            drho_dt = drho_dt + dalpha_rho_dt(i)
                            dgamma_dt = dgamma_dt + dadv_dt(i)*gammas(i)
                            dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_infs(i)
                        end do
                    end if
                    ! ============================================================

                    ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
    !$acc loop seq
                    do i = 1, contxe
                        flux_rsy_vf(-1, k, r, i) = flux_rsy_vf(0, k, r, i) &
                                                     + ds(0)*dalpha_rho_dt(i)
                    end do

    !$acc loop seq
                    do i = momxb, momxe
                        flux_rsy_vf(-1, k, r, i) = flux_rsy_vf(0, k, r, i) &
                                                     + ds(0)*(vel(i - contxe)*drho_dt &
                                                              + rho*dvel_dt(i - contxe))
                    end do

                    flux_rsy_vf(-1, k, r, E_idx) = flux_rsy_vf(0, k, r, E_idx) &
                                                     + ds(0)*(pres*dgamma_dt &
                                                              + gamma*dpres_dt &
                                                              + dpi_inf_dt &
                                                              + rho*vel_dv_dt_sum &
                                                              + 5d-1*drho_dt*vel_K_sum)

                    if (riemann_solver == 1) then
    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsy_vf(-1, k, r, i) = 0d0
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsy_vf(-1, k, r, i) = &
                                1d0/max(abs(vel(dir_idx(1))), sgm_eps) &
                                *sign(1d0, vel(dir_idx(1))) &
                                *(flux_rsy_vf(0, k, r, i) &
                                  + vel(dir_idx(1)) &
                                  *flux_src_rsy_vf(0, k, r, i) &
                                  + ds(0)*dadv_dt(i - E_idx))
                        end do

                    else

    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsy_vf(-1, k, r, i) = flux_rsy_vf(0, k, r, i) - &
                                                         adv(i - E_idx)*flux_src_rsy_vf(0, k, r, i) + &
                                                         ds(0)*dadv_dt(i - E_idx)
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsy_vf(-1, k, r, i) = 0d0
                        end do

                    end if
                    ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============

                end do
            end do

        else
        ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2 ==================
            if (weno_order == 3) then

                call s_convert_primitive_to_flux_variables(q_prim_rsz_vf, &
                                                           F_rsz_vf, &
                                                           F_src_rsz_vf, &
                                                           is1, is2, is3, starty, startx)

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = 1, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            flux_rsz_vf(0, k, r, i) = F_rsz_vf(0, k, r, i) &
                                                        + pi_coef_z(0, 0, cbc_loc)* &
                                                        (F_rsz_vf(1, k, r, i) - &
                                                         F_rsz_vf(0, k, r, i))
                        end do
                    end do
                end do

    !$acc parallel loop collapse(3) gang vector default(present)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                    flux_src_rsz_vf(0, k, r, i) = F_src_rsz_vf(0, k, r, i) + &
                                                    (F_src_rsz_vf(1, k, r, i) - &
                                                     F_src_rsz_vf(0, k, r, i)) &
                                                    *pi_coef_z(0, 0, cbc_loc)
                        end do
                    end do
                end do
                ! ==================================================================

                ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2 =============
            elseif (weno_order == 5) then

                call s_convert_primitive_to_flux_variables(q_prim_rsz_vf, &
                                                           F_rsz_vf, &
                                                           F_src_rsz_vf, &
                                                           is1, is2, is3, starty, startx)


    !$acc parallel loop collapse(4) gang vector default(present)
                do i = 1, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_rsz_vf(j, k, r, i) = F_rsz_vf(j, k, r, i) &
                                                    + pi_coef_z(j, 0, cbc_loc)* &
                                                    (F_rsz_vf(3, k, r, i) - &
                                                     F_rsz_vf(2, k, r, i)) &
                                                    + pi_coef_z(j, 1, cbc_loc)* &
                                                    (F_rsz_vf(2, k, r, i) - &
                                                     F_rsz_vf(1, k, r, i)) &
                                                    + pi_coef_z(j, 2, cbc_loc)* &
                                                    (F_rsz_vf(1, k, r, i) - &
                                                     F_rsz_vf(0, k, r, i))
                            end do
                        end do
                    end do
                end do

    !$acc parallel loop collapse(4) gang vector default(present) 
                do i = advxb, advxe
                    do j = 0, 1
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                        flux_src_rsz_vf(j, k, r, i) = F_src_rsz_vf(j, k, r, i) + &
                                                        (F_src_rsz_vf(3, k, r, i) - &
                                                         F_src_rsz_vf(2, k, r, i)) &
                                                        *pi_coef_z(j, 0, cbc_loc) + &
                                                        (F_src_rsz_vf(2, k, r, i) - &
                                                         F_src_rsz_vf(1, k, r, i)) &
                                                        *pi_coef_z(j, 1, cbc_loc) + &
                                                        (F_src_rsz_vf(1, k, r, i) - &
                                                         F_src_rsz_vf(0, k, r, i)) &
                                                        *pi_coef_z(j, 2, cbc_loc)
                            end do
                        end do
                    end do
                end do

            end if
            ! ==================================================================

            ! FD2 or FD4 of RHS at j = 0 =======================================
    !$acc parallel loop collapse(2) gang vector default(present) private(alpha_rho, vel, adv, mf, dvel_ds, dadv_ds, Re_cbc, dalpha_rho_ds,dvel_dt, dadv_dt, dalpha_rho_dt,L, lambda)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end

                    ! Transferring the Primitive Variables =======================
    !$acc loop seq
                    do i = 1, contxe
                        alpha_rho(i) = q_prim_rsz_vf(0, k, r, i)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        vel(i) = q_prim_rsz_vf(0, k, r, contxe + i)
                    end do

                    vel_K_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_K_sum = vel_K_sum + vel(i)**2d0
                    end do

                    pres = q_prim_rsz_vf(0, k, r, E_idx)

    !$acc loop seq
                    do i = 1, advxe - E_idx
                        adv(i) = q_prim_rsz_vf(0, k, r, E_idx + i)
                    end do

                    if(bubbles) then
                        call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, adv, alpha_rho, 0, k, r)

                    else
                        call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, adv, alpha_rho, Re_cbc, 0, k, r)
                    end if

                    E = gamma*pres + pi_inf + 5d-1*rho*vel_K_sum

                    H = (E + pres)/rho



    !$acc loop seq
                    do i = 1, contxe
                        mf(i)  = alpha_rho(i)/rho
                    end do

                    ! Compute mixture sound speed
                    if (alt_soundspeed) then
                        blkmod1 = ((gammas(1) + 1d0)*pres + &
                                   pi_infs(1))/gammas(1)
                        blkmod2 = ((gammas(2) + 1d0)*pres + &
                                   pi_infs(2))/gammas(2)
                        c = (1d0/(rho*(adv(1)/blkmod1 + adv(2)/blkmod2)))
                    elseif (model_eqns == 3) then
                        c = 0d0
    !$acc loop seq
                        do i = 1, num_fluids
                            c = c + q_prim_rsz_vf(0, k, r, i + advxb - 1)*(1d0/gammas(i) + 1d0)* &
                                (pres + pi_infs(i)/(gammas(i) + 1d0))
                        end do
                        c = c/rho
                    else
                        c = ((H - 5d-1*vel_K_sum)/gamma)
                    end if

                    c = sqrt(c)

    !                  IF (mixture_err .AND. c < 0d0) THEN
    !                    c = sgm_eps
    !                  ELSE
    !                    c = SQRT(c)
    !                  END IF

                    ! ============================================================

                    ! First-Order Spatial Derivatives of Primitive Variables =====

    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = 0d0
                        end do

    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = 0d0
                        end do

                        dpres_ds = 0d0
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = 0d0
                        end do

    !$acc loop seq
                    do j = 0, buff_size

    !$acc loop seq
                        do i = 1, contxe
                            dalpha_rho_ds(i) = q_prim_rsz_vf(j, k, r, i)* &
                                               fd_coef_z(j, cbc_loc) + &
                                               dalpha_rho_ds(i)
                        end do
    !$acc loop seq
                        do i = 1, num_dims
                            dvel_ds(i) = q_prim_rsz_vf(j, k, r, contxe + i)* &
                                         fd_coef_z(j, cbc_loc) + &
                                         dvel_ds(i)
                        end do

                        dpres_ds = q_prim_rsz_vf(j, k, r, E_idx)* &
                                   fd_coef_z(j, cbc_loc) + &
                                   dpres_ds
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = q_prim_rsz_vf(j, k, r, E_idx + i)* &
                                         fd_coef_z(j, cbc_loc) + &
                                         dadv_ds(i)
                        end do

                    end do
                    ! ============================================================

                    ! First-Order Temporal Derivatives of Primitive Variables ====
                    lambda(1) = vel(dir_idx(1)) - c
                    lambda(2) = vel(dir_idx(1))
                    lambda(3) = vel(dir_idx(1)) + c

                   ! call s_compute_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------

                    if((cbc_loc == -1 .and. bczb == -5) .or. (cbc_loc == 1 .and. bcze == -5)) then
                        call s_compute_slip_wall_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -6) .or. (cbc_loc == 1 .and. bcze == -6)) then
                        call s_compute_nonreflecting_subsonic_buffer_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -7) .or. (cbc_loc == 1 .and. bcze == -7)) then
                        call s_compute_nonreflecting_subsonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -8) .or. (cbc_loc == 1 .and. bcze == -8)) then
                        call s_compute_nonreflecting_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -9) .or. (cbc_loc == 1 .and. bcze == -9)) then
                        call s_compute_force_free_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -10) .or. (cbc_loc == 1 .and. bcze == -10)) then
                        call s_compute_constant_pressure_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else if((cbc_loc == -1 .and. bczb == -11) .or. (cbc_loc == 1 .and. bcze == -11)) then  
                        call s_compute_supersonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    else  
                        call s_compute_supersonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
                    end if


                    ! Be careful about the cylindrical coordinate!
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                        dpres_dt = -5d-1*(L(advxe) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                   /y_cc(n)
                    else
                        dpres_dt = -5d-1*(L(advxe) + L(1))
                    end if

    !$acc loop seq
                    do i = 1, contxe
                        dalpha_rho_dt(i) = &
                            -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                    end do

    !$acc loop seq
                    do i = 1, num_dims
                        dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                              (L(1) - L(advxe))/(2d0*rho*c) + &
                                              (dir_flg(dir_idx(i)) - 1d0)* &
                                              L(momxb + i)
                    end do

                    vel_dv_dt_sum = 0d0
    !$acc loop seq
                    do i = 1, num_dims
                        vel_dv_dt_sum = vel_dv_dt_sum +vel(i)*dvel_dt(i)
                    end do 

                    ! The treatment of void fraction source is unclear
                    if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i) !+ adv(i) * vel(dir_idx(1))/y_cc(n)
                        end do
                    else
    !$acc loop seq
                        do i = 1, advxe - E_idx
                            dadv_dt(i) = -L(momxe + i)
                        end do
                    end if

                    drho_dt = 0d0; dgamma_dt = 0d0; dpi_inf_dt = 0d0

                    if (model_eqns == 1) then
                        drho_dt = dalpha_rho_dt(1)
                        dgamma_dt = dadv_dt(1)
                        dpi_inf_dt = dadv_dt(2)
                    else
    !$acc loop seq
                        do i = 1, num_fluids
                            drho_dt = drho_dt + dalpha_rho_dt(i)
                            dgamma_dt = dgamma_dt + dadv_dt(i)*gammas(i)
                            dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_infs(i)
                        end do
                    end if
                    ! ============================================================

                    ! flux_rs_vf and flux_src_rs_vf at j = -1/2 ==================
    !$acc loop seq
                    do i = 1, contxe
                        flux_rsz_vf(-1, k, r, i) = flux_rsz_vf(0, k, r, i) &
                                                     + ds(0)*dalpha_rho_dt(i)
                    end do

    !$acc loop seq
                    do i = momxb, momxe
                        flux_rsz_vf(-1, k, r, i) = flux_rsz_vf(0, k, r, i) &
                                                     + ds(0)*(vel(i - contxe)*drho_dt &
                                                              + rho*dvel_dt(i - contxe))
                    end do

                    flux_rsz_vf(-1, k, r, E_idx) = flux_rsz_vf(0, k, r, E_idx) &
                                                     + ds(0)*(pres*dgamma_dt &
                                                              + gamma*dpres_dt &
                                                              + dpi_inf_dt &
                                                              + rho*vel_dv_dt_sum &
                                                              + 5d-1*drho_dt*vel_K_sum)

                    if (riemann_solver == 1) then
    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsz_vf(-1, k, r, i) = 0d0
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsz_vf(-1, k, r, i) = &
                                1d0/max(abs(vel(dir_idx(1))), sgm_eps) &
                                *sign(1d0, vel(dir_idx(1))) &
                                *(flux_rsz_vf(0, k, r, i) &
                                  + vel(dir_idx(1)) &
                                  *flux_src_rsz_vf(0, k, r, i) &
                                  + ds(0)*dadv_dt(i - E_idx))
                        end do

                    else

    !$acc loop seq
                        do i = advxb, advxe
                            flux_rsz_vf(-1, k, r, i) = flux_rsz_vf(0, k, r, i) - &
                                                         adv(i - E_idx)*flux_src_rsz_vf(0, k, r, i) + &
                                                         ds(0)*dadv_dt(i - E_idx)
                        end do

    !$acc loop seq
                        do i = advxb, advxe
                            flux_src_rsz_vf(-1, k, r, i) = 0d0
                        end do

                    end if
                    ! END: flux_rs_vf and flux_src_rs_vf at j = -1/2 =============

                end do
            end do 
        end if           



        ! END: FD2 or FD4 of RHS at j = 0 ==================================

        ! The reshaping of outputted data and disssociation of the FD and PI
        ! coefficients, or CBC coefficients, respectively, based on selected
        ! CBC coordinate direction.

        
        call s_finalize_cbc(flux_vf, flux_src_vf, &
                            cbc_dir, cbc_loc, &
                            ix, iy, iz)


    end subroutine s_cbc ! -------------------------------------------------

    !>  The L variables for the slip wall CBC, see pg. 451 of
        !!      Thompson (1990). At the slip wall (frictionless wall),
        !!      the normal component of velocity is zero at all times,
        !!      while the transverse velocities may be nonzero.
        !!  @param dflt_int Default null integer
    subroutine s_compute_slip_wall_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, advxe
         L(i) = 0d0
        end do

        L(advxe) = L(1)

    end subroutine s_compute_slip_wall_L ! ---------------------------------

    !>  The L variables for the nonreflecting subsonic buffer CBC
        !!      see pg. 13 of Thompson (1987). The nonreflecting subsonic
        !!      buffer reduces the amplitude of any reflections caused by
        !!      outgoing waves.
        !!  @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_buffer_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i !< Generic loop iterator

        L(1) = (5d-1 - 5d-1*sign(1d0, lambda(1)))*lambda(1) &
               *(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = (5d-1 - 5d-1*sign(1d0, lambda(2)))*lambda(2) &
                   *(dadv_ds(i - momxe))
        end do

        L(advxe) = (5d-1 - 5d-1*sign(1d0, lambda(3)))*lambda(3) &
                         *(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_nonreflecting_subsonic_buffer_L ! -------------

    !>  The L variables for the nonreflecting subsonic inflow CBC
        !!      see pg. 455, Thompson (1990). This nonreflecting subsonic
        !!      CBC assumes an incoming flow and reduces the amplitude of
        !!      any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, advxe
            L(i) = 0d0
        end do

    end subroutine s_compute_nonreflecting_subsonic_inflow_L ! -------------

    !>  The L variables for the nonreflecting subsonic outflow
        !!      CBC see pg. 454 of Thompson (1990). This nonreflecting
        !!      subsonic CBC presumes an outgoing flow and reduces the
        !!      amplitude of any reflections caused by outgoing waves.
        !! @param dflt_int Default null integer
    subroutine s_compute_nonreflecting_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        ! bubble index
        L(advxe) = 0d0

    end subroutine s_compute_nonreflecting_subsonic_outflow_L ! ------------

    !>  The L variables for the force-free subsonic outflow CBC,
        !!      see pg. 454 of Thompson (1990). The force-free subsonic
        !!      outflow sets to zero the sum of all of the forces which
        !!      are acting on a fluid element for the normal coordinate
        !!      direction to the boundary. As a result, a fluid element
        !!      at the boundary is simply advected outward at the fluid
        !!      velocity.
        !! @param dflt_int Default null integer
    subroutine s_compute_force_free_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = L(1) + 2d0*rho*c*lambda(2)*dvel_ds(dir_idx(1))

    end subroutine s_compute_force_free_subsonic_outflow_L ! ---------------

    !>  L variables for the constant pressure subsonic outflow
        !!      CBC see pg. 455 Thompson (1990). The constant pressure
        !!      subsonic outflow maintains a fixed pressure at the CBC
        !!      boundary in absence of any transverse effects.
        !! @param dflt_int Default null integer
    subroutine s_compute_constant_pressure_subsonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i !> Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = -L(1)

    end subroutine s_compute_constant_pressure_subsonic_outflow_L ! --------

    !>  L variables for the supersonic inflow CBC, see pg. 453
        !!      Thompson (1990). The supersonic inflow CBC is a steady
        !!      state, or nearly a steady state, CBC in which only the
        !!      transverse terms may generate a time dependence at the
        !!      inflow boundary.
        !! @param dflt_int Default null integer
    subroutine s_compute_supersonic_inflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i

        do i = 1, advxe
            L(i) = 0d0
        end do

    end subroutine s_compute_supersonic_inflow_L ! -------------------------

    !>  L variables for the supersonic outflow CBC, see pg. 453
        !!      of Thompson (1990). For the supersonic outflow CBC, the
        !!      flow evolution at the boundary is determined completely
        !!      by the interior data.
        !! @param dflt_int Default null integer
    subroutine s_compute_supersonic_outflow_L(dflt_int, lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds) ! --------------
!$acc routine seq
        integer, intent(IN) :: dflt_int
        real(kind(0d0)), dimension(:), intent(IN) :: lambda, mf, dalpha_rho_ds, dvel_ds, dadv_ds
        real(kind(0d0)), intent(IN) :: rho, c, dpres_ds
        real(kind(0d0)), dimension(:), intent(INOUT) :: L

        integer :: i !< Generic loop iterator

        L(1) = lambda(1)*(dpres_ds - rho*c*dvel_ds(dir_idx(1)))

        do i = 2, momxb
            L(i) = lambda(2)*(c*c*dalpha_rho_ds(i - 1) - mf(i - 1)*dpres_ds)
        end do

        do i = momxb + 1, momxe
            L(i) = lambda(2)*(dvel_ds(dir_idx(i - contxe)))
        end do

        do i = E_idx, advxe - 1
            L(i) = lambda(2)*(dadv_ds(i - momxe))
        end do

        L(advxe) = lambda(3)*(dpres_ds + rho*c*dvel_ds(dir_idx(1)))

    end subroutine s_compute_supersonic_outflow_L ! ------------------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      selected CBC.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, & ! ------
                                cbc_dir, cbc_loc, &
                                ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc
        type(bounds_info), intent(IN) :: ix, iy, iz


        integer :: i, j, k, r !< Generic loop iterators

        ! Configuring the coordinate direction indexes and flags


        ! Determining the indicial shift based on CBC location
        


        ! END: Allocation/Association of Primitive and Flux Variables ======

        if (cbc_dir == 1) then
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        elseif (cbc_dir == 2) then
            is1%beg = 0; is1%end = buff_size; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0d0, 1d0, 0d0/)
        else
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0d0, 0d0, 1d0/)
        end if

        dj = max(0, cbc_loc)

        !$acc update device(is1, is2, is3, dir_idx, dir_flg, dj)

        ! Reshaping Inputted Data in x-direction ===========================
        if (cbc_dir == 1) then

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsx_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(dj*(m - 2*j) + j, k, r)
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsx_vf(j, k, r, momxb) = &
                            q_prim_vf(momxb)%sf(dj*(m - 2*j) + j, k, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsx_vf(j, k, r, i) = &
                                flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsx_vf(j, k, r, momxb) = &
                            flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_rsx_vf(j, k, r, advxb) = &
                            flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsx_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)                
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsx_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Inputted Data in x-direction ======================

            ! Reshaping Inputted Data in y-direction ===========================
        elseif (cbc_dir == 2) then

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsy_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(k, dj*(n - 2*j) + j, r)
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsy_vf(j, k, r, momxb + 1) = &
                            q_prim_vf(momxb + 1)%sf(k, dj*(n - 2*j) + j, r)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsy_vf(j, k, r, i) = &
                                flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsy_vf(j, k, r, momxb + 1) = &
                            flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_rsy_vf(j, k, r, advxb) = &
                            flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsy_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsy_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Inputted Data in y-direction ======================

            ! Reshaping Inputted Data in z-direction ===========================
        else

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, sys_size
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = 0, buff_size
                            q_prim_rsz_vf(j, k, r, i) = &
                                q_prim_vf(i)%sf(r, k, dj*(p - 2*j) + j)
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsz_vf(j, k, r, momxe) = &
                            q_prim_vf(momxe)%sf(r, k, dj*(p - 2*j) + j)* &
                            sign(1d0, -real(cbc_loc, kind(0d0)))
                    end do
                end do
            end do

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsz_vf(j, k, r, i) = &
                                flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsz_vf(j, k, r, momxe) = &
                            flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_rsz_vf(j, k, r, advxb) = &
                            flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsz_vf(j, k, r, i) = &
                                    flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsz_vf(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Inputted Data in z-direction ======================

        ! Association of the procedural pointer to the appropriate procedure
        ! that will be utilized in the evaluation of L variables for the CBC

        ! ==================================================================

        ! ==================================================================

    end subroutine s_initialize_cbc ! --------------------------------------

    !>  Deallocation and/or the disassociation procedures that
        !!      are necessary in order to finalize the CBC application
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir CBC coordinate direction
        !!  @param cbc_loc CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_finalize_cbc(flux_vf, flux_src_vf, & ! -------------------
                              cbc_dir, cbc_loc, &
                              ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf

        integer, intent(IN) :: cbc_dir, cbc_loc
        type(bounds_info), intent(IN) :: ix, iy, iz


        integer :: i, j, k, r !< Generic loop iterators

        ! Determining the indicial shift based on CBC location
        dj = max(0, cbc_loc)
        !$acc update device(dj)

        ! Reshaping Outputted Data in x-direction ==========================
        if (cbc_dir == 1) then

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_rsx_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_rsx_vf(j, k, r, momxb)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_src_rsx_vf(j, k, r, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                    flux_src_rsx_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_src_rsx_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in x-direction =====================

            ! Reshaping Outputted Data in y-direction ==========================
        elseif (cbc_dir == 2) then

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_rsy_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_rsy_vf(j, k, r, momxb + 1)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_src_rsy_vf(j, k, r, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                    flux_src_rsy_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_src_rsy_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in y-direction =====================

            ! Reshaping Outputted Data in z-direction ==========================
        else

!$acc parallel loop collapse(4) gang vector default(present)
            do i = 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_rsz_vf(j, k, r, i)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_rsz_vf(j, k, r, momxe)
                    end do
                end do
            end do

!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_src_rsz_vf(j, k, r, advxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
!$acc parallel loop collapse(4) gang vector default(present)
                do i = advxb + 1, advxe
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                    flux_src_rsz_vf(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
!$acc parallel loop collapse(3) gang vector default(present)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_src_rsz_vf(j, k, r, advxb)* &
                                sign(1d0, -real(cbc_loc, kind(0d0)))
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Outputted Data in z-direction =====================

        ! Deallocation/Disassociation of Primitive and Flux Variables ======

        ! ==================================================================

        ! Nullifying procedural pointer used in evaluation of L for the CBC


    end subroutine s_finalize_cbc ! ----------------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_cbc_module() ! -----------------------------------

        if (all((/bc_x%beg, bc_x%end/) > -5) &
            .and. &
            (n > 0 .and. all((/bc_y%beg, bc_y%end/) > -5)) &
            .and. &
            (p > 0 .and. all((/bc_z%beg, bc_z%end/) > -5))) return

        ! Deallocating the cell-average primitive variables
        deallocate (q_prim_rsx_vf)
        if(weno_order > 1) then
            deallocate (F_rsx_vf, F_src_rsx_vf)
        end if
        deallocate (flux_rsx_vf, flux_src_rsx_vf)


        if(n > 0) then
            deallocate (q_prim_rsy_vf)
            if(weno_order > 1) then
                deallocate (F_rsy_vf, F_src_rsy_vf)
            end if
            deallocate (flux_rsy_vf, flux_src_rsy_vf)
        end if
        if(p > 0) then
            deallocate (q_prim_rsz_vf)
            if(weno_order > 1) then
                deallocate (F_rsz_vf, F_src_rsz_vf)
            end if
            deallocate (flux_rsz_vf, flux_src_rsz_vf)
        end if

        ! Deallocating the cell-average partial densities, the velocity, the
        ! advection variables, the mass fractions and also the Weber numbers
        deallocate (alpha_rho, vel, adv, mf)

        ! Deallocating the first-order spatial derivatives, in s-direction,
        ! of the partial densities, the velocity and the advected variables
        deallocate (dalpha_rho_ds, dvel_ds, dadv_ds)

        ! Deallocating L, see Thompson (1987, 1990)
        deallocate (L)

        ! Deallocating the cell-width distribution in the s-direction
        deallocate (ds)

        ! Deallocating CBC Coefficients in x-direction =====================
        if (any((/bc_x%beg, bc_x%end/) <= -5)) then
            deallocate (fd_coef_x); if (weno_order > 1) deallocate (pi_coef_x)
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in y-direction =====================
        if (n > 0 .and. any((/bc_y%beg, bc_y%end/) <= -5)) then
            deallocate (fd_coef_y); if (weno_order > 1) deallocate (pi_coef_y)
        end if
        ! ==================================================================

        ! Deallocating CBC Coefficients in z-direction =====================
        if (p > 0 .and. any((/bc_z%beg, bc_z%end/) <= -5)) then
            deallocate (fd_coef_z); if (weno_order > 1) deallocate (pi_coef_z)
        end if
        ! ==================================================================

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables

    end subroutine s_finalize_cbc_module ! ---------------------------------

end module m_cbc
