!>
!! @file m_cbc.f90
!! @brief Contains module m_cbc

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

#:include 'macros.fpp'

module m_cbc

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use m_compute_cbc

    use m_thermochem, only: &
        get_mixture_energy_mass, get_mixture_specific_heat_cv_mass, &
        get_mixture_specific_heat_cp_mass, gas_constant, &
        get_mixture_molecular_weight, get_species_enthalpies_rt, &
        molecular_weights, get_species_specific_heats_r, &
        get_mole_fractions, get_species_specific_heats_r

    implicit none

    private; public :: s_initialize_cbc_module, s_cbc, s_finalize_cbc_module

    !! The cell-average primitive variables. They are obtained by reshaping (RS)
    !! q_prim_vf in the coordinate direction normal to the domain boundary along
    !! which the CBC is applied.

    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_rsx_vf
    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_rsy_vf
    real(wp), allocatable, dimension(:, :, :, :) :: q_prim_rsz_vf
    $:GPU_DECLARE(create='[q_prim_rsx_vf,q_prim_rsy_vf,q_prim_rsz_vf]')

    !! Cell-average fluxes (src - source). These are directly determined from the
    !! cell-average primitive variables, q_prims_rs_vf, and not a Riemann solver.

    real(wp), allocatable, dimension(:, :, :, :) :: F_rsx_vf, F_src_rsx_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: F_rsy_vf, F_src_rsy_vf !<
    real(wp), allocatable, dimension(:, :, :, :) :: F_rsz_vf, F_src_rsz_vf !<
    $:GPU_DECLARE(create='[F_rsx_vf,F_src_rsx_vf,F_rsy_vf,F_src_rsy_vf,F_rsz_vf,F_src_rsz_vf]')

    !! There is a CCE bug that is causing some subset of these variables to interfere
    !! with variables of the same name in m_riemann_solvers.fpp, and giving this versions
    !! unique "_l" names works around the bug. Other private module allocatable arrays
    !! in `acc declare create` clauses don't have this problem, so we still need to
    !! isolate this bug.

    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsx_vf_l, flux_src_rsx_vf_l !<
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsy_vf_l, flux_src_rsy_vf_l
    real(wp), allocatable, dimension(:, :, :, :) :: flux_rsz_vf_l, flux_src_rsz_vf_l
    $:GPU_DECLARE(create='[flux_rsx_vf_l,flux_src_rsx_vf_l,flux_rsy_vf_l,flux_src_rsy_vf_l,flux_rsz_vf_l,flux_src_rsz_vf_l]')

    real(wp) :: dpres_ds !< Spatial derivatives in s-dir of pressure
    $:GPU_DECLARE(create='[dpres_ds]')

    real(wp), allocatable, dimension(:) :: ds !< Cell-width distribution in the s-direction

    ! CBC Coefficients

    real(wp), allocatable, dimension(:, :) :: fd_coef_x !< Finite diff. coefficients x-dir
    real(wp), allocatable, dimension(:, :) :: fd_coef_y !< Finite diff. coefficients y-dir
    real(wp), allocatable, dimension(:, :) :: fd_coef_z !< Finite diff. coefficients z-dir

    !! The first dimension identifies the location of a coefficient in the FD
    !! formula, while the last dimension denotes the location of the CBC.

    ! Bug with NVHPC when using nullified pointers in a declare create
    !    real(wp), pointer, dimension(:, :) :: fd_coef => null()

    real(wp), allocatable, dimension(:, :, :) :: pi_coef_x !< Polynomial interpolant coefficients in x-dir
    real(wp), allocatable, dimension(:, :, :) :: pi_coef_y !< Polynomial interpolant coefficients in y-dir
    real(wp), allocatable, dimension(:, :, :) :: pi_coef_z !< Polynomial interpolant coefficients in z-dir

    $:GPU_DECLARE(create='[ds,fd_coef_x,fd_coef_y,fd_coef_z,pi_coef_x,pi_coef_y,pi_coef_z]')

    !! The first dimension of the array identifies the polynomial, the
    !! second dimension identifies the position of its coefficients and the last
    !! dimension denotes the location of the CBC.

    type(int_bounds_info) :: is1, is2, is3 !< Indical bounds in the s1-, s2- and s3-directions
    $:GPU_DECLARE(create='[is1,is2,is3]')

    integer :: dj
    integer :: bcxb, bcxe, bcyb, bcye, bczb, bcze
    integer :: cbc_dir, cbc_loc
    integer :: flux_cbc_index
    $:GPU_DECLARE(create='[dj,bcxb,bcxe,bcyb,bcye,bczb,bcze]')
    $:GPU_DECLARE(create='[cbc_dir, cbc_loc,flux_cbc_index]')

    !! GRCBC inputs for subsonic inflow and outflow conditions consisting of
    !! inflow velocities, pressure, density and void fraction as well as
    !! outflow velocities and pressure

    real(wp), allocatable, dimension(:) :: pres_in, pres_out, Del_in, Del_out
    real(wp), allocatable, dimension(:, :) :: vel_in, vel_out
    real(wp), allocatable, dimension(:, :) :: alpha_rho_in, alpha_in
    $:GPU_DECLARE(create='[pres_in,pres_out,Del_in,Del_out]')
    $:GPU_DECLARE(create='[vel_in,vel_out]')
    $:GPU_DECLARE(create='[alpha_rho_in,alpha_in]')

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_cbc_module

        integer :: i
        logical :: is_cbc
        type(int_bounds_info) :: idx1, idx2

        if (chemistry) then
            flux_cbc_index = sys_size
        else
            flux_cbc_index = adv_idx%end
        end if
        $:GPU_UPDATE(device='[flux_cbc_index]')

        call s_any_cbc_boundaries(is_cbc)

        if (is_cbc .eqv. .false.) return

        if (n == 0) then
            is2%beg = 0

        else
            is2%beg = -buff_size
        end if

        is2%end = n - is2%beg

        if (p == 0) then
            is3%beg = 0

        else
            is3%beg = -buff_size
        end if
        is3%end = p - is3%beg

        @:ALLOCATE(q_prim_rsx_vf(0:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:sys_size))

        if (weno_order > 1 .or. muscl_order > 1) then

            @:ALLOCATE(F_rsx_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:flux_cbc_index))

            @:ALLOCATE(F_src_rsx_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        @:ALLOCATE(flux_rsx_vf_l(-1:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, 1:flux_cbc_index))

        @:ALLOCATE(flux_src_rsx_vf_l(-1:buff_size, &
            is2%beg:is2%end, &
            is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        if (n > 0) then

            if (m == 0) then
                is2%beg = 0

            else
                is2%beg = -buff_size
            end if

            is2%end = m - is2%beg

            if (p == 0) then
                is3%beg = 0

            else
                is3%beg = -buff_size
            end if
            is3%end = p - is3%beg

            @:ALLOCATE(q_prim_rsy_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1 .or. muscl_order > 1) then

                @:ALLOCATE(F_rsy_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, 1:flux_cbc_index))

                @:ALLOCATE(F_src_rsy_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, adv_idx%beg:adv_idx%end))

            end if

            @:ALLOCATE(flux_rsy_vf_l(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:flux_cbc_index))

            @:ALLOCATE(flux_src_rsy_vf_l(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        if (p > 0) then

            if (n == 0) then
                is2%beg = 0

            else
                is2%beg = -buff_size
            end if

            is2%end = n - is2%beg

            if (m == 0) then
                is3%beg = 0

            else
                is3%beg = -buff_size
            end if
            is3%end = m - is3%beg

            @:ALLOCATE(q_prim_rsz_vf(0:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:sys_size))

            if (weno_order > 1 .or. muscl_order > 1) then

                @:ALLOCATE(F_rsz_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, 1:flux_cbc_index))

                @:ALLOCATE(F_src_rsz_vf(0:buff_size, &
                    is2%beg:is2%end, &
                    is3%beg:is3%end, adv_idx%beg:adv_idx%end))

            end if

            @:ALLOCATE(flux_rsz_vf_l(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, 1:flux_cbc_index))

            @:ALLOCATE(flux_src_rsz_vf_l(-1:buff_size, &
                is2%beg:is2%end, &
                is3%beg:is3%end, adv_idx%beg:adv_idx%end))

        end if

        ! Allocating the cell-width distribution in the s-direction
        @:ALLOCATE(ds(0:buff_size))

        if (recon_type == WENO_TYPE) then
            idx1%beg = 0
            idx1%end = weno_polyn - 1
            idx2%beg = 0
            idx2%end = weno_order - 3
        else if (recon_type == MUSCL_TYPE) then
            idx1%beg = 0
            idx1%end = muscl_polyn
            idx2%beg = 0
            idx2%end = muscl_order - 1
        end if
        ! Allocating/Computing CBC Coefficients in x-direction
        if (all((/bc_x%beg, bc_x%end/) <= -5) .and. all((/bc_x%beg, bc_x%end/) >= -13)) then

            @:ALLOCATE(fd_coef_x(0:buff_size, -1:1))

            if (weno_order > 1 .or. muscl_order > 1) then
                @:ALLOCATE(pi_coef_x(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:1))
            end if

            call s_compute_cbc_coefficients(1, -1)
            call s_compute_cbc_coefficients(1, 1)

        elseif (bc_x%beg <= -5 .and. bc_x%beg >= -13) then

            @:ALLOCATE(fd_coef_x(0:buff_size, -1:-1))

            if (weno_order > 1 .or. muscl_order > 1) then
                @:ALLOCATE(pi_coef_x(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:-1))
            end if

            call s_compute_cbc_coefficients(1, -1)

        elseif (bc_x%end <= -5 .and. bc_x%end >= -13) then

            @:ALLOCATE(fd_coef_x(0:buff_size, 1:1))

            if (weno_order > 1 .or. muscl_order > 1) then
                @:ALLOCATE(pi_coef_x(idx1%beg:idx1%end, idx2%beg:idx2%end, 1:1))
            end if

            call s_compute_cbc_coefficients(1, 1)

        end if

        ! Allocating/Computing CBC Coefficients in y-direction
        if (n > 0) then

            if (all((/bc_y%beg, bc_y%end/) <= -5) .and. all((/bc_y%beg, bc_y%end/) >= -13)) then

                @:ALLOCATE(fd_coef_y(0:buff_size, -1:1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_y(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:1))
                end if

                call s_compute_cbc_coefficients(2, -1)
                call s_compute_cbc_coefficients(2, 1)

            elseif (bc_y%beg <= -5 .and. bc_y%beg >= -13) then

                @:ALLOCATE(fd_coef_y(0:buff_size, -1:-1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_y(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:-1))
                end if

                call s_compute_cbc_coefficients(2, -1)

            elseif (bc_y%end <= -5 .and. bc_y%end >= -13) then

                @:ALLOCATE(fd_coef_y(0:buff_size, 1:1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_y(idx1%beg:idx1%end, idx2%beg:idx2%end, 1:1))
                end if

                call s_compute_cbc_coefficients(2, 1)

            end if

        end if

        ! Allocating/Computing CBC Coefficients in z-direction
        if (p > 0) then

            if (all((/bc_z%beg, bc_z%end/) <= -5) .and. all((/bc_z%beg, bc_z%end/) >= -13)) then

                @:ALLOCATE(fd_coef_z(0:buff_size, -1:1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_z(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:1))
                end if

                call s_compute_cbc_coefficients(3, -1)
                call s_compute_cbc_coefficients(3, 1)

            elseif (bc_z%beg <= -5 .and. bc_z%beg >= -13) then

                @:ALLOCATE(fd_coef_z(0:buff_size, -1:-1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_z(idx1%beg:idx1%end, idx2%beg:idx2%end, -1:-1))
                end if

                call s_compute_cbc_coefficients(3, -1)

            elseif (bc_z%end <= -5 .and. bc_z%end >= -13) then

                @:ALLOCATE(fd_coef_z(0:buff_size, 1:1))

                if (weno_order > 1 .or. muscl_order > 1) then
                    @:ALLOCATE(pi_coef_z(idx1%beg:idx1%end, idx2%beg:idx2%end, 1:1))
                end if

                call s_compute_cbc_coefficients(3, 1)

            end if

        end if

        $:GPU_UPDATE(device='[fd_coef_x,fd_coef_y,fd_coef_z, &
            & pi_coef_x,pi_coef_y,pi_coef_z]')

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        bcxb = bc_x%beg
        bcxe = bc_x%end

        $:GPU_UPDATE(device='[bcxb, bcxe]')

        if (n > 0) then
            bcyb = bc_y%beg
            bcye = bc_y%end

            $:GPU_UPDATE(device='[bcyb, bcye]')
        end if

        if (p > 0) then
            bczb = bc_z%beg
            bcze = bc_z%end

            $:GPU_UPDATE(device='[bczb, bcze]')
        end if

        ! Allocate GRCBC inputs
        @:ALLOCATE(pres_in(1:num_dims), pres_out(1:num_dims))
        @:ALLOCATE(Del_in(1:num_dims), Del_out(1:num_dims))
        @:ALLOCATE(vel_in(1:num_dims, 1:num_dims), vel_out(1:num_dims, 1:num_dims))
        @:ALLOCATE(alpha_rho_in(1:num_fluids, 1:num_dims), alpha_in(1:num_fluids, 1:num_dims))

        ! Assign and update GRCBC inputs
        #:for CBC_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (${CBC_DIR}$ <= num_dims) then
                vel_in(${CBC_DIR}$, 1) = bc_${XYZ}$%vel_in(1)
                vel_out(${CBC_DIR}$, 1) = bc_${XYZ}$%vel_out(1)
                if (n > 0) then
                    vel_in(${CBC_DIR}$, 2) = bc_${XYZ}$%vel_in(2)
                    vel_out(${CBC_DIR}$, 2) = bc_${XYZ}$%vel_out(2)
                    if (p > 0) then
                        vel_in(${CBC_DIR}$, 3) = bc_${XYZ}$%vel_in(3)
                        vel_out(${CBC_DIR}$, 3) = bc_${XYZ}$%vel_out(3)
                    end if
                end if
                Del_in(${CBC_DIR}$) = maxval(d${XYZ}$)
                Del_out(${CBC_DIR}$) = maxval(d${XYZ}$)
                pres_in(${CBC_DIR}$) = bc_${XYZ}$%pres_in
                pres_out(${CBC_DIR}$) = bc_${XYZ}$%pres_out
                do i = 1, num_fluids
                    alpha_rho_in(i, ${CBC_DIR}$) = bc_${XYZ}$%alpha_rho_in(i)
                    alpha_in(i, ${CBC_DIR}$) = bc_${XYZ}$%alpha_in(i)
                end do
            end if
        #:endfor
        $:GPU_UPDATE(device='[vel_in,vel_out,pres_in,pres_out, &
            & Del_in,Del_out,alpha_rho_in,alpha_in]')

    end subroutine s_initialize_cbc_module

    !>  Compute CBC coefficients
        !!  @param cbc_dir_in CBC coordinate direction
        !!  @param cbc_loc_in CBC coordinate location
    subroutine s_compute_cbc_coefficients(cbc_dir_in, cbc_loc_in)
        ! Description: The purpose of this subroutine is to compute the grid
        !              dependent FD and PI coefficients, or CBC coefficients,
        !              provided the CBC coordinate direction and location.

        ! CBC coordinate direction and location
        integer, intent(in) :: cbc_dir_in, cbc_loc_in

        ! Cell-boundary locations in the s-direction
        real(wp), dimension(0:buff_size + 1) :: s_cb

        ! Generic loop iterator
        integer :: i

        ! Associating CBC coefficients pointers
        call s_associate_cbc_coefficients_pointers(cbc_dir_in, cbc_loc_in)

        ! Determining the cell-boundary locations in the s-direction
        s_cb(0) = 0._wp

        do i = 0, buff_size
            s_cb(i + 1) = s_cb(i) + ds(i)
        end do

        ! Computing CBC1 Coefficients
        #:for CBC_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (cbc_dir_in == ${CBC_DIR}$ .and. recon_type == WENO_TYPE) then
                if (weno_order == 1) then

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0._wp
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -2._wp/(ds(0) + ds(1))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -fd_coef_${XYZ}$ (0, cbc_loc_in)

                    ! Computing CBC2 Coefficients
                elseif (weno_order == 3) then

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0._wp
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -6._wp/(3._wp*ds(0) + 2._wp*ds(1) - ds(2))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -4._wp*fd_coef_${XYZ}$ (0, cbc_loc_in)/3._wp
                    fd_coef_${XYZ}$ (2, cbc_loc_in) = fd_coef_${XYZ}$ (0, cbc_loc_in)/3._wp

                    pi_coef_${XYZ}$ (0, 0, cbc_loc_in) = (s_cb(0) - s_cb(1))/(s_cb(0) - s_cb(2))

                    ! Computing CBC4 Coefficients
                else

                    fd_coef_${XYZ}$ (:, cbc_loc_in) = 0._wp
                    fd_coef_${XYZ}$ (0, cbc_loc_in) = -50._wp/(25._wp*ds(0) + 2._wp*ds(1) &
                                                               - 1.e1_wp*ds(2) + 1.e1_wp*ds(3) &
                                                               - 3._wp*ds(4))
                    fd_coef_${XYZ}$ (1, cbc_loc_in) = -48._wp*fd_coef_${XYZ}$ (0, cbc_loc_in)/25._wp
                    fd_coef_${XYZ}$ (2, cbc_loc_in) = 36._wp*fd_coef_${XYZ}$ (0, cbc_loc_in)/25._wp
                    fd_coef_${XYZ}$ (3, cbc_loc_in) = -16._wp*fd_coef_${XYZ}$ (0, cbc_loc_in)/25._wp
                    fd_coef_${XYZ}$ (4, cbc_loc_in) = 3._wp*fd_coef_${XYZ}$ (0, cbc_loc_in)/25._wp

                    pi_coef_${XYZ}$ (0, 0, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(1))*(s_cb(1) - s_cb(2))* &
                         (s_cb(1) - s_cb(3)))/((s_cb(1) - s_cb(4))* &
                                               (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(2)))
                    pi_coef_${XYZ}$ (0, 1, cbc_loc_in) = &
                        ((s_cb(1) - s_cb(0))*(s_cb(1) - s_cb(2))* &
                         ((s_cb(1) - s_cb(3))*(s_cb(1) - s_cb(3)) - &
                          (s_cb(0) - s_cb(4))*((s_cb(3) - s_cb(1)) + &
                                               (s_cb(4) - s_cb(1)))))/ &
                        ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                         (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                    pi_coef_${XYZ}$ (0, 2, cbc_loc_in) = &
                        (s_cb(1) - s_cb(0))*((s_cb(1) - s_cb(2))* &
                                             (s_cb(1) - s_cb(3)) + ((s_cb(0) - s_cb(2)) + &
                                                                    (s_cb(1) - s_cb(3)))*(s_cb(0) - s_cb(4)))/ &
                        ((s_cb(2) - s_cb(0))*(s_cb(0) - s_cb(3))* &
                         (s_cb(0) - s_cb(4)))
                    pi_coef_${XYZ}$ (1, 0, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(2))*(s_cb(2) - s_cb(1))* &
                         (s_cb(2) - s_cb(3)))/((s_cb(2) - s_cb(4))* &
                                               (s_cb(4) - s_cb(0))*(s_cb(4) - s_cb(1)))
                    pi_coef_${XYZ}$ (1, 1, cbc_loc_in) = &
                        ((s_cb(0) - s_cb(2))*(s_cb(1) - s_cb(2))* &
                         ((s_cb(1) - s_cb(3))*(s_cb(2) - s_cb(3)) + &
                          (s_cb(0) - s_cb(4))*((s_cb(1) - s_cb(3)) + &
                                               (s_cb(2) - s_cb(4)))))/ &
                        ((s_cb(0) - s_cb(3))*(s_cb(1) - s_cb(3))* &
                         (s_cb(0) - s_cb(4))*(s_cb(1) - s_cb(4)))
                    pi_coef_${XYZ}$ (1, 2, cbc_loc_in) = &
                        ((s_cb(1) - s_cb(2))*(s_cb(2) - s_cb(3))* &
                         (s_cb(2) - s_cb(4)))/((s_cb(0) - s_cb(2))* &
                                               (s_cb(0) - s_cb(3))*(s_cb(0) - s_cb(4)))

                end if
            end if
        #:endfor

        ! END: Computing CBC4 Coefficients

        ! Nullifying CBC coefficients

    end subroutine s_compute_cbc_coefficients

    !!  The goal of the procedure is to associate the FD and PI
    !!      coefficients, or CBC coefficients, with the appropriate
    !!      targets, based on the coordinate direction and location
    !!      of the CBC.
    !!  @param cbc_dir_in CBC coordinate direction
    !!  @param cbc_loc_in CBC coordinate location
    subroutine s_associate_cbc_coefficients_pointers(cbc_dir_in, cbc_loc_in)

        integer, intent(in) :: cbc_dir_in, cbc_loc_in

        integer :: i !< Generic loop iterator

        ! Associating CBC Coefficients in x-direction
        if (cbc_dir_in == 1) then

            !fd_coef => fd_coef_x; if (weno_order > 1) pi_coef => pi_coef_x

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dx(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dx(m - i)
                end do
            end if

            ! Associating CBC Coefficients in y-direction
        elseif (cbc_dir_in == 2) then

            !fd_coef => fd_coef_y; if (weno_order > 1) pi_coef => pi_coef_y

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dy(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dy(n - i)
                end do
            end if

            ! Associating CBC Coefficients in z-direction
        else

            !fd_coef => fd_coef_z; if (weno_order > 1) pi_coef => pi_coef_z

            if (cbc_loc_in == -1) then
                do i = 0, buff_size
                    ds(i) = dz(i)
                end do
            else
                do i = 0, buff_size
                    ds(i) = dz(p - i)
                end do
            end if

        end if

        $:GPU_UPDATE(device='[ds]')

    end subroutine s_associate_cbc_coefficients_pointers

    !>  The following is the implementation of the CBC based on
        !!      the work of Thompson (1987, 1990) on hyperbolic systems.
        !!      The CBC is indirectly applied in the computation of the
        !!      right-hand-side (RHS) near the relevant domain boundary
        !!      through the modification of the fluxes.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param cbc_dir_norm CBC coordinate direction
        !!  @param cbc_loc_norm CBC coordinate location
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                     cbc_dir_norm, cbc_loc_norm, &
                     ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf

        integer, intent(in) :: cbc_dir_norm, cbc_loc_norm

        type(int_bounds_info), intent(in) :: ix, iy, iz

        ! First-order time derivatives of the partial densities, density,
        ! velocity, pressure, advection variables, and the specific heat
        ! ratio and liquid stiffness functions

        real(wp), dimension(num_fluids) :: dalpha_rho_dt
        real(wp) :: drho_dt
        real(wp), dimension(num_dims) :: dvel_dt
        real(wp) :: dpres_dt
        real(wp), dimension(num_fluids) :: dadv_dt
        real(wp) :: dgamma_dt
        real(wp) :: dpi_inf_dt
        real(wp) :: dqv_dt
        real(wp), dimension(contxe) :: alpha_rho, dalpha_rho_ds, mf
        real(wp), dimension(2) :: Re_cbc
        real(wp), dimension(num_vels) :: vel, dvel_ds
        real(wp), dimension(num_fluids) :: adv_local, dadv_ds
        real(wp), dimension(sys_size) :: L
        real(wp), dimension(3) :: lambda

        real(wp) :: rho         !< Cell averaged density
        real(wp) :: pres        !< Cell averaged pressure
        real(wp) :: E           !< Cell averaged energy
        real(wp) :: H           !< Cell averaged enthalpy
        real(wp) :: gamma       !< Cell averaged specific heat ratio
        real(wp) :: pi_inf      !< Cell averaged liquid stiffness
        real(wp) :: qv          !< Cell averaged fluid reference energy
        real(wp) :: c
        real(wp) :: Ma
        real(wp) :: T, sum_Enthalpies
        real(wp) :: Cv, Cp, e_mix, Mw, R_gas
        real(wp), dimension(num_species) :: Ys, h_k, dYs_dt, dYs_ds, Xs, Gamma_i, Cp_i

        real(wp) :: vel_K_sum, vel_dv_dt_sum

        integer :: i, j, k, r !< Generic loop iterators

        ! Reshaping of inputted data and association of the FD and PI
        ! coefficients, or CBC coefficients, respectively, hinging on
        ! selected CBC coordinate direction

        cbc_dir = cbc_dir_norm
        cbc_loc = cbc_loc_norm

        $:GPU_UPDATE(device='[cbc_dir, cbc_loc]')

        call s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                              ix, iy, iz)

        call s_associate_cbc_coefficients_pointers(cbc_dir, cbc_loc)

        #:for CBC_DIR, XYZ in [(1, 'x'), (2, 'y'), (3, 'z')]
            if (cbc_dir == ${CBC_DIR}$ .and. recon_type == WENO_TYPE) then

                ! PI2 of flux_rs_vf and flux_src_rs_vf at j = 1/2
                if (weno_order == 3) then

                    call s_convert_primitive_to_flux_variables(q_prim_rs${XYZ}$_vf, &
                                                               F_rs${XYZ}$_vf, &
                                                               F_src_rs${XYZ}$_vf, &
                                                               is1, is2, is3, idwbuff(2)%beg, idwbuff(3)%beg)

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = 1, flux_cbc_index
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                flux_rs${XYZ}$_vf_l(0, k, r, i) = F_rs${XYZ}$_vf(0, k, r, i) &
                                                                  + pi_coef_${XYZ}$ (0, 0, cbc_loc)* &
                                                                  (F_rs${XYZ}$_vf(1, k, r, i) - &
                                                                   F_rs${XYZ}$_vf(0, k, r, i))
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=3)
                    do i = advxb, advxe
                        do r = is3%beg, is3%end
                            do k = is2%beg, is2%end
                                flux_src_rs${XYZ}$_vf_l(0, k, r, i) = F_src_rs${XYZ}$_vf(0, k, r, i) + &
                                                                      (F_src_rs${XYZ}$_vf(1, k, r, i) - &
                                                                       F_src_rs${XYZ}$_vf(0, k, r, i)) &
                                                                      *pi_coef_${XYZ}$ (0, 0, cbc_loc)
                            end do
                        end do
                    end do

                    ! PI4 of flux_rs_vf and flux_src_rs_vf at j = 1/2, 3/2
                else
                    call s_convert_primitive_to_flux_variables(q_prim_rs${XYZ}$_vf, &
                                                               F_rs${XYZ}$_vf, &
                                                               F_src_rs${XYZ}$_vf, &
                                                               is1, is2, is3, idwbuff(2)%beg, idwbuff(3)%beg)

                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = 1, flux_cbc_index
                        do j = 0, 1
                            do r = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    flux_rs${XYZ}$_vf_l(j, k, r, i) = F_rs${XYZ}$_vf(j, k, r, i) &
                                                                      + pi_coef_${XYZ}$ (j, 0, cbc_loc)* &
                                                                      (F_rs${XYZ}$_vf(3, k, r, i) - &
                                                                       F_rs${XYZ}$_vf(2, k, r, i)) &
                                                                      + pi_coef_${XYZ}$ (j, 1, cbc_loc)* &
                                                                      (F_rs${XYZ}$_vf(2, k, r, i) - &
                                                                       F_rs${XYZ}$_vf(1, k, r, i)) &
                                                                      + pi_coef_${XYZ}$ (j, 2, cbc_loc)* &
                                                                      (F_rs${XYZ}$_vf(1, k, r, i) - &
                                                                       F_rs${XYZ}$_vf(0, k, r, i))
                                end do
                            end do
                        end do
                    end do

                    $:GPU_PARALLEL_LOOP(collapse=4)
                    do i = advxb, advxe
                        do j = 0, 1
                            do r = is3%beg, is3%end
                                do k = is2%beg, is2%end
                                    flux_src_rs${XYZ}$_vf_l(j, k, r, i) = F_src_rs${XYZ}$_vf(j, k, r, i) + &
                                                                          (F_src_rs${XYZ}$_vf(3, k, r, i) - &
                                                                           F_src_rs${XYZ}$_vf(2, k, r, i)) &
                                                                          *pi_coef_${XYZ}$ (j, 0, cbc_loc) + &
                                                                          (F_src_rs${XYZ}$_vf(2, k, r, i) - &
                                                                           F_src_rs${XYZ}$_vf(1, k, r, i)) &
                                                                          *pi_coef_${XYZ}$ (j, 1, cbc_loc) + &
                                                                          (F_src_rs${XYZ}$_vf(1, k, r, i) - &
                                                                           F_src_rs${XYZ}$_vf(0, k, r, i)) &
                                                                          *pi_coef_${XYZ}$ (j, 2, cbc_loc)
                                end do
                            end do
                        end do
                    end do

                end if

                ! FD2 or FD4 of RHS at j = 0
                $:GPU_PARALLEL_LOOP(collapse=2, private='[alpha_rho, vel, adv_local, &
                    & mf, dvel_ds, dadv_ds, Re_cbc, dalpha_rho_ds,dvel_dt, &
                    & dadv_dt, dalpha_rho_dt, L, lambda, Ys, dYs_dt, &
                    & dYs_ds, h_k, Cp_i, Gamma_i, Xs]')
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end

                        ! Transferring the Primitive Variables
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            alpha_rho(i) = q_prim_rs${XYZ}$_vf(0, k, r, i)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            vel(i) = q_prim_rs${XYZ}$_vf(0, k, r, contxe + i)
                        end do

                        vel_K_sum = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            vel_K_sum = vel_K_sum + vel(i)**2._wp
                        end do

                        pres = q_prim_rs${XYZ}$_vf(0, k, r, E_idx)

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, advxe - E_idx
                            adv_local(i) = q_prim_rs${XYZ}$_vf(0, k, r, E_idx + i)
                        end do

                        if (bubbles_euler) then
                            call s_convert_species_to_mixture_variables_bubbles_acc(rho, gamma, pi_inf, qv, adv_local, alpha_rho, Re_cbc)
                        else
                            call s_convert_species_to_mixture_variables_acc(rho, gamma, pi_inf, qv, adv_local, alpha_rho, Re_cbc)
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            mf(i) = alpha_rho(i)/rho
                        end do

                        if (chemistry) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = chemxb, chemxe
                                Ys(i - chemxb + 1) = q_prim_rs${XYZ}$_vf(0, k, r, i)
                            end do

                            call get_mixture_molecular_weight(Ys, Mw)
                            R_gas = gas_constant/Mw
                            T = pres/rho/R_gas
                            call get_mixture_specific_heat_cp_mass(T, Ys, Cp)
                            call get_mixture_energy_mass(T, Ys, e_mix)
                            E = rho*e_mix + 5.e-1_wp*rho*vel_K_sum
                            if (chem_params%gamma_method == 1) then
                                !> gamma_method = 1: Ref. Section 2.3.1 Formulation of doi:10.7907/ZKW8-ES97.
                                call get_mole_fractions(Mw, Ys, Xs)
                                call get_species_specific_heats_r(T, Cp_i)
                                Gamma_i = Cp_i/(Cp_i - 1.0_wp)
                                gamma = sum(Xs(:)/(Gamma_i(:) - 1.0_wp))
                            else if (chem_params%gamma_method == 2) then
                                !> gamma_method = 2: c_p / c_v where c_p, c_v are specific heats.
                                call get_mixture_specific_heat_cv_mass(T, Ys, Cv)
                                gamma = 1.0_wp/(Cp/Cv - 1.0_wp)
                            end if
                        else
                            E = gamma*pres + pi_inf + 5.e-1_wp*rho*vel_K_sum
                        end if

                        H = (E + pres)/rho

                        ! Compute mixture sound speed
                        call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv_local, vel_K_sum, 0._wp, c)

                        ! First-Order Spatial Derivatives of Primitive Variables

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            dalpha_rho_ds(i) = 0._wp
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            dvel_ds(i) = 0._wp
                        end do

                        dpres_ds = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, advxe - E_idx
                            dadv_ds(i) = 0._wp
                        end do

                        if (chemistry) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_species
                                dYs_ds(i) = 0._wp
                            end do
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do j = 0, buff_size

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, contxe
                                dalpha_rho_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, i)* &
                                                   fd_coef_${XYZ}$ (j, cbc_loc) + &
                                                   dalpha_rho_ds(i)
                            end do
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_dims
                                dvel_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, contxe + i)* &
                                             fd_coef_${XYZ}$ (j, cbc_loc) + &
                                             dvel_ds(i)
                            end do

                            dpres_ds = q_prim_rs${XYZ}$_vf(j, k, r, E_idx)* &
                                       fd_coef_${XYZ}$ (j, cbc_loc) + &
                                       dpres_ds
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, advxe - E_idx
                                dadv_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, E_idx + i)* &
                                             fd_coef_${XYZ}$ (j, cbc_loc) + &
                                             dadv_ds(i)
                            end do

                            if (chemistry) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 1, num_species
                                    dYs_ds(i) = q_prim_rs${XYZ}$_vf(j, k, r, chemxb - 1 + i)* &
                                                fd_coef_${XYZ}$ (j, cbc_loc) + &
                                                dYs_ds(i)
                                end do
                            end if
                        end do

                        ! First-Order Temporal Derivatives of Primitive Variables
                        lambda(1) = vel(dir_idx(1)) - c
                        lambda(2) = vel(dir_idx(1))
                        lambda(3) = vel(dir_idx(1)) + c

                        Ma = vel(dir_idx(1))/c

                        if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_SLIP_WALL) .or. &
                            (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_SLIP_WALL)) then
                            call s_compute_slip_wall_L(lambda, L, rho, c, dpres_ds, dvel_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_NR_SUB_BUFFER) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_NR_SUB_BUFFER)) then
                            call s_compute_nonreflecting_subsonic_buffer_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_NR_SUB_INFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_NR_SUB_INFLOW)) then
                            call s_compute_nonreflecting_subsonic_inflow_L(lambda, L, rho, c, dpres_ds, dvel_ds)
                            ! Add GRCBC for Subsonic Inflow
                            if (bc_${XYZ}$%grcbc_in) then
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = 2, momxb
                                    L(2) = c**3._wp*Ma*(alpha_rho(i - 1) - alpha_rho_in(i - 1, ${CBC_DIR}$))/Del_in(${CBC_DIR}$) - c*Ma*(pres - pres_in(${CBC_DIR}$))/Del_in(${CBC_DIR}$)
                                end do
                                if (n > 0) then
                                    L(momxb + 1) = c*Ma*(vel(dir_idx(2)) - vel_in(${CBC_DIR}$, dir_idx(2)))/Del_in(${CBC_DIR}$)
                                    if (p > 0) then
                                        L(momxb + 2) = c*Ma*(vel(dir_idx(3)) - vel_in(${CBC_DIR}$, dir_idx(3)))/Del_in(${CBC_DIR}$)
                                    end if
                                end if
                                $:GPU_LOOP(parallelism='[seq]')
                                do i = E_idx, advxe - 1
                                    L(i) = c*Ma*(adv_local(i + 1 - E_idx) - alpha_in(i + 1 - E_idx, ${CBC_DIR}$))/Del_in(${CBC_DIR}$)
                                end do
                                L(advxe) = rho*c**2._wp*(1._wp + Ma)*(vel(dir_idx(1)) + vel_in(${CBC_DIR}$, dir_idx(1))*sign(1, cbc_loc))/Del_in(${CBC_DIR}$) + c*(1._wp + Ma)*(pres - pres_in(${CBC_DIR}$))/Del_in(${CBC_DIR}$)
                            end if
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_NR_SUB_OUTFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_NR_SUB_OUTFLOW)) then
                            call s_compute_nonreflecting_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
                            ! Add GRCBC for Subsonic Outflow (Pressure)
                            if (bc_${XYZ}$%grcbc_out) then
                                L(advxe) = c*(1._wp - Ma)*(pres - pres_out(${CBC_DIR}$))/Del_out(${CBC_DIR}$)

                                ! Add GRCBC for Subsonic Outflow (Normal Velocity)
                                if (bc_${XYZ}$%grcbc_vel_out) then
                                    L(advxe) = L(advxe) + rho*c**2._wp*(1._wp - Ma)*(vel(dir_idx(1)) + vel_out(${CBC_DIR}$, dir_idx(1))*sign(1, cbc_loc))/Del_out(${CBC_DIR}$)
                                end if
                            end if
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_FF_SUB_OUTFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_FF_SUB_OUTFLOW)) then
                            call s_compute_force_free_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_CP_SUB_OUTFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_CP_SUB_OUTFLOW)) then
                            call s_compute_constant_pressure_subsonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_SUP_INFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_SUP_INFLOW)) then
                            call s_compute_supersonic_inflow_L(L)
                        else if ((cbc_loc == -1 .and. bc${XYZ}$b == BC_CHAR_SUP_OUTFLOW) .or. &
                                 (cbc_loc == 1 .and. bc${XYZ}$e == BC_CHAR_SUP_OUTFLOW)) then
                            call s_compute_supersonic_outflow_L(lambda, L, rho, c, mf, dalpha_rho_ds, dpres_ds, dvel_ds, dadv_ds, dYs_ds)
                        end if

                        ! Be careful about the cylindrical coordinate!
                        if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                            dpres_dt = -5.e-1_wp*(L(advxe) + L(1)) + rho*c*c*vel(dir_idx(1)) &
                                       /y_cc(n)
                        else
                            dpres_dt = -5.e-1_wp*(L(advxe) + L(1))
                        end if

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            dalpha_rho_dt(i) = &
                                -(L(i + 1) - mf(i)*dpres_dt)/(c*c)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            dvel_dt(dir_idx(i)) = dir_flg(dir_idx(i))* &
                                                  (L(1) - L(advxe))/(2._wp*rho*c) + &
                                                  (dir_flg(dir_idx(i)) - 1._wp)* &
                                                  L(momxb + i - 1)
                        end do

                        vel_dv_dt_sum = 0._wp
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, num_dims
                            vel_dv_dt_sum = vel_dv_dt_sum + vel(i)*dvel_dt(i)
                        end do

                        if (chemistry) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_species
                                dYs_dt(i) = -1._wp*L(chemxb + i - 1)
                            end do
                        end if

                        ! The treatment of void fraction source is unclear
                        if (cyl_coord .and. cbc_dir == 2 .and. cbc_loc == 1) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, advxe - E_idx
                                dadv_dt(i) = -L(momxe + i) !+ adv_local(i) * vel(dir_idx(1))/y_cc(n)
                            end do
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, advxe - E_idx
                                dadv_dt(i) = -L(momxe + i)
                            end do
                        end if

                        drho_dt = 0._wp; dgamma_dt = 0._wp; dpi_inf_dt = 0._wp; dqv_dt = 0._wp

                        if (model_eqns == 1) then
                            drho_dt = dalpha_rho_dt(1)
                            dgamma_dt = dadv_dt(1)
                            dpi_inf_dt = dadv_dt(2)
                        else
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_fluids
                                drho_dt = drho_dt + dalpha_rho_dt(i)
                                dgamma_dt = dgamma_dt + dadv_dt(i)*gammas(i)
                                dpi_inf_dt = dpi_inf_dt + dadv_dt(i)*pi_infs(i)
                                dqv_dt = dqv_dt + dalpha_rho_dt(i)*qvs(i)
                            end do
                        end if

                        ! flux_rs_vf_l and flux_src_rs_vf_l at j = -1/2
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = 1, contxe
                            flux_rs${XYZ}$_vf_l(-1, k, r, i) = flux_rs${XYZ}$_vf_l(0, k, r, i) &
                                                               + ds(0)*dalpha_rho_dt(i)
                        end do

                        $:GPU_LOOP(parallelism='[seq]')
                        do i = momxb, momxe
                            flux_rs${XYZ}$_vf_l(-1, k, r, i) = flux_rs${XYZ}$_vf_l(0, k, r, i) &
                                                               + ds(0)*(vel(i - contxe)*drho_dt &
                                                                        + rho*dvel_dt(i - contxe))
                        end do

                        if (chemistry) then
                            ! Evolution of LODI equation of energy for real gases adjusted to perfect gas, doi:10.1006/jcph.2002.6990
                            call get_species_enthalpies_rt(T, h_k)
                            sum_Enthalpies = 0._wp
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_species
                                h_k(i) = h_k(i)*gas_constant/molecular_weights(i)*T
                                sum_Enthalpies = sum_Enthalpies + (rho*h_k(i) - pres*Mw/molecular_weights(i)*Cp/R_gas)*dYs_dt(i)
                            end do
                            flux_rs${XYZ}$_vf_l(-1, k, r, E_idx) = flux_rs${XYZ}$_vf_l(0, k, r, E_idx) &
                                                                   + ds(0)*((E/rho + pres/rho)*drho_dt + rho*vel_dv_dt_sum + Cp*T*L(2)/(c*c) + sum_Enthalpies)
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = 1, num_species
                                flux_rs${XYZ}$_vf_l(-1, k, r, i - 1 + chemxb) = flux_rs${XYZ}$_vf_l(0, k, r, chemxb + i - 1) &
                                                                                + ds(0)*(drho_dt*Ys(i) + rho*dYs_dt(i))
                            end do
                        else
                            flux_rs${XYZ}$_vf_l(-1, k, r, E_idx) = flux_rs${XYZ}$_vf_l(0, k, r, E_idx) &
                                                                   + ds(0)*(pres*dgamma_dt &
                                                                            + gamma*dpres_dt &
                                                                            + dpi_inf_dt &
                                                                            + dqv_dt &
                                                                            + rho*vel_dv_dt_sum &
                                                                            + 5.e-1_wp*drho_dt*vel_K_sum)
                        end if

                        if (riemann_solver == 1) then
                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf_l(-1, k, r, i) = 0._wp
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_src_rs${XYZ}$_vf_l(-1, k, r, i) = &
                                    1._wp/max(abs(vel(dir_idx(1))), sgm_eps) &
                                    *sign(1._wp, vel(dir_idx(1))) &
                                    *(flux_rs${XYZ}$_vf_l(0, k, r, i) &
                                      + vel(dir_idx(1)) &
                                      *flux_src_rs${XYZ}$_vf_l(0, k, r, i) &
                                      + ds(0)*dadv_dt(i - E_idx))
                            end do

                        else

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_rs${XYZ}$_vf_l(-1, k, r, i) = flux_rs${XYZ}$_vf_l(0, k, r, i) + &
                                                                   ds(0)*dadv_dt(i - E_idx)
                            end do

                            $:GPU_LOOP(parallelism='[seq]')
                            do i = advxb, advxe
                                flux_src_rs${XYZ}$_vf_l(-1, k, r, i) = flux_src_rs${XYZ}$_vf_l(0, k, r, i)
                            end do

                        end if
                        ! END: flux_rs_vf_l and flux_src_rs_vf_l at j = -1/2

                    end do
                end do
            end if
        #:endfor

        ! END: FD2 or FD4 of RHS at j = 0

        ! The reshaping of outputted data and disssociation of the FD and PI
        ! coefficients, or CBC coefficients, respectively, based on selected
        ! CBC coordinate direction.
        call s_finalize_cbc(flux_vf, flux_src_vf)
    end subroutine s_cbc

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are required for the setup of the
        !!      selected CBC.
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
        !!  @param ix Index bound in the first coordinate direction
        !!  @param iy Index bound in the second coordinate direction
        !!  @param iz Index bound in the third coordinate direction
    subroutine s_initialize_cbc(q_prim_vf, flux_vf, flux_src_vf, &
                                ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: flux_vf, flux_src_vf

        type(int_bounds_info), intent(in) :: ix, iy, iz

        integer :: i, j, k, r !< Generic loop iterators

        ! Configuring the coordinate direction indexes and flags

        ! Determining the indicial shift based on CBC location

        ! END: Allocation/Association of Primitive and Flux Variables

        if (cbc_dir == 1) then
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1._wp, 0._wp, 0._wp/)
        elseif (cbc_dir == 2) then
            is1%beg = 0; is1%end = buff_size; is2 = ix; is3 = iz
            dir_idx = (/2, 1, 3/); dir_flg = (/0._wp, 1._wp, 0._wp/)
        else
            is1%beg = 0; is1%end = buff_size; is2 = iy; is3 = ix
            dir_idx = (/3, 1, 2/); dir_flg = (/0._wp, 0._wp, 1._wp/)
        end if

        dj = max(0, cbc_loc)
        $:GPU_UPDATE(device='[is1,is2,is3,dj]')
        $:GPU_UPDATE(device='[dir_idx,dir_flg]')

        ! Reshaping Inputted Data in x-direction
        if (cbc_dir == 1) then

            $:GPU_PARALLEL_LOOP(collapse=4)
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

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsx_vf(j, k, r, momxb) = &
                            q_prim_vf(momxb)%sf(dj*(m - 2*j) + j, k, r)* &
                            sign(1._wp, -1._wp*cbc_loc)
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsx_vf_l(j, k, r, i) = &
                                flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsx_vf_l(j, k, r, momxb) = &
                            flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsx_vf_l(j, k, r, i) = &
                                    flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsx_vf_l(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Inputted Data in x-direction

            ! Reshaping Inputted Data in y-direction
        elseif (cbc_dir == 2) then

            $:GPU_PARALLEL_LOOP(collapse=4)
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

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsy_vf(j, k, r, momxb + 1) = &
                            q_prim_vf(momxb + 1)%sf(k, dj*(n - 2*j) + j, r)* &
                            sign(1._wp, -1._wp*cbc_loc)
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsy_vf_l(j, k, r, i) = &
                                flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsy_vf_l(j, k, r, momxb + 1) = &
                            flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsy_vf_l(j, k, r, i) = &
                                    flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsy_vf_l(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Inputted Data in y-direction

            ! Reshaping Inputted Data in z-direction
        else

            $:GPU_PARALLEL_LOOP(collapse=4)
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

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = 0, buff_size
                        q_prim_rsz_vf(j, k, r, momxe) = &
                            q_prim_vf(momxe)%sf(r, k, dj*(p - 2*j) + j)* &
                            sign(1._wp, -1._wp*cbc_loc)
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_rsz_vf_l(j, k, r, i) = &
                                flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_rsz_vf_l(j, k, r, momxe) = &
                            flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_rsz_vf_l(j, k, r, i) = &
                                    flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_rsz_vf_l(j, k, r, advxb) = &
                                flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Inputted Data in z-direction

        ! Association of the procedural pointer to the appropriate procedure
        ! that will be utilized in the evaluation of L variables for the CBC

    end subroutine s_initialize_cbc

    !>  Deallocation and/or the disassociation procedures that
        !!      are necessary in order to finalize the CBC application
        !!  @param flux_vf Cell-boundary-average fluxes
        !!  @param flux_src_vf Cell-boundary-average flux sources
    subroutine s_finalize_cbc(flux_vf, flux_src_vf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(inout) :: flux_vf, flux_src_vf

        integer :: i, j, k, r !< Generic loop iterators

        ! Determining the indicial shift based on CBC location
        dj = max(0, cbc_loc)
        $:GPU_UPDATE(device='[dj]')

        ! Reshaping Outputted Data in x-direction
        if (cbc_dir == 1) then

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_rsx_vf_l(j, k, r, i)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do
            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                            flux_rsx_vf_l(j, k, r, momxb)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                    flux_src_rsx_vf_l(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(dj*((m - 1) - 2*j) + j, k, r) = &
                                flux_src_rsx_vf_l(j, k, r, advxb)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if
            ! END: Reshaping Outputted Data in x-direction

            ! Reshaping Outputted Data in y-direction
        elseif (cbc_dir == 2) then

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_rsy_vf_l(j, k, r, i)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxb + 1)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                            flux_rsy_vf_l(j, k, r, momxb + 1)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                    flux_src_rsy_vf_l(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(k, dj*((n - 1) - 2*j) + j, r) = &
                                flux_src_rsy_vf_l(j, k, r, advxb)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if

            ! END: Reshaping Outputted Data in y-direction

            ! Reshaping Outputted Data in z-direction
        else

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, flux_cbc_index
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_rsz_vf_l(j, k, r, i)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end do

            $:GPU_PARALLEL_LOOP(collapse=3)
            do r = is3%beg, is3%end
                do k = is2%beg, is2%end
                    do j = -1, buff_size
                        flux_vf(momxe)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                            flux_rsz_vf_l(j, k, r, momxe)
                    end do
                end do
            end do

            if (riemann_solver == 1) then
                $:GPU_PARALLEL_LOOP(collapse=4)
                do i = advxb, advxe
                    do r = is3%beg, is3%end
                        do k = is2%beg, is2%end
                            do j = -1, buff_size
                                flux_src_vf(i)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                    flux_src_rsz_vf_l(j, k, r, i)
                            end do
                        end do
                    end do
                end do
            else
                $:GPU_PARALLEL_LOOP(collapse=3)
                do r = is3%beg, is3%end
                    do k = is2%beg, is2%end
                        do j = -1, buff_size
                            flux_src_vf(advxb)%sf(r, k, dj*((p - 1) - 2*j) + j) = &
                                flux_src_rsz_vf_l(j, k, r, advxb)* &
                                sign(1._wp, -1._wp*cbc_loc)
                        end do
                    end do
                end do
            end if

        end if
        ! END: Reshaping Outputted Data in z-direction

    end subroutine s_finalize_cbc

    ! Detext if the problem has any characteristic boundary conditions
    pure elemental subroutine s_any_cbc_boundaries(toggle)

        logical, intent(inout) :: toggle

        toggle = .false.

        #:for BC in {-5, -6, -7, -8, -9, -10, -11, -12, -13}
            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == ${BC}$)) then
                toggle = .true.
            end if
        #:endfor

    end subroutine s_any_cbc_boundaries

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_cbc_module

        logical :: is_cbc

        call s_any_cbc_boundaries(is_cbc)

        if (is_cbc .eqv. .false.) return

        ! Deallocating the cell-average primitive variables
        @:DEALLOCATE(q_prim_rsx_vf)
        if (weno_order > 1 .or. muscl_order > 1) then
            @:DEALLOCATE(F_rsx_vf, F_src_rsx_vf)
        end if
        @:DEALLOCATE(flux_rsx_vf_l, flux_src_rsx_vf_l)

        if (n > 0) then
            @:DEALLOCATE(q_prim_rsy_vf)
            if (weno_order > 1 .or. muscl_order > 1) then
                @:DEALLOCATE(F_rsy_vf, F_src_rsy_vf)
            end if
            @:DEALLOCATE(flux_rsy_vf_l, flux_src_rsy_vf_l)
        end if
        if (p > 0) then
            @:DEALLOCATE(q_prim_rsz_vf)
            if (weno_order > 1 .or. muscl_order > 1) then
                @:DEALLOCATE(F_rsz_vf, F_src_rsz_vf)
            end if
            @:DEALLOCATE(flux_rsz_vf_l, flux_src_rsz_vf_l)
        end if

        ! Deallocating the cell-width distribution in the s-direction
        @:DEALLOCATE(ds)

        ! Deallocating GRCBC inputs
        @:DEALLOCATE(vel_in, vel_out, pres_in, pres_out, Del_in, Del_out, alpha_rho_in, alpha_in)

        ! Deallocating CBC Coefficients in x-direction
        if (any((/bc_x%beg, bc_x%end/) <= -5) .and. any((/bc_x%beg, bc_x%end/) >= -13)) then
            @:DEALLOCATE(fd_coef_x)
            if (weno_order > 1 .or. muscl_order > 1) then
                @:DEALLOCATE(pi_coef_x)
            end if
        end if

        ! Deallocating CBC Coefficients in y-direction
        if (n > 0 .and. any((/bc_y%beg, bc_y%end/) <= -5) .and. &
            any((/bc_y%beg, bc_y%end/) >= -13 .and. bc_y%beg /= -14)) then
            @:DEALLOCATE(fd_coef_y)
            if (weno_order > 1 .or. muscl_order > 1) then
                @:DEALLOCATE(pi_coef_y)
            end if
        end if

        ! Deallocating CBC Coefficients in z-direction
        if (p > 0 .and. any((/bc_z%beg, bc_z%end/) <= -5) .and. any((/bc_z%beg, bc_z%end/) >= -13)) then
            @:DEALLOCATE(fd_coef_z)
            if (weno_order > 1 .or. muscl_order > 1) then
                @:DEALLOCATE(pi_coef_z)
            end if
        end if

    end subroutine s_finalize_cbc_module

end module m_cbc
