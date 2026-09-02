!>
!! @file
!! @brief Contains module m_global_parameters_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared global parameters and equation-index setup for all three executables. Each per-target m_global_parameters uses this
!! module (default-public, so all symbols re-export to downstream consumers via two-hop use-association).
module m_global_parameters_common

#ifdef MFC_MPI
    use mpi
#endif

    use m_derived_types
    use m_thermochem, only: num_species
    use m_constants, only: model_eqns_gamma_law, model_eqns_5eq, model_eqns_6eq, recon_type_weno, recon_type_muscl, name_len, &
        & dflt_int, dflt_real, eos_stiffened_gas, eos_ideal_gas

    implicit none

    ! All namelist-bound scalar and array declarations (per-target, regenerated at build time by the ninja custom command)
    #:include 'generated_decls.fpp'

    ! Case-optimization declarations: parameters (under MFC_CASE_OPTIMIZATION) or plain variables for
    ! num_dims, num_vels, weno_polyn, muscl_polyn, weno_num_stencils, wenojs, igr, etc.
    #:include 'generated_case_opt_decls.fpp'

    !> @name Annotations of the structure of the state and flux vectors in terms of the size and configuration of the system of
    !! equations
    !> @{
    integer            :: sys_size  !< Number of unknowns in system of equations
    type(eqn_idx_info) :: eqn_idx   !< All conserved-variable equation index ranges and scalars
    !> @}

    !> @name Chemistry modeling (Fypp compile-time constant; same value in all targets)
    !> @{
    logical, parameter :: chemistry = .${chemistry}$.
    !> @}

    !> @name Hypoelastic shear stress state (identical across all three executables)
    !> @{
    integer                  :: shear_num              !< Number of shear stress components
    integer, dimension(3)    :: shear_indices          !< Indices of the stress components that represent shear stress
    integer                  :: shear_BC_flip_num      !< Number of shear stress components to reflect for boundary conditions
    integer, dimension(3, 2) :: shear_BC_flip_indices  !< Shear stress BC reflection indices (1:3, 1:shear_BC_flip_num)
    !> @}

    !> @name Material properties derived from fluid_pp
    !> @{ One declaration is shared by all executables and initialized by m_variables_conversion after the case parameters have been
    !! read.
    !> gammas is the stored form 1/(gamma - 1), not the ratio of specific heats; isentrope_n and isentrope_B are the same EOS
    !! written as p + B = const*rho**n.
    real(wp), allocatable, dimension(:) :: gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps
    $:GPU_DECLARE(create='[gammas, isentrope_n, pi_infs, isentrope_B, cvs, qvs, qvps]')
    !> @}

    !> @name Fluids participating in shear and bulk viscosity
    !> @{
    integer, dimension(2)                :: Re_size = 0
    integer                              :: Re_size_max = 0
    integer, allocatable, dimension(:,:) :: Re_idx
    $:GPU_DECLARE(create='[Re_size, Re_size_max, Re_idx]')
    !> @}

    !> Minimum cell widths. These are distinct from the per-cell width arrays named dx, dy, and dz in simulation and post-process.
    real(wp) :: dx_min, dy_min, dz_min

    $:GPU_DECLARE(create='[sys_size, eqn_idx]')
    $:GPU_DECLARE(create='[shear_num, shear_indices, shear_BC_flip_num, shear_BC_flip_indices]')

    !> @name Processor coordinates and parallel-IO addressing (identical declaration across all three targets)
    !> @{
    integer, allocatable, dimension(:) :: proc_coords      !< Processor coordinates in MPI_CART_COMM
    integer, allocatable, dimension(:) :: start_idx        !< Starting cell-center index of local processor in global grid
    integer                            :: num_procs_x = 1  !< Number of MPI ranks in x-direction
    integer                            :: num_procs_y = 1  !< Number of MPI ranks in y-direction
    integer                            :: num_procs_z = 1  !< Number of MPI ranks in z-direction
    !> @}

    !> @name MPI info for parallel IO with Lustre file systems (identical across all three targets)
    !> @{
    character(len=name_len) :: mpiiofs
    integer                 :: mpi_info_int
    !> @}

contains

    !> Initialize equation-index state (eqn_idx and sys_size) from the namelist parameters. This is the shared skeleton: it covers
    !! the model_eqns dispatch, all eqn_idx field assignments, and the hypoelastic/surface-tension/chemistry extensions.
    !!
    !! @param nmom_in  Number of carried moments per R0 location (per-target: pre/post pass an
    !!   integer variable; sim passes its integer parameter nmom = 6).  Used only in the 5eq
    !!   qbmm bubble-index calculation (eqn_idx%bub%end = eqn_idx%adv%end + nb_in*nmom_in).
    !!
    !! Per-target callers are responsible for the following after this call:
    !!   - qbmm_idx allocations and fills (diverge between pre vs sim/post)
    !!   - sim-only: gam = bub_pp%gam_g, nmomsp/nmomtot, Re_idx allocation, GPU_UPDATE calls
    !!   - post-only: beta_idx increment (bubbles_lagrange path), offset/grid allocations
    impure subroutine s_initialize_eqn_idx(nmom_in, nb_in, six_eqn_alf_is_advected)

        integer, intent(in) :: nmom_in
        integer, intent(in) :: nb_in
        logical, intent(in) :: six_eqn_alf_is_advected

        ! Gamma/Pi_inf Model

        if (model_eqns == model_eqns_gamma_law) then
            ! Annotating structure of the state and flux vectors belonging to the system of
            ! equations defined by the selected number of spatial dimensions and the gamma/pi_inf model
            eqn_idx%cont%beg = 1
            eqn_idx%cont%end = eqn_idx%cont%beg
            eqn_idx%mom%beg = eqn_idx%cont%end + 1
            eqn_idx%mom%end = eqn_idx%cont%end + num_vels
            eqn_idx%E = eqn_idx%mom%end + 1
            eqn_idx%adv%beg = eqn_idx%E + 1
            eqn_idx%adv%end = eqn_idx%adv%beg + 1
            eqn_idx%gamma = eqn_idx%adv%beg
            eqn_idx%pi_inf = eqn_idx%adv%end
            sys_size = eqn_idx%adv%end

            ! Volume Fraction Model (5-equation model)
        else if (model_eqns == model_eqns_5eq) then
            ! Annotating structure of the state and flux vectors belonging to the system of
            ! equations defined by the selected number of spatial dimensions and the volume fraction model
            eqn_idx%cont%beg = 1
            eqn_idx%cont%end = num_fluids
            eqn_idx%mom%beg = eqn_idx%cont%end + 1
            eqn_idx%mom%end = eqn_idx%cont%end + num_vels
            eqn_idx%E = eqn_idx%mom%end + 1

            if (igr) then
                ! IGR: volume fractions after energy (N-1 for N fluids; skipped when num_fluids=1)
                eqn_idx%adv%beg = eqn_idx%E + 1
                eqn_idx%adv%end = eqn_idx%E + num_fluids - 1
            else
                ! WENO/MUSCL + Riemann tracks a total of (N) volume fractions for N fluids
                eqn_idx%adv%beg = eqn_idx%E + 1
                eqn_idx%adv%end = eqn_idx%E + num_fluids
            end if

            sys_size = eqn_idx%adv%end

            if (bubbles_euler) then
                eqn_idx%alf = eqn_idx%adv%end
            else
                eqn_idx%alf = 1
            end if

            if (bubbles_euler) then
                eqn_idx%bub%beg = sys_size + 1
                if (qbmm) then
                    eqn_idx%bub%end = eqn_idx%adv%end + nb_in*nmom_in
                else
                    if (.not. polytropic) then
                        eqn_idx%bub%end = sys_size + 4*nb_in
                    else
                        eqn_idx%bub%end = sys_size + 2*nb_in
                    end if
                end if
                sys_size = eqn_idx%bub%end

                if (adv_n) then
                    eqn_idx%n = eqn_idx%bub%end + 1
                    sys_size = eqn_idx%n
                end if
            end if

            if (mhd) then
                eqn_idx%B%beg = sys_size + 1
                if (n == 0) then
                    eqn_idx%B%end = sys_size + 2  ! 1D: By, Bz
                else
                    eqn_idx%B%end = sys_size + 3  ! 2D/3D: Bx, By, Bz
                end if
                sys_size = eqn_idx%B%end
            end if

            ! Volume Fraction Model (6-equation model)
        else if (model_eqns == model_eqns_6eq) then
            ! Annotating structure of the state and flux vectors belonging to the system of
            ! equations defined by the selected number of spatial dimensions and the volume fraction model
            eqn_idx%cont%beg = 1
            eqn_idx%cont%end = num_fluids
            eqn_idx%mom%beg = eqn_idx%cont%end + 1
            eqn_idx%mom%end = eqn_idx%cont%end + num_vels
            eqn_idx%E = eqn_idx%mom%end + 1
            eqn_idx%adv%beg = eqn_idx%E + 1
            eqn_idx%adv%end = eqn_idx%E + num_fluids
            if (six_eqn_alf_is_advected) eqn_idx%alf = eqn_idx%adv%end
            eqn_idx%int_en%beg = eqn_idx%adv%end + 1
            eqn_idx%int_en%end = eqn_idx%adv%end + num_fluids
            sys_size = eqn_idx%int_en%end
        end if

        if (model_eqns == model_eqns_5eq .or. model_eqns == model_eqns_6eq) then
            if (hypoelasticity) then
                eqn_idx%stress%beg = sys_size + 1
                eqn_idx%stress%end = sys_size + (num_dims*(num_dims + 1))/2
                if (cyl_coord) eqn_idx%stress%end = eqn_idx%stress%end + 1
                ! number of stresses is 1 in 1D, 3 in 2D, 4 in 2D-Axisym, 6 in 3D
                sys_size = eqn_idx%stress%end

                ! shear stress index is 2 for 2D and 2,4,5 for 3D. Readers test the whole array
                ! rather than the first shear_num entries, so unused slots must not be garbage.
                shear_indices = 0
                if (num_dims == 1) then
                    shear_num = 0
                else if (num_dims == 2) then
                    shear_num = 1
                    shear_indices(1) = eqn_idx%stress%beg - 1 + 2
                    shear_BC_flip_num = 1
                    shear_BC_flip_indices(1:2,1) = shear_indices(1)
                    ! Both x-dir and y-dir: flip tau_xy only
                else if (num_dims == 3) then
                    shear_num = 3
                    shear_indices(1:3) = eqn_idx%stress%beg - 1 + (/2, 4, 5/)
                    shear_BC_flip_num = 2
                    shear_BC_flip_indices(1,1:2) = shear_indices((/1, 2/))
                    shear_BC_flip_indices(2,1:2) = shear_indices((/1, 3/))
                    shear_BC_flip_indices(3,1:2) = shear_indices((/2, 3/))
                    ! x-dir: flip tau_xy and tau_xz; y-dir: flip tau_xy and tau_yz; z-dir: flip tau_xz and tau_yz
                end if
            end if

            if (surface_tension) then
                eqn_idx%c = sys_size + 1
                sys_size = eqn_idx%c
            end if

            if (cont_damage) then
                eqn_idx%damage = sys_size + 1
                sys_size = eqn_idx%damage
            end if

            if (hyper_cleaning) then
                eqn_idx%psi = sys_size + 1
                sys_size = eqn_idx%psi
            end if
        end if

        if (chemistry) then
            eqn_idx%species%beg = sys_size + 1
            eqn_idx%species%end = sys_size + num_species
            sys_size = eqn_idx%species%end
        end if

    end subroutine s_initialize_eqn_idx

    !> Configure MPI parallel I/O settings and allocate processor coordinate arrays. Shared across all three executables;
    !! num_dims/num_vels are computed here for pre/post unconditionally and for sim only when not case-optimized (in which case they
    !! are compile-time parameters). Callers must have already populated n and p (grid dimensions).
    impure subroutine s_initialize_parallel_io_common

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
#endif

        ! Under case-optimization, num_dims and num_vels are compile-time parameters; skip assignment.
        #:if not MFC_CASE_OPTIMIZATION
            num_dims = 1 + min(1, n) + min(1, p)

            if (mhd) then
                num_vels = 3
            else
                num_vels = num_dims
            end if
        #:endif

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

#ifdef MFC_MPI
        ! Option for Lustre file system (Darter/Comet/Stampede)
        write (mpiiofs, '(A)') '/lustre_'
        mpiiofs = trim(mpiiofs)
        call MPI_INFO_CREATE(mpi_info_int, ierr)
        call MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

        ! Option for UNIX file system (Hooke/Thomson) WRITE(mpiiofs, '(A)') '/ufs_' mpiiofs = TRIM(mpiiofs) mpi_info_int =
        ! MPI_INFO_NULL

        allocate (start_idx(1:num_dims))
#endif

    end subroutine s_initialize_parallel_io_common

    !> Shared finalize core: deallocate proc_coords and start_idx. Per-target finalize routines call this first, then handle their
    !! own extras (qbmm_idx, grid arrays, MPI_IO_DATA null/dealloc - those reference per-target typed variables and stay
    !! per-target).
    impure subroutine s_finalize_global_parameters_common

        deallocate (proc_coords)

#ifdef MFC_MPI
        if (parallel_io) then
            deallocate (start_idx)
        end if
#endif

    end subroutine s_finalize_global_parameters_common

    !> Assign default values to the user-input parameters that are shared across all three executables (pre_process, simulation,
    !! post_process). Per-target defaults (bc_io, cfl_dt, precision, nb, chem_params, bub_pp scalar companions, fluid_pp loop, patch
    !! arrays, output flags) remain in the per-target s_assign_default_values_to_user_inputs routines, which call this subroutine
    !! first, then call s_update_cell_bounds(cells_bounds, m, n, p) (cells_bounds is per-target), then apply their own assignments.
    impure subroutine s_assign_common_defaults

        ! Logistics
        case_dir = '.'

        ! Computational domain parameters (m/n/p set here; caller must call s_update_cell_bounds after)
        m = dflt_int; n = 0; p = 0

        cyl_coord = .false.

        ! CFL adaptive time-stepping flags
        cfl_adap_dt = .false.
        cfl_const_dt = .false.

        ! Time-stepping bookkeeping
        n_start = dflt_int
        t_step_start = dflt_int

        ! Simulation algorithm
        model_eqns = dflt_int
        relax = .false.
        relax_model = dflt_int
        hypoelasticity = .false.
        cont_damage = .false.
        hyper_cleaning = .false.

        ! Condensed-phase reactive burn
        reactive_burn = .false.
        rburn%k = dflt_real
        rburn%pign = dflt_real
        rburn%pref = dflt_real
        rburn%n = dflt_real
        rburn%ta = 0._wp

        ! Case-optimization params: under case-opt these are compile-time constants in sim (skip assignment); in pre/post
        ! MFC_CASE_OPTIMIZATION is always False so the block always executes there.
        #:if not MFC_CASE_OPTIMIZATION
            recon_type = recon_type_weno
            weno_order = dflt_int
            muscl_order = dflt_int
            num_fluids = dflt_int
            igr = .false.
            igr_order = dflt_int
            mhd = .false.
            relativity = .false.
            viscous = .false.
        #:endif

        ! Not a case-optimization parameter, so it stays a runtime variable in every build.
        riemann_solver = dflt_int

        ! Tait EOS

        ! Bubble modeling flags and parameters
        R0ref = dflt_real
        bubbles_euler = .false.
        polydisperse = .false.
        poly_sigma = dflt_real
        qbmm = .false.
        surface_tension = .false.
        adv_n = .false.
        sigma = dflt_real
        bubbles_lagrange = .false.

        ! Immersed boundaries
        ib = .false.
        num_ibs = dflt_int

        ! MHD (background field)
        Bx0 = dflt_real

        ! Output and I/O options
        parallel_io = .false.
        file_per_process = .false.
        down_sample = .false.
        fft_wrt = .false.

        ! Mixture conversion and sound-speed behavior
        avg_state = dflt_int
        alt_soundspeed = .false.
        mixture_err = .false.
        sigR = dflt_real
        dx_min = dflt_real
        dy_min = dflt_real
        dz_min = dflt_real

    end subroutine s_assign_common_defaults

end module m_global_parameters_common
