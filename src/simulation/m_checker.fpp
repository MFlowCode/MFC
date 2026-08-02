!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic
    use m_constants, only: model_eqns_5eq, riemann_solver_hll, riemann_solver_hllc, riemann_solver_hlld, recon_type_weno, &
        & recon_type_muscl, muscl_order_first_order, riemann_solver_lax_friedrichs, wave_speeds_pressure

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        call s_check_inputs_time_stepping

        call s_check_inputs_hll_non_conservative

        call s_check_inputs_hypo_branch

        @:PROHIBIT(chemistry .and. chem_params%reaction_substeps < 0, &
                   & "chem_params%reaction_substeps must be >= 0 (0 = reaction source in the flow RHS; > 0 = operator-split sub-stepping)")

        @:PROHIBIT(chemistry .and. igr .and. chem_params%reaction_substeps > 0, &
                   & "operator-split reaction sub-stepping (reaction_substeps > 0) is not supported with igr: the reactor reads the post-flow (rho, e, T) state, which the IGR update path does not guarantee")

        @:PROHIBIT(chemistry .and. chem_params%adap_substeps .and. chem_params%reaction_substeps < 1, &
                   & "chem_params%adap_substeps requires reaction_substeps >= 1 (the operator-split floor)")

        @:PROHIBIT(chemistry .and. chem_params%adap_substeps &
                   & .and. chem_params%reaction_substeps_max < chem_params%reaction_substeps, &
                   & "chem_params%reaction_substeps_max must be >= reaction_substeps when adap_substeps = T")

        ! Chemistry with Euler bubbles is not currently supported: the IBM image-point
        ! interpolation selects the bubbles/QBMM branch before the chemistry branch, so the
        ! species state (Ys_IP) is not carried when both are enabled. Disallow until implemented.
        @:PROHIBIT(chemistry .and. (bubbles_euler .or. qbmm), &
                   & "chemistry is not currently supported with Euler bubbles (bubbles_euler/qbmm)")

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")
        @:PROHIBIT(many_ib_patch_parallelism .and. .not. ib, "many_ib_patch_parallelism requires ib to be enabled")

        @:PROHIBIT(bf_spatial_support .and. (n == 0 .or. p /= 0), &
                   & "bf_spatial_support is implemented for 2D only (it forces mom%beg and mom%beg+1)")

        ! Condensed-phase reactive burn assumes exactly two fluids (reactant=1, product=2) that share the
        ! stiffened-gas EOS and differ only in qv; violating these silently corrupts the mass/energy balance.
        @:PROHIBIT(reactive_burn .and. num_fluids /= 2, "reactive_burn requires num_fluids = 2 (reactant then product)")
        @:PROHIBIT(reactive_burn .and. .not. f_approx_equal(fluid_pp(1)%gamma, fluid_pp(2)%gamma), &
                   & "reactive_burn requires fluid_pp(1)%gamma == fluid_pp(2)%gamma (reactant and product share the EOS)")
        @:PROHIBIT(reactive_burn .and. .not. f_approx_equal(fluid_pp(1)%pi_inf, fluid_pp(2)%pi_inf), &
                   & "reactive_burn requires fluid_pp(1)%pi_inf == fluid_pp(2)%pi_inf (reactant and product share the EOS)")
        @:PROHIBIT(reactive_burn .and. fluid_pp(1)%qv <= fluid_pp(2)%qv, &
                   & "reactive_burn requires fluid_pp(1)%qv > fluid_pp(2)%qv (reactant releases energy on conversion to product)")
        @:PROHIBIT(reactive_burn .and. rburn%pref <= 0._wp, &
                   & "reactive_burn requires rburn%pref > 0 (it normalizes the pressure drive (p - rburn%pign)/rburn%pref and is used as a divisor)")
        ! The rate uses rburn%k, rburn%pign, rburn%n directly; each defaults to the sentinel dflt_real,
        ! so an unset value silently produces spurious ignition (pign), NaN via drive**n (n), or a
        ! backward reaction (k). Require each to be set to a physical value.
        @:PROHIBIT(reactive_burn .and. rburn%k <= 0._wp, &
                   & "reactive_burn requires rburn%k > 0 (rate coefficient [1/s]; unset defaults to a negative sentinel that runs the reaction backward)")
        @:PROHIBIT(reactive_burn .and. f_is_default(rburn%pign), &
                   & "reactive_burn requires rburn%pign to be set (ignition pressure threshold [Pa]; unset defaults to a negative sentinel, so the reactant ignites everywhere from t = 0)")
        @:PROHIBIT(reactive_burn .and. rburn%n < 0._wp, &
                   & "reactive_burn requires rburn%n >= 0 (pressure-drive exponent; unset defaults to a negative sentinel, so drive**n overflows to Inf and the field goes NaN)")
        @:PROHIBIT(reactive_burn .and. model_eqns /= 2 .and. model_eqns /= 3, &
                   & "reactive_burn requires model_eqns = 2 or 3 (the 5-equation pressure-equilibrium or 6-equation multi-fluid model)")
        @:PROHIBIT(reactive_burn .and. rburn%ta < 0._wp, &
                   & "reactive_burn requires rburn%ta >= 0 (activation temperature [K]; 0 disables the Arrhenius factor)")
        @:PROHIBIT(reactive_burn .and. rburn%ta > 0._wp .and. fluid_pp(1)%cv <= 0._wp, &
                   & "reactive_burn with rburn%ta > 0 requires fluid_pp(1)%cv > 0 (the reactant temperature T = (p + pi_inf)/((gamma - 1) cv rho) needs a physical heat capacity; cv = 0 silently disables the Arrhenius factor)")

        if (ib .and. chemistry) then
            call s_check_inputs_ib_injection
        end if

        if (num_particle_clouds > 0) then
            call s_check_inputs_particle_clouds
        end if

        if (synthetic_turbulence) then
            call s_check_inputs_synthetic_turbulence
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(muscl_order == muscl_order_first_order .and. int_comp > 0, &
                   & "int_comp requires muscl_order >= 2 (muscl_order=1 leaves the reconstruction workspace uninitialised)")

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    impure subroutine s_check_inputs_hll_non_conservative

        @:PROHIBIT((riemann_solver == riemann_solver_hll) .and. hll_u_interface .and. cyl_coord .and. p > 0, &
                   & "HLL Method 2 is not supported for 3D cylindrical geometry")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == riemann_solver_hll .and. (.not. hll_u_interface) .and. cyl_coord &
                   & .and. p == 0, "alt_soundspeed with HLL Method 1 is not supported for 2D axisymmetric geometry")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == riemann_solver_hll .and. cyl_coord .and. p > 0, &
                   & "alt_soundspeed with HLL is not currently supported for 3D cylindrical geometry")
        @:PROHIBIT(hll_u_interface .and. riemann_solver /= riemann_solver_hll, &
                   & "hll_u_interface requires the HLL Riemann solver (riemann_solver = 1)")
        @:PROHIBIT(hll_u_interface .and. mhd, &
                   & "HLL Method 2 does not support MHD (the MHD path zeroes the shared interface-velocity trace)")
        @:PROHIBIT(surface_tension .and. riemann_solver == riemann_solver_hll .and. (.not. hll_u_interface), &
                   & "surface_tension requires a shared interface-velocity representation (HLL Method 2 or HLLC)")
        @:PROHIBIT(surface_tension .and. riemann_solver == riemann_solver_lax_friedrichs, &
                   & "surface_tension requires a shared interface-velocity representation (HLL Method 2 or HLLC)")

    end subroutine s_check_inputs_hll_non_conservative

    impure subroutine s_check_inputs_hypo_branch

        @:PROHIBIT(hypoelasticity .and. cyl_coord .and. p > 0, "3D cylindrical hypoelasticity is not supported")

        ! Hypoelasticity solver restrictions
        @:PROHIBIT(hypoelasticity .and. model_eqns /= model_eqns_5eq, "hypoelasticity requires model_eqns = 2")
        @:PROHIBIT(hypoelasticity .and. riemann_solver /= riemann_solver_hll .and. riemann_solver /= riemann_solver_hllc &
                   & .and. riemann_solver /= riemann_solver_hlld, &
                   & "hypoelasticity requires HLL (1), HLLC (2), or HLLD (4) Riemann solver")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. n == 0, &
                   & "HLLD hypoelasticity requires at least 2D (n must be > 0)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. num_fluids /= 2, &
                   & "HLLD hypoelasticity currently requires exactly 2 fluid components")
        @:PROHIBIT(hypoelasticity .and. mhd, "MHD and hypoelasticity cannot be enabled together")
        @:PROHIBIT(hypoelasticity .and. bubbles_euler, "Hypoelasticity does not support Euler-Euler bubbles")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. viscous, &
                   & "HLLD hypoelasticity does not support viscous effects (the dual-pass omits the viscous source term)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. surface_tension, &
                   & "HLLD hypoelasticity does not support surface tension (the dual-pass omits the surface-tension source term)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. cont_damage, &
                   & "HLLD hypoelasticity does not support continuum damage (the dual-pass does not damage-scale the shear modulus)")
        @:PROHIBIT(hypoelasticity .and. riemann_solver == riemann_solver_hlld .and. chemistry, &
                   & "HLLD hypoelasticity does not support chemistry")
        @:PROHIBIT(riemann_solver == riemann_solver_hlld .and. (.not. mhd) .and. (.not. hypoelasticity), &
                   & "HLLD is only available for MHD or hypoelasticity")

        ! Feature flag prerequisites
        @:PROHIBIT(riemann_hypo_ADC .and. .not. hypoelasticity, "riemann_hypo_ADC requires hypoelasticity = T")
        @:PROHIBIT(riemann_hypo_ADC .and. riemann_solver /= riemann_solver_hllc .and. riemann_solver /= riemann_solver_hlld, &
                   & "riemann_hypo_ADC only applies to hypo HLLC/HLLD")
        @:PROHIBIT(riemann_hypo_ADC .and. (bubbles_euler .or. surface_tension .or. chemistry .or. cont_damage), &
                   & "riemann_hypo_ADC does not support bubbles, surface tension, chemistry, or continuum damage (the ADC HLL blend omits their flux components)")
        @:PROHIBIT(hypo_hll_interface_rhs .and. .not. hypoelasticity, "hypo_hll_interface_rhs requires hypoelasticity = T")
        @:PROHIBIT(hypo_hll_interface_rhs .and. riemann_solver /= riemann_solver_hll, &
                   & "hypo_hll_interface_rhs requires HLL Riemann solver (riemann_solver = 1)")
        @:PROHIBIT(alt_soundspeed .and. riemann_solver == riemann_solver_hlld .and. .not. hypoelasticity, &
                   & "alt_soundspeed with HLLD requires hypoelasticity = T")
        @:PROHIBIT(alt_soundspeed .and. num_fluids /= 2, &
                   & "alt_soundspeed requires exactly 2 fluid components (the Kapila K coefficient is a two-fluid closure)")
        @:PROHIBIT(cont_damage .and. alt_soundspeed, "Continuum damage does not support alt_soundspeed")
        @:PROHIBIT(hypoelasticity .and. igr, "Hypoelasticity is not compatible with IGR")
        @:PROHIBIT(hypoelasticity .and. fd_order /= 1 .and. fd_order /= 2 .and. fd_order /= 4, &
                   & "hypoelasticity requires fd_order to be set to 1, 2, or 4 (the finite-difference coefficients are initialized unconditionally)")
        @:PROHIBIT(hypoelasticity .and. wave_speeds == wave_speeds_pressure, &
                   & "Pressure-based wave speeds (wave_speeds = 2) omit the elastic longitudinal speed and are not supported with hypoelasticity")
        @:PROHIBIT(riemann_solver == riemann_solver_hlld .and. wave_speeds == wave_speeds_pressure, &
                   & "HLLD uses its own direct wave-speed estimates; wave_speeds = 2 is not supported")
        @:PROHIBIT(riemann_solver == riemann_solver_hlld .and. low_Mach > 0, "low_Mach corrections are not implemented for HLLD")
        @:PROHIBIT(riemann_hypo_ADC .and. low_Mach > 0, &
                   & "riemann_hypo_ADC does not support low_Mach (the ADC HLL reference flux uses pre-correction velocities)")
        @:PROHIBIT(riemann_hypo_ADC .and. ADC_kappa <= 0._wp, "ADC_kappa must be positive")
        @:PROHIBIT(hypoelasticity .and. (cfl_const_dt .or. cfl_adap_dt), &
                   & "Automatic CFL time stepping uses the acoustic sound speed only and does not bound the elastic characteristic speed; set dt explicitly for hypoelastic runs")

    end subroutine s_check_inputs_hypo_branch

    !> Validates that each burning immersed-boundary patch injects a species index within the mechanism. inj_species indexes the
    !! image-point mass-fraction array Ys_IP(1:num_species) in m_ibm; an out-of-range value is an out-of-bounds write (silent
    !! corruption). Only reachable with chemistry.
    impure subroutine s_check_inputs_ib_injection

        integer :: i

        do i = 1, num_ibs
            @:PROHIBIT(patch_ib(i)%inj_species > num_species, &
                       & "patch_ib inj_species must be <= num_species (it indexes the image-point species mass fractions; an out-of-range value writes out of bounds)")
        end do

    end subroutine s_check_inputs_ib_injection

    !> Checks that each active particle cloud has a valid packing_method specified
    impure subroutine s_check_inputs_particle_clouds

        integer          :: i
        character(len=5) :: idxStr

        do i = 1, num_particle_clouds
            call s_int_to_str(i, idxStr)
            @:PROHIBIT(particle_cloud(i)%packing_method == dflt_int, &
                       & "particle_cloud("//trim(idxStr) &
                       & //")%packing_method must be specified (1 = rejection sampling, 2 = lattice)")
            @:PROHIBIT(particle_cloud(i)%packing_method /= 1 .and. particle_cloud(i)%packing_method /= 2, &
                       & "particle_cloud("//trim(idxStr) //")%packing_method must be 1 (rejection sampling) or 2 (lattice)")
        end do

    end subroutine s_check_inputs_particle_clouds

    !> Checks that each active synthetic-turbulence forcing zone has a fully specified position and a positive size in every active
    !! dimension
    impure subroutine s_check_inputs_synthetic_turbulence

        integer          :: i, d
        character(len=5) :: idxStr

        @:PROHIBIT(num_turbulent_sources <= 0, "num_turbulent_sources must be > 0 when synthetic_turbulence is enabled")

        do i = 1, num_turbulent_sources
            call s_int_to_str(i, idxStr)
            do d = 1, num_dims
                @:PROHIBIT(f_is_default(turb_pos(i, d)), &
                           & "turb_pos("//trim(idxStr) &
                           & //",:) must be specified for all num_dims when synthetic_turbulence is enabled")
                @:PROHIBIT(f_is_default(synth_L(i, d)) .or. synth_L(i, d) <= 0._wp, &
                           & "synth_L("//trim(idxStr)//",:) must be positive for all num_dims when synthetic_turbulence is enabled")
            end do
        end do

    end subroutine s_check_inputs_synthetic_turbulence

end module m_checker
