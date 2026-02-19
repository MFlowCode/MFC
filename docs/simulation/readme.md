@mainpage MFC Simulation

The simulation component is the core flow solver. It advances the governing equations in time using high-order finite-volume methods on structured grids with GPU acceleration via OpenACC/OpenMP offloading.

## Modules

### Simulation

- @ref m_acoustic_src "m_acoustic_src"
- @ref m_body_forces "m_body_forces"
- @ref m_bubbles "m_bubbles"
- @ref m_bubbles_ee "m_bubbles_EE"
- @ref m_bubbles_el "m_bubbles_EL"
- @ref m_bubbles_el_kernels "m_bubbles_EL_kernels"
- @ref m_cbc "m_cbc"
- @ref m_checker "m_checker"
- @ref m_compute_cbc "m_compute_cbc"
- @ref m_compute_levelset "m_compute_levelset"
- @ref m_data_output "m_data_output"
- @ref m_derived_variables "m_derived_variables"
- @ref m_fftw "m_fftw"
- @ref m_global_parameters "m_global_parameters"
- @ref m_hyperelastic "m_hyperelastic"
- @ref m_hypoelastic "m_hypoelastic"
- @ref m_ib_patches "m_ib_patches"
- @ref m_ibm "m_ibm"
- @ref m_igr "m_igr"
- @ref m_mpi_proxy "m_mpi_proxy"
- @ref m_muscl "m_muscl"
- @ref m_pressure_relaxation "m_pressure_relaxation"
- @ref m_qbmm "m_qbmm"
- @ref m_rhs "m_rhs"
- @ref m_riemann_solvers "m_riemann_solvers"
- @ref m_sim_helpers "m_sim_helpers"
- @ref m_start_up "m_start_up"
- @ref m_surface_tension "m_surface_tension"
- @ref m_time_steppers "m_time_steppers"
- @ref m_viscous "m_viscous"
- @ref m_weno "m_weno"

### Common (shared)

- @ref m_boundary_common "m_boundary_common"
- @ref m_checker_common "m_checker_common"
- @ref m_chemistry "m_chemistry"
- @ref m_constants "m_constants"
- @ref m_derived_types "m_derived_types"
- @ref m_finite_differences "m_finite_differences"
- @ref m_helper "m_helper"
- @ref m_helper_basic "m_helper_basic"
- @ref m_model "m_model"
- @ref m_mpi_common "m_mpi_common"
- @ref m_phase_change "m_phase_change"
- @ref m_variables_conversion "m_variables_conversion"

## See Also

- [Home & User Guide](../documentation/index.html)
- [Pre-Process API](../pre_process/index.html)
- [Post-Process API](../post_process/index.html)
