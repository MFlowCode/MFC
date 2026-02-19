@mainpage MFC Pre-Process

The pre-process component generates initial conditions and computational meshes for MFC simulations. It supports patch-based geometry construction, multi-component material initialization, and immersed boundary geometry.

## Modules

### Pre-Process

| Module | Description |
|--------|-------------|
| @ref m_assign_variables "m_assign_variables" | Assigns initial primitive variables to computational cells based on patch geometry |
| @ref m_boundary_conditions "m_boundary_conditions" | Applies spatially varying boundary condition patches along domain edges and faces |
| @ref m_check_ib_patches "m_check_ib_patches" | Validates geometry parameters and constraints for immersed boundary patches |
| @ref m_check_patches "m_check_patches" | Validates geometry parameters and constraints for initial condition patches |
| @ref m_checker "m_checker" | Checks pre-process input file parameters for compatibility and correctness |
| @ref m_data_output "m_data_output" | Writes grid and initial condition data to serial or parallel output files |
| @ref m_global_parameters "m_global_parameters" | Defines global parameters for the computational domain, simulation algorithm, and initial conditions |
| @ref m_grid "m_grid" | Generates uniform or stretched rectilinear grids with hyperbolic-tangent spacing |
| @ref m_icpp_patches "m_icpp_patches" | Constructs initial condition patch geometries (lines, circles, rectangles, spheres, etc.) on the grid |
| @ref m_initial_condition "m_initial_condition" | Assembles initial conditions by layering prioritized patches via constructive solid geometry |
| @ref m_mpi_proxy "m_mpi_proxy" | Broadcasts user inputs and decomposes the domain across MPI ranks for pre-processing |
| @ref m_perturbation "m_perturbation" | Perturbs initial mean flow fields with random noise, mixing-layer instabilities, or simplex noise |
| @ref m_simplex_noise "m_simplex_noise" | 2D and 3D simplex noise generation for procedural initial condition perturbations |
| @ref m_start_up "m_start_up" | Reads and validates user inputs, loads existing grid/IC data, and initializes pre-process modules |

### Common (shared)

| Module | Description |
|--------|-------------|
| @ref m_boundary_common "m_boundary_common" | Noncharacteristic and processor boundary condition application for ghost cells and buffer regions |
| @ref m_checker_common "m_checker_common" | Shared input validation checks for grid dimensions and AMD GPU compiler limits |
| @ref m_chemistry "m_chemistry" | Multi-species chemistry interface for thermodynamic properties, reaction rates, and transport coefficients |
| @ref m_compile_specific "m_compile_specific" | Platform-specific file and directory operations: create, delete, inquire, getcwd, and basename |
| @ref m_constants "m_constants" | Compile-time constant parameters: default values, tolerances, and physical constants |
| @ref m_delay_file_access "m_delay_file_access" | Rank-staggered file access delays to prevent I/O contention on parallel file systems |
| @ref m_derived_types "m_derived_types" | Shared derived types for field data, patch geometry, bubble dynamics, and MPI I/O structures |
| @ref m_finite_differences "m_finite_differences" | Finite difference operators for computing divergence of velocity fields |
| @ref m_helper "m_helper" | Utility routines for bubble model setup, coordinate transforms, array sampling, and special functions |
| @ref m_helper_basic "m_helper_basic" | Basic floating-point utilities: approximate equality, default detection, and coordinate bounds |
| @ref m_model "m_model" | Binary STL file reader and processor for immersed boundary geometry |
| @ref m_mpi_common "m_mpi_common" | MPI communication layer: domain decomposition, halo exchange, reductions, and parallel I/O setup |
| @ref m_nvtx "m_nvtx" | NVIDIA NVTX profiling API bindings for GPU performance instrumentation |
| @ref m_phase_change "m_phase_change" | Phase transition relaxation solvers for liquid-vapor flows with cavitation and boiling |
| @ref m_precision_select "m_precision_select" | Working-precision kind selection (half/single/double) and corresponding MPI datatype parameters |
| @ref m_variables_conversion "m_variables_conversion" | Conservative-to-primitive variable conversion, mixture property evaluation, and pressure computation |

## See Also

- [Home & User Guide](../documentation/index.html)
- [Simulation API](../simulation/index.html)
- [Post-Process API](../post_process/index.html)
