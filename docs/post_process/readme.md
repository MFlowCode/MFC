@mainpage MFC Post-Process

The post-process component reads raw simulation output and computes derived quantities for visualization. It produces silo/HDF5 files compatible with VisIt, ParaView, and other visualization tools.

## Modules

### Post-Process

| Module | Description |
|--------|-------------|
| @ref m_checker "m_checker" | Validates post-process input parameters and output format consistency |
| @ref m_data_input "m_data_input" | Reads raw simulation grid and conservative-variable data for a given time-step and fills buffer regions |
| @ref m_data_output "m_data_output" | Writes post-processed grid and flow-variable data to Silo-HDF5 or binary database files |
| @ref m_derived_variables "m_derived_variables" | Computes derived flow quantities (sound speed, vorticity, Schlieren, etc.) from conservative and primitive variables |
| @ref m_global_parameters "m_global_parameters" | Global parameters for the post-process: domain geometry, equation of state, and output database settings |
| @ref m_mpi_proxy "m_mpi_proxy" | MPI gather and scatter operations for distributing post-process grid and flow-variable data |
| @ref m_start_up "m_start_up" | Reads and validates user inputs, allocates variables, and configures MPI decomposition and I/O for post-processing |

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
- [Pre-Process API](../pre_process/index.html)
- [Simulation API](../simulation/index.html)
