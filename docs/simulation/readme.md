@mainpage MFC Simulation

The simulation component is the core flow solver. It advances the governing equations in time using high-order finite-volume methods on structured grids with GPU acceleration via OpenACC/OpenMP offloading.

## Modules

### Simulation

| Module | Description |
|--------|-------------|
| @ref m_acoustic_src "m_acoustic_src" | Applies acoustic pressure source terms including focused, planar, and broadband transducers |
| @ref m_body_forces "m_body_forces" | Computes gravitational and user-defined body force source terms for the momentum equations |
| @ref m_bubbles "m_bubbles" | Shared bubble-dynamics procedures (radial acceleration, wall pressure, sound speed) for ensemble- and volume-averaged models |
| @ref m_bubbles_ee "m_bubbles_EE" | Computes ensemble-averaged (Euler--Euler) bubble source terms for radius, velocity, pressure, and mass transfer |
| @ref m_bubbles_el "m_bubbles_EL" | Tracks Lagrangian bubbles and couples their dynamics to the Eulerian flow via volume averaging |
| @ref m_bubbles_el_kernels "m_bubbles_EL_kernels" | Kernel functions (Gaussian, delta) that smear Lagrangian bubble effects onto the Eulerian grid |
| @ref m_cbc "m_cbc" | Characteristic boundary conditions (CBC) for slip walls, non-reflecting subsonic inflow/outflow, and supersonic boundaries |
| @ref m_checker "m_checker" | Validates simulation input parameters for consistency and supported configurations |
| @ref m_compute_cbc "m_compute_cbc" | Characteristic boundary condition (CBC) computations for subsonic inflow, outflow, and slip walls |
| @ref m_compute_levelset "m_compute_levelset" | Computes signed-distance level-set fields and surface normals for immersed-boundary patch geometries |
| @ref m_data_output "m_data_output" | Writes solution data, run-time stability diagnostics (ICFL, VCFL, CCFL, Rc), and probe/center-of-mass files |
| @ref m_derived_variables "m_derived_variables" | Derives diagnostic flow quantities (vorticity, speed of sound, numerical Schlieren, etc.) from conservative and primitive variables |
| @ref m_fftw "m_fftw" | Forward and inverse FFT wrappers (FFTW/cuFFT/hipFFT) for azimuthal Fourier filtering in cylindrical geometries |
| @ref m_global_parameters "m_global_parameters" | Global parameters for the computational domain, fluid properties, and simulation algorithm configuration |
| @ref m_hyperelastic "m_hyperelastic" | Computes the left Cauchy--Green deformation tensor and hyperelastic stress source terms |
| @ref m_hypoelastic "m_hypoelastic" | Computes hypoelastic stress-rate source terms and damage-state evolution |
| @ref m_ib_patches "m_ib_patches" | Immersed boundary patch geometry constructors for 2D and 3D shapes |
| @ref m_ibm "m_ibm" | Ghost-node immersed boundary method: locates ghost/image points, computes interpolation coefficients, and corrects the flow state |
| @ref m_igr "m_igr" | Iterative ghost rasterization (IGR) for sharp immersed boundary treatment |
| @ref m_mpi_proxy "m_mpi_proxy" | MPI halo exchange, domain decomposition, and buffer packing/unpacking for the simulation solver |
| @ref m_muscl "m_muscl" | MUSCL reconstruction with interface sharpening for contact-preserving advection |
| @ref m_pressure_relaxation "m_pressure_relaxation" | Pressure relaxation for the six-equation multi-component model via Newton--Raphson equilibration and volume-fraction correction |
| @ref m_qbmm "m_qbmm" | Quadrature-based moment methods (QBMM) for polydisperse bubble moment inversion and transport |
| @ref m_rhs "m_rhs" | Assembles the right-hand side of the governing equations using finite-volume flux differencing, Riemann solvers, and physical source terms |
| @ref m_riemann_solvers "m_riemann_solvers" | Approximate and exact Riemann solvers (HLL, HLLC, HLLD, exact) for the multicomponent Navier--Stokes equations |
| @ref m_sim_helpers "m_sim_helpers" | Simulation helper routines for enthalpy computation, CFL calculation, and stability checks |
| @ref m_start_up "m_start_up" | Reads input files, loads initial conditions and grid data, and orchestrates solver initialization and finalization |
| @ref m_surface_tension "m_surface_tension" | Computes capillary source fluxes and color-function gradients for the diffuse-interface surface tension model |
| @ref m_time_steppers "m_time_steppers" | Total-variation-diminishing (TVD) Runge--Kutta time integrators (1st-, 2nd-, and 3rd-order SSP) |
| @ref m_viscous "m_viscous" | Computes viscous stress tensors and diffusive flux contributions for the Navier--Stokes equations |
| @ref m_weno "m_weno" | WENO/WENO-Z/TENO reconstruction with optional monotonicity-preserving bounds and mapped weights |

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
- [Post-Process API](../post_process/index.html)
