@mainpage API Documentation

MFC's source code is organized into three components that form a complete simulation pipeline. Each component has full module-level API documentation.

### [Pre-Process](../pre_process/index.html)

The pre-process component generates initial conditions and computational meshes for MFC simulations. It supports patch-based geometry construction, multi-component material initialization, and immersed boundary geometry.

### [Simulation](../simulation/index.html)

The simulation component is the core flow solver. It advances the governing equations in time using high-order finite-volume methods on structured grids with GPU acceleration via OpenACC/OpenMP offloading.

### [Post-Process](../post_process/index.html)

The post-process component reads raw simulation output and computes derived quantities for visualization. It produces silo/HDF5 files compatible with VisIt, ParaView, and other visualization tools.

All three components share a set of **common modules** for MPI communication, variable conversion, derived types, and utility functions. These are documented within each component's API reference.
