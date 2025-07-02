# Temporally evolving turbulent mixing layer (3D)

Reference:
> Bell, J. H., & Mehta, R. D. (1990). Development of a two-stream mixing layer from tripped and untripped boundary layers. AIAA Journal, 28(12), 2034-2042.
> Rogers, M. M., & Moser, R. D. (1994). Direct simulation of a self‐similar turbulent mixing layer. Physics of Fluids, 6(2), 903-923.
> Pantano, C., & Sarkar, S. (2002). A study of compressibility effects in the high-speed turbulent shear layer using direct simulation. Journal of Fluid Mechanics, 451, 329-371.
> Vaghefi, S. N. S. (2014). Simulation and modeling of compressible turbulent mixing layer. State University of New York at Buffalo.
> Wang, X., Wang, J., & Chen, S. (2022). Compressibility effects on statistics and coherent structures of compressible turbulent mixing layers. Journal of Fluid Mechanics, 947, A38.

## Description
This directory (`turbulence_stat`) contains sub-directories and Matlab scripts for computing turbulence statistics of 3D temporally evolving turbulent mixing layer as described below:

Files:
- `set_user_inputs.m`: User input parameters are defined in this file.
- `run_turbulence.m`: Turbulence statistics (growth rate, Reynolds stress, TKE budget) are computed and plots are generated in this file.
- `average_tke_over_self_similar.m`: Averaging over self-similar period of TKE budget is performed in this file.
- `submit_batch_job_*.sh`: Example scripts for batch job submission on DoD Carpenter (PBS system) and NCSA Delta (Slurm system)

Directories:
- `log`: Log files from batch job are stored in this directory.
- `reference_data`: Data from reference papers are provided in this directory.
- `variables`: User inputs from `set_user_inputs.m` are saved as Matlab data (`user_inputs.mat`) in this directory.
- `results`: Outputs are stored in sub-directories in the directory.
