## The Main Idea behind the implemented Out-of-Core Strategy for Grace-Hopper

To run MFC out-of-core on Grace-Hopper using Unified Memory we implement a zero-copy strategy.

We start by setting preferred location CPU for all buffers by hooking into the allocate macro and setting `NVIDIA_ALLOC_MODE=2`.
This way we disable access counter based migrations and keep everything on the CPU memory, freeing up as much GPU memory as possible.

Then, for the "important" buffers that are frequently accessed from the GPU, we reset preferred location to GPU in order to place them (and directly populate them) in GPU memory.
This is done by the `PREFER_GPU` macro that has been manually placed in the code right after the allocations of the "important" buffers.
To activate these hints we export `NVIDIA_MANUAL_GPU_HINTS=1`.

To allow fine grained control and be able to simulate larger sizes, we also use the following environment variables:
- With `NVIDIA_IGR_TEMPS_ON_GPU` we control how many temporaries from the IGR module are to be placed in GPU memory.
- With `NVIDIA_VARS_ON_GPU` we control how many of the `q_cons_ts(1)%vf(j)%sf` arrays we place in GPU memory.

It is important to note that we have rearranged the timestep updates in the 3rd order TVD Runge Kutta scheme in a way that allows us to pass only `q_cons_ts(1)` to the `compute_rhs` routines.
This way, in order to keep the computation of `compute_rhs` (mostly) on GPU data, we only need to store `q_cons_ts(1)` (fully or even partially) in GPU memory.
Thus, we choose to keep `q_cons_ts(2)` in CPU memory for the full lifetime of the simulation, freeing up space in GPU memory that allows for bumping up the size of the simulation, without sacrificing performance.
In the timestep updates between the `compute_rhs` calls, we access both `q_cons_ts(1)` and `q_cons_ts(2)` directly from the physical location where they reside (zero-copy), simultaneously pulling data from GPU memory and CPU memory (through C2C), making good use of Grace-Hopper.

Note: This rearrangement most likely "breaks" the timestepper for different physics cases, but we can easily fix it in a later step.

## Example Workflow for Out-of-Core Strategy based on Unified Memory

```shell
# Allocate a node
salloc -A g183 --partition normal -t 02:00:00 -N 1 -n 4 --cpus-per-task=71

# Start uenv
uenv start --view=modules icon/25.2:v1

# cd to root directory of MFC
cd MFC-Wilfong

# Load modules
. ./mfc.sh load -c san -m g

# Build
export MFC_CUDA_CC=90
./mfc.sh build --gpu -j 71 --single --unified --verbose

# Run pre_process and simulation binaries with case optimization (in an interactive job)
./mfc.sh run examples/3D_IGR_perf_test/case.py --case-optimization -t pre_process simulation --gpu -N 1 -n 4 -j 71 -c santis

# Run pre_process and simulation binaries with case optimization (in an batch job)
./mfc.sh run examples/3D_IGR_perf_test/case.py --case-optimization -t pre_process simulation --gpu -N 1 -n 4 -j 71 -c santis -e batch -p normal -a g183 -w 00:15:00
```
The environment variables `NVIDIA_ALLOC_MODE`, `NVIDIA_MANUAL_GPU_HINTS`, `NVIDIA_VARS_ON_GPU`, and `NVIDIA_IGR_TEMPS_ON_GPU`, can be set appropriately in `toolchain/templates/santis.mako`, to configure a run with ALL buffers either in GPU or in CPU memory, or a run with SOME buffers in GPU memory and the rest in CPU memory.
