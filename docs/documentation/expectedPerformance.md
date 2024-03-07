# Performance Results

MFC has been benchmarked on several CPUs and GPU devices.
This page shows a summary of these results.

## Expected time-steps/hour

The following table outlines observed performance as nanoseconds per grid point (ns/GP) per right-hand side evaluation (lower is better).
A 3D, inviscid, 5-equation problem with two advected species is solved for various problem sizes (grid cells per CPU die or GPU device) and hardware.

| Hardware             | # Cores | 1M GPs      | 4M GPs      | 8M GPs | Compiler    | Computer      |
| ---:                 | :----:  |    :----:      |  :---:         | :---:        | :----:      | :---          |
| NVIDIA V100          | 1 (device)       | 95.6         | 104          | 104        | NVHPC 22.11 | PACE Phoenix  |
| NVIDIA V100          | 1 (device)      | 101         |  104         | 104        | NVHPC 22.11 | OLCF Summit   |
| NVIDIA A100          | 1 (device)      | 71         | 56          | 59        | NVHPC 23.5  | Wingtip       |
| AMD MI250X           | 1 (GCD)      | 108          | 90       | 96      | CCE 16.0.1  | OLCF Frontier |
| Intel Xeon Gold 6226 | 12 (cores)     | 1963           | 1688           | 1686         | GNU 10.3.0  | PACE Phoenix  |
| Apple M2     | 6 (cores)      | 2919           | 245          | 4500        | GNU 13.2.0  | N/A           |

__All results are in nanoseconds per grid point per RHS evaluation.__

## Weak scaling

Weak scaling results are obtained by increasing the problem size with the number of processes so that work per process remains constant.

### AMD MI250X GPU

MFC weask scales to (at least) 65,536 AMD MI250X GPUs on OLCF Frontier with 96% efficiency.
This corresponds to 87% of the entire machine.

<img src="../res/weakScaling/frontier.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

### NVIDIA V100 GPU

MFC weak scales to (at least) 13,824 V100 NVIDIA V100 GPUs on OLCF Summit with 97% efficiency.
This corresponds to 50% of the entire machine.

<img src="../res/weakScaling/summit.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

### IBM Power9 CPU
MFC Weak scales to 13,824 Power9 CPU cores on OLCF Summit to within 1% of ideal scaling.

<img src="../res/weakScaling/cpuScaling.svg" style="height: 50%; width:50%; border-radius: 10pt"/>

## Strong scaling

Strong scaling results are obtained by keeping the problem size constant and increasing the number of processes so that work per process decreases.

### NVIDIA V100 GPU

For these tests, the base case utilizes 8 GPUs with one MPI process per GPU.
The performance is analyzed at two different problem sizes of 16M and 64M grid points, with the base case using 2M and 8M grid points per process.

#### 16M Grid Points

<img src="../res/strongScaling/strongScaling16.svg" style="width: 50%; border-radius: 10pt"/>

#### 64M Grid Points
<img src="../res/strongScaling/strongScaling64.svg" style="width: 50%; border-radius: 10pt"/>

### IBM Power9 CPU

CPU strong scaling tests are done with problem sizes of 16, 32, and 64M grid points, with the base case using 2, 4, and 8M cells per process.

<img src="../res/strongScaling/cpuStrongScaling.svg" style="width: 50%; border-radius: 10pt"/>
