# Performance Results

MFC has been benchmarked on several CPUs and GPU devices.
This page shows a summary of these results.

## Expected time-steps/hour

The following table outlines expected performance in terms of the number of time steps per hour, rounded to the nearest hundred (higher is better).
A 3D inviscid, 6-equation problem is solved for various problem sizes (grid cells) and hardware.
A 3rd order (3-stage) Runge-Kutta time-stepper is used.
CPU results utilize an entire processor die.

| Hardware             | # Cores | Steps/Hr (1M pts)      | Steps/Hr (4M pts)      | Steps/Hr  (8M pts)   | Compiler    | Computer      |
| ---:                 | :----:  |    :----:      |  :---:         | :---:        | :----:      | :---          |
| NVIDIA V100          | 1 (device)      | 88.5k          | 18.7k          | N/A          | NVHPC 22.11 | PACE Phoenix  |
| NVIDIA V100          | 1 (device)      | 78.8k          | 18.8k          | N/A          | NVHPC 22.11 | OLCF Summit   |
| NVIDIA A100          | 1 (device)      | 114.4k         | 34.6k          | 16.5k        | NVHPC 23.5  | Wingtip       |
| AMD MI250X           | 1 (GCD)      | 77.5k          | 22.3k          | 11.2k        | CCE 16.0.1  | OLCF Frontier |
| Intel Xeon Gold 6226 | 12 (cores)     | 2.5k           | 0.7k           | 0.4k         | GNU 10.3.0  | PACE Phoenix  |
| Apple Silicon M2     | 6 (cores)       | 2.8k           | 0.6k           | 0.2k         | GNU 13.2.0  | N/A           |

We also show the expected performance of MFC for the same problem as above, except for the 5-equation model used, in the table below.
It is presented in the same manner as the one above.

| Hardware             | # Cores | Steps/Hr (1M pts)      | Steps/Hr (4M pts)      | Steps/Hr  (8M pts)   | Compiler    | Computer      |
| ---:                 | :----:  |    :----:      |  :---:         | :---:        | :----:      | :---          |
| NVIDIA V100          | 1 (device)       | 113.4k         | 26.2k          | 13.0k        | NVHPC 22.11 | PACE Phoenix  |
| NVIDIA V100          | 1 (device)      | 107.7k         | 26.3k          | 13.1k        | NVHPC 22.11 | OLCF Summit   |
| NVIDIA A100          | 1 (device)      | 153.5k         | 48.0k          | 22.5k        | NVHPC 23.5  | Wingtip       |
| AMD MI250X           | 1 (GCD)      | 104.2k         | 31.0k          | 14.8k        | CCE 16.0.1  | OLCF Frontier |
| Intel Xeon Gold 6226 | 12 (cores)     | 5.4k           | 1.6k           | 0.8k         | GNU 10.3.0  | PACE Phoenix  |
| Apple Silicon M2     | 6 (cores)      | 3.7k           | 11.0k          | 0.3k         | GNU 13.2.0  | N/A           |

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
