# Docker

## Navigating Docker Desktop/CLI

- Install Docker on [Mac](https://docs.docker.com/desktop/setup/install/mac-install/), [Windows](https://docs.docker.com/desktop/setup/install/windows-install/), or [Linux](https://docs.docker.com/desktop/setup/install/linux/).

Via Docker Desktop GUI,
- Search for [sbryngelson/mfc](https://hub.docker.com/r/sbryngelson/mfc) repository where all the MFC images are stored then pull a release tag (e.g., `latest-cpu`).

    Read through the **Tag Details** below to distinguish between them. Docker Desktop's left sidebar has two key tabs: **Images** stores your program copies, while **Containers** shows instances of those images. You can launch multiple containers from a single image.

- Start a container by clicking the Run button in the Images tab.

    Use the *Exec* section to interact with MFC directly via terminal, the *Files* section to transfer files between your device and container, and the *Stats* section to display resource usage.

Or via Docker CLI,

- Pull from [sbryngelson/mfc](https://hub.docker.com/r/sbryngelson/mfc) repository and run the latest MFC container.

```bash
docker run -it --rm --entrypoint bash sbryngelson/mfc:latest-cpu
```
<br>

**Selecting OS/ARCH:** 

Docker by default selects the compatible architecture when pulling and running a container. However, you can manually specify your platform (i.e., `linux/amd64` for most devices or `linux/arm64` for Apple Silicon).
```bash
docker run -it --rm --entrypoint bash --platform linux/amd64 sbryngelson/mfc:latest-cpu
```
<br>

**What is Next:** 

After starting a container, the primary working directory is `/opt/MFC`, where all necessary files are located. Read through [Example Cases](https://mflowcode.github.io/documentation/md_examples.html) to get familiar with running cases, then review [Case Files](https://mflowcode.github.io/documentation/md_case.html) to write a custom case file.


<br>

## Running Containers

Start a CPU container.
```bash
docker run -it --rm --entrypoint bash sbryngelson/mfc:latest-cpu
```
Start a GPU container.
```bash
docker run -it --rm --gpus all --entrypoint bash sbryngelson/mfc:latest-gpu
```
**Note:** `--gpus all` exposes the container to available GPUs, and only NVIDIA GPUs are currently supported. Make sure your device's CUDA version is at least 12.3 to avoid backward compatibility issues.

<br>
<br>

**Mounting Directory:** 

Mount a directory to `mnt` inside the container to easily transfer files between the host and the container, e.g. `cp -r <source> /mnt/destination>`.
```bash
docker run -it --rm --entrypoint bash -v "$PWD":/mnt sbryngelson/mfc:latest-cpu
```
<br>

**Shared Memory:** 

When encountering MPI memory binding errors resulting in failed tests and cases, increase the shared memory size or disable MPI inside the container (`./mfc.sh --no-mpi`).
```bash
docker run -it --rm --entrypoint bash --shm-size=<e.g., 4gb> sbryngelson/mfc:latest-cpu
```




### **For Portability,**

On the source machine,
- Pull and save the image.
```bash
docker pull sbryngelson/mfc:latest-cpu
docker save -o mfc:latest-cpu.tar sbryngelson/mfc:latest-cpu
```
On the target machine,
- Load and run the image.
```bash
docker load -i mfc:latest-cpu.tar
docker run -it --rm mfc:latest-cpu
```

<br>

## HPC Cluster Usage (Apptainer/Singularity)

### **Interactive Shell**
```bash
apptainer shell --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt  docker://sbryngelson/mfc:latest-gpu
Apptainer>cd /opt/MFC
```
or
```bash
apptainer exec --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt  docker://sbryngelson/mfc:latest-gpu bash -c "cd /opt/MFC && bash"
```
**Note:** To run MFC on CPU cores, omit `--nv` and use `mfc:latest-cpu` container image.

### **For Portability,**
On the source machine,
- Pull and translate the image into `.sif` format.
```bash
apptainer build mfc:latest-gpu.sif docker://sbryngelson/mfc:latest-gpu
```
On the target machine,
- Load and start an interactive shell.
```bash
apptainer shell --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt mfc:latest-gpu.sif
```



### Slurm Job
Below is a slurm job template. Refer to your HPC user guide for instructions on properly loading and using apptainer.
```bash
#!/bin/bash
#SBATCH --job-name=mfc-sim
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --time=06:00:00
#SBATCH --partition=
#SBATCH --output=mfc-sim-%j.out
#SBATCH --error=mfc-sim-%j.err

cd $SLURM_SUBMIT_DIR

# Define container image
CONTAINER="mfc:latest-gpu.sif"

apptainer exec --nv --fakeroot --writable-tmpfs \
--bind "$PWD":/mnt \
  $CONTAINER \
  bash -c "cd /opt/MFC && ./mfc.sh run sim/case.py -- -c <computer>"
```
Where,

`/sim` directory should all the simulation files, including the case setup (`case.py`).
 
`--nv --fakeroot --writable-tmpfs`, these flags are critical to:
- Grant access to the host system's NVIDIA GPUs and its CUDA libraries.
- Enable root-like permissions inside the container without actual root access.
- Allow temporary write access to the container filesystem.



## Tag Details
### Base Images
- CPU images (v4.3.0-latest releases) are built on **Ubuntu 22.04**.
- GPU images (v4.3.0-latest releases) are built on **NVHPC SDK 23.11 (CUDA 12.3) & Ubuntu 22.04**.

### Tag Structure
- **`vx.x.x`** - Official MFC release versions (recommended: use `latest` release)
- **`cpu/gpu`** - Build configurations for CPU or GPU acceleration.
- **`ubuntu-xx.xx`** - Base Ubuntu version (standard = `amd64`, `-arm` = `arm64`)

### Available Tags
```
mfc:latest-xxx                   # Latest version (amd64 & arm64)
mfc:vx.x.x-cpu                   # CPU version    (amd64 & arm64)
mfc:vx.x.x-gpu                   # GPU version    (amd64 & arm64)
mfc:vx.x.x-xxx-ubuntu-xx.xx      # amd64 natively-supported version
mfc:vx.x.x-xxx-ubuntu-xx.xx-arm  # arm64 natively-supported version
```
### **Architecture Support**
You can specify the desired architecture with `--platform <os>/<arch>` - either `linux/amd64` or `linux/arm64`. If unsure, Docker automatically selects the compatible image with your system architecture. If native support isn't available, QEMU emulation is enabled for the following architectures albeit with degraded performance.
```
linux/amd64
linux/amd64/v2
linux/amd64/v3
linux/arm64
linux/riscv64
linux/ppc64le
linux/s390x
linux/386
linux/mips64le
linux/mips64
linux/loong64
linux/arm/v7
linux/arm/v6
```


<br>
<br>
