# Containers

## Navigating Docker

Install Docker on 
* [MacOS](https://docs.docker.com/desktop/setup/install/mac-install/)
* [Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
* [Linux](https://docs.docker.com/desktop/setup/install/linux/)

### Docker Desktop GUI

- Search for the [sbryngelson/mfc](https://hub.docker.com/r/sbryngelson/mfc) repository. All the MFC containers are stored here.

- Find and pull a release tag (e.g., `latest-cpu`).

    * Read through the **Tag Details** below to distinguish between them. Docker Desktop's left sidebar has two key tabs: **Images** stores your program copies, and **Containers** shows instances of those images. You can launch multiple containers from a single image.

- Start a container by navigating to the `Images` tab and clicking the `Run` button.

    * Use the *Exec* section to interact with MFC directly via terminal, the *Files* section to transfer files between your device and container, and the *Stats* section to display resource usage.

### Docker CLI

You can navigate Docker entirely from the command line.
From a bash-like shell, pull from the [sbryngelson/mfc](https://hub.docker.com/r/sbryngelson/mfc) repository and run the latest MFC container:
```bash
docker run -it --rm --entrypoint bash sbryngelson/mfc:latest-cpu
```

**Selecting OS/ARCH:**  Docker selects the compatible architecture by default when pulling and running a container.
You can manually specify your platform if something seems to go wrong, as Docker may suggest doing so.
For example, `linux/amd64` handles many *nix-based x86 architectures, and `linux/arm64` handles Apple Silicon and Arm-based *nix devices.
You can specify it like this: 
```bash
docker run -it --rm --entrypoint bash --platform linux/amd64 sbryngelson/mfc:latest-cpu
```

**What's Next?** 

Once a container has started, the primary working directory is `/opt/MFC`, and all necessary files are located there.
You can check out the usual MFC documentation, such as the [Example Cases](https://mflowcode.github.io/documentation/md_examples.html), to get familiar with running cases.
Then, review the [Case Files](https://mflowcode.github.io/documentation/md_case.html) to write a custom case file.

## Details on Running Containers

Let's take a closer look at running MFC within a container.
Kick off a CPU container:
```bash
docker run -it --rm --entrypoint bash sbryngelson/mfc:latest-cpu
```
Or, start a GPU container.
```bash
docker run -it --rm --gpus all --entrypoint bash sbryngelson/mfc:latest-gpu
```

**Note:** `--gpus all` exposes the container to available GPUs, and _only NVIDIA GPUs are currently supported_.
[Ensure your device's CUDA version is at least 12.3](https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version) to avoid backward compatibility issues.

**Mounting Directory**

Mount a directory to `mnt` inside the container to easily transfer files between your host computer and the separate container.
For example, `cp -r <source> /mnt/destination>` moves something from your source computer to the container (reverse the order for the reverse to happen!).
```bash
docker run -it --rm --entrypoint bash -v "$PWD":/mnt sbryngelson/mfc:latest-cpu
```

**Shared Memory**

If you run a job with multiple MPI ranks, you could run into _MPI memory binding errors_.
This can manifest as a failed test (launched via `./mfc.sh test`) and running cases with `./mfc.sh run -n X <path/to/case.py>` where `X > 1`.
To avoid this issue, you can increase the shared memory size (to keep MPI working):
```bash
docker run -it --rm --entrypoint bash --shm-size=<e.g., 4gb> sbryngelson/mfc:latest-cpu
```
or avoid MPI altogether via `./mfc.sh <your commands> --no-mpi`.


### Portability

On the source machine, pull and save the image:
```bash
docker pull sbryngelson/mfc:latest-cpu
docker save -o mfc:latest-cpu.tar sbryngelson/mfc:latest-cpu
```
On the target machine, load and run the image:
```bash
docker load -i mfc:latest-cpu.tar
docker run -it --rm mfc:latest-cpu
```

## Using Supercomputers/Clusters via Apptainer/Singularity

### Interactive Shell

```bash
apptainer shell --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt  docker://sbryngelson/mfc:latest-gpu
Apptainer>cd /opt/MFC
```
or
```bash
apptainer exec --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt  docker://sbryngelson/mfc:latest-gpu bash -c "cd /opt/MFC && bash"
```
To run MFC on CPUs, omit `--nv` and use the `mfc:latest-cpu` container image.

### For Portability

On the source machine, pull and translate the image into `.sif` format:
```bash
apptainer build mfc:latest-gpu.sif docker://sbryngelson/mfc:latest-gpu
```
On the target machine, load and start an interactive shell:
```bash
apptainer shell --nv --fakeroot --writable-tmpfs --bind "$PWD":/mnt mfc:latest-gpu.sif
```

### Slurm Job

Below is an example Slurm batch job script.
Refer to your machine's user guide for instructions on properly loading and using Apptainer.
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

In the above,
* The `/sim` directory should have all the simulation files, including the case setup (`case.py`).
* The `--nv --fakeroot --writable-tmpfs` set of flags are needed to
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

### Example Tags

```shell
mfc:latest-xxx                   # Latest version (amd64 & arm64)
mfc:vx.x.x-cpu                   # CPU version    (amd64 & arm64)
mfc:vx.x.x-gpu                   # GPU version    (amd64 & arm64)
mfc:vx.x.x-xxx-ubuntu-xx.xx      # amd64 natively-supported version
mfc:vx.x.x-xxx-ubuntu-xx.xx-arm  # arm64 natively-supported version
```

### Architecture Support

You can specify your architecture with `--platform <os>/<arch>`, typically either `linux/amd64` or `linux/arm64`.
If you are unsure, Docker automatically selects the compatible image with your system architecture.
If native support isn't available, QEMU emulation is enabled for the following architectures, albeit with degraded performance.
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
