# Building

## Build Environment

To fetch, build, and run MFC and its dependencies on a UNIX-like system, you
must have installed common utilities such as GNU's Make, Python3, its development
headers and libraries, a C/C++ compiler
(GCC, NVHPC, etc., but *not Clang*), and an MPI wrapper (like Open MPI). 
Below are some commands for popular operating systems and package managers.

[Anaconda](https://www.anaconda.com/) may interfere with the building process.
If an issue arises, you can either uninstall the affected Anaconda packages,
change the ordering of directory paths in your `$PATH`, or make aliases to the
correct binaries.

### *nix

- **Via [Aptitude](https://wiki.debian.org/Aptitude):**

```console
$ sudo apt update
$ sudo apt upgrade
$ sudo apt install tar wget make cmake gcc g++ \
                   python3 python3-dev         \
                   "openmpi-*" libopenmpi-dev
```

- **Via [Pacman](https://wiki.archlinux.org/title/pacman):**

```console
$ sudo pacman -Syu
$ sudo pacman -S base-devel coreutils  \
                 git ninja gcc-fortran \
                 cmake openmpi python3 \
                 python-pip openssh    \
                 python-virtualenv vim \
                 wget tree
```

If you wish to build MFC using [NVidia's NVHPC SDK](https://developer.nvidia.com/hpc-sdk), follow the instructions [here](https://developer.nvidia.com/nvidia-hpc-sdk-downloads).

### Windows

On Windows, you can either use Intel Compilers with the standard Microsoft toolchain, [Docker](https://docs.docker.com/get-docker/) or
the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/) for a Linux experience.

#### Windows + Intel (Native)

Install the latest version of:
- [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/)
- [Intel® oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
- [Intel® oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)

Then, in order to initialize your development environment, open a terminal window and run:
```console
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```

To follow this guide, please replace `./mfc.sh` with `mfc.bat` when running any
commands. `./mfc.sh` is intended Unix-like systems. You will also have access to the `.sln`
Microsoft Visual Studio solution files for an IDE (Integrated Development 
Environment).

#### Windows + Docker

See the instructions in the "Docker (Cross-Platform)" section of this document.

#### Windows + WSL

Install the latest version of the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/)
as well as a distribution such as Ubuntu which can be found [here](https://apps.microsoft.com/store/detail/ubuntu/9PDXGNCFSCZV). Acquiring an   interactive session is as simple as typing `wsl` in your
command prompt, or alternatively, selecting the distribution from the dropdown menu
available in the [Microsoft Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701).

You can now follow the appropriate instructions for your distribution.

### MacOS (x86 and Apple Silicon)

**Note:** macOS remains the most difficult platform to consistently compile MFC on.
If you run into issues, we suggest you try using Docker (instructions above).

  - **MacOS v10.15 (Catalina) or newer [ZSH]** (Verify with `echo $SHELL`)

```console
$ touch ~/.zshrc
$ open ~/.zshrc
```

  - **Older than MacOS v10.15 (Catalina) [BASH]** (Verify with `echo $SHELL`)
  
```console
$ touch ~/.bash_profile
$ open ~/.bash_profile
```
  
An editor should open. Please paste the following lines into it before saving the file. If you wish to use a version of GNU's GCC other than 11, modify the first assignment. These lines ensure that LLVM's Clang, and Apple's modified version of GCC, won't be used to compile MFC. Further reading on `open-mpi` incompatibility with `clang`-based `gcc` on macOS: [here](https://stackoverflow.com/questions/27930481/how-to-build-openmpi-with-homebrew-and-gcc-4-9). We do *not* support `clang` due to conflicts with our Silo dependency.

```console
# === MFC MPI Installation ===
export MFC_GCC_VER=11
export OMPI_MPICC=gcc-$MFC_GCC_VER
export OMPI_CXX=g++-$MFC_GCC_VER
export OMPI_FC=gfortran-$MFC_GCC_VER
export CC=gcc-$MFC_GCC_VER
export CXX=g++-$MFC_GCC_VER
export FC=gfortran-$MFC_GCC_VER
# === MFC MPI Installation ===
```

**Close the open editor and terminal window**. Open a **new terminal** window before executing the commands bellow.

```console
$ brew install wget make python make cmake coreutils gcc@$MFC_GCC_VER
$ HOMEBREW_MAKE_JOBS=$(nproc) brew install --cc=gcc-$MFC_GCC_VER --verbose --build-from-source open-mpi
```

They will download the dependencies MFC requires to build itself. `open-mpi` will be compiled from source, using the version of GCC we specified above with the environment variables `HOMEBREW_CC` and `HOMEBREW_CXX`. Building this package might take a while.

### Docker (Cross-Platform)

Docker is a lightweight, cross-platform, and performant alternative to Virtual Machines (VMs).
We build a Docker Image that contains the packages required to build and run MFC on your local machine.
  
First install Docker and Git:
- Windows: [Docker](https://docs.docker.com/get-docker/) + [Git](https://git-scm.com/downloads).
- macOS: `brew install git docker` (requires [Homebrew](https://brew.sh/)).
- Other systems:
```console
$ sudo apt install git docker # Debian / Ubuntu via Aptitude
$ sudo pacman -S git docker   # Arch Linux via Pacman
```

Once Docker and Git are installed on your system, clone MFC with

```console
$ git clone https://github.com/MFlowCode/MFC
$ cd MFC 
```

To fetch the prebuilt Docker image and enter an interactive bash session with the
recommended settings applied, run

```console
$ ./mfc.sh  docker # If on \*nix/macOS
  .\mfc.bat docker # If on Windows
```

We automatically mount and configure the proper permissions in order for you to
access your local copy of MFC, available at `~/MFC`. You will be logged-in as the
`me` user with root permissions.

:warning: The state of your container is entirely transient, except for the MFC mount.
Thus, any modification outside of `~/MFC` should be considered as permanently lost upon
session exit.

## Fetch & Build

MFC can be built without a helper script by directly running CMake but they
(`mfc.sh` on Linux and `mfc.bat` on Windows) offer convenience features that
assist in building, testing, optimization, as well as interactive and batch execution. To
build MFC without a helper script, consult the [CMakeLists.txt](CMakeLists.txt)
file for a full list of options, as well as [toolchain/dependencies/CMakeLists.txt](toolchain/dependencies/CMakeLists.txt)
for a CMake superbuild file that fetches and builds MFC's main dependencies with
supported versions.

+ **Fetch MFC:**

```console
$ git clone https://github.com/MFlowCode/MFC
$ cd MFC
```

+ **(Optional) Configure MFC defaults in [defaults.yaml](defaults.yaml):**

If you wish, you can override MFC's default build parameters in [defaults.yaml](defaults.yaml), a file intended for user customization. This can greatly reduce the number of command-line arguments you have to pass to [mfc.sh](mfc.sh)` in the following sections. You can do this at any time.

+ **Build MFC's codes in `release-cpu` mode with 8 threads:**

```console
$ ./mfc.sh build -t pre_process simulation post_process -j 8
```

To build MFC in different configurations (herein, *modes*), the `-m <mode>` option
can be specified to each call to `mfc.sh`. A full list of modes is located in
[defaults.yaml](defaults.yaml). It can be modified to work with system, and additional
modes can be created at your discretion. The default mode is `release-cpu` but
you can use others such as `release-gpu`.

**IMPORTANT NOTE**: This same mode will be used for any future commands such as `./mfc.sh test` and `./mfc.sh run` until you specify `-m` again (in any of these commands).

+ **Run MFC's tests in `release-cpu` mode with 8 threads:**

```console
$ ./mfc.sh test -j 8
```

Please refer to the [Testing](#testing-mfc) section of this document for more information. 
