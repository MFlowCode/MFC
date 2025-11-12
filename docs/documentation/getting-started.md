# Getting Started

## Fetching MFC

You can either download MFC's [latest release from GitHub](https://github.com/MFlowCode/MFC/releases/latest) or clone the repository:

```shell
git clone https://github.com/MFlowCode/MFC.git
cd MFC
```

## Install via Homebrew (macOS)

On macOS, install prebuilt MFC via Homebrew:

```bash
brew install mflowcode/mfc/mfc
```

Run a quick example:

```bash
mkdir -p ~/mfc_quickstart && cd ~/mfc_quickstart
cp $(brew --prefix mfc)/examples/1D_sodshocktube/case.py .
# Use -n X to choose the number of MPI processes
mfc case.py -n 2
```

Notes:
- The Homebrew package uses a simplified syntax: just `mfc <case.py>` to run cases.
- Developer commands like `build`, `test`, `clean` are available when you clone the repo and use `./mfc.sh`.
- The package bundles a Python venv and prebuilt binaries; no additional setup is required.
- Examples are installed at `$(brew --prefix mfc)/examples/`.

## Build Environment

MFC can be built in multiple ways on various operating systems.
Please select your desired configuration from the list below:

  <summary><h2>*nix</h2></summary>

- **On supported clusters:** Load environment modules

```shell
source ./mfc.sh load
```

- **Via Aptitude:**

```shell
sudo apt update
sudo apt upgrade
sudo apt install tar wget make cmake gcc g++      \
                 python3 python3-dev python3-venv \
                 openmpi-bin libopenmpi-dev       \
                 libhdf5-dev libfftw3-dev         \
                 libblas-dev liblapack-dev
```

If you wish to build MFC using [NVIDIA's NVHPC SDK](https://developer.nvidia.com/hpc-sdk), first follow the instructions [here](https://developer.nvidia.com/nvidia-hpc-sdk-downloads).


  <summary><h2>Windows</h2></summary>

On Windows, you can either use Intel Compilers with the standard Microsoft toolchain,
or the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/)
for a Linux experience.

 <details>

   <summary><h3>Windows + WSL (Recommended)</h3></summary>

Install [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/) on Windows 11:
Either
1. Open a terminal with administrator privileges and run the following command:
```shell
wsl --install
```
Or
1. Open the Start menu, search for "Windows Features", and select "Turn Windows features on or off". Enable "Windows Subsystem for Linux" by checking the corresponding box.
2. Open the Microsoft Store, search for "Linux", and install your preferred distribution (e.g., [Ubuntu](https://apps.microsoft.com/store/detail/ubuntu/9PDXGNCFSCZV))

Useful software to install for using WSL on Windows:
- [Windows Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701)
- [Visual Studio Code](https://code.visualstudio.com/) and the [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) extension

Once you have WSL installed, you can follow the instructions for *nix systems above (for Ubuntu, see `Via Aptitude` section).

  </details>

  <details>

   <summary><h3>Native Windows (Intel)</h3></summary>

Install the latest version of:
- [Microsoft Visual Studio Community](https://visualstudio.microsoft.com/)
- Intel® oneAPI Base Toolkit
- Intel® oneAPI HPC Toolkit
- [Strawberry Perl](https://strawberryperl.com/) (Install and add `C:\strawberry\perl\bin\perl.exe` or your installation path to your [PATH](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/))
Please note that Visual Studio must be installed first, and the oneAPI Toolkits need to be configured with the installed Visual Studio, even if you plan to use a different IDE.

Then, to initialize your development environment, run the following command (or your installation path) in the command prompt:
```shell
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```
Alternatively, you can run the following command in Powershell:
```shell
cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
```
You could verify the initialization by typing `where mpiexec` in the command prompt terminal (does not work in Powershell), which should return the path to the Intel MPI executable.
To continue following this guide, please stay in the initialized terminal window. Replace `./mfc.sh` with `.\mfc.bat` for all commands.

If `.\mfc.bat build` produces errors, please run the command again. Repeating this process three times should resolve all errors (once each for pre_process, simulation, and post_process). If the same error persists after each attempt, please verify that you have installed all required software and properly initialized the development environment. If uncertain, you could try deleting the build directory and starting over.

You will also have access to the `.sln` Microsoft Visual Studio solution files for an IDE (Integrated Development Environment).

  </details>

  <summary><h3>MacOS</h3></summary>

Using [Homebrew](https://brew.sh/) you can install the necessary dependencies
before configuring your environment:

```shell
brew install coreutils python cmake fftw hdf5 gcc boost open-mpi lapack
echo -e "export BOOST_INCLUDE='$(brew --prefix --installed boost)/include'" | tee -a ~/.bash_profile ~/.zshrc
. ~/.bash_profile 2>/dev/null || . ~/.zshrc 2>/dev/null
! [ -z "${BOOST_INCLUDE+x}" ] && echo 'Environment is ready!' || echo 'Error: $BOOST_INCLUDE is unset. Please adjust the previous commands to fit with your environment.'
```

They will download the dependencies MFC requires to build itself.

</details>

## Building MFC

MFC can be built with support for various (compile-time) features:

| Feature            | Enable      | Disable        | Default | Description                                                     |
| :----------------: | :---------: | :------------: | :-----: | --------------------------------------------------------------- |
| **MPI**            | `--mpi`     | `--no-mpi`     | On      | Allows MFC to run on multiple processors (and nodes). |
| **GPU**            | `--gpu`     | `--no-gpu`     | Off     | Enables GPU acceleration via OpenACC.                           |
| **Debug**          | `--debug`   | `--no-debug`   | Off     | Requests the compiler build MFC in debug mode.                  |
| **GCov**           | `--gcov`    | `--no-gcov`    | Off     | Build MFC with coverage flags on.                              |
| **Unified Memory** | `--unified` | `--no-unified` | Off     | Build MFC with unified CPU/GPU memory (GH200 superchip only)  |
| **Single**         | `--single`  | `--no-single`  | Off     | Build MFC in single precision     

_⚠️ The `--gpu` option requires that your compiler supports OpenACC for Fortran for your target GPU architecture._

When these options are given to `mfc.sh`, they will be remembered when you issue future commands.
You can enable and disable features anytime by passing any of the arguments above.
For example, if you previously built MFC with MPI support and no longer wish to run using MPI, you can pass `--no-mpi` once, making the change permanent.

MFC comprises three codes, each being a separate _target_.
By default, all targets (`pre_process`, `simulation`, and `post_process`) are selected.
To only select a subset, use the `-t` (i.e., `--targets`) argument.
For a detailed list of options, arguments, and features, please refer to `./mfc.sh build --help`.

Most first-time users will want to build MFC using 8 threads (or more!) with MPI support:
```shell
./mfc.sh build -j 8
```

Some examples:
- Build MFC using 8 threads with MPI and GPU acceleration: `./mfc.sh build --gpu -j 8`.
- Build MFC using a single thread without MPI, GPU, and Debug support: `./mfc.sh build --no-mpi`.
- Build MFC's `simulation` code in Debug mode with MPI and GPU support: `./mfc.sh build --debug --gpu -t simulation`.

## Using Containers

Instead of building MFC from scratch, you can use containers to quickly access a pre-built version of MFC and its dependencies.
In brief, you can run the latest MFC container:
```bash
docker run -it --rm --entrypoint bash sbryngelson/mfc:latest-cpu
```
Please refer to the [Docker](https://mflowcode.github.io/documentation/docker.html) document for more information.

## Running the Test Suite

Run MFC's test suite with 8 threads:
```shell
./mfc.sh test -j 8
```

Please refer to the [Testing](testing.md) document for more information.

## Running an Example Case

MFC has example cases in the `examples` folder. You can run such a case interactively using 2 tasks by typing:

```shell
./mfc.sh run examples/2D_shockbubble/case.py -n 2
```

Please refer to the [Running](running.md) document for more information on `case.py` files and how to run them.
