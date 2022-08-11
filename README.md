# <image src="doc/MFC.png" />

<p align="center">
  <a href="http://dx.doi.org/10.1016/j.cpc.2020.107396">
    <img src="https://zenodo.org/badge/doi/10.1016/j.cpc.2020.107396.svg" />
  </a>
  <a href="https://github.com/ComputationalFlowPhysics/MFC-develop/actions">
    <img src="https://github.com/ComputationalFlowPhysics/MFC-develop/workflows/CI/badge.svg" />
  </a>
  <a href="https://lbesson.mit-license.org/">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" />
  </a>
  <a href="https://github.com/MFlowCode/MFC-develop/commit/">
    <img src="https://badgen.net/github/last-commit/MFlowCode/MFC-develop" />
  </a>
  <a href="https://hub.docker.com/repository/docker/henryleberre/mfc">
    <img src="https://shields.io/docker/image-size/henryleberre/mfc" />
  </a>
</p>

Welcome to MFC! 
The MFC is a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows.

<p align="center">
 <a href="#authors">Authors</a> | 
 <a href="#publications">Publications</a> | 
 <a href="#build-environment">Build Environment</a> | 
 <a href="#running-mfc">Running</a> | 
 <a href="#testing-mfc">Testing</a> | 
 <a href="#development">Development</a> | 
 <a href="#useful-scripts">Useful Scripts</a> | 
 <a href="https://github.com/MFlowCode/MFC/raw/master/doc/MFC_user_guide.pdf">User's Guide</a> | 
 <a href="https://mflowcode.github.io/">Documentation</a>
</p>

## Authors

<p align="justify">
  This is the documentation for the MFC (Multicomponent Flow Code).
  The MFC is a simulation software for multi-component, multi-phase, and bubbly flows. 
  MFC was first developed by the Colonius research group at Caltech.
  Now it is developed and maintained by the groups of Professors <a href="https://comp-physics.group">Spencer Bryngelson</a>, <a href="https://colonius.caltech.edu/">Tim Colonius</a>, and <a href="https://vivo.brown.edu/display/mrodri97">Mauro Rodriguez</a> (alphabetical).
  We try to maintain a list of current and past developers in the <a href="AUTHORS">AUTHORS</a> file!
</p>
 
## Publications
 
### Primary Paper
 
  The paper that describes the MFC's capabilities:
* <a href="https://doi.org/10.1016/j.cpc.2020.107396">
        S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396
        </a>

<details>
  <summary><h3 style="margin-vertical: 0;">Related publications</h3></summary>
 
  Several publications have used the MFC in various stages of its 
  development. A partial list is included here.
 
  Journal papers:
  * <a href="https://arxiv.org/abs/2112.14172">
      S. H. Bryngelson, R. O. Fox, T. Colonius (2021) arXiv: 2112.14172.
    </a>
  * <a href="https://asa.scitation.org/doi/full/10.1121/10.0000746">
      S. H. Bryngelson and T. Colonius (2020) Journal of the Acoustical Society of America, Vol. 147, pp. 1126-1135
    </a>
  * <a href="https://www.sciencedirect.com/science/article/pii/S0021999119307855">
      K. Schmidmayer, S. H. Bryngelson, T. Colonius (2020) Journal of Computational Physics, Vol. 402, 109080
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/BryngelsonSchmidmayerColonius2019.pdf">
      S. H. Bryngelson, K. Schmidmayer, T. Colonius (2019) International Journal of Multiphase Flow, Vol. 115, pp. 137-143  
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/MaedaColonius2019.pdf">
      K. Maeda and T. Colonius (2019) Journal of Fluid Mechanics, Vol. 862, pp. 1105-1134 
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/MaedaColonius2018c.pdf">
      K. Maeda and T. Colonius (2018) Journal of Computational Physics, Vol. 371, pp. 994-1017 
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/MengColonius2018.pdf">
      J. C. Meng and T. Colonius (2018) Journal of Fluid Mechanics,  Vol. 835, pp. 1108-1135 
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/MaedaColonius2017.pdf">
      K. Maeda and T. Colonius (2017) Wave Motion, Vol. 75, pp. 36-49 
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/MengColonius2015.pdf">
      J. C. Meng and T. Colonius (2015) Shock Waves, Vol. 25(4), pp. 399-414 
    </a>
  * <a href="http://colonius.caltech.edu/pdfs/CoralicColonius2014.pdf">
      V. Coralic and T. Colonius (2014) Journal of Computational Physics, Vol. 274, pp. 95-121 
    </a>
   
  Ph.D. Disserations:
  * <a href="https://thesis.library.caltech.edu/11395/">
      J.-C. Veilleux (2019) Ph.D. thesis, California Institute of Technology 
    </a>
  * <a href="https://thesis.library.caltech.edu/11007/">
      K. Maeda (2018) Ph.D. thesis, California Institute of Technology 
    </a>
  * <a href="https://thesis.library.caltech.edu/9764/">
      J. Meng (2016) Ph.D. thesis, California Institute of Technology
    </a>
  * <a href="https://thesis.library.caltech.edu/8758/">
      V. Coralic (2014) Ph.D. thesis, California Institute of Technology
    </a>

</details>
 
## Build Environment

<p align="justify">
  To fetch, build, and run MFC and its dependencies on a UNIX-like system, you
  must have installed common utilities such as GNU's Make, Python3, its development
  headers and libraries, a C/C++ compiler
  (GCC, NVHPC, etc., but *not Clang*), and an MPI wrapper (like Open MPI). 
  Below are some commands for popular operating systems and package managers.
<p>

[Anaconda](https://www.anaconda.com/) may interfere with the building process.
If an issue arises, you can either uninstall the affected Anaconda packages,
change the ordering of directory paths in your `$PATH`, or make aliases to the
correct binaries.

<details>
  <summary><h3>*nix</h3></summary>

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
  
</details>

<details>
  <summary><h3>Windows</h3></summary>

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
  
</details>

<details>
  <summary><h3>MacOS (x86 and Apple Silicon)</h3></summary>
 
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

</details>

<details>
  <summary><h3>Docker (Cross-Platform)</h3></summary>

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
  $ .\mfc.bat docker # If on Windows
  ```
  
  We automatically mount and configure the proper permissions in order for you to
  access your local copy of MFC, available at `~/MFC`. You will be logged-in as the
  `me` user with root permissions.
  
  :warning: The state of your container is entirely transient, except for the MFC mount.
  Thus, any modification outside of `~/MFC` should be considered as permanently lost upon
  session exit.

</details>
 
## Fetch, Configure, Build, and Test MFC

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

### User Configuration (`defaults.yaml`)

The `mfc.sh` script used in the previous section is configured through the file named `defaults.yaml`.

## Running MFC

MFC can be run using `mfc.sh`'s `run` command. It supports both interactive and
batch execution, the latter being designed for multi-socket systems, namely supercomputers,
equipped with a scheduler such as PBS, SLURM, and LSF. A full (and updated) list
of available arguments can be acquired with `./mfc.sh run -h`. Example Python input
files can be found in the [samples/](samples/) directory. They print a Python dictionary containing input parameters for MFC. Their contents, and a guide to filling them out, are documented
in the user manual. A commented, tutorial script
can also be found in [samples/3d_sphbubcollapse/](samples/3D_sphbubcollapse/case.py).

The skeleton for an input file may look like the following:

```python
#!/usr/bin/env python3

import json

# Configuring case dictionary
print(json.dumps({
  # Insert case parameters here
  ...
}))
```

Thus, you can run your case file with Python to view the computed case dictionary
that will be processed by MFC when you run. This is particularly useful when
computations are done in Python to generate the case.

<details>
  <summary><h3>Interactive Execution</h3></summary>

  To run all stages of MFC, that is [pre_process](src/pre_process_code/), [simulation](src/simulation_code/), and [post_process](src/post_process_code/)   on the sample case [2D_shockbubble](samples/2D_shockbubble/),
  
  ```console
  $ ./mfc.sh run samples/2D_shockbubble/case.py
  ```
  
  If you want to run a subset of the available stages, you can use the `-t` argument.
  To use multiple threads, use the `-n` option along with the number of threads you wish to use.
  If a (re)build is required, it will be done automatically, with the number of threads
  specified with the `-j` option.
  
  For example,
  
  - Running [pre_process](src/pre_process_code/) with 2 cores:
  
  ```console
  $ ./mfc.sh run samples/2D_shockbubble/case.py -t pre_process -n 2
  ```
  
  - Running [simulation](src/simulation_code/) and [post_process](src/post_process_code/)
  using 4 cores:
  
  ```console
  $ ./mfc.sh run samples/2D_shockbubble/case.py -t simulation post_process -n 4
  ```
  
  Most parameters have sensible defaults which can be overridden in [defaults.yaml](defaults.yaml):
  
  https://github.com/MFlowCode/MFC-develop/blob/d74e714b08562a9f8f815112e05df54c99c8c18f/defaults.yaml#L12-L21
  
  On some computer clusters, MFC might select the wrong MPI program to execute your application
  because it uses a general heuristic for its selection. Notably, `srun` is known to fail on some SLURM
  systems when using GPUs or MPI implementations from different vendors, whereas `mpirun` functions properly. To override and manually specify which
  MPI program you wish to run your application with, please use the `-b <program name>` option (i.e `--binary`).
  
  Additional flags can be appended to the MPI executable call using the `-f` (i.e `--flags`) option.
  
  Please refer to `./mfc.sh run -h` for a complete list of arguments and options, along with their defaults.

</details>

<details>
  <summary><h3>Batch Execution</h3></summary>

  The MFC detects which scheduler your system is using and handles the creation and
  execution of batch scripts. The batch engine is requested with the `-e batch` option.
  Whereas the interactive engine can execute all of MFC's codes in succession, the batch engine
  requires you to only specify one target with the `-t` option. The number of nodes and GPUs can, 
  respectively be specified with the `-N` (i.e `--nodes`) and `-g` (i.e `--gpus-per-node`) options.
  
  ```console
  $ ./mfc.sh run samples/2D_shockbubble/case.py -e batch -N 2 -n 4 -g 4 -t simulation
  ```
  
  Other useful arguments include:
  
  - `-# <job name>` to name your job. (i.e `--name`)
  - `-@ sample@example.com` to receive emails from the scheduler. (i.e `--email`)
  - `-w hh:mm:ss` to specify the job's maximum allowed walltime. (i.e `--walltime`)
  - `-a <account name>` to identify the account to be charged for the job. (i.e `--account`)
  - `-p <partition name>` to select the job's partition. (i.e `--partition`)
  
  Since some schedulers don't have a standardized syntax to request certain resources, MFC can only
  provide support for a restricted subset of common configuration options. If MFC fails
  to execute on your system, or if you wish to adjust how the program runs and resources
  are requested to be allocated, you are invited to modify the template batch script for your queue system.
  Upon execution of `./mfc.sh run`, MFC fills in the template with runtime parameters, to
  generate the batch file it will submit. These files are located in the [templates](templates/)
  directory. To request GPUs, modification of the template will be required on most systems.
  
  - Lines that begin with `#>` are ignored and won't figure in the final batch
  script, not even as a comment.
  
  - Statements of the form `${expression}` are string-replaced to provide
  runtime parameters, most notably execution options. They reference the variables in the
  same format as those under the "run" section of [defaults.yaml](defaults.yaml), replacing
  `-` for `_`. You can perform therein any Python operation recognized by the built-in `expr()` function.
  
  As an example, one might request GPUs on a SLURM system using the following:
  
  ```
  #SBATCH --gpus=v100-32:{gpus_per_node*nodes}
  ```
  
  - Statements of the form `{MFC::expression}` tell MFC where to place the common code,
  across all batch files, that is required for proper execution. They are not intended to be
  modified by users.
  
  **Disclaimer**: IBM's JSRUN on LSF-managed computers does not use the traditional node-based approach to
  allocate resources. Therefore, the MFC constructs equivalent resource-sets in task and GPU count.

</details>

### Example Runs

- Oak Ridge National Laboratory's [Summit](https://www.olcf.ornl.gov/summit/):

```console
$ ./mfc.sh run samples/2D_shockbubble/case.py -e batch    \
               -N 2 -n 4 -g 4 -t simulation -a <redacted>
```

- University of California, San Diego's [Expanse](https://www.sdsc.edu/services/hpc/expanse/):

```console
$ ./mfc.sh run samples/2D_shockbubble/case.py -e batch -p GPU -t simulation \
               -N 2 -n 8 -g 8 -f="--gpus=v100-32:16" -b mpirun –w 00:30:00
```

## Testing MFC
 
To run MFC's test suite, run
```console
$ ./mfc.sh test -j <thread count>
```

It will generate and run test cases, comparing their output to that of previous runs from versions of MFC considered to be accurate. *golden files*, stored in the `tests/` directory contain this data, by aggregating `.dat` files generated when running MFC. A test is considered passing when our error tolerances are met, in order to maintain a high level of stability and accuracy. Run `./mfc.sh test -h` for a full list of accepted arguments.

Most notably, you can consult the full list of tests by running
```
$ ./mfc.sh test -l
```

To restrict to a given range, use the `--from` (`-f`) and `--to` (`-t`) options. To run a 
(non-contiguous) subset of tests, use the `--only` (`-o`) option instead.

<details>
  <summary><h3>Creating Tests</h3></summary>

  To (re)generate *golden files*, append the `-g` (i.e `--generate`) option:
  ```console
  $ ./mfc.sh test -g -j 8
  ```
  
  Adding a new test case can be done by modifying [cases.py](toolchain/mfc/tests/cases.py). The function `generate_cases` is responsible for generating the list of test cases. Loops and conditionals are used to vary parameters, whose defaults can be found in the `BASE_CFG` case object within [case.py](toolchain/mfc/tests/case.py). The function operates on two variables:
  
  - `stack`: A stack that holds the variations to the default case parameters. By pushing and popping the stack inside loops and conditionals, it is easier to nest test case descriptions, as it holds the variations that are common to all future test cases within the same indentation level (in most scenarios).
  
  - `cases`: A list that holds fully-formed `Case` objects, that will be returned at the end of the function. 
  
  Internally a test case is described as:
  ```python
  @dataclasses.dataclass(init=False)
  class Case:
      trace:  str
      params: dict
      ppn:    int
  ```
  
  where:
  - The `trace` is a string that contains a human-readable description of what parameters were varied, or more generally what the case is meant to test. **Each `trace` must be distinct.**
  - `params` is the fully resolved case dictionary, as would appear in a Python case input file.
  - `ppn` is the number of processes per node to use when running the case.
  
  To illustrate, consider the following excerpt from `generate_cases`:
  
  ```python
  for weno_order in [3, 5]:
    stack.push(f"weno_order={weno_order}", {'weno_order': weno_order})

    for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
        stack.push(f"(mapped_weno={mapped_weno},mp_weno={mp_weno})", {
            'mapped_weno': mapped_weno,
            'mp_weno':     mp_weno
        })

        if not (mp_weno == 'T' and weno_order != 5):
            cases.append(create_case(stack, '', {}))

        stack.pop()

    stack.pop()
  ```
  
  When pushing to the stack, or creating a new case with the `create_case` function, you must specify:
  - `stack`: The current stack.
  - `trace`: A human-readable string describing what you are currently varying.
  - `variations`: A Python dictionary with case parameter variations.
  - (Optional) `ppn`: The number of processes per node to use (default is 1).

  If a trace is empty (that is, the empty string `""`), it will not appear in the final trace, but any case parameter variations associated with it will still be applied.

  Finally, the case is appended to the `cases` list, which will be returned by the `generate_cases` function.
  
</details>

## Development

### Fypp

MFC uses [Fypp](https://github.com/aradi/fypp), a Python-based Fortran preprocessor to reduce code duplication. `.fpp` files are converted into regular `.f90` files as part of the build process. Documentation for Fypp can be found [here](https://fypp.readthedocs.io/en/stable/). 

You can inspect the generated `.f90` files in `src/<code name>/autogen/`.

## Useful Scripts

### Loading Modules

On computer systems that run using environment modules, with implementations like [TACC's Lmod](https://github.com/TACC/Lmod), MFC provides a script that can load modules that have been used by contributors.

```console
$ . ./mfc.sh load
``` 

The list of modules offered by a system is subject to change. The aforementioned script serves as a convenient way to load modules that should work for most users of MFC. 

### OpenACC Memory Profiling

You can append `-DMFC_MEMORY_DUMP` to `release-gpu`'s Fortran compiler options in [defaults.yaml](defaults.yaml) to make the [simulation code](src/simulation_code/) call `acc_present_dump()` at various stages of program execution to obtain a printout of on-device memory usage. The [mem_parse.sh](misc/mem_parse.sh) script can be given as an argument the path to a file containing MFC's output, in order to aggregate the data and produce tables with formatted output.

## License
 
Copyright 2022.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

## Acknowledgements
 
<p align="justify">
  The development of the MFC was supported in part by multiple current and past grants from the US Office of Naval Research (ONR), the US National Institute of Health (NIH), and the US National Science Foundation (NSF).
  MFC computations utilize the Extreme Science and Engineering Discovery Environment (XSEDE), under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson) and ORNL Summit under allocation CFD154 (PI Bryngelson).
</p>
