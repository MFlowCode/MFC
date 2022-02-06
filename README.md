# Multi-component Flow Code (MFC)

[![DOI](https://zenodo.org/badge/doi/10.1016/j.cpc.2020.107396.svg)](http://dx.doi.org/10.1016/j.cpc.2020.107396)
[![YourActionName Actions Status](https://github.com/ComputationalFlowPhysics/MFC-develop/workflows/CI/badge.svg)](https://github.com/ComputationalFlowPhysics/MFC-develop/actions)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![GitHub latest commit](https://badgen.net/github/last-commit/MFlowCode/MFC-develop)](https://github.com/MFlowCode/MFC-develop/commit/)


Welcome to MFC! 
The MFC is a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows.

# Authors

This is the documentation for the MFC (Multicomponent Flow Code).
The MFC is a simulation software for multi-component, multi-phase, and bubbly flows. 
MFC was first developed by the Colonius research group at Caltech.
Now it is developed and maintained by the groups of Professors <a href="https://colonius.caltech.edu/">Tim Colonius</a>, <a href="https://comp-physics.group">Spencer Bryngelson</a>, and <a href="https://vivo.brown.edu/display/mrodri97">Mauro Rodriguez</a>.
We try to maintain a list of current and past developers in the `AUTHORS` file!

# Documentation
 
  The following codes are documented, please follow the links to see their Doxygen:
* <a href="https://mflowcode.github.io/pre_process/namespaces.html">Pre_process</a> 
* <a href="https://mflowcode.github.io/simulation/namespaces.html">Simulation</a> 
* <a href="https://mflowcode.github.io/post_process/namespaces.html">Post_process</a>
 

## User's guide
 
  A user's guide is included 
  <a href="https://github.com/MFlowCode/MFC/raw/master/doc/MFC_user_guide.pdf">here.</a>
 
## MFC paper
 
  The paper that describes the MFC's capabilities:
* <a href="https://doi.org/10.1016/j.cpc.2020.107396">
        S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396
        </a>
  
## Related publications
 
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

## Installing MFC

To fetch, build, and run MFC and its dependencies on a UNIX-like system, you must have installed common utilities such as GNU's Make, Python3, its developement headers and libraries, a C/C++ compiler
(GCC, NVHPC, etc., but *not Clang*), and an MPI wrapper (like Open MPI). 
Below are some commands for popular operating systems and package managers.

### \*nix

* Via [Aptitude](https://wiki.debian.org/Aptitude)
```
sudo apt install tar wget make cmake gcc g++ python3 "openmpi-*" python python-dev python3-dev python3-venv libopenmpi-dev
```

### MacOS (including x86 and M1/Apple Silicon)

* Via [Homebrew](https://brew.sh/)

You can modify the assignment on the first line to have the GCC major version you wish to have installed and use.

```
USE_GCC_VERSION=11
brew install wget make python make cmake gcc@$USE_GCC_VERSION
HOMEBREW_CC=gcc-$USE_GCC_VERSION; HOMEBREW_CXX=g++-$USE_GCC_VERSION; brew install open-mpi
```

Further reading on `open-mpi` incompatibility with `clang`-based `gcc` on macOS: [here](https://stackoverflow.com/questions/27930481/how-to-build-openmpi-with-homebrew-and-gcc-4-9). We do *not* support `clang` due to conflicts with our Silo dependency.

### Fetch and build MFC

The following commands fetch and build MFC and its required dependencies. 
The dependencies are built to the `.mfc/` directory within your MFC installation. 
This should have no impact on your local installation(s) of these packages.

+ Fetch MFC

```
git clone --recursive https://github.com/MFlowCode/MFC
cd MFC
```

+ Build MFC and its dependencies with 8 threads in `release-cpu` mode:

| Argument     | Flag | Default Value | Possible values |
| ------------ | ---- | ------------- | --------------- |
| Threads      | `-j` | `1` | From `1` to the maximum logical thread count of your processor. |
| Target(s)    | `-t` | `MFC` | `MFC_PreProcess` / `MFC_Simulation` / `MFC_PostProcess` / `FFTW3` / `HDF5` / `SILO` |
| Configuration | `-cc` | `release-cpu` | `release-cpu` / `release-gpu` / `debug-cpu` / `debug-gpu` |

```
chmod +x ./mfc.sh
./mfc.sh -j 8 -t release-cpu
```

+ Run MFC's tests to make sure it was correctly built and your environment is adequate

```
./mfc.sh --test
```

More information about `mfc.sh` can be found in its help page (`./mfc.sh -h`).

## User Configuration (`mfc.user.yaml`)

The `mfc.sh` script used in the previous section is configured through the file named `mfc.user.yaml`.

# Running

The MFC can be run by changing into
a case directory and executing the appropriate Python input file.
Example Python input files can be found in the 
`example_cases` directories, and they are called `input.py`.
Their contents, and a guide to filling them out, are documented
in the user manual. A commented, tutorial script
can also be found in [example_cases/3d_sphbubcollapse/](example_cases/3D_sphbubcollapse/)
MFC can be executed as  

```
./input.py MFC_PreProcess
```

which will generate the restart and grid files that will be read 
by the simulation code. Then  

```
./input.py MFC_Simulation
```

will execute the flow solver. The last (optional) step
is to post treat the data files and output HDF5 databases
for the flow variables via  

```
./input.py MFC_PostProcess
```

# Development

## Fypp

MFC uses [Fypp](https://github.com/aradi/fypp), a Python-based Fortran preprocessor to reduce code duplication. `.fpp` files are converted into regular `.f90` files as part of the build process. Documentation for Fypp can be found [here](https://fypp.readthedocs.io/en/stable/). 

You can inspect the generated `.f90` files located in `.mfc/___current___/src/<name of target>/src`.

# License
 
Copyright 2022.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

# Acknowledgements
 
The development of the MFC  was supported in part by multiple current and past grants from the US Office of Naval Research (ONR), the US National Institute of Health (NIH), and the US National Science Foundation (NSF).
MFC computations utilize the Extreme Science and Engineering Discovery Environment (XSEDE), under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson) and ORNL Summit under allocation CFD154 (PI Bryngelson).
