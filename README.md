# Multi-component Flow Code (MFC)

[![DOI](https://zenodo.org/badge/doi/10.1016/j.cpc.2020.107396.svg)](http://dx.doi.org/10.1016/j.cpc.2020.107396)
[![YourActionName Actions Status](https://github.com/ComputationalFlowPhysics/MFC-develop/workflows/CI/badge.svg)](https://github.com/ComputationalFlowPhysics/MFC-develop/actions)

Welcome to the MFC! 
The MFC is a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows.

# Authors

This is the documentation for the MFC (Multicomponent Flow Code).
The MFC is a simulation software for multi-component, multi-phase, and bubbly flows. 
MFC was first developed by the Colonius research group at Caltech.
Now it is developed and maintained by the groups of Professors <a href="https://colonius.caltech.edu/">Tim Colonius</a>, <a href="https://comp-physics.group">Spencer Bryngelson</a>, and <a href="https://vivo.brown.edu/display/mrodri97">Mauro Rodriguez</a>.

# Documentation
 
  The following codes are documented, please follow the links to see their Doxygen:
* <a href="https://mflowcode.github.io/pre_process/namespaces.html">Pre_process</a> 
* <a href="https://mflowcode.github.io/simulation/namespaces.html">Simulation</a> 
* <a href="https://mflowcode.github.io/post_process/namespaces.html">Post_process</a>
 

## User's guide
 
  A user's guide is included 
  <a href="https://github.com/ComputationalFlowPhysics/MFC/raw/master/doc/MFC_user_guide.pdf">here.</a>
 
## MFC paper
 
  The paper that describes the MFC's capabilities:
* <a href="https://doi.org/10.1016/j.cpc.2020.107396">
        S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2020) Computer Physics Communications 4655, 107396
        </a>
  
## Related publications
 
  Several publications have used the MFC in various stages of its 
  development. A partial list is included here.
 
  Refereed journal publications:
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



# Simple Installation

To get MFC running as fast as possible without having to configure the
dependencies yourself, you can follow the following steps on most UNIX-like
systems. This method is best suited for development and Continous Integration
(CI) workflows.

To fetch, build, and run MFC and its dependencies on a UNIX-like system, you
must have installed common utilities (Tar, Wget, GNU Make, ...), Python2,
Python3, Python's developement headers and libraries, a C/C++ compiler
((GCC and G++) or Clang), an MPI Fortran compiler
(like [Open MPI](https://www.open-mpi.org/)). Here are some commands for popular
Operating Systems and package managers, 

Package Manager                                            | Suggested Command 
---                                                        | ---
[Aptitude](https://wiki.debian.org/Aptitude) (Debian-like) | `sudo apt install tar wget make gcc g++ python3 openmpi-* python python-dev python3-dev`
[Homebrew](https://brew.sh/) (macOS)                       | `brew install wget make python open-mpi make gcc`

The following commands fetch and build the required dependencies to to the
`dependencies/` folder within your MFC installation. This should have no impact on your
local installations of these packages.

```
cd dependencies
sh ./install.sh -j <num_threads>
cd ..
```

You can now run `make -j <num_threads>` on MFC to build it, and check your
installation by running tests with `make tests`.

```
make -j <num_threads>
make tests
```

# Advanced Installation

MFC is split into 3 components with (mostly) layered dependencies. The following
table can be used to resolve the minimal amount of packages you are required to
install and configure depending on your specific use case.

Codes             | Description                             | Required Packages
----------------- | --------------------------------------- | ---
Pre process Code  | "The pre-processor generates initial conditions and spatial grids from the physical patches specified in the Python input file andexports them as binary files to be read by the simulator." [1] | <ul><li>UNIX Utilities</li><li>GNU Make</li><li>C/C++ Compiler</li><li>Fortran MPI Compiler</li><li>[Python](https://www.python.org/)</li>
Simulation Code   | "The simulator, given the initial-condition files generated bythe pre-processor, solves the corresponding governing flowequations with the specified boundary conditions using ourinterface-capturingnumericalmethod." [1] | <ul><li>[FFTW](https://www.fftw.org/)</li><li>[Lapack](http://www.netlib.org/lapack/)</li></ul>
Post process Code | "The post-processor reads simulation data and exports HDF5/Silo databases that include variables and derived variables, asspecifiedintheinputfile." [1] | <ul><li>[HDF5](https://www.hdfgroup.org/solutions/hdf5/)</li><li>[Silo](https://wci.llnl.gov/simulation/computer-codes/silo)*</li></ul>

<sub>[1] <a href="https://doi.org/10.1016/j.cpc.2020.107396">
        S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2020) Computer Physics Communications 4655, 107396
</a></sub>

<sup>*Please note that <a href="https://wci.llnl.gov/simulation/computer-codes/silo">Silo</a> depends on <a href="https://www.hdfgroup.org/solutions/hdf5/">HDF5</a>.

Once you have decided on which dependencies you want to install, you can:

- Build them from source
- Obtain them through your package manager
- Modify and run the [install.sh](dependencies/install.sh) script above

You can now modify the [Makefile.user](Makefile.user) file at the root of MFC's 
source tree to specify the paths to those packages' binaries and header files, 
as well as to supply additional command line arguments.

Variable           | Description
------------------ | ---
`FC`               | Binary to execute to compile and preprocess Fortran files.
`FFLAGS`           | Command-line arguments to supply to your Fortran compiler of choice.
`lapack_lib_dir`   | Path to the [Lapack](http://www.netlib.org/lapack/)'s install location.
`fftw_lib_dir`     | Path to the parent folder containing [FFTW](https://www.fftw.org/)'s binaries (/lib).
`fftw_include_dir` | Path to the parent folder containing [FFTW](https://www.fftw.org/)'s header files (/include).
`silo_lib_dir`     | Path to the parent folder containing [Silo](https://wci.llnl.gov/simulation/computer-codes/silo)'s binaries (/lib).
`silo_include_dir` | Path to the parent folder containing [Silo](https://wci.llnl.gov/simulation/computer-codes/silo)'s header files (/include).

You can now run `make -j <num_threads>` on MFC to build it, and check your installation by running tests with `make tests`.

```
make -j <num_threads>
make tests
```

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
python pre_process
```

which will generate the restart and grid files that will be read 
by the simulation code. Then  

```
python simulation
```

will execute the flow solver. The last (optional) step
is to post treat the data files and output HDF5 databases
for the flow variables via  

```
python post_process
```

# License
 
Copyright 2021.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

# Acknowledgements
 
The development of the MFC  was supported in part by multiple current and past grants from the US Office of Naval Research (ONR), the US National Institute of Health (NIH), and the US National Science Foundation (NSF).
MFC computations utilize the Extreme Science and Engineering Discovery Environment (XSEDE), under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson) and ORNL Summit under allocation CFD154 (PI Bryngelson).
