# Documentation
 
  The following codes are documented, please follow the links to see their Doxygen:
* <a href="https://mfc-caltech.github.io/pre_process/index.html">Pre_process</a> 
* <a href="https://mfc-caltech.github.io/simulation/index.html">Simulation</a> 
* <a href="https://mfc-caltech.github.io/post_process/index.html">Post_process</a>
    
# Authors
 
  This is the documentation for the MFC.
  The MFC is a simulation software for multi-component, multi-phase,
  and bubbly flows. MFC was developed at Caltech by a group
  of post-doctoral scientists and graduate research students
  under the supervision of Professor Tim Colonius. These contributors 
  include:
* Dr. Spencer Bryngelson
* Dr. Kevin Schmidmayer
* Dr. Vedran Coralic
* Dr. Jomela Meng
* Dr. Kazuki Maeda  

  and their contact information is located in the `AUTHORS` file in the source code.
 
# Installation
 
  The documents that describe how to configure and install the MFC are located in the 
  source code as `CONFIGURE` and `INSTALL`. They are also described here.
 
## Step 1: Configure and ensure dependencies can be located
 
 
### Main dependencies: MPI and Python 
  If you do not have Python, it can be installed via
  Homebrew on OSX (https://brew.sh/) as:  
`brew install python`
 
  or compiled via your favorite package manager on Unix systems.
 
  An MPI fortran compiler is required for all systems.
  If you do not have one, Homebrew can take care of this
  on OSX:  
`brew install open-mpi`  
 
  or compiled via another package manager on Unix systems.
 
### Simulation code dependency: FFTW 

If you already have FFTW compiled:
* Specify the location of your FFTW library and
      include files in Makefile.user (`fftw_lib_dir` and
      `fftw_include_dir`)  


If you do not have FFTW compiler, the library and
  installer are included in this package. Just:  
`cd installers`  
`./install_fftw.sh`  
 
### Post process code dependency: Silo/HDF5
 
  Post-processing of parallel data files is not required,
  but can indeed be handled with the MFC. For this, HDF5
  and Silo must be installed
 
  On OSX, a custom Homebrew tap for Silo is included in the installers
  directory. You can use it via  
`cd installers`  
`brew install silo.rb`  
 
  This will install silo and its dependences (including HDF5)
  in their usual locations (`/usr/local/lib` and
  `/usr/local/include`)
 
  On Unix systems, you can install via a package manager or
  from source. On CentOS (also Windows 7), HDF5
  binaries can be found at
      https://support.hdfgroup.org/ftp/HDF5/current18/bin/
  Untar this archive in your intended location via  
`tar -zxf [your HDF5 archive]`  
  
  Silo should be downloaded at
      https://wci.llnl.gov/simulation/computer-codes/silo/downloads
  then  
`tar -zxf [your Silo archive]`  
`cd [your Silo archive]`  
`./configure --prefix=[target installation directory] --enable-pythonmodule --enable-optimization --disable-hzip --disable-fpzip --enableportable-binary FC=mpif90 F77=mpif77 -with-hdf5=[your hdf5 directory]/include,/[your hdf5 directory]/lib --disable-silex`  
`make`  
`make install`  
 
  Add the following line to your `~/.bash_profile`:  
  `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/[your silo directory]/lib:/[your hdf5 directory]/lib`
 
  Finally:  
`source ~/.bash_profile`  
  
  You will then need to modify `silo_lib_dir` and `silo_include_dir` in
  `Makefile.user` to point to your silo directory.
 
## Step 2: Build and test
 
  Once all dependencies have been installed, the MFC can be built via
`make`
 
  from the MFC directory. This will build all MFC components. Individual
  components can be built via  
`make [component]`  
 
  where `[component]` is one of `pre_process`, `simulation`, or `post_process`.
 
  Once this is completed, you can ensure that the software is working
  as intended by  
`make test`  
 
 
# License
 
  MFC is free software: you can redistribute it and/or modify it under 
  the terms of the GNU General Public License, either version 3 
  of the License, or any later version. 
  A copy of the GNU General Public License is included with the MFC, and is
  also located at http://www.gnu.org/licenses.
 
  The MFC is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 
  Copyright 2019 California Institute of Technology
 
# Useful documents
  
## User's guide
 
  A user's guide is included here:
 
## MFC paper
 
  The paper that describes the MFC's capabilities is located here:
  
## Related publications
 
  Several publications have used the MFC in various stages of its 
  development. A partial list is included here.
 
  Refereed journal publications:
* <a href="https://arxiv.org/pdf/1903.08242.pdf">
        K. Schmidmayer, S. H. Bryngelson, T. Colonius (2019) under review, arXiv: 1903.08242
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
* <a href="http://colonius.caltech.edu/pdfs/CoralicColonius2013.pdf">
        V. Coralic and T. Colonius (2013) European Journal of Mechanics B-Fluids, Vol. 40, pp. 64-74 
        </a>
 
 
Ph.D. Disserations:
* <a href="https://thesis.library.caltech.edu/11007/">
        K. Maeda (2018) Ph.D. thesis, California Institute of Technology 
        </a>
* <a href="https://thesis.library.caltech.edu/9764/">
        J. Meng (2016) Ph.D. thesis, California Institute of Technology
        </a>
* <a href="https://thesis.library.caltech.edu/8758/">
        V. Coralic (2014) Ph.D. thesis, California Institute of Technology
        </a>

# Acknowledgements
 
The development of the MFC  was supported in part by multiple past grants from the US Office of 
Naval Research (ONR), the US National Institute of 
Health (NIH), and the US National Science Foundation (NSF), as well as current ONR grant numbers 
N0014-17-1-2676 and N0014-18-1-2625 and NIH grant number 2P01-DK043881.
The computations presented here utilized the Extreme Science
and Engineering Discovery Environment, which is supported under NSF
grant number CTS120005. K.M. acknowledges support from the Funai Foundation
for Information Technology via the Overseas Scholarship.
