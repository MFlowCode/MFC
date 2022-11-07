# Welcome to MFC!

The MFC is a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows. The source code of MFC is available [on GitHub](https://github.com/MFlowCode/MFC) and is open to contributions!

## Code Documentation

- [Pre Process](../pre_process/)
- [Simulation](../simulation/)
- [Post Process](../post_process/)
 
## Publications
 
### Primary Paper
 
The paper that describes the MFC's capabilities:
* <a href="https://doi.org/10.1016/j.cpc.2020.107396">
    S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396
</a>

```
@article{Bryngelson_2021,
	title = {MFC: An open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver},
	author = {Spencer H. Bryngelson and Kevin Schmidmayer and Vedran Coralic and Jomela C. Meng and Kazuki Maeda and Tim Colonius},
	journal = {Computer Physics Communications},
	doi = {10.1016/j.cpc.2020.107396},
	year = 2021,
	month = {may},
	publisher = {Elsevier {BV}},
	pages = {107396},
}
```

### Related publications
 
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

### User Configuration

The [mfc.sh](https://github.com/MFlowCode/MFC/blob/master/mfc.sh) script used in the previous section is configured through the [defaults.yaml](https://github.com/MFlowCode/MFC/blob/master/defaults.yaml) file.

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

You can append `-DMFC_MEMORY_DUMP` to `release-gpu`'s Fortran compiler options in [defaults.yaml](defaults.yaml) to make the [simulation code](src/simulation/) call `acc_present_dump()` at various stages of program execution to obtain a printout of on-device memory usage. The [mem_parse.sh](misc/mem_parse.sh) script can be given as an argument the path to a file containing MFC's output, in order to aggregate the data and produce tables with formatted output.

## License
 
Copyright 2022.
MFC is under the MIT license (see [LICENSE](https://github.com/MFlowCode/MFC/blob/master/LICENSE) file for full text).

## Acknowledgements
 
The development of the MFC was supported in part by multiple current and past grants from the US Office of Naval Research (ONR), the US National Institute of Health (NIH), and the US National Science Foundation (NSF).
MFC computations utilize the Extreme Science and Engineering Discovery Environment (XSEDE), under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson) and ORNL Summit under allocation CFD154 (PI Bryngelson).
