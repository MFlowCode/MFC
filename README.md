<p align="center">
  <a href="http://mflowcode.github.io/">
    <img src="docs/res/readme.png" alt="MFC Banner" width="500"/>
  </a>
</p>

<p align="center">
  <a href="http://dx.doi.org/10.1016/j.cpc.2020.107396">
    <img src="https://zenodo.org/badge/doi/10.1016/j.cpc.2020.107396.svg" />
  </a>
  <a href="https://github.com/MFlowCode/MFC/actions">
    <img src="https://github.com/MFlowCode/MFC/actions/workflows/test.yml/badge.svg" />
  </a>
  <a href="https://lbesson.mit-license.org/">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" />
  </a>
</p>

<p align="justify">
  Welcome to the home of MFC!
  MFC simulates compressible multi-component and multi-phase flows, amongst other things. 
  It scales <b>ideally to exascale</b>; tens of thousands of GPUs on NVIDIA- and AMD-GPU Machines, like Oak Ridge Summit and Frontier.
  MFC is written in Fortran and makes use of metaprogramming to keep the code short (about 20K lines).
  Get in touch with the maintainers, like <a href="mailto:shb@gatech.edu">Spencer</a>, if you have questions!
  We have an active Slack channel and development team.

  MFC has high-level documentation, visualizations, and more on [its website](https://mflowcode.github.io/).
</p>

## Getting started

You can navigate [to this webpage](https://mflowcode.github.io/documentation/md_getting-started.html) to get started using MFC!
It's rather straightforward.
On MacOS, you can use [Homebrew](https://brew.sh) to install MFC's modest set of dependencies.
From a command line, issue
```console
brew install wget make python make cmake coreutils gcc openmpi
```
Now you're ready to build and test MFC!
Clone it to a convenient directory via
```console
git clone https://github.com/mflowcode/MFC.git
cd MFC
```
then build and test!
```console
./mfc.sh build -j 8
./mfc.sh test -j 8
```
And you're done!
You can learn more about MFC's capabilities [via its documentation](https://mflowcode.github.io/documentation/index.html) or play with the examples located in the `examples/` directory (some are [shown here](https://mflowcode.github.io/documentation/md_examples.html))!
Also, see below.

## Run an example

Examples are in the `examples/` directory.
For example, MFC can execute high-fidelity simulations of shock-droplet interactions.
Try
```console
./mfc.sh run ./examples/3d_shockdroplet/case.py -n 8
```
where `8` is the number of cores the example will run on.
You can visualize the output data, located in `examples/3d_shockdroplet/silo_hdf5`, via Paraview, Visit, or your other favorite software.
The result looks something like the below!

<p align="center">
    <img src="docs/res/shockdrop.png" alt="Shock Droplet Example" width="700"/>
</p>

## Citation

If you use MFC, consider citing it:

<p align="center">
  <a href="https://doi.org/10.1016/j.cpc.2020.107396">
    S. H. Bryngelson, K. Schmidmayer, V. Coralic, K. Maeda, J. Meng, T. Colonius (2021) Computer Physics Communications 4655, 107396
  </a>
</p>

```bibtex
@article{Bryngelson_2021,
  title = {{MFC: A}n open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver},
  author = {Spencer H. Bryngelson and Kevin Schmidmayer and Vedran Coralic and Jomela C. Meng and Kazuki Maeda and Tim Colonius},
  journal = {Computer Physics Communications},
  doi = {10.1016/j.cpc.2020.107396},
  year = {2021},
  pages = {107396},
}
```

## License
 
Copyright 2021-2024.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

## Acknowledgements
 
<p align="justify">
  MFC development was supported by multiple current and past grants from the US Department of Defense, National Institute of Health (NIH), Department of Energy (DOE), and the National Science Foundation (NSF).
  MFC computations use OLCF Frontier, Summit, and Wombat under allocation CFD154 (PI Bryngelson) and ACCESS-CI under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson).
</p>
