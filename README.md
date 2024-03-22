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
  <a href="https://join.slack.com/t/mflowcode/shared_invite/zt-y75wibvk-g~zztjknjYkK1hFgCuJxVw">
    <img src="https://img.shields.io/badge/slack-MFC-purple.svg?logo=slack" />
  </a>
  <a href="https://lbesson.mit-license.org/">
    <img src="https://img.shields.io/badge/License-MIT-blue.svg" />
  </a>
</p>

Welcome to the home of MFC!
MFC simulates compressible multi-component and multi-phase flows, [amongst other things](#what-else-can-this-thing-do). 
It scales <b>ideally to exascale</b>; [tens of thousands of GPUs on NVIDIA- and AMD-GPU machines](#is-this-really-exascale), like Oak Ridge Summit and Frontier.
MFC is written in Fortran and makes use of metaprogramming to keep the code short (about 20K lines).
  
Get in touch with the maintainers, like <a href="mailto:shb@gatech.edu">Spencer</a>, if you have questions!
We have an [active Slack channel](https://join.slack.com/t/mflowcode/shared_invite/zt-y75wibvk-g~zztjknjYkK1hFgCuJxVw) and development team.
MFC has high-level documentation, visualizations, and more on [its website](https://mflowcode.github.io/).


## An example

We keep many examples.
Here's one!
MFC can execute high-fidelity simulations of shock-droplet interaction (see `examples/3d_shockdroplet`)

<p align="center">
    <img src="docs/res/shockdrop.png" alt="Shock Droplet Example" width="700"/>
</p>

Another example is the high-Mach flow over an airfoil, shown below.

<p align="center">
    <img src="docs/res/airfoil.png" alt="Airfoil Example" width="700"/><br/>
</p>


## Getting started

You can navigate [to this webpage](https://mflowcode.github.io/documentation/md_getting-started.html) to get started using MFC!
It's rather straightforward.
We'll give a brief intro. here for MacOS.
Using [brew](https://brew.sh), install MFC's modest set of dependencies:
```console
brew install wget python cmake gcc@13 mpich
```
You're now ready to build and test MFC!
Put it to a convenient directory via
```console
git clone https://github.com/mflowcode/MFC.git
cd MFC
```
and make sure MFC knows what compilers to use by putting the following in your `~/.profile`
```console
export CC=gcc-13
export CXX=g++-13
export FC=gfortran-13
```
and source that file, build, and test!
```console
source ~/.profile
./mfc.sh build -j 8
./mfc.sh test -j 8
```
And... you're done!

You can learn more about MFC's capabilities [via its documentation](https://mflowcode.github.io/documentation/index.html) or play with the examples located in the `examples/` directory (some are [shown here](https://mflowcode.github.io/documentation/md_examples.html))!

The shock-droplet interaction case above was run via
```console
./mfc.sh run ./examples/3d_shockdroplet/case.py -n 8
```
where `8` is the number of cores the example will run on.
You can visualize the output data, located in `examples/3d_shockdroplet/silo_hdf5`, via Paraview, Visit, or your other favorite software.

## Is this really exascale

[OLCF Frontier](https://www.olcf.ornl.gov/frontier/) is the first exascale supercomputer.
The weak scaling of MFC on this machine is below, showing near-ideal utilization. 

<p align="center">
    <img src="docs/res/scaling.png" alt="Scaling" width="400"/>
</p>


## What else can this thing do

MFC has many features.
They are organized below, just click the drop-downs!

<details>
<summary>Physics</summary>

* 1-3D
* Compressible
* Multi- and single-component
	* 4, 5, and 6 equation models for multi-component/phase features
* Multi- and single-phase 
	* Phase change via p, pT, and pTg schemes
* Grids
	* 1-3D Cartesian, cylindrical, axi-symmetric. 
	* Arbitrary grid stretching for multiple domain regions available.
	* Complex/arbitrary geometries via immersed boundary methods 
	* STL geometry files supported
* Sub-grid Euler-Euler multiphase models for bubble dynamics and similar
* Viscous effects (high-order accurate representations)
* Ideal and stiffened gas equations of state
* Acoustic wave generation (one- and two-way sound sources)
</details>

<details>
<summary>Numerics</summary>

* Shock and interface capturing schemes
	* First-order upwinding, WENO3 and 5. 
	* Reliable handling of high density ratios.
* Exact and approximate (e.g., HLL, HLLC) Riemann solvers
* Boundary conditions: Periodic, reflective, extrapolation/Neumann, slip/no-slip, non-reflecting characteristic buffers, inflows, outflows, and more.
* Runge-Kutta orders 1-3 (SSP TVD)
* Interface sharpening (THINC-like)
</details>

<details>
<summary>Large-scale and accelerated simulation</summary>

* GPU compatible on NVIDIA (P/V/A/H100, etc.) and AMD (MI200+) hardware
* Ideal weak scaling to 100% of leadership class machines
	* \>10K GPUs on [OLCF Summit](https://www.olcf.ornl.gov/summit/) (V100-based)
	* \>60K GPUs on world's first exascale computer, [OLCF Frontier](https://www.olcf.ornl.gov/frontier/) (MI250X-based)
* Near roofline behavior
</details>

<details>
<summary>Software robustness and other features</summary>

* [Fypp](https://fypp.readthedocs.io/en/stable/fypp.html) metaprogramming for code readability, performance, and portability
* Continuous Integration (CI)
	* Regression test cases on CPU and GPU hardware with each PR. Performed with GNU, Intel, and NVIDIA compilers.
	* Benchmarking to avoid performance regressions and identify speed-ups
* Continuous Deployment (CD) of [website](https://mflowcode.github.io) and [API documentation](https://mflowcode.github.io/documentation/index.html)
</details>


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
 
Copyright 2021-2024 Spencer Bryngelson and Tim Colonius.
MFC is under the MIT license (see [LICENSE](LICENSE) file for full text).

## Acknowledgements

Multiple federal sponsors have supported MFC development, including the US Department of Defense (DOD), National Institutes of Health (NIH), Department of Energy (DOE), and National Science Foundation (NSF).
MFC computations use OLCF Frontier, Summit, and Wombat under allocation CFD154 (PI Bryngelson) and ACCESS-CI under allocations TG-CTS120005 (PI Colonius) and TG-PHY210084 (PI Bryngelson).
