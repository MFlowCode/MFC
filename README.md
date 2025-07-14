<p align="center">
  <a href="http://mflowcode.github.io/">
    <img src="docs/res/readme.png" alt="MFC Banner" width="500"/>
  </a>
</p>

<p align="center">
  <a href="http://dx.doi.org/10.1016/j.cpc.2020.107396" target="_blank">
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
  <a href="https://codecov.io/github/MFlowCode/MFC" target="_blank">
    <img src="https://codecov.io/github/MFlowCode/MFC/graph/badge.svg?token=8SY043QND4">
  </a>
  <a href="https://github.com/MFlowCode/MFC/blob/master/.github/CONTRIBUTING.md" target="_blank">
    <img src="https://img.shields.io/badge/Contributing-Guide-orange?style=flat">
  </a>
</p>

Welcome to the home of MFC!
MFC simulates compressible multi-component and multi-phase flows, [amongst other things](#what-else-can-this-thing-do). 
MFC is written in Fortran and uses metaprogramming to keep the code short (about 20K lines).

MFC is used on the latest leadership-class supercomputers.
It scales <b>ideally to exascale</b>; [tens of thousands of GPUs on NVIDIA- and AMD-GPU machines](#is-this-really-exascale) on Oak Ridge Summit and Frontier.
MFC is a SPEChpc benchmark candidate, part of the JSC JUPITER Early Access Program, and used OLCF Frontier and LLNL El Capitan early access systems.
  
Get in touch with <a href="mailto:shb@gatech.edu">Spencer</a> if you have questions!
We have an [active Slack channel](https://join.slack.com/t/mflowcode/shared_invite/zt-y75wibvk-g~zztjknjYkK1hFgCuJxVw) and development team.
MFC has high- and low-level documentation, visualizations, and more on [its website](https://mflowcode.github.io/).

## An example

We keep many examples.
Here are some of them!
MFC can execute high-fidelity simulations of shock-droplet interaction (see `examples/3d_shockdroplet`)

<p align="center">
    <img src="docs/res/shockdrop.png" alt="Shock Droplet Example" width="700"/>
</p>

This one simulates high-Mach flow over an airfoil:

<p align="center">
    <img src="docs/res/airfoil.png" alt="Airfoil Example" width="700"/><br/>
</p>

And here is a high amplitude acoustic wave reflecting and emerging through a circular orifice:

<p align="center">
    <img src="docs/res/orifice.png" alt="Orifice Example" width="700"/><br/>
</p>


## Getting started

You can navigate [to this webpage](https://mflowcode.github.io/documentation/md_getting-started.html) to get started using MFC!
It's rather straightforward.
We'll give a brief intro. here for MacOS.
Using [brew](https://brew.sh), install MFC's dependencies:
```shell
brew install coreutils python cmake fftw hdf5 gcc boost open-mpi
```
You're now ready to build and test MFC!
Put it to a convenient directory via
```shell
git clone https://github.com/MFlowCode/MFC
cd MFC
```
and be sure MFC knows where to find Boost by appending to your dotfiles and sourcing them again
```shell
echo -e "export BOOST_INCLUDE='$(brew --prefix --installed boost)/include'" | tee -a ~/.bash_profile ~/.zshrc
. ~/.bash_profile 2>/dev/null || . ~/.zshrc 2>/dev/null
! [ -z "${BOOST_INCLUDE+x}" ] && echo 'Environment is ready!' || echo 'Error: $BOOST_INCLUDE is unset. Please adjust the previous commands to fit with your environment.'
```
then you can build MFC and run the test suite!
```shell
./mfc.sh build -j $(nproc)
./mfc.sh test -j $(nproc)
```
And... you're done!

You can learn more about MFC's capabilities [via its documentation](https://mflowcode.github.io/documentation/index.html) or play with the examples located in the `examples/` directory (some are [shown here](https://mflowcode.github.io/documentation/md_examples.html))!

The shock-droplet interaction case above was run via
```shell
./mfc.sh run -n $(nproc) ./examples/3d_shockdroplet/case.py 
```
where `$(nproc)` is the number of cores the example will run on (and the number of physical cores on your CPU device).
You can visualize the output data in `examples/3d_shockdroplet/silo_hdf5` via Paraview, Visit, or your favorite software.

## Is this _really_ exascale?

[OLCF Frontier](https://www.olcf.ornl.gov/frontier/) is the first exascale supercomputer.
The weak scaling of MFC on this machine shows near-ideal utilization. 

<p align="center">
    <img src="docs/res/scaling.png" alt="Scaling" width="400"/>
</p>


## What else can this thing do

MFC has many features.
They are organized below.

### Physics

* 1-3D
* Compressible
	* Low Mach number treatment available
* Multi- and single-component
	* 4, 5, and 6 equation models for multi-component/phase features
   	* Kapila and Allaire models 5-equation models
* Multi- and single-phase 
	* Phase change via p, pT, and pTg schemes
* Grids
	* 1-3D Cartesian, cylindrical, axisymmetric. 
	* Arbitrary grid stretching for multiple domain regions.
	* Complex/arbitrary geometries via immersed boundary method
	* STL geometry files supported
* Surface tension for multiphase cases
* Sub-grid dynamics
	* Euler-Euler particle models for bubble dynamics and similar
	* Euler-Lagrange bubble dynamics
* Viscous effects (high-order accurate representations)
* Ideal and stiffened gas equations of state
* Body forces
* Acoustic wave generation (one- and two-way sound sources)
* Magnetohydrodynamics (MHD)
* Relativistic Magnetohydrodynamics (RMHD)

### Numerics

* Shock and interface capturing schemes
	* First-order upwinding
 	* WENO reconstructions of order 3, 5, and 7
  	* WENO variants: WENO-JS, WENO-M, WENO-Z, TENO
   	* Monotonicity-preserving reconstructions
	* Reliable handling of large density ratios
* Exact and approximate (e.g., HLL, HLLC, HLLD) Riemann solvers
* Boundary conditions
	* Periodic, reflective, extrapolation/Neumann
	* Slip and no-slip
 	* Thompson-based characteristic BCs: non-reflecting sub/supersonic buffers, inflows, outflows
	* Generalized characteristic relaxation boundary conditions
* Runge-Kutta orders 1-3 (SSP TVD), adaptive time stepping
* RK4-5 operator splitting for Euler-Lagrange modeling
* Interface sharpening (THINC-like)

### Large-scale and accelerated simulation

* GPU compatible on NVIDIA ([P/V/A/H]100, GH200, etc.) and AMD (MI[1/2/3]00+) GPU and APU hardware
* Ideal weak scaling to 100% of the largest GPU and superchip supercomputers
 	* \>36K AMD APUs (MI300A) on [LLNL El Capitan](https://hpc.llnl.gov/hardware/compute-platforms/el-capitan)
   	* \>3K AMD APUs (MI300A) on [LLNL Tuolumne](https://hpc.llnl.gov/hardware/compute-platforms/tuolumne)
	* \>33K AMD GPUs (MI250X) on the first exascale computer, [OLCF Frontier](https://www.olcf.ornl.gov/frontier/) 
	* \>10K NVIDIA GPUs (V100) on [OLCF Summit](https://www.olcf.ornl.gov/summit/) 
* Near compute roofline behavior
* RDMA (remote data memory access; GPU-GPU direct communication) via GPU-aware MPI on NVIDIA (CUDA-aware MPI) and AMD GPU systems
* Optional single-precision computation and storage

### Software robustness and other features

* [Fypp](https://fypp.readthedocs.io/en/stable/fypp.html) metaprogramming for code readability, performance, and portability
* Continuous Integration (CI)
	* \>300 Regression tests with each PR.
 		* Performed with GNU (GCC), Intel (oneAPI), Cray (CCE), and NVIDIA (NVHPC) compilers on NVIDIA and AMD GPUs.
		* Line-level test coverage reports via [Codecov](https://app.codecov.io/gh/MFlowCode/MFC) and `gcov`
	* Benchmarking to avoid performance regressions and identify speed-ups
* Continuous Deployment (CD) of [website](https://mflowcode.github.io) and [API documentation](https://mflowcode.github.io/documentation/index.html)

## Citation

If you use MFC, consider citing it as below.
Ref. 1 includes all modern MFC features, including GPU acceleration and many new physics features.
If referencing MFC's (GPU) performance, consider citing ref. 1 and 2, which describe the solver and how it was crafted.
The original open-source release of MFC is ref. 3, which should be cited for provenance as appropriate.

```bibtex
@article{Wilfong_2025,
  author = {Wilfong, Benjamin and {Le Berre}, Henry and Radhakrishnan, Anand and Gupta, Ansh and Vaca-Revelo, Diego and Adam, Dimitrios and Yu, Haocheng and Lee, Hyeoksu and Chreim, Jose Rodolfo and {Carcana Barbosa}, Mirelys and Zhang, Yanjun and Cisneros-Garibay, Esteban and Gnanaskandan, Aswin and {Rodriguez Jr.}, Mauro and Budiardja, Reuben D. and Abbott, Stephen and Colonius, Tim and Bryngelson, Spencer H.},
  title = {{MFC 5.0: A}n exascale many-physics flow solver},
  journal = {arXiv preprint arXiv:2503.07953},
  year = {2025},
  doi = {10.48550/arXiv.2503.07953}
}

@article{Radhakrishnan_2024,
  title   = {Method for portable, scalable, and performant {GPU}-accelerated simulation of multiphase compressible flow},
  author  = {A. Radhakrishnan and H. {Le Berre} and B. Wilfong and J.-S. Spratt and M. {Rodriguez Jr.} and T. Colonius and S. H. Bryngelson},
  journal = {Computer Physics Communications},
  year    = {2024},
  volume  = {302},
  pages   = {109238},
  doi     = {10.1016/j.cpc.2024.109238}
}

@article{Bryngelson_2021,
  title   = {{MFC: A}n open-source high-order multi-component, multi-phase, and multi-scale compressible flow solver},
  author  = {S. H. Bryngelson and K. Schmidmayer and V. Coralic and J. C. Meng and K. Maeda and T. Colonius},
  journal = {Computer Physics Communications},
  year    = {2021},
  volume  = {266},
  pages   = {107396},
  doi     = {10.1016/j.cpc.2020.107396}
}
```

## License
 
Copyright 2021 Spencer Bryngelson and Tim Colonius.
MFC is under the MIT license (see [LICENSE](LICENSE) for full text).

## Acknowledgements

Federal sponsors have supported MFC development, including the US Department of Defense (DOD), the National Institutes of Health (NIH), the Department of Energy (DOE), and the National Science Foundation (NSF).

MFC computations have used many supercomputing systems. A partial list is below
  * OLCF Frontier and Summit, and testbeds Wombat, Crusher, and Spock (allocation CFD154, PI Bryngelson)
  * LLNL Tuolumne and Lassen, El Capitan early access system Tioga
  * PSC Bridges(1/2), NCSA Delta, SDSC Comet and Expanse, Purdue Anvil, TACC Stampede(1-3), and TAMU ACES via ACCESS-CI allocations from Bryngelson, Colonius, Rodriguez, and more.
  * DOD systems Onyx, Carpenter, Nautilus, and Narwhal via the DOD HPCMP program
  * Sandia National Labs systems Doom and Attaway and testbed systems Weaver and Vortex


## Contributors ✨

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/sbryngelson"><img src="https://avatars.githubusercontent.com/u/39740593?v=4?s=100" width="100px;" alt="Spencer Bryngelson"/><br /><sub><b>Spencer Bryngelson</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=sbryngelson" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/henryleberre"><img src="https://avatars.githubusercontent.com/u/24208000?v=4?s=100" width="100px;" alt="Henry Le Berre"/><br /><sub><b>Henry Le Berre</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=henryleberre" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/wilfonba"><img src="https://avatars.githubusercontent.com/u/48168887?v=4?s=100" width="100px;" alt="Ben Wilfong"/><br /><sub><b>Ben Wilfong</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=wilfonba" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/haochey"><img src="https://avatars.githubusercontent.com/u/98496194?v=4?s=100" width="100px;" alt="Haocheng Yu"/><br /><sub><b>Haocheng Yu</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=haochey" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/hyeoksu-lee"><img src="https://avatars.githubusercontent.com/u/121892187?v=4?s=100" width="100px;" alt="Hyeoksu Lee"/><br /><sub><b>Hyeoksu Lee</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=hyeoksu-lee" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/anandrdbz"><img src="https://avatars.githubusercontent.com/u/62814442?v=4?s=100" width="100px;" alt="anandrdbz"/><br /><sub><b>anandrdbz</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=anandrdbz" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/anshgupta1234"><img src="https://avatars.githubusercontent.com/u/38015438?v=4?s=100" width="100px;" alt="Ansh Gupta"/><br /><sub><b>Ansh Gupta</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=anshgupta1234" title="Code">💻</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/ChrisZYJ"><img src="https://avatars.githubusercontent.com/u/120074479?v=4?s=100" width="100px;" alt="Yanjun Zhang"/><br /><sub><b>Yanjun Zhang</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=ChrisZYJ" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/okBrian"><img src="https://avatars.githubusercontent.com/u/63064763?v=4?s=100" width="100px;" alt="Brian Ok"/><br /><sub><b>Brian Ok</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=okBrian" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/RasmitDevkota"><img src="https://avatars.githubusercontent.com/u/20760979?v=4?s=100" width="100px;" alt="Rasmit Devkota"/><br /><sub><b>Rasmit Devkota</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=RasmitDevkota" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/js-spratt"><img src="https://avatars.githubusercontent.com/u/37048159?v=4?s=100" width="100px;" alt="Jean-Sebastien Spratt"/><br /><sub><b>Jean-Sebastien Spratt</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=js-spratt" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://www.linkedin.com/in/malmahrouqi3"><img src="https://avatars.githubusercontent.com/u/145478595?v=4?s=100" width="100px;" alt="Mohammed S. Al-Mahrouqi"/><br /><sub><b>Mohammed S. Al-Mahrouqi</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=Malmahrouqi3" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/AiredaleDev"><img src="https://avatars.githubusercontent.com/u/71239694?v=4?s=100" width="100px;" alt="Cameron"/><br /><sub><b>Cameron</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=AiredaleDev" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/XZTian64"><img src="https://avatars.githubusercontent.com/u/195029395?v=4?s=100" width="100px;" alt="Xuzheng Tian"/><br /><sub><b>Xuzheng Tian</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=XZTian64" title="Code">💻</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/DimAdam-01"><img src="https://avatars.githubusercontent.com/u/153962156?v=4?s=100" width="100px;" alt="Dimitrios Adam"/><br /><sub><b>Dimitrios Adam</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=DimAdam-01" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/JRChreim"><img src="https://avatars.githubusercontent.com/u/22774282?v=4?s=100" width="100px;" alt="JRChreim"/><br /><sub><b>JRChreim</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=JRChreim" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/max-Hawkins"><img src="https://avatars.githubusercontent.com/u/37495064?v=4?s=100" width="100px;" alt="Max Hawkins"/><br /><sub><b>Max Hawkins</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=max-Hawkins" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/prathi-wind"><img src="https://avatars.githubusercontent.com/u/128338570?v=4?s=100" width="100px;" alt="Tanush Prathi"/><br /><sub><b>Tanush Prathi</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=prathi-wind" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/dgvacarevelo"><img src="https://avatars.githubusercontent.com/u/143528660?v=4?s=100" width="100px;" alt="Diego Vaca Revelo"/><br /><sub><b>Diego Vaca Revelo</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=dgvacarevelo" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/belericant"><img src="https://avatars.githubusercontent.com/u/17969913?v=4?s=100" width="100px;" alt="Eric Dong"/><br /><sub><b>Eric Dong</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=belericant" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://vivo.brown.edu/display/mrodri97"><img src="https://avatars.githubusercontent.com/u/19245739?v=4?s=100" width="100px;" alt="Mauro Rodriguez Jr"/><br /><sub><b>Mauro Rodriguez Jr</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=mrodrig6" title="Code">💻</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/aricer123"><img src="https://avatars.githubusercontent.com/u/107273558?v=4?s=100" width="100px;" alt="aricer123"/><br /><sub><b>aricer123</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=aricer123" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Sam-Briney"><img src="https://avatars.githubusercontent.com/u/46232496?v=4?s=100" width="100px;" alt="Sam-Briney"/><br /><sub><b>Sam-Briney</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=Sam-Briney" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/abbotts"><img src="https://avatars.githubusercontent.com/u/13371405?v=4?s=100" width="100px;" alt="Steve Abbott"/><br /><sub><b>Steve Abbott</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=abbotts" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/arjunj05"><img src="https://avatars.githubusercontent.com/u/88004187?v=4?s=100" width="100px;" alt="arjunj05"/><br /><sub><b>arjunj05</b></sub></a><br /><a href="https://github.com/MFlowCode/MFC/commits?author=arjunj05" title="Code">💻</a></td>
    </tr>
  </tbody>
  <tfoot>
    <tr>
      <td align="center" size="13px" colspan="7">
        <img src="https://raw.githubusercontent.com/all-contributors/all-contributors-cli/1b8533af435da9854653492b1327a23a4dbd0a10/assets/logo-small.svg">
          <a href="https://all-contributors.js.org/docs/en/bot/usage">Add your contributions</a>
        </img>
      </td>
    </tr>
  </tfoot>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

Thanks goes to these wonderful people ✨
