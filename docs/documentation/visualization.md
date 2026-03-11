@page visualization Flow visualization

# Flow visualization

After running `post_process` on a simulation (see @ref running "Running"), MFC produces output in either Silo-HDF5 format (`format=1`) or binary format (`format=2`).
These can be visualized using MFC's built-in CLI tool or external tools like ParaView and VisIt.

---

## Quick visualization with `./mfc.sh viz`

MFC includes a built-in visualization command that renders images and videos directly from post-processed output — no external GUI tools needed.

### Basic usage

```bash
# Launch the terminal UI (default mode)
./mfc.sh viz case_dir/

# Launch with a specific variable pre-selected
./mfc.sh viz case_dir/ --var pres
```

The command auto-detects the output format (binary or Silo-HDF5) and dimensionality (1D, 2D, or 3D).
By default it launches an interactive terminal UI that works over SSH.
Use `--interactive` for a browser-based UI (supports 3D), `--png` to save images, or `--mp4` for video.

### Exploring available data

Before plotting, you can inspect what data is available:

```bash
# List all available timesteps
./mfc.sh viz case_dir/ --list-steps

# List all available variables at a given timestep
./mfc.sh viz case_dir/ --list-vars --step 0
```

### Timestep selection

The `--step` argument accepts several formats:

| Format | Example | Description |
|--------|---------|-------------|
| Single | `--step 1000` | One timestep |
| Range | `--step 0:10000:500` | Start:end:stride (inclusive) |
| Last | `--step last` | Most recent available timestep |
| All | `--step all` | Every available timestep |

### Rendering options

Customize the appearance of plots:

```bash
# Custom colormap and color range
./mfc.sh viz case_dir/ --var rho --step 1000 --png --cmap RdBu --vmin 0.5 --vmax 2.0

# Higher resolution
./mfc.sh viz case_dir/ --var pres --step 500 --png --dpi 300

# Logarithmic color scale
./mfc.sh viz case_dir/ --var pres --step 1000 --png --log-scale
```

| Option | Description | Default |
|--------|-------------|---------|
| `--cmap` | Matplotlib colormap name | `viridis` |
| `--vmin` | Minimum color scale value | auto |
| `--vmax` | Maximum color scale value | auto |
| `--dpi` | Image resolution (dots per inch) | 150 |
| `--log-scale` | Use logarithmic color scale | off |
| `--output` | Output directory for images | `case_dir/viz/` |

### 3D slicing

For 3D simulations, `viz` extracts a 2D slice for plotting.
By default, it slices at the midplane along the z-axis.

> [!NOTE]
> To limit memory use, 3D batch rendering is capped at 500 timesteps and
> `--interactive` mode at 50. Use `--step start:end:stride` to stay within
> these limits when processing many steps.

```bash
# Default z-midplane slice
./mfc.sh viz case_dir/ --var pres --step 500 --png

# Slice along the x-axis at x=0.25
./mfc.sh viz case_dir/ --var pres --step 500 --png --slice-axis x --slice-value 0.25

# Slice by array index
./mfc.sh viz case_dir/ --var pres --step 500 --png --slice-axis y --slice-index 50
```

### Video generation

Generate MP4 videos from a range of timesteps:

```bash
# Basic video (10 fps)
./mfc.sh viz case_dir/ --var pres --step 0:10000:100 --mp4

# Custom frame rate
./mfc.sh viz case_dir/ --var rho --step all --mp4 --fps 24

# Video with fixed color range
./mfc.sh viz case_dir/ --var rho --step 0:5000:50 --mp4 --vmin 0.1 --vmax 1.0
```

Videos are saved as `case_dir/viz/<varname>.mp4`.
The color range is automatically computed from the first, middle, and last frames unless `--vmin`/`--vmax` are specified.

### Tiled 1D rendering

For 1D cases, omitting `--var` (or passing `--var all`) renders all variables in a single tiled figure:

```bash
# Tiled plot of all variables at the last timestep
./mfc.sh viz case_dir/ --step last --png

# Equivalent explicit form
./mfc.sh viz case_dir/ --var all --step last --png
```

Each variable gets its own subplot with automatic LaTeX-style axis labels.
Tiled mode is available for 1D and 2D data. For 3D data, omitting `--var` auto-selects the first variable.

### Interactive mode

Launch a browser-based interactive viewer with `--interactive`:

```bash
./mfc.sh viz case_dir/ --interactive

# Custom port
./mfc.sh viz case_dir/ --interactive --port 9000
```

The interactive viewer provides a Dash web UI with:
- Variable and timestep selection
- Live plot updates
- Pan, zoom, and hover inspection

> [!NOTE]
> Interactive mode requires the `dash` Python package.

### Terminal UI (TUI)

The default mode launches a live terminal UI that works over SSH with no browser or port-forwarding needed:

```bash
./mfc.sh viz case_dir/

# Start with a specific variable pre-selected
./mfc.sh viz case_dir/ --var pres
```

The TUI loads all timesteps and renders plots directly in the terminal using Unicode block characters.
It supports 1D and 2D data only (use `--interactive` for 3D).

**Keyboard shortcuts:**

| Key | Action |
|-----|--------|
| `.` / `→` | Next timestep |
| `,` / `←` | Previous timestep |
| `Space` | Toggle autoplay |
| `l` | Toggle logarithmic scale |
| `f` | Freeze / unfreeze color range |
| `↑` / `↓` | Select variable (in sidebar) |
| `c` | Cycle colormap |
| `q` | Quit |

### Plot styling

Axis labels use LaTeX-style math notation — for example, `pres` is labeled as \f$p\f$, `vel1` as \f$u\f$, and `alpha1` as \f$\alpha_1\f$.
Plots use serif fonts and the Computer Modern math font for consistency with publication figures.

### Format selection

The output format is auto-detected from the case directory.
To override:

```bash
./mfc.sh viz case_dir/ --var pres --step 0 --format binary
./mfc.sh viz case_dir/ --var pres --step 0 --format silo
```

> [!NOTE]
> Reading Silo-HDF5 files requires the `h5py` Python package.
> If it is not installed, you will see a clear error message with installation instructions.
> Alternatively, use `format=2` (binary) in your case file to produce binary output, which has no extra dependencies.

### Complete option reference

Run `./mfc.sh viz -h` for a full list of options.

---

## Visualizing with ParaView

ParaView is an open-source interactive parallel visualization and graphical analysis tool for viewing scientific data.
Post-processed data in Silo-HDF5 format (`format=1`) can be opened directly in ParaView.
Paraview 5.11.0 has been confirmed to work with the MFC databases for some parallel environments.
Nevertheless, the installation and configuration of Paraview can be environment-dependent and are left to the user.

The user can launch Paraview and open the index files under `/silo_hdf5/root`.
Once the database is loaded, flow field variables contained in the database can be added to the render view.
Further information on Paraview can be found in its [documentation](https://docs.paraview.org/en/latest/).
The figure below shows the iso-contour of the liquid void fraction (`alpha1`) in the database generated by the example case `3D_sphbubcollapse`.

![](../res/paraview.png)

### Visualizing data in cylindrical coordinates

Visualizing data in cylindrical coordinates requires a coordinate transformation of the raw data in the database file.
In Paraview, this coordinate transformation can be accomplished with the following steps:

1. Apply a `clean to grid` filter to the raw data

2. Apply a `calculator` filter to the cleaned data
    - Set the calculator `attribute type` to point data
    - Check the box for `Coordinate Results`
    - Enter the formula `coordsX*cos(coordsY)*iHat + coordsX*sin(coordsY)*jHat + coordsZ*kHat`
    - click apply

These steps will transform the raw data into cylindrical coordinates.
For many cases, this step will require resizing the render view window.

## Visualizing with VisIt

VisIt is an alternative open-source interactive parallel visualization and graphical analysis tool for viewing scientific data.
Versions of VisIt after 2.6.0 have been confirmed to work with the MFC databases for some parallel environments.
Nevertheless, installation and configuration of VisIt can be environment-dependent and are left to the user.
Further remarks on parallel flow visualization, analysis, and processing of the MFC database using VisIt can also be found in \cite Coralic15 and \cite Meng16.

The user can launch VisIt and open the index files under `/silo_hdf5/root`.
Once the database is loaded, flow field variables contained in the database can be added to the plot.
The figure below shows the iso-contour of the liquid void fraction (`alpha1`) in the database generated by the example case `3D_sphbubcollapse`.
For analysis and processing of the database using VisIt's capability, the user is encouraged to address [VisIt user manual](https://wci.llnl.gov/simulation/computer-codes/visit/manuals).

![](../res/visit.png)

*Iso-contour of the liquid void fraction (`alpha1`) in the database generated by example case `3D_sphbubcollapse`*

## Serial data output

If ``parallel_io = 'F'``, MFC will output the conservative variables to a directory `D/`. 
If multiple cores are used (\f$\mathtt{ppn > 1}\f$), then a separate file is created for each core.
If only one coordinate dimension (`n = 0` and `p = 0`) exists, the primitive variables will also be written to `D/`.
The file names correspond to the variables associated with each equation solved by MFC.
They are written at every `t_step_save` time step.
The conservative variables are

\f[ (\rho \alpha)_{1}, \dots, (\rho\alpha)_{N_c}, \rho u_{1}, \dots, \rho u_{N_d}, E, \alpha_1, \dots, \alpha_{N_c} \f]

and the primitive variables are

\f[ (\rho \alpha)_1, \dots, (\rho\alpha)_{N_c}, u_1, \dots, u_{N_d}, p, \alpha_1, \dots, \alpha_{N_c} \f]

where $N_c$ are the number of components `num_fluids` and $N_d$ is the number of spatial dimensions. 
There are exceptions: if `model_eqns = 3`, then the six-equation model appends these variables with the internal energies of each component.
If there are sub-grid bubbles `bubbles = T`, then the bubble variables are also written. 
These depend on the bubble dynamics model used.
If ``polytropic = 'T'``, then the conservative variables are appended by 

\f[ n_b R_1, n_b \dot{R}_1, \dots, n_b R_{N_b}, n_b \dot{R}_{N_b} \f]

where $n_B$ is the bubble number density, and $N_b$ is the number of bubble sizes (see the matching variable in the input file, `Nb`).
The primitive bubble variables do not include $n_B$:

\f[ R_1, \dot{R}_1, \dots, R_{N_b}, \dot{R}_{N_b} \f]

## Remote Visualization on PACE Phoenix

### Step 1: Setting up your Environment

Begin by downloading the `.zip` file [here](https://www.dropbox.com/scl/fi/bdk8702oas8zqu0mk24vx/paceParaview.zip?rlkey=bov1s6lra0z7dhhrh6etniucx&st=2m9xvls4&dl=0).
This file contains two things:

- A bash script to automate the job submission process and provide instructions for remote connection
- A prebuilt Paraview 5.11 binary

Place the file `paceParview.zip` in your scratch direction on Phoenix and unzip it using `unzip paceParaview.zip`.
Enter the new directory `paceParaview` and run `tar -xvf ParaView-5.11.0-egl-MPI-Linux-Python3.9-x86_64.tar.gz` to decompress the compiled binary.
Now that you have the binary on Phoenix, you must download Paraview 5.11 on your local machine.
Paraview binaries can be downloaded [here](https://www.paraview.org/download/).
Select `v5.11` from the version drop-down bar and install a `5.11.0` version of Paraview.

### Step 2: Customizing the script

While all of the bash script's options could be passed as command-line arguments, hardcoding certain unlikely-to-change options saves time.
The following is a list of required and suggested updates to make to `pace-paraview-server`.

- (Optional) Update line 4 to customize the job name that will show up in the scheduler
- (Required) Update line 6 to point towards the location of your Paraview bin directory
- (Optional) Update line 51 to reflect the default account you'll use to run Paraview jobs
- (Optional) Update line 52 to reflect the default wall time requested for your job

### Step 3: Running pace-paraview-server

Before running `pace-paraview-server` for the first time, you must update its permissions by running `chmod u+x pace-paraview-server` in your command line.
Once this has been done, you can run `./pace-paraview-server` with the following options:

- `--account` specifies the charge account for the job.
If you updated line 51 of `pace-paraview-server` to reflect a default account, this option is optional; otherwise, it is required.
- `--nodes` specifies the number of nodes to request (default 1)
- `--mem` specifies the memory per node to request (default is to request all memory)
- `--gres` specifies the GPU resources to request.
These instructions have only been verified to work with `--gres gpu:V100:2`.
- `--time` specifies the time limit for your job.
This is optional and defaults to whatever is set on line 52 of `pace-paraview-server`.

Once you run `./pace-paraview-server <options>`, it'll take a bit to start up. 
In the meantime, you'll see the below message:

```shell
Submitted batch job <job #>
Waiting for ParaView server to start. This may take several minutes  ...
```

When initializing is done, you should see a dialogue with some recommended next steps, numbered 1-4. 
Below is a slightly altered version of that dialogue:

1) Create the appropriate port forwarding for your local ParaView session to connect with.
* On your local machine, run the following from a terminal: `nodeIdentifier` is the remote node running the ParaView server process (given in the output of the batch script), and `paceSystemIdentifier` is the name of the PACE system (however, this is configured with your `.ssh/config`).
* This terminal session must not be killed for your ParaView session as it maintains the port forwarding.
    * `ssh -L 8722:<nodeIdentifier>:53723 <paceSystemIdentifier>`

2) Once you have `Paraview5.11.0` on your machine, select `File -> Connect..` to open the remote connection dialogue box.
* Double-click the existing configuration if you've already set up the PACE connection.
* Click `Add Server` If you have not set up the PACE connection.
This will create a new dialogue box where you can specify a configuration name and set the `Port` to `8722`.
Once this is done, click `configure` and then `save` on the next dialogue box.


<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-02-13</div>
