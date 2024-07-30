## How to Profile MFC using `omniperf`

* Start an interactive session with the desired number of nodes and total tasks using `salloc -A [account] -J interactive -t 2:00:00 -p batch -N [nnodes] -n [total tasks]`

* Generate MFC input files by running `./mfc.sh run [path to casefile] -N [nnodes] -n  [total tasks] --gpu -t pre_process simulation --case-optimization`

* Move to the simulation directory using `cd [path to casefile]`

* `module load` the following modules in the order listed
    - rocm/5.5.1
    - cray-python
    - omniperf

* `omniperf profile -n [profile name] -- [path to MFC beginning with /]/build/install/bin/simulation`

* `omniperf analyze --gui -p [path to casefile]/workloads/[profile name]/mi200`

* Determine what login node you're on, call it [node name]

* Open a new terminal window and log into Frontier using `ssh -L8050:localhost:8050 username@[node name].frontier.olcf.ornl.gov`

* Open a web browser and navigate to `http://localhost:8050/`
