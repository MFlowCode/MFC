This example case contains an automated convergence test using a 1D, two-component advection case.
The case can be run by executing the bash script `./submitJobs.sh` in a terminal after enabling
execution permissions with `chmod +x ./submitJobs.sh` and setting the `ROOT_DIR` and `MFC_DIR`
variables. By default the script runs the case for 6 different grid resolutions with 1st,
3rd, and 5th, order spatial reconstructions. These settings can be modified by editing the variables
at the top of the script. You can also run different model equations by setting the `ME` variable
and different Riemann solvers by setting the `RS` variable.

Once the simulations have been run, you can generate convergence plots with matplotlib by running
`python3 plot.py` in a terminal. This will generate plots of the L1, L2, and Linf error norms and
save the results to `errors.csv`.
