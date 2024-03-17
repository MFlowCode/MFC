# This script is ran from the 3D_weak_scaling case directory after running
# MFC with --omni -n <name>. To analyze, run chmod u+x ./analyze.sh followed
# by ./analyze.sh <name>

omniperf analyze -p workloads/$1/mi200 --metric 0 7.1.5 7.1.6 7.1.7 7.1.8 7.1.9 16.3.1 16.3.2 16.3.7 17.3.2 17.3.3 17.3.8
