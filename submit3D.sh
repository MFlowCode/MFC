#!/bin/bash


#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbl_base/input.py -e batch -p wholenode -w 01:00:00 -# 3Dmfc -N 6 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbl_base/input.py -e batch -p wholenode -w 06:00:00 -# 3Dmfc -N 6 -n 128 -t simulation -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbl_base/input.py -e batch -p wholenode -w 00:30:00 -# 3Dmfc -N 6 -n 128 -t post_process -b mpirun

#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbg_base/input.py -e batch -p wholenode -w 01:00:00 -# 3Dmfc -N 6 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbg_base/input.py -e batch -p wholenode -w 06:00:00 -# 3Dmfc -N 6 -n 128 -t simulation -b mpirun
./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/3Gbg_base/input.py -e batch -p wholenode -w 00:30:00 -# 3Dmfc -N 6 -n 128 -t post_process -b mpirun


