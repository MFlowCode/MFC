#!/bin/bash

#./mfc.sh run /anvil/projects/x-mch220010/undex/input.py -e batch -p wholenode -w 01:00:00 -# undex -N 2 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/rc/input.py -e batch -p wholenode -w 01:00:00 -# rc -N 2 -n 128 -t pre_process -b mpirun


./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/test/input.py -e batch -p wholenode -w 01:00:00 -# undex -N 2 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/test/input.py -e batch -p wholenode -w 01:00:00 -# undex -N 2 -n 128 -t simulation -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/mfcruns/test/input.py -e batch -p wholenode -w 01:00:00 -# undex -N 2 -n 128 -t post_process -b mpirun

