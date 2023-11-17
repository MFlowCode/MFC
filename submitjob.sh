#!/bin/bash

#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/base/input.py -e batch -p wholenode -w 01:00:00 -# basepre -N 8 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s1/input.py -e batch -p wholenode -w 01:00:00 -# s1pre -N 8 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s2/input.py -e batch -p wholenode -w 01:00:00 -# s2pre -N 8 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s3/input.py -e batch -p wholenode -w 01:00:00 -# s3pre -N 8 -n 128 -t pre_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s4/input.py -e batch -p wholenode -w 01:00:00 -# s4pre -N 8 -n 128 -t pre_process -b mpirun


#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/base/input.py -e batch -p wholenode -w 24:00:00 -# basesim -N 8 -n 128 -t simulation -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s1/input.py -e batch -p wholenode -w 24:00:00 -# s1sim -N 8 -n 128 -t simulation -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s2/input.py -e batch -p wholenode -w 24:00:00 -# s2sim -N 8 -n 128 -t simulation -b mpirun
./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s3/input.py -e batch -p wholenode -w 8:00:00 -# s3sim -N 8 -n 128 -t simulation -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s4/input.py -e batch -p wholenode -w 8:00:00 -# s4sim -N 8 -n 128 -t simulation -b mpirun


#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/base/input.py -e batch -p wholenode -w 01:00:00 -# basepost -N 8 -n 128 -t post_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s1/input.py -e batch -p wholenode -w 01:00:00 -# s1post -N 8 -n 128 -t post_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s2/input.py -e batch -p wholenode -w 01:00:00 -# s2post -N 8 -n 128 -t post_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s3/input.py -e batch -p wholenode -w 02:00:00 -# s3post -N 8 -n 128 -t post_process -b mpirun
#./mfc.sh run /anvil/projects/x-mch220010/mrdz/dfd2023/s4/input.py -e batch -p wholenode -w 02:00:00 -# s4post -N 8 -n 128 -t post_process -b mpirun

