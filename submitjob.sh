#!/bin/bash

#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/base/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# basepre  -t pre_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s1/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# s1pre  -t pre_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s2/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# s2pre  -t pre_process -b mpirun
./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s3/input.py -e batch -p GPU-shared -N 1 -n 2 -g 1 -w 01:00:00 -# s3pre  -t pre_process #-b mpirun
./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s4/input.py -e batch -p GPU-shared -N 1 -n 2 -g 1 -w 01:00:00 -# s4pre  -t pre_process #-b mpirun


#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/base/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 24:00:00 -# basesim  -t simulation -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s1/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 24:00:00 -# s1sim  -t simulation -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s2/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 24:00:00 -# s2sim  -t simulation -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s3/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 8:00:00 -# s3sim  -t simulation -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s4/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 8:00:00 -# s4sim  -t simulation -b mpirun


#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/base/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# basepost  -t post_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s1/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# s1post  -t post_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s2/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 01:00:00 -# s2post  -t post_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s3/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 02:00:00 -# s3post  -t post_process -b mpirun
#./mfc.sh run /ocean/projects/mch220006p/mrdz/dfd2023/s4/input.py -e batch -p GPU-shared -N 1 -n 4 -g 4 -w 02:00:00 -# s4post  -t post_process -b mpirun

