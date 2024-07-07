#!/bin/bash

#./mfc.sh run ./examples/3D_phasechange_bubble/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 8 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./examples/3D_phasechange_bubble/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/mcarcanabarbosa/testingpc4f/4testing/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta


