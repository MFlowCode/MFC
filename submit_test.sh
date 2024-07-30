#!/bin/bash

#./mfc.sh run ./examples/3D_phasechange_bubble/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 8 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./examples/3D_phasechange_bubble/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run ./examples/3D_ctr_test/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/2speed/3dpc-noel-ptg.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/2speed/3dpc-noel-ptg.py -p batch -N 1 -n 4 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/2speed/3dpc-noel-ptg.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/4speed/3dpc-noel-ptg.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 00:10:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/4speed/3dpc-noel-ptg.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 04:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/testingpc4f/4speed/3dpc-noel-ptg.py -e batch -p gpuA100x4 -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t post_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.3/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.4/case.py -e batch -p gpuA100x4 -N 1 -n 4 -g 4 -w 01:00:00 -# pre_bubingel -t pre_process -a bciv-delta-gpu -c delta 
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.4/case.py -e batch -p gpuA100x4 -N 1 -n 4 -g 4 -w 06:00:00 -# sim_bubingel -t simulation -a bciv-delta-gpu -c delta 

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.5/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.5/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test6.5/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -e batch -p gpuA100x4 -N 1 -n 4 -g 4 -w 06:00:00 -# pre_bubinwater -t pre_process -a bciv-delta-gpu -c delta 
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4/case.py -e batch -p gpuA100x4 -N 1 -n 4 -g 4 -w 06:00:00 -# sim_bubinwater -t simulation -a bciv-delta-gpu -c delta 

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4.1/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4.1/case.py -p batch -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t simulation -c delta
./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test4.1/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta


#./mfc.sh run /users/mrodri97/scratch/ctr2024/test6/case.py -p batch -N 1 -n 16 -g 0 -w 01:00:00 -# test1 -t pre_process -c oscar
#./mfc.sh run /users/mrodri97/scratch/ctr2024/test6/case.py -p batch -N 1 -n 16 -g 0 -w 02:00:00 -# test1 -t simulation -c oscar
#./mfc.sh run /users/mrodri97/scratch/ctr2024/test6/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c oscar


#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/elcom/bubliq/25wv/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta
#./mfc.sh run /scratch/bciv/mcarcanabarbosa/ctr2024/finalruns/elcom/bubliq/75wv/case.py -p batch -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -c delta

