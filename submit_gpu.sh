#!/bin/bash

### A100s
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwex_pre -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim2/input.py -e batch -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# bwim_pre -t pre_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwex_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwim_sim -t simulation -a bciv-delta-gpu -c delta

### A40s
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwex_pre -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwim_pre -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/test/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 01:00:00 -# test -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim2/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwim_pre -t pre_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwex_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwim_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwex_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwim_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim2/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 12:00:00 -# bwim_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/test/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 03:00:00 -# test_sim -t simulation -a bciv-delta-gpu -c delta


#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterex/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 00:30:00 -# bwex_post -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 00:30:00 -# bwim_post -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim2/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwim_post -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/mancia2024/bubwaterim2/input.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 01:00:00 -# bwim_post -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/test/input.py -e batch -p gpuA40x4 -N 5 -n 4 -g 1 -w 00:30:00 -# test_post -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/3Dsph_hyper_prestress_input.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/3Dsph_hyper_prestress_input.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/3Dsph_hyper_prestress_input.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hyper/hyper_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hypoe/hypo_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hypoe/hypo_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/cav2024/hypoe/hypo_gel.py -p gpuA40x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test61/case.py -e batch -p gpuA100x4 -N 1 -n 2 -g 0 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test61/case.py -e batch -p gpuA100x4 -N 1 -n 1 -g 0 -w 02:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta 
#./mfc.sh run /scratch/bciv/rodrigu1/ctr2024/test61/case.py -e batch -p gpuA100x4 -N 1 -n 1 -g 0 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

./mfc.sh run /scratch/bciv/rodrigu1/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/rodrigu1/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/rodrigu1/3D_bubble_channel/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/3D_bubble_channel/case.py -p gpuA100x4 -N 1 -n 2 -g 1 -w 01:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/rodrigu1/3D_lungwave/case.py -p gpuA100x4 -N 1 -n 1 -g 1 -w 01:00:00 -# test1 -t post_process -a bciv-delta-gpu -c delta

