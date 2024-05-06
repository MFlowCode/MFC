#./mfc.sh run /scratch/bciv/sremillard/sph_col_32/3Dshinput_new.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 1:00:00 -# test1 -t simulation -a bciv-delta-gpu -c delta 

## Pre-process
./mfc.sh run /scratch/bciv/sremillard/sph_col_2atm/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# sph2 -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/sremillard/sph_col_5atm/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# sph5 -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/sremillard/pert_col_2atm_tenth/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# base -t pre_process -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/sremillard/sph_col_1atm/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# sph_1 -t pre_process -a bciv-delta-gpu -c delta
