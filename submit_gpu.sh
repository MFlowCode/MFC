
## Pre-process
#./mfc.sh run /scratch/bciv/sremillard/sph_col_2atm/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# sph2 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_50atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 1:00:00 -# sph50 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_20atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 1:00:00 -# sph20 -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_10atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 1:00:00 -# sph10 -t pre_process -a bciv-delta-gpu -c delta

#./mfc.sh run /scratch/bciv/sremillard/pert_col_2atm_tenth/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# base -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_1atm/3Dshinput_new.py -e batch -p gpuA100x4 -N 8 -n 4 -g 4 -w 1:00:00 -# sph_1 -t pre_process -a bciv-delta-gpu -c delta


## simulation
./mfc.sh run /scratch/bciv/sremillard/sph_col_50atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 5:30:00 -# sph50 -t simulation -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/sremillard/sph_col_20atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 5:30:00 -# sph20 -t simulation -a bciv-delta-gpu -c delta
./mfc.sh run /scratch/bciv/sremillard/sph_col_10atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 4 -w 5:30:00 -# sph10 -t simulation -a bciv-delta-gpu -c delta


