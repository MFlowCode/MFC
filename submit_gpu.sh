
## Pre-process
#./mfc.sh run /scratch/bciv/sremillard/sph_col_20atm_BD/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# sph20prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_35atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# sph30prp -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 10q_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_hf/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20fh_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20q_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_half/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20h_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_50atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 50q_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20p1_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20p1_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_35atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 35p1_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_10atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# sph10_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_hf/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 10fh_prpN -t pre_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_p1e/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 0:15:00 -# 20p1e_prpN -t pre_process -a bciv-delta-gpu -c delta





## simulation
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_hf/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 08:00:00 -# 10fh_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/sph_col_10atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 10:00:00 -# sph10_sim -t simulation -a bciv-delta-gpu -c delta




#./mfc.sh run /scratch/bciv/sremillard/sph_col_35atm/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 09:00:00 -# sph35sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 10:00:00 -# 20p1_sim -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 09:00:00 -# 10p1_sim -t simulation -a bciv-delta-gpu -c delta


#./mfc.sh run /scratch/bciv/sremillard/sph_col_20atm_BD/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 04:30:00 -# sph20 -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 01:00:00 -# 10q -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_hf/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 07:30:00 -# 20fh -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 07:30:00 -# 20q -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_half/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 10:00:00 -# 20h -t simulation -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_50atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 5 -n 4 -g 1 -w 10:00:00 -# 50q -t simulation -a bciv-delta-gpu -c delta

## post_process - Check if sim_data is TRUE!!!!!

#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_p1/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 03:00:00 -# 20p1_sim -t post_process -a bciv-delta-gpu -c delta


./mfc.sh run /scratch/bciv/sremillard/sph_col_20atm_BD/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 01:30:00 -# sph20 -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 01:30:00 -# 10q -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_10atm_hf/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 02:30:00 -# 10fh -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /projects/bciv/sremillard/pert_20q/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 03:00:00 -# 20q -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_20atm_half/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 01:30:00 -# 20h -t post_process -a bciv-delta-gpu -c delta
#./mfc.sh run /scratch/bciv/sremillard/pert_col_50atm_quart/3Dinput.py -e batch -p gpuA100x4 -N 1 -n 4 -g 1 -w 01:30:00 -# 50q -t post_process -a bciv-delta-gpu -c delta



