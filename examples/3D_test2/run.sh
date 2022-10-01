#!/bin/bash
#BSUB -P CFD154
#BSUB -W 15
#BSUB -nnodes 1




export PGI_ACC_TIME=0
export PGI_ACC_NOTIFY=
export PGI_ACC_DEBUG=0

module restore MFC-GPU
for proc in 1
do
	jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1  ../../.mfc/___current___/build/bin/pre_process 
	time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../.mfc/___current___/build/bin/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 ../src/simulation/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../src/simulation/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1341_gpu -f -o profile_anand_weno ../../src/simulation/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1565_gpu -f -o profile_anand_mp_weno ../../src/simulation/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_riemann_solvers_s_hllc_riemann_solver_acc_2539_gpu -f -o profile_anand_riemann ../../.mfc/___current___/build/bin/simulation
done
