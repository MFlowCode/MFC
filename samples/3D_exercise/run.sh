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
	time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 ../../.mfc/___current___/build/bin/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 ../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_902_gpu -f -o profile_anand_alt ../../src/simulation_code/simulation
done
