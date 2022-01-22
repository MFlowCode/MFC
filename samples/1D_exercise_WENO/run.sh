#!/bin/bash
# module purge
# module load lsf-tools/2.0 DefApps nvhpc/21.9 spectrum-mpi/10.3.1.2-20200121 cuda/11.2.0




export PGI_ACC_TIME=0
export PGI_ACC_NOTIFY=
export PGI_ACC_DEBUG=0

for proc in 1
do
	#size=$((100 * proc))
	#sed -i 's/m = .*/m = '"$size"'/' simulation.inp
	#sed -i 's/m = .*/m = '"$size"'/' pre_process.inp
	jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1  ../../.mfc/___current___/build/bin/MFC_PreProcess 
	time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1  ../../.mfc/___current___/build/bin/MFC_Simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_902_gpu -f -o profile_anand_alt ../../src/simulation_code/simulation
done
