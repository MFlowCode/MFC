#!/bin/bash
#BSUB -P CFD154
#BSUB -W 20
#BSUB -nnodes 1




export PGI_ACC_TIME=0
export PGI_ACC_NOTIFY=
export PGI_ACC_DEBUG=0

module restore MFC-GPU
for proc in 1
do
	sizex=$((100*1))
	sizey=$((100*1))
	sizez=$((100*1))	
	sed -i 's/m = .*/m = '"$sizex"'/' simulation.inp
	sed -i 's/n = .*/n = '"$sizey"'/' simulation.inp
	sed -i 's/p = .*/p = '"$sizez"'/' simulation.inp
	sed -i 's/t_step_stop = .*/t_step_stop = 10/' simulation.inp
	sed -i 's/precision = .*/precision = 2/' simulation.inp
	sed -i 's/m = .*/m = '"$sizex"'/' pre_process.inp
	sed -i 's/n = .*/n = '"$sizey"'/' pre_process.inp
	sed -i 's/p = .*/p = '"$sizez"'/' pre_process.inp
	jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 ../../src/pre_process_code/pre_process
	time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1  ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1341_gpu -f -o profile_anand_weno ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1565_gpu -f -o profile_anand_mp_weno ../../src/simulation_code/simulation
	#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_riemann_solvers_s_hllc_riemann_solver_acc_4212_gpu -f -o profile_anand_riemann ../../src/simulation_code/simulation
done
