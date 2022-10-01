#!/bin/bash




export PGI_ACC_TIME=0
export PGI_ACC_NOTIFY=
export PGI_ACC_DEBUG=0

srun -N 1 -n 1 -G 1 ../../build/___current___/build/bin/pre_process 
time srun -N 1 -n 1 -G 1 ../../build/___current___/build/bin/simulation
#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 ../src/simulation/simulation
#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nsys profile --stats=true ../../src/simulation/simulation
#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1341_gpu -f -o profile_anand_weno ../../src/simulation/simulation
#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_weno_s_weno_alt_1565_gpu -f -o profile_anand_mp_weno ../../src/simulation/simulation
#time jsrun --smpiargs="-gpu" -r$proc -a1 -c1 -g1 nv-nsight-cu-cli --set full -k m_riemann_solvers_s_hllc_riemann_solver_acc_2539_gpu -f -o profile_anand_riemann ../../.mfc/___current___/build/bin/simulation
