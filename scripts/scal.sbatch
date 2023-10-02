#!/bin/bash -l
#SBATCH --partition=ccm
#SBATCH --constraint=skylake&opa
#SBATCH --time=10000

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=8

source ~/.bashrc.local
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

WORK_DIR=./
OUT_DIR=./vis/

cd ${WORK_DIR}
mkdir -p ${OUT_DIR}
env | grep -i slurm

# full run
#time mpirun --report-bindings -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}./bin/mobility -Nobj 512 -loop_rad 0.45 -ts_tol 1e-7 -geom_tol 1e-8 -gmres_tol 1e-9 -quad_tol 1e-9 -omp ${OMP_NUM_THREADS} -out_path "${OUT_DIR}" -precond ""

# one time-step
time mpirun --report-bindings -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}./bin/mobility -geom "data/rigidgeom-sed512" -dt 0.15 -T 0.15 -loop_rad 0.45 -ts_tol 1e-7 -geom_tol 1e-8 -gmres_tol 1e-9 -quad_tol 1e-9 -omp ${OMP_NUM_THREADS} -out_path "${OUT_DIR}" -precond ""

