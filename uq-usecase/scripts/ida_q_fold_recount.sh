#!/bin/bash
# Fold recount on the Q-limited dispatch at 200 and 400 MW (ida_q_fold_recount.py).
# Fixed-state fold scan at alphaP 0.8, X = 0.01, all buses. usage: sbatch ida_q_fold_recount.sh
#SBATCH --job-name=ida_q_fold
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/q_limit/fold_recount_slurm.log

set -u
# One BLAS thread per worker; avoids oversubscription with 16 Pool workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?}"
$PY "$S/ida_q_fold_recount.py"
echo "done"
