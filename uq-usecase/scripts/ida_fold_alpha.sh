#!/bin/bash
# Fold location vs alphaP (ida_fold_alpha.py): 200 buses x 3 alphaP, fixed-state continuation.
# usage: sbatch ida_fold_alpha.sh
#SBATCH --job-name=ida_fold_alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/fold_alpha_slurm.log

set -u
# One BLAS thread per worker; avoids oversubscription with 16 Pool workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?}"
$PY "$S/ida_fold_alpha.py"
echo "done"
