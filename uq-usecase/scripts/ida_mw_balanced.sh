#!/bin/bash
# Balanced datacenter-MW sweep at bus 148 (ida_mw_balanced.py): PF re-solve, case.json per MW,
# fixed-state fold scan at alphaP 0.8. usage: sbatch ida_mw_balanced.sh
#SBATCH --job-name=ida_mw_balanced
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/mw_balanced_slurm.log

set -u
# One BLAS thread per worker; avoids oversubscription with 16 Pool workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?}"
$PY "$S/ida_mw_balanced.py"
echo "done"
