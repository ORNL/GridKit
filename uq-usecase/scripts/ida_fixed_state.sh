#!/bin/bash
# Fixed-state algebraic IC test (ida_fixed_state_test.py), dc200 mode, all 200 fault buses.
# usage: sbatch ida_fixed_state.sh
#SBATCH --job-name=ida_fixed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state_slurm.log

set -u
# One BLAS thread per worker; avoids oversubscription with 16 Pool workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?}"
$PY "$S/ida_fixed_state_test.py"
echo "done"
