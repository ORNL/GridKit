#!/bin/bash
# Hunt-matrix fold scans (ida_hunt_matrix.py). Forwards args to the driver.
# usage: sbatch [--output=LOG] ida_hunt_matrix.sh --dispatch h1879 --modes qlim --mw 100 150 ... --tag NAME
#SBATCH --job-name=ida_hunt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/hunt/slurm_%j.log

set -u
# One BLAS thread per worker; avoids oversubscription with 16 Pool workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?} args=$*"
$PY -u "$S/ida_hunt_matrix.py" "$@"
echo "done"
