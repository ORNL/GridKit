#!/bin/bash
# Fault-on nose trace at alphaP 0.8 (ida_nose_trace.py): 2 buses, fixed-state continuation.
# usage: sbatch ida_nose_trace.sh
#SBATCH --job-name=ida_nose_trace
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/nose_trace_slurm.log

set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "host=$(hostname)"
$PY "$S/ida_nose_trace.py"
echo "done"
