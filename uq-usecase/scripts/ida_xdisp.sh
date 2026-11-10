#!/bin/bash
# Cross-dispatch IDA bus sweep at alphaP 0.8, 200 MW balanced cases (hunt/<d>/cases/unlim_mw200.json).
# Tests whether IDA IC failures match fold sets on dispatches other than h1879.
# Runs ida_zoom.py --stage bus per dispatch. Compare to hunt/<d>/<d>.csv fold rows (mode unlim, 200 MW).
#SBATCH --job-name=ida_xdisp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=03:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-campaigns/ida_zoom/xdisp/slurm_%j.log

set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
Z=/scratch/isatkaus/gridkit-campaigns/ida_zoom

echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?}"
for d in h1879 h2151 h4770; do
  src="$Z/xdisp/$d"
  mkdir -p "$src"
  cp "$Z/hunt/$d/cases/unlim_mw200.json" "$src/case.json"
  echo "== $d"
  $PY -u "$S/ida_zoom.py" --src "$src" --name "xdisp_$d" --stage bus --alphas 0.8 --fault-bus 1
done
echo "done"
