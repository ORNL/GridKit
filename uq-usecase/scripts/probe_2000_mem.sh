#!/bin/bash
# Single-process memory probe at ACTIVSg2000: one dynamic run, then one stage-A batch.
# usage: cd /scratch/$USER/scale-probe && sbatch /home/isatkaus/gridkit/uq-usecase/scripts/probe_2000_mem.sh
#SBATCH --job-name=probe2k
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=probe2k_%j.out

PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts

echo "== dynamic run"
/usr/bin/time -v "$PY" "$S/probe_2000_mem.py" dyn 2>&1 | grep -E "^\{|Maximum resident|Elapsed|Exit status|Traceback|Error"
echo "== stage A (dc-opf, ac)"
/usr/bin/time -v "$PY" "$S/probe_2000_mem.py" stagea 2>&1 | grep -E "^rc=|Maximum resident|Elapsed|Exit status|Traceback|Error"
