#!/bin/bash
#SBATCH --job-name=afront1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/kfs2/projects/scidac/scidac-data/gridkit-campaigns/afront1/slurm/afront1.log

set -u
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
TOTAL=272
W=16
PER=$(( (TOTAL + W - 1) / W ))

echo "host=$(hostname) workers=$W per=$PER total=$TOTAL"
for i in $(seq 0 $((W-1))); do
  START=$(( i * PER ))
  [ "$START" -ge "$TOTAL" ] && break
  $PY "$S/alphap_frontier.py" --name afront1 --start "$START" --count "$PER" \
      > "/kfs2/projects/scidac/scidac-data/gridkit-campaigns/afront1/slurm/w${i}.log" 2>&1 &
done
wait
echo "all workers done"
