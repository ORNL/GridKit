#!/bin/bash
#SBATCH --job-name=afront2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/kfs2/projects/scidac/scidac-data/gridkit-campaigns/afront2/slurm/afront2.log

set -u
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
OUT=/kfs2/projects/scidac/scidac-data/gridkit-campaigns/afront2/slurm
TOTAL=136
W=16
PER=$(( (TOTAL + W - 1) / W ))

echo "host=$(hostname) workers=$W per=$PER total=$TOTAL"
for i in $(seq 0 $((W-1))); do
  START=$(( i * PER ))
  [ "$START" -ge "$TOTAL" ] && break
  $PY "$S/alphap_frontier.py" --name afront2 \
      --hours 6164,7549 --fault-buses 66,108,148,177 \
      --start "$START" --count "$PER" > "$OUT/w${i}.log" 2>&1 &
done
wait
echo "all workers done"
