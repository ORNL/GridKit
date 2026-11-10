#!/bin/bash
# CCT at alphaP 0.2 vs 1.0 over the alphap_frontier design (272 dispatches).
# usage: sbatch cctcheck.sh <campaign-name> <fault-buses>
#   e.g. sbatch cctcheck.sh cct181 181
#SBATCH --job-name=cctcheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short

set -u
NAME=$1; FB=$2
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
OUT=/kfs2/projects/scidac/scidac-data/gridkit-campaigns/$NAME/slurm
TOTAL=272   # alphap_frontier design: 17 buses x 4 hours x 2 fracs x 2 MW
W=16
PER=$(( (TOTAL + W - 1) / W ))
mkdir -p "$OUT"

echo "host=$(hostname) name=$NAME fault=$FB workers=$W per=$PER total=$TOTAL"
for i in $(seq 0 $((W-1))); do
  START=$(( i * PER ))
  [ "$START" -ge "$TOTAL" ] && break
  $PY "$S/cct_check.py" --name "$NAME" --fault-buses "$FB" --alphas 0.2,1.0 \
      --start "$START" --count "$PER" > "$OUT/w${i}.log" 2>&1 &
done
wait
echo "all workers done"
