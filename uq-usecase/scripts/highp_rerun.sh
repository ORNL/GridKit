#!/bin/bash
# Re-run the datacenter (sited) arm of an existing campaign at another dc alphaP.
# Same conditions, same dispatches; baseline rows (no datacenter) are dropped
# because they do not depend on alphaP.
# usage: sbatch highp_rerun.sh <src-campaign> <dst-campaign> <alphaP> <workers>
#   e.g. sbatch highp_rerun.sh faultscreen1 faultscreen1_ap1 1.0 16
#SBATCH --job-name=highp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --account=scidac
#SBATCH --partition=short

set -u
SRC=$1; DST=$2; ALPHA=$3; W=$4
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
ROOT=/kfs2/projects/scidac/scidac-data/gridkit-campaigns
OUT=$ROOT/$DST/slurm
mkdir -p "$OUT"

# sited rows only: dc_bus is column 6
awk -F, 'NR==1 || $6!=""' "$ROOT/$SRC/conditions.csv" > "$ROOT/$DST/conditions.csv"
N=$(( $(wc -l < "$ROOT/$DST/conditions.csv") - 1 ))
PER=$(( (N + W - 1) / W ))
echo "host=$(hostname) src=$SRC dst=$DST alphaP=$ALPHA sited=$N workers=$W per=$PER"

for i in $(seq 0 $((W-1))); do
  START=$(( i * PER ))
  [ "$START" -ge "$N" ] && break
  $PY "$S/sweep_worker.py" --conditions "$ROOT/$DST/conditions.csv" --name "$DST" \
      --start "$START" --count "$PER" --dc-alpha-p "$ALPHA" > "$OUT/w${i}.log" 2>&1 &
done
wait
echo "all workers done"
