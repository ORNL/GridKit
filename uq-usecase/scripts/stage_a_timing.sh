#!/bin/bash
# Stage-A timing: is the per-worker fixed cost in the campaigns Julia launch
# overhead, or contention from 16 workers starting at once?
# Runs four experiments on one node: one process alone at 4 and 52 dispatches,
# then 16 concurrent processes at 4 and 52 dispatches each. Each process
# prints its launches (count, cases, wall) to a JSON file.
# usage: sbatch stage_a_timing.sh
#SBATCH --job-name=stageA_timing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --account=scidac
#SBATCH --partition=short
#SBATCH --output=/scratch/isatkaus/gridkit-timing/slurm-%j.out

set -u
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
S=/home/isatkaus/gridkit/uq-usecase/scripts
DATA=/kfs2/projects/scidac/scidac-data/gridkit-campaigns
ROOT=/scratch/isatkaus/gridkit-timing
OUT=$ROOT/results
mkdir -p "$OUT"
FAULT=$DATA/faultscreen1_ap08_qlim/conditions.csv
ALEA=$DATA/aleatoric1_ap08_qlim/conditions.csv
echo "host=$(hostname) start=$(date -Is)"

# 1. one process alone, fault-sweep shape (4 dispatches)
$PY "$S/stage_a_timing.py" --conditions "$FAULT" --start 0 --n 4 --out "$OUT/alone_fault_n4.json"
# 2. one process alone, aleatoric shape (52 dispatches)
$PY "$S/stage_a_timing.py" --conditions "$ALEA" --start 0 --n 52 --out "$OUT/alone_alea_n52.json"

# 3. 16 concurrent processes, 4 dispatches each (fault-sweep shape)
for i in $(seq 0 15); do
  $PY "$S/stage_a_timing.py" --conditions "$FAULT" --start $((i * 4)) --n 4 \
      --out "$OUT/conc_fault_n4_w$i.json" > "$OUT/conc_fault_n4_w$i.log" 2>&1 &
done
wait

# 4. 16 concurrent processes, 52 dispatches each (aleatoric shape)
for i in $(seq 0 15); do
  $PY "$S/stage_a_timing.py" --conditions "$ALEA" --start $((i * 52)) --n 52 \
      --out "$OUT/conc_alea_n52_w$i.json" > "$OUT/conc_alea_n52_w$i.log" 2>&1 &
done
wait
echo "done $(date -Is)"
